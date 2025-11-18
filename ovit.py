#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
OVITOで4H-SiC構造（LAMMPS data 等）を読み込み、
・Si/Cのタイプ名と色を設定
・適切なカットオフでSi–C結合線を自動生成
・球表示/きれいなライティング（AO）
・高解像度PNG出力
を行うスクリプト。

★ 使い方（例）
    python3 render_4HSiC.py \
        --input /mnt/data/1_sic_vac.data \
        --output ./render_4HSiC.png \
        --cutoff 2.1 \
        --radius 0.35 \
        --dpi 300

※ OVITO 3.x系のPython APIを想定。バージョン差による属性名の違いにも可能な範囲で対応。
"""

import argparse
import os
from math import radians
from typing import Dict, Optional


def _vector3(value: str):
    """Parse comma-separated vector strings passed via CLI."""
    parts = value.split(",")
    if len(parts) != 3:
        raise argparse.ArgumentTypeError("Vector must have three comma-separated values (e.g. 1,0,0).")
    try:
        return tuple(float(p) for p in parts)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"Invalid float in vector: {value}") from exc

# OVITO API
from ovito.io import import_file
from ovito.modifiers import CreateBondsModifier, PythonScriptModifier
from ovito.vis import Viewport, OpenGLRenderer, OSPRayRenderer, ParticlesVis


def _get_particle_types_container(data):
    """Return the particle types container, handling API naming variants."""
    particles = getattr(data, "particles_", None)
    if particles is None:
        particles = getattr(data, "particles", None)
    if particles is None:
        return None
    types = getattr(particles, "particle_types_", None)
    if types is None:
        types = getattr(particles, "particle_types", None)
    return types


def _type_by_id(types, tid: int):
    """Access particle type by ID across OVITO API variants."""
    if types is None:
        return None
    for accessor in ("type_by_id_", "type_by_id", "type"):
        func = getattr(types, accessor, None)
        if func is None:
            continue
        try:
            return func(tid)
        except Exception:
            continue
    return None


def _type_label_modifier_factory(radius_map: Dict[int, Optional[float]]):
    """Create a PythonScriptModifier function bound to the requested radii."""
    clean_radius = {tid: r for tid, r in (radius_map or {}).items() if r is not None}

    def _apply(frame, data):
        """粒子タイプに名前/色/半径を設定。"""
        types = _get_particle_types_container(data)
        if types is None:
            return

        def _apply_to_type(tid: int, name: str, color):
            try:
                t = _type_by_id(types, tid)
                if t is None:
                    return
                t.name = name
                t.color = color
                desired_radius = clean_radius.get(tid)
                if desired_radius is not None:
                    try:
                        t.radius = desired_radius
                    except Exception:
                        pass
            except Exception:
                pass

        _apply_to_type(1, "Si", (0.0, 0.65, 1.0))  # ブルー
        _apply_to_type(2, "C",  (0.50, 1.00, 0.83))   # オレンジ


    return _apply


def build_pipeline(input_path: str, cutoff: float, radius_map: Dict[int, Optional[float]]):
    """OVITOのPipelineを構築して返す。"""
    # LAMMPS data 読み込み（SiCは charge スタイル）
    node = import_file(input_path, atom_style='charge', multiple_frames=False)

    # 粒子タイプ名/色/半径を付与
    label_modifier = PythonScriptModifier(function=_type_label_modifier_factory(radius_map))
    node.modifiers.append(label_modifier)

    # --- Si–C 結合生成（タイプごとのカットオフ設定）---
    bond_mod = CreateBondsModifier()
    try:
        mode_value = CreateBondsModifier.Mode.Pairwise  # 新しめの名称
    except AttributeError:
        mode_value = CreateBondsModifier.Mode.Pairing   # 旧名称にフォールバック
    bond_mod.mode = mode_value

    # 1=Si, 2=C の想定
    bond_mod.set_pairwise_cutoff(1, 2, cutoff)  # Si–C のみ結合
    bond_mod.set_pairwise_cutoff(1, 1, 0.0)     # Si–Si 無効
    bond_mod.set_pairwise_cutoff(2, 2, 0.0)     # C–C 無効

    node.modifiers.append(bond_mod)

    return node


def _set_bond_vis_width(bonds_vis, bond_width: float) -> bool:
    """Helper to set bond width on various OVITO bond vis objects."""
    if bonds_vis is None:
        return False
    try:
        bonds_vis.enabled = True
    except Exception:
        pass
    for attr in ("width", "line_width"):
        if hasattr(bonds_vis, attr):
            try:
                setattr(bonds_vis, attr, bond_width)
                try:
                    bonds_vis.transparency = 0.0
                except Exception:
                    pass
                return True
            except Exception:
                continue
    return False


def _configure_bond_visual(node, scene_obj, bond_width: float):
    """Try multiple ways to apply the requested bond width."""
    # 1) Scene display object (typical case after add_to_scene)
    try:
        if _set_bond_vis_width(scene_obj.display.bonds, bond_width):
            return
    except Exception:
        pass

    # 2) Directに pipeline.source.data.particles.bonds.vis を触る（ユーザー例）
    try:
        source_particles = getattr(node.source.data, "particles", None)
        bonds_vis = getattr(getattr(source_particles, "bonds", None), "vis", None)
        if _set_bond_vis_width(bonds_vis, bond_width):
            return
    except Exception:
        pass

    # 3) モディファイア適用後のデータを一度評価し、そのvisを設定
    try:
        evaluated = node.compute()
        particles = getattr(evaluated, "particles", None)
        bonds_vis = getattr(getattr(particles, "bonds", None), "vis", None)
        _set_bond_vis_width(bonds_vis, bond_width)
    except Exception:
        pass


def setup_visual(node, bond_width: float, specular: float):
    """見栄えの設定（球サイズ/形状/ AO / 背景 等）。"""
    scene_obj = node.add_to_scene()
   
    # 粒子の見た目
    try:
        pvis: ParticlesVis = scene_obj.display
    except Exception:
        pvis = None

    if pvis is not None:
        try:
            pvis.shape = ParticlesVis.Shape.Sphere
        except Exception:
            pass
        try:
            pvis.use_particle_types = True
        except Exception:
            pass
        try:
            pvis.use_particle_radii = True
        except Exception:
            pass
        try:
            pvis.radius = 0.0
        except Exception:
            pass
        try:
            pvis.shading = ParticlesVis.Shading.Smooth
        except Exception:
            try:
                pvis.shading = ParticlesVis.Shading.Phong
            except Exception:
                pass
        try:
            material = pvis.material
            if material is not None:
                try:
                    material.specular_intensity = specular
                except Exception:
                    pass
        except Exception:
            pass
        try:
            pvis.outline_width = 0.0
        except Exception:
            pass

    # Bonds設定
    _configure_bond_visual(node, scene_obj, bond_width)

    # ======================================================
    # ★ シミュレーションセルの黒枠を完全に消す処理を追加
    # ======================================================
    try:
        cell_vis = node.source.data.cell.vis
        cell_vis.enabled = False
    except Exception:
        pass

    return scene_obj



def _clamp_specular_brightness(value: float) -> float:
    return max(0.0, min(1.0, value))


def render_image(node, output_path: str, width: int, height: int, dpi: int,
                camera_dir, camera_pos, zoom_factor: float, fov_deg: float,
                specular_brightness: float):
    vp = Viewport(type=Viewport.Type.Perspective)
    vp.zoom_all()
    if zoom_factor and abs(zoom_factor - 1.0) > 1e-4:
        try:
            vp.zoom(zoom_factor)
        except Exception:
            pass

    # --- ここを追加 ---
    vp.camera_dir = camera_dir  # 方向（例: a軸側）
    vp.camera_pos = camera_pos
    vp.fov = radians(fov_deg)
    # -----------------

    # レンダラ（OSPRay優先、無ければOpenGL）
    renderer = None
    try:
        renderer = OSPRayRenderer()
    except Exception:
        renderer = OpenGLRenderer()
    try:
        renderer.antialiasing_level = 2
    except Exception:
        pass
    try:
        renderer.shadows_enabled = False
    except Exception:
        pass
    try:
        renderer.ambient_occlusion_enabled = True
        renderer.ambient_occlusion_intensity = 0.6
        renderer.ambient_occlusion_samples = 64
    except Exception:
        pass

    # OSPRay固有のプリンシプルスペキュラー設定
    try:
        renderer.principled_specular_brightness = _clamp_specular_brightness(specular_brightness)
    except Exception:
        pass

    # 背景は白
    try:
        renderer.background_color = (1.0, 1.0, 1.0)
    except Exception:
        pass

    # 画像サイズ（物理解像度を意識するならdpiも付けて保存）
    # OVITOはdpiをメタデータ保存に使うだけなので、ピクセル数も適切に。
    vp.render_image(
        filename=output_path,
        size=(width, height),
        renderer=renderer,
    )

    # PNGメタとしてdpiを反映したい場合（Pillowがあれば上書き保存）
    try:
        from PIL import Image  # type: ignore
        im = Image.open(output_path)
        im.save(output_path, dpi=(dpi, dpi))
    except Exception:
        pass


def main():
    parser = argparse.ArgumentParser(description="Render 4H-SiC with OVITO")
    parser.add_argument("--input", type=str, default='/Users/kyou/Library/CloudStorage/Box-Box/output/4HSiC/vacancy/vashishta-k/disp/4hsic_q/2.sic',
                        help="入力構造ファイル（LAMMPS data など）")
    parser.add_argument("--output", type=str, default="./render_4HSiC.png",
                        help="出力PNGパス")
    parser.add_argument("--cutoff", type=float, default=0.001,
                        help="結合生成のカットオフ [Å]")
    parser.add_argument("--radius", type=float, default=None,
                        help="粒子半径（タイプ未指定の場合のデフォルト）")
    parser.add_argument("--radius-si", type=float, default=0.3,
                        help="Si の粒子半径（未指定なら --radius を利用）")
    parser.add_argument("--radius-c", type=float, default=0.1,
                        help="C の粒子半径（未指定なら --radius を利用）")
    parser.add_argument("--bondwidth", type=float, default=0.11,
                        help="結合線の太さ")
    parser.add_argument("--width", type=int, default=2400,
                        help="画像幅[px]")
    parser.add_argument("--height", type=int, default=1800,
                        help="画像高さ[px]")
    parser.add_argument("--dpi", type=int, default=300,
                        help="PNGメタデータとして埋め込むDPI")
    parser.add_argument("--zoom", type=float, default=1.0,
                        help="zoom_all 後に適用する倍率 (<1でズームイン, >1でズームアウト)")
    parser.add_argument("--camera-dir", type=_vector3, default=(0, 0, -1.0),
                        help="視線方向（\"x,y,z\"）")
    parser.add_argument("--camera-pos", type=_vector3, default=(5, 5, 20.0),
                        help="カメラ位置（\"x,y,z\"）")
    parser.add_argument("--fov", type=float, default=35.0,
                        help="視野角[deg]")
    parser.add_argument("--specular", type=float, default=0.80,
                        help="球のハイライト強度")
    parser.add_argument("--shininess", type=float, default=0.8,
                        help="OSPRay principled_specular_brightness（0〜1を推奨）")

    args = parser.parse_args()

    # 粒子タイプ別半径を設定（指定がなければデフォルト値）
    default_radius = args.radius
    radius_map = {
        1: args.radius_si if args.radius_si is not None else default_radius,
        2: args.radius_c if args.radius_c is not None else default_radius,
    }

    # Pipeline構築
    node = build_pipeline(args.input, args.cutoff, radius_map)
    setup_visual(node, bond_width=args.bondwidth, specular=args.specular)

    # レンダリング
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    render_image(node, args.output, args.width, args.height, args.dpi,
                 camera_dir=args.camera_dir, camera_pos=args.camera_pos,
                 zoom_factor=args.zoom, fov_deg=args.fov,
                 specular_brightness=args.shininess)

    print(f"[OK] Rendered → {args.output}")


if __name__ == "__main__":
    main()
