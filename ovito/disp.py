#!/usr/bin/env python3
"""
Render displacement vectors stored in an XYZ file with OVITO's Python API.

Example
-------
ovitos ovito/disp.py \
    --input /Users/kyou/Library/CloudStorage/Box-Box/output/Al/vacancy/eam_alloy2004/disp_weight_0/dump/displacement_ovit/displacement_ovit0.xyz \
    --output ./displacement_ovit0.png
"""

from __future__ import annotations

import argparse
import os
from math import radians
from typing import Sequence, Tuple

from ovito.io import import_file
from ovito.modifiers import ColorCodingModifier, ComputePropertyModifier
from ovito.vis import OpenGLRenderer, OSPRayRenderer, ParticlesVis, Viewport


def _vector(value: str) -> Tuple[float, float, float]:
    """Parse comma separated xyz arguments."""
    try:
        parts = [float(v) for v in value.split(",")]
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"Vector must be 'x,y,z': {value}") from exc
    if len(parts) != 3:
        raise argparse.ArgumentTypeError(f"Vector must have 3 components: {value}")
    return tuple(parts)  # type: ignore


def _color(value: str) -> Tuple[float, float, float]:
    vec = _vector(value)
    for c in vec:
        if c < 0.0 or c > 1.0:
            raise argparse.ArgumentTypeError(f"Color components must be in [0,1]: {value}")
    return vec


def _build_pipeline(input_path: str):
    """Load the displacement file and add helpful modifiers."""
    node = import_file(input_path, multiple_frames=False)
    node.modifiers.append(
        ComputePropertyModifier(
            expressions=["sqrt(disp.X*disp.X + disp.Y*disp.Y + disp.Z*disp.Z)"],
            output_property="disp_mag",
        )
    )
    color_mod = ColorCodingModifier(property="disp_mag", gradient=ColorCodingModifier.Rainbow())
    node.modifiers.append(color_mod)
    return node


def _setup_display(node, sphere_radius: float, vector_color, vector_scale: float, vector_width: float):
    """Configure particle and vector visualization."""
    scene_obj = node.add_to_scene()
    try:
        pvis: ParticlesVis = scene_obj.display
    except Exception:
        return scene_obj

    try:
        pvis.shape = ParticlesVis.Shape.Sphere
    except Exception:
        pass
    try:
        pvis.radius = sphere_radius
    except Exception:
        pass
    try:
        pvis.use_particle_radii = False
    except Exception:
        pass
    try:
        pvis.vector_display.enabled = True
        pvis.vector_display.vector_property = "disp"
        pvis.vector_display.color = vector_color
        pvis.vector_display.scale = vector_scale
        pvis.vector_display.line_width = vector_width
        pvis.vector_display.arrow_size = vector_width * 3.0
    except Exception:
        pass
    return scene_obj


def _render(node, output_path: str, size: Sequence[int], camera_pos, camera_dir, fov_deg: float, zoom: float):
    """Render a still image."""
    vp = Viewport(type=Viewport.Type.Perspective)
    try:
        vp.zoom_all()
    except Exception:
        pass
    try:
        vp.camera_pos = camera_pos
        vp.camera_dir = camera_dir
        vp.fov = radians(fov_deg)
    except Exception:
        pass
    try:
        vp.zoom(zoom)
    except Exception:
        pass

    try:
        renderer = OSPRayRenderer()
    except Exception:
        renderer = OpenGLRenderer()
    try:
        renderer.background_color = (1.0, 1.0, 1.0)
        renderer.antialiasing_level = 2
        renderer.ambient_occlusion_enabled = True
        renderer.ambient_occlusion_intensity = 0.6
    except Exception:
        pass

    os.makedirs(os.path.dirname(os.path.abspath(output_path)) or ".", exist_ok=True)
    vp.render_image(filename=output_path, size=size, renderer=renderer)


def main():
    parser = argparse.ArgumentParser(description="Visualize displacement vectors with OVITO.")
    default_input = "/Users/kyou/Library/CloudStorage/Box-Box/output/Al/vacancy/eam_alloy2004/disp_weight_0/dump/displacement_ovit/displacement_ovit0.xyz"
    parser.add_argument("--input", type=str, default=default_input, help="XYZ file with disp property.")
    parser.add_argument("--output", type=str, default="./displacement.png", help="Rendered PNG path.")
    parser.add_argument("--width", type=int, default=1600, help="Image width in pixels.")
    parser.add_argument("--height", type=int, default=1200, help="Image height in pixels.")
    parser.add_argument("--sphere-radius", type=float, default=0.6, help="Particle sphere radius.")
    parser.add_argument("--vector-color", type=_color, default=(0.1, 0.1, 0.1), help="Vector RGB color (0-1).")
    parser.add_argument("--vector-scale", type=float, default=8.0, help="Multiplicative arrow scale.")
    parser.add_argument("--vector-width", type=float, default=0.06, help="Vector line width.")
    parser.add_argument("--camera-pos", type=_vector, default=(16.0, 14.0, 20.0), help="Camera position 'x,y,z'.")
    parser.add_argument("--camera-dir", type=_vector, default=(-0.6, -0.4, -0.7), help="Camera direction unit vector.")
    parser.add_argument("--fov", type=float, default=32.0, help="Field of view in degrees.")
    parser.add_argument("--zoom", type=float, default=1.0, help="Extra zoom factor after zoom_all.")

    args = parser.parse_args()

    node = _build_pipeline(args.input)
    _setup_display(node, args.sphere_radius, args.vector_color, args.vector_scale, args.vector_width)
    _render(node, args.output, (args.width, args.height), args.camera_pos, args.camera_dir, args.fov, args.zoom)
    print(f"[OK] Rendered {args.output}")


if __name__ == "__main__":
    main()
