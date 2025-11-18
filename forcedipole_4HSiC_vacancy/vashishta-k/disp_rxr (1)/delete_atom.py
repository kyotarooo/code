#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
SiCのLAMMPS dataファイルから、中心に最も近い原子を1つ削除して空孔を作るスクリプト。
削除する原子タイプ (1=Si, 2=C) を指定可能。
Atomsの形式: id type x y z q
"""

import math
import os

######## 環境変数の取得 ########
output_dir = os.environ.get("OUTPUT_PATH")
if output_dir is None:
    raise ValueError("環境変数 OUTPUT_PATH error")

############ Load input file count ############
with open(f"{output_dir}/4hsic_q/filename.txt", 'r') as f_num:
    lines = f_num.readlines()
    file_num = len(lines)
print(file_num)
# ======== 設定 ========
input_file = "in.sicq"        # 入力ファイル名
output_file = "sic_vac.data"   # 出力ファイル名
target_type = 2                # 削除したい原子タイプ (1=Si, 2=C)
# ========================
with open(f"{output_dir}/4hsic_vacancy/filename.txt" , "w") as f_f:
    for k in range(file_num):
        # ----------------------------------------------------------------------
        # ファイル読込
        # ----------------------------------------------------------------------
        with open(f"{output_dir}/4hsic_q/{k+1}_{input_file}", "r") as f:
            lines = f.readlines()

        # Atoms セクションを見つける
        atoms_start_idx = None
        for i, line in enumerate(lines):
            if line.strip().lower().startswith("atoms"):
                atoms_start_idx = i
                break

        if atoms_start_idx is None:
            raise RuntimeError("Atoms セクションが見つかりません")

        # データ開始位置（空行をスキップ）
        data_start = atoms_start_idx + 1
        while data_start < len(lines) and lines[data_start].strip() == "":
            data_start += 1

        # ----------------------------------------------------------------------
        # 境界を取得
        # ----------------------------------------------------------------------
        xlo = xhi = ylo = yhi = zlo = zhi = None
        for line in lines:
            parts = line.strip().split()
            if len(parts) == 4:
                if parts[2] == "xlo" and parts[3] == "xhi":
                    xlo, xhi = float(parts[0]), float(parts[1])
                elif parts[2] == "ylo" and parts[3] == "yhi":
                    ylo, yhi = float(parts[0]), float(parts[1])
                elif parts[2] == "zlo" and parts[3] == "zhi":
                    zlo, zhi = float(parts[0]), float(parts[1])

        if None in (xlo, xhi, ylo, yhi, zlo, zhi):
            raise RuntimeError("ボックス境界が読み取れません")

        xc, yc, zc = 0.5 * (xlo + xhi), 0.5 * (ylo + yhi), 0.5 * (zlo + zhi)

        # ----------------------------------------------------------------------
        # Atoms の行をパース
        # ----------------------------------------------------------------------
        atom_lines = []
        for i in range(data_start, len(lines)):
            line = lines[i]
            if line.strip() == "":
                break
            atom_lines.append(line)

        # 指定タイプの中で中心に最も近い原子を探す
        closest_idx = None
        closest_dist2 = 1000000

        for idx, line in enumerate(atom_lines):
            parts = line.strip().split()
            if len(parts) < 6:
                continue  # 不完全な行はスキップ

            atom_id = int(parts[0])
            atom_type = int(parts[1])
            x, y, z = float(parts[3]), float(parts[4]), float(parts[5])

            if atom_type != target_type:
                continue

            d2 = (x - xc)**2 + (y - yc)**2 + (z - zc)**2
            
            
            if  d2 < closest_dist2:
                closest_dist2 = d2
                closest_idx = idx
                

        if closest_idx is None:
            raise RuntimeError(f"タイプ {target_type} の原子が見つかりません")

        deleted_line = atom_lines[closest_idx]
        #print("削除する原子:", deleted_line.strip())

        # ----------------------------------------------------------------------
        # 新しいAtomsセクションを作成
        # ----------------------------------------------------------------------
        new_atom_lines = [l for i, l in enumerate(atom_lines) if i != closest_idx]

        # 原子数を1減らす
        new_lines = []
        for line in lines[:atoms_start_idx]:
            if line.strip().endswith("atoms"):
                parts = line.strip().split()
                if parts[0].isdigit():
                    n_old = int(parts[0])
                    n_new = n_old - 1
                    line = line.replace(str(n_old), str(n_new), 1)
            new_lines.append(line)
            
            
        # Atomsセクションを再構築
        new_lines.append(lines[atoms_start_idx])  # "Atoms" の行

        # "Atoms" の直後の空行を追加（フォーマット維持）
        idx_after_atoms = atoms_start_idx + 1
        while idx_after_atoms < len(lines) and lines[idx_after_atoms].strip() == "":
            new_lines.append(lines[idx_after_atoms])
            idx_after_atoms += 1

        # ---- IDを詰めながら原子リストを出力 ----
        new_id = 1
        for i, line in enumerate(new_atom_lines):
            parts = line.strip().split()
            parts[0] = str(new_id)  # 新しいIDを割り当て
            new_id += 1
            new_line = " ".join(parts) + "\n"
            new_lines.append(new_line)
        
        f_f.write(f"{i+1}_{output_file}\n")

        # 出力
        with open(f"{output_dir}/4hsic_vacancy/{k+1}_{output_file}", "w") as f:
            f.writelines(new_lines)

        print(f"\n中心近傍のタイプ {target_type} 原子を1個削除しました。")
        print(f"出力ファイル: {i+1}_{output_file}")
