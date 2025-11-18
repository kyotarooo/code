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

######## 削除する原子の指定 ########  どの位置の原子を削除するか ←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←
atom_type_to_delete = int(os.environ.get("ATOM"))
if atom_type_to_delete is None:
    raise ValueError("環境変数 ATOM error")

######## 繰り返し単位 ########
a1 = 3.07659964
a2 = 5.32882689
a3 = 2.51203309 * 4

######## ファイル名 ########
input_file = "in.sicq"        # 入力ファイル名
output_file = "sic_vac.data"   # 出力ファイル名

####### 単位胞のデータ(4H-SiC) ########
base_atoms = [
    [1, 0.000000, 0.000000, 0.000000, 0.000000],
    [1, 0.000000, 1.538300, 2.664413, 0.000000],
    [1, 0.000000, 0.000000, 1.776276, 2.512033],
    [1, 0.000000, 1.538300, 4.440689, 2.512033],
    [1, 0.000000, 0.000000, 0.000000, 5.024066],
    [1, 0.000000, 1.538300, 2.664413, 5.024066],
    [1, 0.000000, 0.000000, 3.552551, 7.536099],
    [1, 0.000000, 1.538300, 0.888138, 7.536099],
    [2, 0.000000, 0.000000, 1.776276, 0.628008],
    [2, 0.000000, 1.538300, 4.440689, 0.628008],
    [2, 0.000000, 0.000000, 0.000000, 3.140041],
    [2, 0.000000, 1.538300, 2.664413, 3.140041],
    [2, 0.000000, 0.000000, 3.552551, 5.652074],
    [2, 0.000000, 1.538300, 0.888138, 5.652074],
    [2, 0.000000, 0.000000, 0.000000, 8.164108],
    [2, 0.000000, 1.538300, 2.664413, 8.164108],
]

with open(f"{output_dir}/4hsic_vacancy/atom_type_delete_{atom_type_to_delete}/filename.txt" , "w") as f_f:
    for m in range(8):
        m = m + 1
        ####### 繰り返し個数の取得 ########
        with open(f"{output_dir}/include/{m}_include_n", 'r') as f_include:
            n = {}
            lines = f_include.readlines()
            for i, line in enumerate(lines):
                n[i] = int(line.strip().split()[3])
                
        ######## 計算領域 ######## 
        sca = {}
        sca[0] = a1 * int(n[0])
        sca[1] = a2 * int(n[1])
        sca[2] = a3 * int(n[2])
        
        with open(f"{output_dir}/supercell.txt" , "a") as f_supercell:
            f_supercell.write(f"{sca[0]} {sca[1]} {sca[2]}\n")
            
        # 偶数と奇数で場合分け        
        for i in range(len(n)):
            if n[i] % 2 != 0:
                n[i] = n[i] - 1
                
        ######## 基準点　########
        ca = {}
        ca[0] = a1 * int(n[0] / 2)
        ca[1] = a2 * int(n[1] / 2)
        ca[2] = a3 * int(n[2] / 2)
        

        ######## 座標データファイルの読み込み ########
        with open(f"{output_dir}/4hsic_q/{m}_{input_file}", "r") as f:
            lines = f.readlines()

            # Atoms セクションを見つける
            atoms_start_idx = None
            for i, line in enumerate(lines):
                if line.strip().lower().startswith("atoms"):
                    atoms_start_idx = i
                    break

            # データ開始位置（空行をスキップ）
            data_start = atoms_start_idx + 1
            while data_start < len(lines) and lines[data_start].strip() == "":
                data_start += 1
                
            # 削除する原子の座標
            xc, yc, zc = ca[0] + base_atoms[atom_type_to_delete][2], ca[1] + base_atoms[atom_type_to_delete][3], ca[2] + base_atoms[atom_type_to_delete][4]

            # 座標データの読み込み開始 
            atom_lines = []
            for i in range(data_start, len(lines)):
                line = lines[i]
                if line.strip() == "":
                    break
                atom_lines.append(line)

            # 指定した座標に最も近い原子を探す
            closest_idx = None
            closest_dist2 = 1000000

            for idx, line in enumerate(atom_lines):
                parts = line.strip().split()
                atom_id = int(parts[0])
                atom_type = int(parts[1])
                x, y, z = float(parts[3]), float(parts[4]), float(parts[5])
                d2 = (x - xc)**2 + (y - yc)**2 + (z - zc)**2
                
                if  d2 < closest_dist2:
                    closest_dist2 = d2
                    closest_idx = idx

            deleted_line = atom_lines[closest_idx]
            print(closest_dist2)
            print("削除する原子:", deleted_line.strip())
            # with open(f"{output_dir}/dump/deleted_atom/deleted_atom_{atom_type_to_delete}/deleted_atom{m}", "w") as f_delete:
            #     f_delete.write(deleted_line.strip())
            
            # 新しい座標データの作成開始
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
            
            f_f.write(f"{n[0]}x{n[1]}x{n[2]}_{output_file} ")

            # 出力
            with open(f"{output_dir}/4hsic_vacancy/atom_type_delete_{atom_type_to_delete}/{m}_{output_file}", "w") as f:
                f.writelines(new_lines)

            print(f"出力ファイル: {m}.{atom_type_to_delete}_{output_file}")
