#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
目的: Vsi-4Hの構造作成,qeで残留応力法によってelastic dipoleの算出

SiCのLAMMPS dataファイル(in.sicq)を読み込み、
V_Si-4H (シリコン空孔 + 4つの水素終端) 構造を作成して、
Quantum Espresso (QE) 用の入力ファイルを出力するスクリプト。
"""

import math
import os
import sys
import numpy as np

######## 環境変数の取得 ########
output_dir = os.environ.get("OUTPUT_PATH")
if output_dir is None:
    raise ValueError("環境変数 OUTPUT_PATH error")

######## 削除する原子の指定 ######## 
atom_type_to_delete = int(os.environ.get("ATOM"))
if atom_type_to_delete is None:
    raise ValueError("環境変数 ATOM error")

######## 定数 ########
C_H_BOND_LENGTH = 1.09  # C-H 結合長 [A]
target_element_type = 1

# SiC格子定数 (qeで最小原子数で格子定数を求めたものを採用)
a1 = 3.093600904
a2 = 5.356555526
a3 = 10.124684058

####### 単位胞のデータ(4H-SiC)(qeで最小原子数で緩和後) ########
base_atoms = [                                  #(六方晶site or 立方晶site) : target atom
    [1, 0.0000000000, -0.0000994450, 0.0015372068], # (k):0  Si
    [1, 1.5468004521,  2.6781782946, 0.0015371710], # (k):1  Si
    [1, 0.0000000000,  1.7856874137, 2.5321748613], # (h):2  Si
    [1, 1.5468004521,  4.4639651533, 2.5321748506], # (h):3  Si
    [1, 0.0000000000,  0.0000992163, 5.0638790676], # (k):4  Si
    [1, 1.5468004521,  2.6783769557, 5.0638790359], # (k):5  Si
    [1, 0.0000000000,  3.5708678565, 7.5945169444], # (h):6  Si
    [1, 1.5468004521,  0.8925900885, 7.5945170057], # (h):7  Si
    [2, 0.0000000000,  1.7857354326, 0.6282891090], # (k):8  C
    [2, 1.5468004521,  4.4640131655, 0.6282890842], # (k):9  C
    [2, 0.0000000000, -0.0000857488, 3.1659263841], # (h):10 C
    [2, 1.5468004521,  2.6781919900, 3.1659263928], # (h):11 C
    [2, 0.0000000000,  3.5708198446, 5.6906310383], # (k):12 C
    [2, 1.5468004521,  0.8925420940, 5.6906311183], # (k):13 C
    [2, 0.0000000000,  0.0000854903, 8.2282685407], # (h):14 C
    [2, 1.5468004521,  2.6783632474, 8.2282685188], # (h):15 C
]
#1: Si
#2: C

# 最近接距離の2乗 (Si-C間距離)
nnd_sq = (0.00-0.00)**2 + (1.776276-0.00)**2 + (0.628008 - 0.00)**2 # nearest neibor distance ** 2

######### includeファイルからスーパーセルの繰り返し数(nx, ny, nz)を取得 ########
def get_supercell_size(m):
    n = [0, 0, 0]
    include_file = f"{output_dir}/include/{m}_include_n"
    with open(include_file, 'r') as f:
        for i, line in enumerate(f):
            parts = line.strip().split()
            if len(parts) >= 4:
                n[i] = int(parts[3])
    return n

######## LAMMPS dataファイルを読み込み、Boxサイズと原子リストを返す ########　← sic.pyで生成されるのがlammpsのファイル形式だから
def read_lammps_data(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()

    atoms = []
    box = {}
    reading_atoms = False
    
    for line in lines:
        parts = line.strip().split()
        if not parts: continue
        
        if "xlo xhi" in line:
            box['x'] = float(parts[1]) - float(parts[0])
        elif "ylo yhi" in line:
            box['y'] = float(parts[1]) - float(parts[0])
        elif "zlo zhi" in line:
            box['z'] = float(parts[1]) - float(parts[0])
        elif "Atoms" in line:
            reading_atoms = True
            continue
        
        if reading_atoms:
            # format: id type q x y z
            if len(parts) >= 6 and parts[0].isdigit():
                atoms.append({
                    'id': int(parts[0]),
                    'type': int(parts[1]),
                    'q': float(parts[2]),
                    'x': float(parts[3]),
                    'y': float(parts[4]),
                    'z': float(parts[5])
                })
    return box, atoms

######## Quantum Espresso入力ファイルを作成 ########
def write_qe_input(filename, box, atoms, id_str, calculation='relax'):
    
    # 擬ポテンシャルファイル名
    pseudo_si = "Si.pbe-n-kjpaw_psl.1.0.0.UPF"
    pseudo_c = "C.pbe-n-kjpaw_psl.1.0.0.UPF"
    pseudo_h = "H.pbe-kjpaw_psl.1.0.0.UPF"

    with open(filename, 'w') as f:
        # --- CONTROL ---
        f.write("&CONTROL\n")
        f.write(f"  calculation = '{calculation}',\n") # 構造緩和
        f.write(f"  prefix = 'sic_vsi4h_{id_str}',\n")
        f.write("  outdir = './out/',\n")
        #f.write("  pseudo_dir = '/Users/kyou/q-e/pseudo/',\n")
        f.write("  pseudo_dir = '/home/kyou/q-e/pseudo/',\n")
        f.write("  tstress = .true.,\n")  # 重要: 残留応力法のため応力を計算
        f.write("  tprnfor = .true.,\n")
        f.write("/\n")
        
        # --- SYSTEM ---
        f.write("&SYSTEM\n")
        f.write("  ibrav = 0,\n") # 格子ベクトルを明示的に指定
        f.write(f"  nat = {len(atoms)}, ntyp = 3,\n") # Si, C, H
        f.write("  ecutwfc = 40.0, ecutrho = 320.0,\n")
        f.write("  occupations = 'smearing', smearing = 'gaussian', degauss = 0.01,\n")
        f.write("  nosym = .true.,\n")
        f.write("  noinv = .true.,\n")
        f.write("/\n")
        
        # --- ELECTRONS ---
        f.write("&ELECTRONS\n")
        f.write("  conv_thr = 1.0d-6,\n")
        f.write("  mixing_beta = 0.3,\n")
        f.write("/\n")
        
        # --- IONS ---
        f.write("&IONS\n")
        f.write("  ion_dynamics = 'bfgs',\n")
        f.write("/\n")
        
        # --- ATOMIC_SPECIES ---
        f.write("ATOMIC_SPECIES\n")
        f.write(f"  Si  28.0855  {pseudo_si}\n")
        f.write(f"  C   12.0107  {pseudo_c}\n")
        f.write(f"  H   1.00784  {pseudo_h}\n")
        
        # --- CELL_PARAMETERS ---
        f.write("CELL_PARAMETERS (angstrom)\n")
        f.write(f"  {box['x']:.8f} 0.00000000 0.00000000\n")
        f.write(f"  0.00000000 {box['y']:.8f} 0.00000000\n")
        f.write(f"  0.00000000 0.00000000 {box['z']:.8f}\n")
        
        # --- ATOMIC_POSITIONS ---
        f.write("ATOMIC_POSITIONS (angstrom)\n")
        for atom in atoms:
            # Type 1:Si, 2:C, 3:H
            element = "Si" if atom['type'] == 1 else ("C" if atom['type'] == 2 else "H")
            f.write(f"  {element:<3} {atom['x']:.8f} {atom['y']:.8f} {atom['z']:.8f}\n")
        
        # --- K_POINTS ---
        f.write("K_POINTS automatic\n")
        f.write("  1 1 1 0 0 0\n") # スーパーセルサイズに応じて調整推奨

if __name__ == "__main__":
    b_atom = base_atoms[atom_type_to_delete] 
    
    for m in range(2):
        n = get_supercell_size(m)
        if n is None: continue
        
        # 基準点 (削除中心) の計算
        # 偶数・奇数の補正
        nx, ny, nz = n[0], n[1], n[2]
        if nx % 2 != 0: nx -= 1
        if ny % 2 != 0: ny -= 1
        if nz % 2 != 0: nz -= 1
        
        ca = [a1 * (nx / 2), a2 * (ny / 2), a3 * (nz / 2)]
        
        # 削除ターゲット座標 
        xc = ca[0] + b_atom[1]
        yc = ca[1] + b_atom[2]
        zc = ca[2] + b_atom[3] 
        
        # 入力ファイルの読み込み
        input_filepath = f"{output_dir}/4hsic_q/{m}_in.sicq"
            
        box, atoms = read_lammps_data(input_filepath)
        
        # 1. 最も近いSi原子(削除対象)を探す
        closest_idx = -1
        min_dist2 = float('inf')
        
        for i, atom in enumerate(atoms):
            if atom['type'] != target_element_type: continue # ターゲットと違う原子以外はスキップ
            
            dx = atom['x'] - xc
            dy = atom['y'] - yc
            dz = atom['z'] - zc
            d2 = dx**2 + dy**2 + dz**2
            
            if d2 < min_dist2:
                min_dist2 = d2
                closest_idx = i
                
        if closest_idx == -1:
            print(f"Error: Target Si atom not found in {m}")
            continue
            
        target_si = atoms[closest_idx]
        print(f"[{m}] Removing Si at: {target_si['x']:.3f}, {target_si['y']:.3f}, {target_si['z']:.3f}")
        
        # 2. 空孔周辺のC原子(Type 2)を探し、Hを追加する
        new_atoms = [] # 新しい原子リスト
        h_count = 0
        
        # Siの座標 (空孔中心)
        v_x, v_y, v_z = target_si['x'], target_si['y'], target_si['z']
        
        for i, atom in enumerate(atoms):
            # 削除対象のSiはリストに加えない
            if i == closest_idx: continue
            
            # 元の原子を追加
            new_atoms.append(atom)
            
            # C原子の場合、隣接チェック
            if atom['type'] == 2:
                dx = atom['x'] - v_x
                dy = atom['y'] - v_y
                dz = atom['z'] - v_z
                
                if dx > box['x'] / 2: dx -= box['x']
                elif dx < -box['x'] / 2: dx += box['x']
                
                if dy > box['y'] / 2: dy -= box['y']
                elif dy < -box['y'] / 2: dy += box['y']
                
                if dz > box['z'] / 2: dz -= box['z']
                elif dz < -box['z'] / 2: dz += box['z']
                d2 = dx**2 + dy**2 + dz**2
                
                # 第一近接距離にあるか判定 (許容誤差 0.1程度)
                if abs(d2 - nnd_sq) < 0.2: 
                    dist = math.sqrt(d2)
                    
                    hx = atom['x'] + (-dx / dist) * C_H_BOND_LENGTH
                    hy = atom['y'] + (-dy / dist) * C_H_BOND_LENGTH
                    hz = atom['z'] + (-dz / dist) * C_H_BOND_LENGTH
                    
                    # 水素原子を追加 (Type 3)
                    new_atoms.append({
                        'id': -1, # 後で採番
                        'type': 3,
                        'q': 0.0,
                        'x': hx, 'y': hy, 'z': hz
                    })
                    h_count += 1
        
        print(f"    Added {h_count} H atoms.")
        print(f" total atoms {len(new_atoms)}")
        with open("./total_atoms", "a") as f:
            f.write(f"total atoms {len(new_atoms)}\n")
        
        # 3. QE入力ファイルの出力
        output_filename = f"{output_dir}/qe_vsi4h_{m}.in"
        write_qe_input(output_filename, box, new_atoms, id_str=str(m))
        print(f"    Generated: {output_filename}")