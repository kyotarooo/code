import numpy as np
import sys

# ================= 設定エリア =================
input_filename = "/Users/kyou/output/4HSiC/vacancy/vashishta-k/residual_stress/dump/defect_lattice/defect_lattice_0/stress_viz_8.dump"
output_filename = "stress_sic_viz.xyz"  # 拡張子を .xyz にします

# 単位変換係数 (bar -> GPa)
conversion_factor = 0.0001
# ==============================================

def process_to_xyz(in_file, out_file):
    print(f"読み込み中: {in_file}")
    
    with open(in_file, 'r') as f:
        lines = f.readlines()

    # --- 1. ボックスサイズ(Lattice)の取得 ---
    # LAMMPS dumpの BOX BOUNDS を探して格子ベクトルを作ります
    box_bounds = []
    start_idx = 0
    col_map = {}
    
    for i, line in enumerate(lines):
        if line.startswith("ITEM: BOX BOUNDS"):
            # 次の3行が x, y, z の範囲
            box_lines = lines[i+1 : i+4]
            for bl in box_lines:
                # "xlo xhi" を取得
                parts = bl.strip().split()
                # 周期境界のフラグ(pp pp pp)等を無視して数値だけ取る
                vals = [float(x) for x in parts if x.replace('.','',1).replace('e','',1).replace('+','',1).replace('-','',1).isdigit()]
                box_bounds.append(vals)
        
        if line.startswith("ITEM: ATOMS"):
            headers = line.strip().split()[2:]
            for idx, name in enumerate(headers):
                col_map[name] = idx
            start_idx = i + 1
            break
            
    # 直交箱(Orthogonal box)と仮定してLattice文字列を作成
    # Lattice="ax 0 0 0 ay 0 0 0 az"
    Lx = box_bounds[0][1] - box_bounds[0][0]
    Ly = box_bounds[1][1] - box_bounds[1][0]
    Lz = box_bounds[2][1] - box_bounds[2][0]
    lattice_str = f'Lattice="{Lx} 0.0 0.0 0.0 {Ly} 0.0 0.0 0.0 {Lz}"'

    # --- 2. データ読み込み ---
    data_block = lines[start_idx:]
    data_array = np.loadtxt(data_block)
    n_atoms = data_array.shape[0]
    print(f"原子数: {n_atoms}")

    # --- 3. 必要なデータの抽出 ---
    # 座標 (x, y, z)
    try:
        idx_x = col_map['x']
        idx_y = col_map['y']
        idx_z = col_map['z']
        pos = data_array[:, [idx_x, idx_y, idx_z]]
        
        # タイプ (type) -> Si, C への変換用
        idx_type = col_map['type']
        atom_types = data_array[:, idx_type].astype(int)
        
        # 応力 (c_st)
        st_cols = []
        for k in range(1, 7):
            st_cols.append(col_map[f"c_st[{k}]"])
        st_raw = data_array[:, st_cols]

        # 体積 (c_voro)
        vol_idx = -1
        possible_names = ["c_voro[1]", "voro[1]", "c_voro"]
        for name in possible_names:
            if name in col_map:
                vol_idx = col_map[name]
                break
        
        if vol_idx == -1:
            raise ValueError("体積データ(c_voro)が見つかりません")
            
        vol = data_array[:, vol_idx]
        
    except Exception as e:
        print(f"❌ データ抽出エラー: {e}")
        return

    # --- 4. 応力計算 (GPa) ---
    print("応力計算中...")
    with np.errstate(divide='ignore', invalid='ignore'):
        stress_gpa = (st_raw / vol[:, np.newaxis]) * conversion_factor
        stress_gpa[~np.isfinite(stress_gpa)] = 0.0

    # Von Mises
    s_xx, s_yy, s_zz = stress_gpa[:, 0], stress_gpa[:, 1], stress_gpa[:, 2]
    s_xy, s_xz, s_yz = stress_gpa[:, 3], stress_gpa[:, 4], stress_gpa[:, 5]
    von_mises = np.sqrt(0.5 * ((s_xx-s_yy)**2 + (s_yy-s_zz)**2 + (s_zz-s_xx)**2 + 6*(s_xy**2 + s_yz**2 + s_xz**2)))

    # --- 5. Extended XYZ形式での書き出し ---
    print(f"書き出し中: {out_file}")
    
    with open(out_file, 'w') as f:
        # 1行目: 原子数
        f.write(f"{n_atoms}\n")
        
        # 2行目: コメント行 (Properties定義 + Lattice)
        # 形式: species:S:1:pos:R:3:von_mises:R:1:stress_tensor:R:6
        props = "Properties=species:S:1:pos:R:3:VonMises_GPa:R:1:Stress_Tensor_GPa:R:6"
        f.write(f"{props} {lattice_str} Origin=\"{box_bounds[0][0]} {box_bounds[1][0]} {box_bounds[2][0]}\"\n")
        
        # 3行目以降: データ
        for i in range(n_atoms):
            # 原子タイプの判定 (LAMMPS: 1=Si, 2=C と仮定)
            species = "Si" if atom_types[i] == 1 else "C"
            if atom_types[i] > 2: species = "X" # 予期せぬタイプ
            
            # 書き込み
            # species x y z vm sxx syy szz sxy sxz syz
            p = pos[i]
            vm = von_mises[i]
            s = stress_gpa[i]
            
            f.write(f"{species} {p[0]:.5f} {p[1]:.5f} {p[2]:.5f} {vm:.5f} "
                    f"{s[0]:.5e} {s[1]:.5e} {s[2]:.5e} {s[3]:.5e} {s[4]:.5e} {s[5]:.5e}\n")

    print("--------------------------------------------------")
    print("✅ 完了しました！")
    print(f"出力ファイル: {output_filename}")
    print("OVITOでこのファイルを開くと、自動的に Si/C として認識されます。")
    print("--------------------------------------------------")

if __name__ == "__main__":
    process_to_xyz(input_filename, output_filename)