"""
目的: 最小構造(計算コストを抑えるため)を読み込み、
      格子定数最適化用のQE入力ファイルを作成する。
"""
import os

# -------- 設定 --------
# sic.pyで出力したディレクトリとID
output_dir = os.environ.get("OUTPUT_PATH", "./sim_data") # 環境変数がない場合のデフォルト
lammps_file = f"{output_dir}/4hsic_q/minimum.sicq" # [1,1,1]で生成されたファイル
qe_input_file = "sic_unitcell_relax.in"

# 擬ポテンシャル (vsi4hgenと同じものを使用)
pseudo_si = "Si.pbe-n-kjpaw_psl.1.0.0.UPF"
pseudo_c = "C.pbe-n-kjpaw_psl.1.0.0.UPF"
# ---------------------

def read_lammps_box_atoms(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()
    atoms = []
    box = {}
    reading_atoms = False
    for line in lines:
        parts = line.strip().split()
        if not parts: continue
        if "xlo xhi" in line: box['x'] = float(parts[1]) - float(parts[0])
        elif "ylo yhi" in line: box['y'] = float(parts[1]) - float(parts[0])
        elif "zlo zhi" in line: box['z'] = float(parts[1]) - float(parts[0])
        elif "Atoms" in line: reading_atoms = True; continue
        if reading_atoms and len(parts) >= 6 and parts[0].isdigit():
            atoms.append({
                'type': int(parts[1]),
                'x': float(parts[3]), 'y': float(parts[4]), 'z': float(parts[5])
            })
    return box, atoms

def write_vcrelax_input(filename, box, atoms):
    with open(filename, 'w') as f:
        f.write("&CONTROL\n")
        f.write("  calculation = 'vc-relax',\n")  # セルサイズも最適化
        f.write("  prefix = 'sic_bulk',\n")
        f.write("  outdir = './tmp/',\n")
        f.write("  pseudo_dir = './pseudo/',\n")
        f.write("  tstress = .true.,\n")
        f.write("  tprnfor = .true.,\n")
        f.write("  disk_io = 'low',\n")
        f.write("/\n")
        
        f.write("&SYSTEM\n")
        f.write("  ibrav = 0,\n")
        f.write(f"  nat = {len(atoms)}, ntyp = 2,\n")
        f.write("  ecutwfc = 40.0, ecutrho = 320.0,\n") # SSSP Efficiency
        f.write("  occupations = 'smearing', smearing = 'gaussian', degauss = 0.01,\n")
        f.write("/\n")
        
        f.write("&ELECTRONS\n")
        f.write("  conv_thr = 1.0d-8,\n") # 厳しめに設定
        f.write("  mixing_beta = 0.3,\n")
        f.write("/\n")
        
        f.write("&IONS\n")
        f.write("  ion_dynamics = 'bfgs',\n")
        f.write("/\n")
        
        f.write("&CELL\n")
        f.write("  cell_dynamics = 'bfgs',\n") # セルの動き方
        f.write("  press_conv_thr = 0.1,\n")   # 圧力の収束条件 (Kbar)
        f.write("  cell_dofree = 'all',\n")    # 全方向に緩和 (直交性を保ちたい場合は 'xyz' など検討)
        f.write("/\n")
        
        f.write("ATOMIC_SPECIES\n")
        f.write(f"  Si  28.0855  {pseudo_si}\n")
        f.write(f"  C   12.0107  {pseudo_c}\n")
        
        f.write("CELL_PARAMETERS (angstrom)\n")
        f.write(f"  {box['x']:.8f} 0.00000000 0.00000000\n")
        f.write(f"  0.00000000 {box['y']:.8f} 0.00000000\n")
        f.write(f"  0.00000000 0.00000000 {box['z']:.8f}\n")
        
        f.write("ATOMIC_POSITIONS (angstrom)\n")
        for atom in atoms:
            elem = "Si" if atom['type'] == 1 else "C"
            f.write(f"  {elem:<3} {atom['x']:.8f} {atom['y']:.8f} {atom['z']:.8f}\n")
            
        f.write("K_POINTS automatic\n")
        # ユニットセルなのでK点は密にする必要があります
        f.write("  6 6 4 0 0 0\n") 

if __name__ == "__main__":
    if not os.path.exists(lammps_file):
        print(f"Error: {lammps_file} が見つかりません.")
    else:
        box, atoms = read_lammps_box_atoms(lammps_file)
        write_vcrelax_input(qe_input_file, box, atoms)
        print(f"Generated: {qe_input_file}")