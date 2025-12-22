import glob
import numpy as np
import os
from ase.io import read
from ase.io.espresso import write_espresso_in
from ase.build import make_supercell

# ================= 設定 =================
# 1. 基準となる「16原子ユニットセル」のファイル
#    (計算済みの sic_bulk_ref.out または .in を指定)
UNIT_CELL_FILE = 'sic_bulk_ref.out' 

# 2. 参照する「欠陥計算」の出力ファイル一覧
DEFECT_FILES = sorted(glob.glob("[0-8].out"))  # 0.out 〜 8.out にマッチ

# 3. QE入力ファイルの共通設定
#    (40Ryのrelax計算用)
input_data = {
    'control': {
        'calculation': 'relax',  # 構造最適化
        'prefix': 'sic_bulk',
        'outdir': './out_bulk/',
        'pseudo_dir': '/home/kyou/q-e/pseudo/', # ← 環境に合わせて確認！
        'tstress': True,
        'tprnfor': True,
    },
    'system': {
        'ibrav': 0,
        'ecutwfc': 40.0,
        'ecutrho': 320.0,
        'occupations': 'smearing',
        'smearing': 'gaussian',
        'degauss': 0.01,
        'nosym': True,   # 対称性を切って安全に
        'noinv': True,
    },
    'electrons': {
        'conv_thr': 1.0e-8,
        'mixing_beta': 0.3,
    },
    'ions': {
        'ion_dynamics': 'bfgs',
    }
}

# 擬ポテンシャル (ファイル名を確認してください)
pseudopotentials = {
    'Si': 'Si.pbe-n-kjpaw_psl.1.0.0.UPF',
    'C':  'C.pbe-n-kjpaw_psl.1.0.0.UPF',
    'H':  'H.pbe-kjpaw_psl.1.0.0.UPF' # Hが含まれない場合は無視されます
}
# ========================================

def main():
    # 1. ユニットセルの読み込み
    try:
        # index=-1 で最終構造を取得
        unit_atoms = read(UNIT_CELL_FILE, index=-1)
        print(f"基準セルを読み込みました: {UNIT_CELL_FILE} ({len(unit_atoms)} atoms)")
    except Exception as e:
        print(f"エラー: 基準ファイル {UNIT_CELL_FILE} が読めません。")
        print(e)
        return

    job_list = []

    # 2. 各サイズのファイルをループ
    for f_path in DEFECT_FILES:
        try:
            # 欠陥セルの読み込み（セルサイズを知るため）
            defect_atoms = read(f_path, index=-1)
            
            # スーパーセルの倍率を計算 (Defect Cell / Unit Cell)
            # 対角成分の比率を取って四捨五入 (例: 6.15 / 3.07 = 2.0 -> 2倍)
            M = np.round(defect_atoms.cell.lengths() / unit_atoms.cell.lengths()).astype(int)
            
            # 倍率行列を作成 (対角行列)
            P = [[M[0], 0, 0], [0, M[1], 0], [0, 0, M[2]]]
            
            # 完全結晶スーパーセルを作成
            bulk_supercell = make_supercell(unit_atoms, P)
            
            # ファイル名生成 (例: 8.out -> bulk_8.in)
            base_name = os.path.splitext(os.path.basename(f_path))[0]
            out_name = f"bulk_{base_name}.in"
            
            # 入力ファイル書き出し
            write_espresso_in(
                out_name, 
                bulk_supercell, 
                input_data=input_data, 
                pseudopotentials=pseudopotentials,
                kpts=(1, 1, 1), # ガンマ点
                koffset=(0, 0, 0)
            )
            
            print(f"作成: {out_name} (Size: {M}, Atoms: {len(bulk_supercell)})")
            job_list.append(out_name)

        except Exception as e:
            print(f"スキップ: {f_path} の処理中にエラー ({e})")

    # 3. 実行用シェルスクリプトの作成
    with open("run_all_bulk.sh", "w") as f:
        f.write("#!/bin/zsh\n")
        f.write("#PBS -N BulkAll\n")
        f.write("#PBS -l amdmi02=1:ppn=32\n")
        f.write("#PBS -e stderr_bulk_all.txt\n")
        f.write("#PBS -o stdout_bulk_all.txt\n\n")
        f.write("cd $PBS_O_WORKDIR\n")
        f.write("source /home/common/intel/oneapi/setvars.sh > /dev/null 2>&1\n")
        f.write("QE_BIN=~/q-e/bin/pw.x\n")
        f.write("export OMP_NUM_THREADS=32\n\n")
        
        for job in job_list:
            outfile = job.replace(".in", ".out")
            f.write(f"echo 'Running {job}...'\n")
            f.write(f"mpirun -np 1 $QE_BIN -in {job} > {outfile}\n")
            f.write(f"echo 'Done.'\n\n")

    print("\n完了: 実行スクリプト 'run_all_bulk.sh' を作成しました。")

if __name__ == "__main__":
    main()