from lammps import lammps
import os

# 実行ファイルの絶対パス
PINN_LAMMPS_PATH = '/Users/kyou/LAMMPS-USER-PINN/src/lmp_mpi' 

# 'name' キーワード引数で実行ファイルを指定
lmp = lammps(name=PINN_LAMMPS_PATH, cmdargs=['-log', 'none']) 
                                       
# 利用可能なスタイルをチェック
available_styles = lmp.available_styles('pair_style') 

print("--- 利用可能なポテンシャルスタイル ---")

if 'pinn' in available_styles:
    print("✅ 成功: 'pair_style pinn' が見つかりました。PINNは利用可能です！")
else:
    print("❌ 失敗: 'pair_style pinn' が見つかりません。")
    
lmp.close()