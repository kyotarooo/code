import matplotlib
# GUI環境に合わせて変更してください (Macなら 'MacOSX', Linuxなら 'TkAgg' など)
matplotlib.use("TkAgg") 

import matplotlib.pyplot as plt
import pandas as pd
import sys
import os

# ================= 設定 =================
DATA_FILE = "cutoff_result.dat"

# ================= データの読み込み =================
if not os.path.exists(DATA_FILE):
    print(f"エラー: {DATA_FILE} が見つかりません。")
    sys.exit(1)

try:
    col_names = ["Ecut (Ry)", "Total Energy (Ry)", "Stress (kbar)", "Time (s)"]
    
    # 【修正1】 sep=r'\s+' とすることで SyntaxWarning を回避
    df = pd.read_csv(DATA_FILE, sep=r'\s+', comment='#', header=None, names=col_names)
    
    # 数値変換とクリーニング
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    df = df.dropna().sort_values(by="Ecut (Ry)")

    # 【修正2】 kbar -> GPa への単位変換 (1 kbar = 0.1 GPa)
    df["Stress (GPa)"] = df["Stress (kbar)"] * 0.1 

    print(df[["Ecut (Ry)", "Total Energy (Ry)", "Stress (GPa)"]])

except Exception as e:
    print(f"エラー: {e}")
    sys.exit(1)


# ================= グラフ化 =================
# エネルギー差分の計算
e_final = df["Total Energy (Ry)"].iloc[-1]
df["Energy Diff (eV)"] = (df["Total Energy (Ry)"] - e_final) * 13.6057

fig, axes = plt.subplots(1, 3, figsize=(24, 6), constrained_layout=True)

# 1. Total Energy
axes[0].plot(df["Ecut (Ry)"], df["Total Energy (Ry)"], 'o-', color='tab:red', markersize=8)
axes[0].set_ylabel('Total Energy (Ry)', fontsize=16)
axes[0].set_xlabel('Cutoff (Ry)', fontsize=16)
axes[0].set_title('Energy Convergence', fontsize=18)
axes[0].grid(True)

# エネルギー差分（右軸）
ax0_twin = axes[0].twinx()
ax0_twin.plot(df["Ecut (Ry)"], df["Energy Diff (eV)"], 'o--', color='tab:red', alpha=0.5)
ax0_twin.set_ylabel('Diff from Final (eV)', color='tab:red', fontsize=14)
ax0_twin.tick_params(axis='y', labelcolor='tab:red')

#収束値を表示
axes[0].axhline(y = e_final, color = "red", linestyle = "--", alpha = 0.5)

# 2. Stress (GPaでプロット)
axes[1].plot(df["Ecut (Ry)"], df["Stress (GPa)"], 's-', color='tab:green', markersize=8, linewidth=2)

# ラベルも変更
axes[1].set_ylabel('Total Stress (xx component) [GPa]', fontsize=16)
axes[1].set_xlabel('Cutoff (Ry)', fontsize=16)
axes[1].set_title('Stress Convergence (GPa)', fontsize=18)
axes[1].grid(True)

# 最後の値を表示 (GPa)
last_stress_gpa = df["Stress (GPa)"].iloc[-1]
axes[1].axhline(y=last_stress_gpa, color='red', linestyle='--', alpha=0.5)

# 3. Time
axes[2].plot(df["Ecut (Ry)"], df["Time (s)"], '^-', color='tab:purple', markersize=8)
axes[2].set_ylabel('Time (s)', fontsize=16)
axes[2].set_xlabel('Cutoff (Ry)', fontsize=16)
axes[2].set_title('Computational Cost', fontsize=18)
axes[2].grid(True)

# フォントサイズ調整
for ax in axes:
    ax.tick_params(axis='both', labelsize=14)

output_img = 'cutoff_stress_convergence_gpa.png'
plt.savefig(output_img)
print(f"\nグラフを '{output_img}' に保存しました。")
plt.show()