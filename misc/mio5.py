import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# ======== スタイル設定 ========
sns.set_theme(style="ticks", context="talk", font="Meiryo")
plt.rcParams["axes.linewidth"] = 1.2  # 軸の太さを少し太く
plt.rcParams["axes.labelweight"] = "bold"  # ラベルを太字に
plt.rcParams["font.size"] = 14  # 全体フォントサイズ

# ======== データ ========
data = {
    "数学": [100, 85, 90, 95, 80, 80, 75, 65, 65, 60, 55, 45, 45],
    "理科": [94, 90, 95, 90, 85, 80, 75, 70, 60, 60, 50, 50, 48],
    "社会": [80, 88, 70, 62, 86, 70, 79, 65, 75, 67, 75, 68, 60],
}
df = pd.DataFrame(data)

# ======== 図設定 ========
fig, ax = plt.subplots(figsize=(6.5, 5.5))

# 散布図＋回帰線
sns.regplot(
    data=df, 
    x="数学", 
    y="理科", 
    scatter_kws={"s": 70, "alpha": 0.8, "color": "#007acc", "edgecolor": "black"},
    line_kws={"color": "crimson", "linewidth": 2.2},
    ax=ax
)

# ======== 軸・タイトル・注記 ========
ax.set_title("数学と理科の得点の相関", fontsize=18, fontweight="bold", pad=15)
ax.set_xlabel("数学の得点", fontsize=16)
ax.set_ylabel("理科の得点", fontsize=16)

# 相関係数を図中に表示
corr = df["数学"].corr(df["理科"])
ax.text(0.05, 0.9, f"相関係数 r = {corr:.2f}", transform=ax.transAxes,
        fontsize=14, bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="gray"))

# 枠・余白などの整理
sns.despine(trim=True)
plt.tight_layout()

# ======== 出力 ========
plt.show()
