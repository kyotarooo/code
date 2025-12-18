import matplotlib
matplotlib.use("Agg") # GUIなしで実行（エラー防止）

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
from pathlib import Path

# ==========================================
# 設定エリア
# ==========================================
FILENAME = "./out/mechanical_behavior.out"  # 読み込むファイル
OUTPUT_VIDEO = "/Users/kyou/Desktop/movie/strain_stress_10s.mp4" # 保存する動画名
DURATION_SEC = 10.0    # 動画の長さ（秒）★ここをParaViewと合わせる！
FPS = 30               # フレームレート（1秒間のコマ数）
# ==========================================

def load_data(path: Path):
    # mechanical_behavior.out の列構成に合わせて調整してください
    # 通常: Col 0=Time, 1=Strain, 2=Plastic, 3=Stress
    try:
        data = np.loadtxt(path, comments="#")
        strain = data[:, 1] # 2列目
        stress = data[:, 3] # 4列目
        return strain, stress
    except Exception as e:
        print(f"Error loading data: {e}")
        return np.array([]), np.array([])

def main():
    path = Path(FILENAME)
    if not path.exists():
        print(f"ファイルが見つかりません: {FILENAME}")
        return

    # データ読み込み
    raw_strain, raw_stress = load_data(path)
    if raw_strain.size == 0:
        return

    stress_mpa = raw_stress * 1e-6 # Pa -> MPa

    # --- フレーム数の計算とデータの間引き（重要） ---
    # 10秒 x 30fps = 300フレーム になるようにデータを間引く・補間する
    total_frames = int(DURATION_SEC * FPS)
    
    # 元データのインデックス（0, 1, 2... N）
    original_indices = np.arange(len(raw_strain))
    # 動画用のインデックス（0, ..., N を 300等分したもの）
    resampled_indices = np.linspace(0, len(raw_strain) - 1, total_frames)

    # 線形補間でデータをサンプリング
    strain_frames = np.interp(resampled_indices, original_indices, raw_strain)
    stress_frames = np.interp(resampled_indices, original_indices, stress_mpa)

    # --- グラフ描画 ---
    fig, ax = plt.subplots(figsize=(6, 4))
    
    # 背景に全データを薄く描く（軌跡として）
    ax.plot(raw_strain * 1e3, stress_mpa, color="lightgray", lw=1.5)
    
    # 動く点（赤）
    red_point, = ax.plot([], [], "o", color="red", markersize=8, zorder=5)
    
    # 現在地までの軌跡（黒）
    current_line, = ax.plot([], [], "-", color="black", lw=2, zorder=4)

    ax.set_xlabel(r"Strain [$10^{-3}$]", fontsize=12)
    ax.set_ylabel("Stress [MPa]", fontsize=12)
    ax.set_title(f"S-S Curve", fontsize=10)
    ax.grid(True, linestyle=":", alpha=0.6)
   

    # 軸の範囲固定
    ax.set_xlim(raw_strain.min()*1e3, raw_strain.max()*1e3 * 1.05)
    ax.set_ylim(stress_mpa.min()*1.5, stress_mpa.max() * 1.1)

    fig.tight_layout()

    # --- アニメーション更新関数 ---
    def update(frame_idx):
        # 現在のデータ
        x = strain_frames[frame_idx] * 1e3
        y = stress_frames[frame_idx]
        
        # 点を移動
        red_point.set_data([x], [y])
        
        # 過去の軌跡を描画（最初から今のフレームまで）
        # ※重くなる場合は current_line の描画を消してもOK
        past_x = strain_frames[:frame_idx+1] * 1e3
        past_y = stress_frames[:frame_idx+1]
        current_line.set_data(past_x, past_y)
        
        return red_point, current_line

    print(f"動画作成中... ({total_frames} frames)")
    
    ani = FuncAnimation(
        fig, update, frames=total_frames, blit=True
    )

    writer = FFMpegWriter(fps=FPS)
    ani.save(OUTPUT_VIDEO, writer=writer, dpi=200)
    print(f"完了: {OUTPUT_VIDEO}")

if __name__ == "__main__":
    main()