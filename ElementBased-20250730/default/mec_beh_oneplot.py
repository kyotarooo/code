import os
import pandas as pd
import matplotlib.pyplot as plt

def plot_stress_strain_in_current_dir(base_dir):
    """
    base_dir 内の 'mechanicalbehavior' ファイルを探して
    応力-ひずみ曲線をグラフ化し、PNGとして保存する
    """
    target_filename = "mechanical_behavior.out"
    file_path = os.path.join(base_dir, target_filename)

    if not os.path.isfile(file_path):
        print(f"エラー: {file_path} が見つかりません。")
        return

    try:
        # データ読み込み（空白区切り、ヘッダーなし）
        df = pd.read_csv(file_path, sep=r'\s+', header=None)

        if df.shape[1] < 4:
            print(f"警告: {file_path} は必要な4列を持っていません。スキップします。")
            return

        df.columns = [f'Column_{i+1}' for i in range(df.shape[1])]

        # プロット作成
        plt.figure(figsize=(10, 6))
        plt.plot(df['Column_2'], df['Column_4'])
        plt.xlabel('Strain (Column 2)', fontsize=14)
        plt.ylabel('Stress (Column 4) [Pa]', fontsize=14)
        plt.title(f'Stress-Strain Curve: {target_filename}', fontsize=16)
        plt.grid(True)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)

        # 出力ファイル名
        output_filename = f"{target_filename}.stress-strain.png"
        output_path = os.path.join(base_dir, output_filename)

        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()

        print(f'✓ グラフを保存しました: {output_path}')

    except Exception as e:
        print(f'✗ エラー: {file_path} の処理中にエラーが発生しました → {e}')

if __name__ == "__main__":
    current_dir = os.path.dirname(os.path.abspath(__file__))
    plot_stress_strain_in_current_dir(current_dir)
