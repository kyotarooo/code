from ovito.io import import_file, export_file
from ovito.data import *
import numpy as np

# ================= 設定エリア =================
# 読み込むファイル名 (※必ずVoronoi体積が含まれているファイル指定してください)
input_filename = "/Users/kyou/output/4HSiC/vacancy/vashishta-k/residual_stress/dump/defect_lattice/defect_lattice_0/stress_viz_8.dump" 

# 出力するファイル名
output_filename = "stress_converted_GPa.dump"

# 単位変換係数 (bar -> GPa)
conversion_factor = 0.0001
# ==============================================

def modify_pipeline(frame, data):
    print(f"--- Frame {frame} 解析開始 ---")
    
    # -------------------------------------------------------
    # 1. 応力データの取得 (バラバラの成分を統合する処理)
    # -------------------------------------------------------
    stress_raw = None
    
    # パターンA: まとまった 'c_st' がある場合
    if 'c_st' in data.particles:
        print("  -> 'c_st' 配列を発見しました。")
        stress_raw = data.particles['c_st']
        
    # パターンB: バラバラの 'c_st[1]'...'c_st[6]' がある場合
    elif 'c_st[1]' in data.particles:
        print("  -> 個別の 'c_st[1]...[6]' を発見。統合します。")
        try:
            s1 = data.particles['c_st[1]'] # xx
            s2 = data.particles['c_st[2]'] # yy
            s3 = data.particles['c_st[3]'] # zz
            s4 = data.particles['c_st[4]'] # xy
            s5 = data.particles['c_st[5]'] # xz
            s6 = data.particles['c_st[6]'] # yz
            # 6成分を (N, 6) の行列にスタックする
            stress_raw = np.column_stack((s1, s2, s3, s4, s5, s6))
        except KeyError as e:
            print(f"❌ エラー: 応力成分の一部が欠けています: {e}")
            return
    else:
        # デバッグ用: 何のデータがあるか表示
        print("❌ エラー: 応力データが見つかりません。")
        print("   現在利用可能なプロパティ一覧:", data.particles.keys())
        return

    # -------------------------------------------------------
    # 2. 体積データの取得
    # -------------------------------------------------------
    vol = None
    if 'c_voro[1]' in data.particles:
         vol = data.particles['c_voro[1]']
    elif 'c_voro' in data.particles:
         c_voro_prop = data.particles['c_voro']
         if c_voro_prop.ndim > 1:
             vol = c_voro_prop[:, 0]
         else:
             vol = c_voro_prop[...]
    
    if vol is None:
        print("❌ エラー: 原子体積 (c_voro または c_voro[1]) が見つかりません。")
        print("   ※LAMMPSで 'compute voronoi/atom' を行い、dumpに含めましたか？")
        print("   現在利用可能なプロパティ一覧:", data.particles.keys())
        return

    # -------------------------------------------------------
    # 3. 計算実行 (bar*A^3 -> GPa)
    # -------------------------------------------------------
    print("  -> GPa単位への変換計算中...")
    with np.errstate(divide='ignore', invalid='ignore'):
        # (N, 6) / (N, 1)
        stress_gpa = (stress_raw / vol[:, np.newaxis]) * conversion_factor
        stress_gpa[~np.isfinite(stress_gpa)] = 0.0

    # フォン・ミーゼス応力
    s_xx = stress_gpa[:, 0]
    s_yy = stress_gpa[:, 1]
    s_zz = stress_gpa[:, 2]
    s_xy = stress_gpa[:, 3]
    s_xz = stress_gpa[:, 4]
    s_yz = stress_gpa[:, 5]

    von_mises = np.sqrt(0.5 * (
        (s_xx - s_yy)**2 + (s_yy - s_zz)**2 + (s_zz - s_xx)**2 + 
        6 * (s_xy**2 + s_yz**2 + s_xz**2)
    ))

    # -------------------------------------------------------
    # 4. 書き込み (Mutable Access)
    # -------------------------------------------------------
    data.particles_.create_property('Stress Tensor (GPa)', data=stress_gpa)
    data.particles_.create_property('Von Mises Stress (GPa)', data=von_mises)
    print("  -> 計算完了。")


# ================= メイン実行部 =================
if __name__ == "__main__":
    print("処理を開始します...")

    try:
        pipeline = import_file(input_filename)
        pipeline.modifiers.append(modify_pipeline)

        print(f"保存中: {output_filename}")
        
        # 出力項目の指定
        # ※もし 'c_voro[1]' が見つからないエラーが出る場合は、下の "c_voro[1]" を削除してください
        export_file(
            pipeline, 
            output_filename, 
            "lammps/dump",
            columns=[
                "Particle Identifier", 
                "Particle Type",       
                "Position.X", "Position.Y", "Position.Z",
                "Stress Tensor (GPa)",      
                "Von Mises Stress (GPa)",
                "c_voro[1]" 
            ]
        )
        print("--------------------------------------------------")
        print("✅ 全処理が正常に完了しました。")
        print(f"出力: {output_filename}")
        print("OVITOで開き、'Von Mises Stress (GPa)' を可視化してください。")
        print("--------------------------------------------------")

    except Exception as e:
        print("\n❌ 予期せぬエラーが発生しました:")
        print(e)