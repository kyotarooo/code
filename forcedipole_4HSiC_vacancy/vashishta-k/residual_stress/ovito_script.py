from ovito.io import import_file, export_file
from ovito.data import *
import numpy as np

# ================= 設定エリア =================
# 読み込むファイル名（ご自身の環境に合わせて書き換えてください）
input_filename = "/Users/kyou/output/4HSiC/vacancy/vashishta-k/residual_stress/dump/defect_lattice/defect_lattice_0/stress_viz_8.dump" 

# 出力するファイル名
output_filename = "stress_converted_GPa.dump"

# LAMMPSの単位系変換係数 (metal units)
# input: bar * A^3
# output: GPa
# 1 bar = 0.0001 GPa
conversion_factor = 0.0001
# ==============================================

def modify_pipeline(frame, data):
    """
    OVITOのパイプライン処理関数
    ここで単位変換と原子タイプの設定を行います
    """
    
    # -------------------------------------------------------
    # 1. 4H-SiCとしての原子タイプ設定 (1=Si, 2=C)
    # -------------------------------------------------------
    # Particle Typeプロパティを取得
    ptypes = data.particles.particle_types
    
    # タイプ1があれば Si に設定 (半径と色を自動設定)
    type1 = ptypes.type_by_id(1)
    if type1:
        type1.name = "Si"
        type1.radius = 1.11  # 共有結合半径(目安)
        type1.color = (0.8, 0.8, 0.8) # 薄いグレー

    # タイプ2があれば C に設定
    type2 = ptypes.type_by_id(2)
    if type2:
        type2.name = "C"
        type2.radius = 0.77
        type2.color = (0.2, 0.2, 0.2) # 濃いグレー

    # -------------------------------------------------------
    # 2. 応力の単位変換 (bar*A^3 -> GPa)
    # -------------------------------------------------------
    
    # LAMMPSから読み込んだデータの取得
    # c_st (6成分), c_voro (体積) があることを前提とします
    try:
        # 配列として取得 (N原子 x 6成分)
        stress_raw = data.particles['c_st'] 
        # 配列として取得 (N原子) 
        # c_voroが配列の場合とタプルの場合があるので安全に取得
        if 'c_voro' in data.particles:
             # c_voro[1]相当を取得 (通常LAMMPS出力は配列の1列目に入ることが多い)
             vol = data.particles['c_voro'][:, 0] if data.particles['c_voro'].ndim > 1 else data.particles['c_voro']
        elif 'c_voro[1]' in data.particles:
             vol = data.particles['c_voro[1]']
        else:
            print("警告: 原子体積のデータ(c_voro)が見つかりません。")
            return

    except KeyError as e:
        print(f"エラー: データ列が見つかりません: {e}")
        return

    # --- 計算処理 (NumPyで高速計算) ---
    
    # 体積が0の原子によるゼロ除算を防ぐ
    with np.errstate(divide='ignore', invalid='ignore'):
        # 応力テンソル (GPa) = (生データ / 体積) * 0.0001
        # vol[:, np.newaxis] は次元を合わせて割り算するため
        stress_gpa = (stress_raw / vol[:, np.newaxis]) * conversion_factor
        
        # 無限大やNaNを0に置換
        stress_gpa[~np.isfinite(stress_gpa)] = 0.0

    # -------------------------------------------------------
    # 3. フォン・ミーゼス応力の計算 (Von Mises Stress)
    # -------------------------------------------------------
    # tensor components: 0:xx, 1:yy, 2:zz, 3:xy, 4:xz, 5:yz
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
    # 4. 新しいプロパティとしてデータに追加
    # -------------------------------------------------------
    # GPa単位のテンソルを書き込み
    data.particles_.create_property('Stress Tensor (GPa)', data=stress_gpa)
    
    # Von Mises応力を書き込み (これをColor Codingに使うと便利)
    data.particles_.create_property('Von Mises Stress (GPa)', data=von_mises)

    print(f"Frame {frame}: 単位変換完了 (Bar*A^3 -> GPa)")


# ================= メイン実行部 =================
print("処理を開始します...")

# 1. ファイル読み込み
pipeline = import_file(input_filename)

# 2. 上記の変換関数をパイプラインに適用
pipeline.modifiers.append(modify_pipeline)

# 3. 計算を実行して保存
# columns引数で出力したい項目を指定できます
print(f"保存中: {output_filename}")
export_file(
    pipeline, 
    output_filename, 
    "lammps/dump",
    columns=[
        "Particle Identifier", # id
        "Particle Type",       # type (Si, Cの名前は数字に戻りますが属性は保持されます)
        "Position.X", "Position.Y", "Position.Z",
        "Stress Tensor (GPa)",      # 計算したGPaテンソル
        "Von Mises Stress (GPa)",   # 計算したMises応力
        "Atomic Volume"             # 元の体積(c_voro)
    ]
)

print("完了しました。")
print("OVITOで出力ファイルを開き、'Von Mises Stress (GPa)' でColor Codingしてください。")