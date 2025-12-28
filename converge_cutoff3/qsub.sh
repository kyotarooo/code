#!/bin/zsh
#PBS -N CutoffKtest
#PBS -l nodes=amdro01:ppn=64
#PBS -e stderr_conv.txt
#PBS -o stdout_conv.txt

cd $PBS_O_WORKDIR

# --- 環境設定 ---
source /home/common/intel/oneapi/setvars.sh > /dev/null 2>&1
QE_BIN=~/q-e/bin/pw.x
export OMP_NUM_THREADS=2

# --- 設定 ---
TEMPLATE="base.in"
RESULT_FILE="convergence_result.dat"

# 保存用ディレクトリ
mkdir -p in out

# --- パラメータ定義 ---
# カットオフエネルギーのリスト (Ry)
ECUT_LIST=(30 40 50 60 70 80 90 100)

# K点のリスト ("kx ky kz")
# SiCのセル形状(a=6.15, c=10.05)から、c軸はa軸より長いため、
# 逆格子空間では kz は kx,ky より小さくて良い (例: 4 4 2)
K_LIST=(
  "2 2 1"
  "3 3 2"
  "4 4 2"
  "5 5 3"
)

# ヘッダー出力
echo "# Ecut(Ry)  K-points    Total_Energy(Ry)     Stress_xx(kbar)    Time(sec)" > $RESULT_FILE
echo "-------------------------------------------------------------------------------"
echo "Job started on $(hostname) at $(date)"

# --- 2重ループ開始 ---
for k_mesh in "${K_LIST[@]}"; do
    # K点の文字列から空白を除去してファイル名用にする (例: "4 4 2" -> "442")
    k_name=$(echo $k_mesh | tr -d ' ')

    for ecut in "${ECUT_LIST[@]}"; do
        
        # PAWポテンシャル用に ecutrho は ecutwfc の 8倍-10倍程度 (ここでは8倍)
        rho=$((ecut * 8))

        # ファイル名設定 (EとKの両方を含める)
        label="E${ecut}_K${k_name}"
        input_file="in/test_${label}.in"
        output_file="out/test_${label}.out"

        # sedで入力ファイルを作成
        # calculationを 'scf' に変更し、ecut と K_POINTS を書き換える
        sed \
            -e "s/calculation.*=.*'relax'/calculation = 'relax'/" \
            -e "s/ecutwfc.*=.*,/ecutwfc = ${ecut}.0,/" \
            -e "s/ecutrho.*=.*,/ecutrho = ${rho}.0,/" \
            -e "/K_POINTS/,/^[[:space:]]*$/c\\
K_POINTS automatic\\
  ${k_mesh} 0 0 0" \
            "$TEMPLATE" > "$input_file"

        echo "Running ${label} (Ecut=${ecut}, K=${k_mesh}) ..."
        start_time=$(date +%s)

        # --- 計算実行 ---
        mpirun -n 32 ${QE_BIN} -in ${input_file} > ${output_file}

        end_time=$(date +%s)
        elapsed=$(( end_time - start_time ))

        # --- 結果抽出 ---
        # エネルギー
        ene=$(grep "!    total energy" ${output_file} | tail -1 | awk '{print $5}')
        # ストレス (xx成分)
        strs=$(grep -A 1 "total   stress" ${output_file} | tail -1 | awk '{print $4}')

        # エラー回避
        if [ -z "$ene" ]; then ene="NaN"; fi
        if [ -z "$strs" ]; then strs="NaN"; fi

        # 画面とファイルに出力
        # 見やすくフォーマット (K点は引用符なしで表示)
        printf "Ecut: %-4s K: %-8s E: %-14s S: %-9s T: %ss\n" "$ecut" "$k_name" "$ene" "$strs" "$elapsed"
        
        # データファイルには解析しやすい形式で書き込み
        echo "$ecut  $k_mesh  $ene  $strs  $elapsed" >> $RESULT_FILE
    done
    # K点が変わるごとに区切り線を入れる（ファイル側）
    echo "" >> $RESULT_FILE
done

echo "-------------------------------------------------------------------------------"
echo "All calculations finished. Check $RESULT_FILE"