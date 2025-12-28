#!/bin/zsh
#PBS -N CutoffTest
#PBS -l nodes=intca01:ppn=36
#PBS -e stderr_cutoff.txt
#PBS -o stdout_cutoff.txt

cd $PBS_O_WORKDIR

# --- 環境設定 ---
source /home/common/intel/oneapi/setvars.sh > /dev/null 2>&1
QE_BIN=~/q-e/bin/pw.x
# スレッド並列数の指定
export OMP_NUM_THREADS=32

# --- 設定 ---
TEMPLATE="base.in"            # 35原子のベースファイル
RESULT_FILE="cutoff_result.dat" # 結果をまとめるファイル

# 保存用ディレクトリの作成（念のため）
mkdir -p in out

# ヘッダー出力
echo "# Ecut(Ry)   Total_Energy(Ry)     Stress_xx(kbar)    Time(sec)" > $RESULT_FILE
echo "---------------------------------------------------------------"
echo "Job started on $(hostname) at $(date)"
echo "Using executable: ${QE_BIN}"

# --- テストループ (10Ry ～ 100Ry) ---
# MEMO: 元のリストは 20 30 ... 100 でした
for ecut in 20 30 40 50 60 70 80 90 100
do
    # PAWポテンシャル用に ecutrho は 8倍
    rho=$((ecut * 8))
    
    # ファイル名
    input_file="in/test_${ecut}.in"
    output_file="out/test_${ecut}.out"

    # sedで置換して新しい入力ファイルを作成
    sed \
  -e "s/calculation.*=.*'relax'/calculation = 'scf'/" \
  -e "s/ecutwfc.*=.*,/ecutwfc = ${ecut}.0,/" \
  -e "s/ecutrho.*=.*,/ecutrho = ${rho}.0,/" \
  -e "/K_POINTS/,/^[[:space:]]*$/c\\
K_POINTS automatic\\
  4 4 2 0 0 0" \
  "$TEMPLATE" > "$input_file"

    # --- 削除しました: $TEMPLATE > $input_file ---

    echo "Running ecut=${ecut} Ry ..."
    start_time=$(date +%s)

    # --- 計算実行 ---
    # pw.x の入力として正しく作成された input_file を指定
    mpirun -n 1 ${QE_BIN} -in ${input_file} > ${output_file}
    
    end_time=$(date +%s)
    elapsed=$(( end_time - start_time ))

    # --- 結果抽出 ---
    # エネルギー
    ene=$(grep "!    total energy" ${output_file} | tail -1 | awk '{print $5}')
    strs=$(grep -A 1 "total   stress" ${output_file} | tail -1 | awk '{print $4}')

    # エラー回避 (空ならNaN)
    if [ -z "$ene" ]; then ene="NaN"; fi
    if [ -z "$strs" ]; then strs="NaN"; fi

    # 画面とファイルに出力
    printf "Ecut: %-4s Energy: %-15s Stress: %-10s Time: %ss\n" "$ecut" "$ene" "$strs" "$elapsed"
    echo "$ecut  $ene  $strs  $elapsed" >> $RESULT_FILE
done

echo "---------------------------------------------------------------"
echo "All calculations finished. Check $RESULT_FILE"