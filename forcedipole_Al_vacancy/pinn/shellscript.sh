#!/bin/zsh

# 実行ファイルのパス
LAMMPS=~/LAMMPS-USER-PINN/src/lmp_mpi

# 入力ファイルリスト
LIST=/Users/kyou/Library/CloudStorage/Box-Box/lammps/pointdefects_Al/displacement/infilename.txt

# 入力ファイルがあるディレクトリ
INPUT_DIR=/Users/kyou/Library/CloudStorage/Box-Box/lammps/pointdefects_Al/displacement/in_file

# リストを1行ずつ読み取って実行
for infile in $(cat "$LIST"); do
    echo "=== Running $infile ==="
    mpirun -np 1 "$LAMMPS" -in "$INPUT_DIR/$infile"
done

## ＜：　入力リダイレクト　　LISTファイルを、infileに読み込んで処理する

