import numpy as np
import sys
import subprocess
import os

class Movie:
    def __init__(self):
        self.n = 1
        self.start_coord = np.zeros(3)
        self.end_coord = np.zeros(3)

class Section:
    def __init__(self):
        self.x = np.zeros(3)
        self.n = np.zeros(3)
        self.s = np.zeros(3)
        self.m = np.zeros(3)

class Volume:
    def __init__(self):
        self.section = Section()
        self.size = np.zeros(3)
        self.nx = np.zeros(3, dtype=int)

def read_movie(movie):
    path = os.path.join(sys.argv[1], "movie.inp")
    with open(path, "r") as f:
        line = f.readline()
        movie.n = int(line.split()[0])
        line = f.readline()
        # 浮動小数点として読み込む
        movie.start_coord = np.array([float(val) for val in line.split()])
        line = f.readline()
        movie.end_coord = np.array([float(val) for val in line.split()])

def read_section(volume):
    # 初期設定を読み込む（テンプレートとしてメモリに保持）
    path = os.path.join(sys.argv[1], "section.inp")
    with open(path, "r") as f:
        lines = f.readlines()
        volume.size = np.array([float(val) for val in lines[0].split()])
        volume.section.x = np.array([float(val) for val in lines[1].split()])
        volume.section.n = np.array([float(val) for val in lines[2].split()])
        volume.section.s = np.array([float(val) for val in lines[3].split()])
        volume.section.m = np.array([float(val) for val in lines[4].split()])
        volume.nx = np.array([int(val) for val in lines[5].split()])

def section_position(i, movie, volume):
    if movie.n > 1:
        dcoord = (movie.end_coord - movie.start_coord) / (movie.n - 1.0)
        volume.section.x = movie.start_coord + i * dcoord
    else:
        volume.section.x = movie.start_coord

def write_section(volume):
    # 【修正箇所】
    # C++プログラムが読み込む "本元の" section.inp を上書きする
    path = os.path.join(sys.argv[1], "section.inp")
    
    size = volume.size
    x = volume.section.x
    n = volume.section.n
    s = volume.section.s
    m = volume.section.m
    nx = volume.nx
    
    with open(path, "w") as f:
        f.write(f"{size[0]} {size[1]} {size[2]}\n")
        f.write(f"{x[0]} {x[1]} {x[2]}\n")
        f.write(f"{n[0]} {n[1]} {n[2]}\n")
        f.write(f"{s[0]} {s[1]} {s[2]}\n")
        f.write(f"{m[0]} {m[1]} {m[2]}\n")
        f.write(f"{nx[0]} {nx[1]} {nx[2]}\n")

def take_sceen(i):
    base_dir = sys.argv[1]
    
    # C++を実行（書き換わった section.inp を読み込んで計算する）
    command = f"./section_vtk {base_dir}"
    subprocess.run(command, shell=True)

    # 生成されたファイルを movie フォルダへ移動
    src_file = os.path.join(base_dir, "section.vtk")
    dst_file = os.path.join(base_dir, "movie", f"section_{i + 1:04d}.vtk")
    
    if os.path.exists(src_file):
        os.rename(src_file, dst_file)
    else:
        print(f"Error: Frame {i+1} failed. {src_file} not found.")

def main():
    if len(sys.argv) < 2:
        print("Usage: python make_movie.py <target_directory>")
        sys.exit(1)

    target_dir = sys.argv[1]
    
    # movieフォルダ作成
    os.makedirs(os.path.join(target_dir, "movie"), exist_ok=True)

    movie = Movie()
    volume = Volume()

    read_movie(movie)
    read_section(volume) # ここで元の設定をメモリに保存

    print(f"Generating {movie.n} frames...")
    for i in range(movie.n):
        section_position(i, movie, volume)
        write_section(volume) # ここでファイルを更新
        take_sceen(i)
        print(f"\rProcessed {i+1}/{movie.n}", end="")
    
    print("\nDone.")

if __name__ == "__main__":
    main()