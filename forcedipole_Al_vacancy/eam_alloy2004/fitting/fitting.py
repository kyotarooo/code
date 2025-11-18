import numpy as np # type: ignore
import sys

conv_eV_to_J = 1.602176634E-19
conv_J_to_eV = 1/conv_eV_to_J

def kelvin_solution(g, v, rx, green):
    r = np.linalg.norm(rx)#ベクトルのノルムを計算、rx=(x1,x2,x3)
    rx /= r#rxについて更新、微分した形になっている。
    c = 1.0 / (16.0 * np.pi * g * (1.0 - v) * r * r)#gとは？\nuと何が違う？
    d = np.zeros((3, 3))#3×3の行列を作成、要素は全て０
    for i in range(3):
        d[i][i] = 1.0#対角成分のみ１にする

    for i in range(3):
        for j in range(3):
            for k in range(3):
                green[i][j][k] = c * (
                    -1.0 * rx[i] * d[j][k]
                    - rx[j] * d[i][k]
                    + (3.0 - 4.0 * v) * rx[k] * d[i][j]
                    + 3.0 * rx[i] * rx[j] * rx[k]
                )

def main():
    a = np.zeros((9, 9))
    b = np.zeros(9)

    iidx = np.array([0, 0, 0, 1, 1, 1, 2, 2, 2])
    jidx = np.array([0, 1, 2, 0, 1, 2, 0, 1, 2])

    with open(f"{sys.argv[1]}/displacement9.inp", "r") as f:#data.inpというファイルをsys.argv[1]で指定されたフォルダから読み込む、python script.py foldernameの場合、sys.argv[1]にfoldernameが入る 
        lines = f.readlines()#ファイルの内容を行ごとに読み込んで、各行を要素とするリストを返す

        g = float((lines[0].split())[0])#わかりやすいように、（）をつけた、ただのインデック操作
        v = float(lines[0].split()[1])
        xc = np.array([float(lines[1].split()[i]) for i in range(3)])#xcは欠陥の座標
        n = int(lines[2].split()[0])
        idx = 3
        for _ in range(n):#一時的にアンダーバーに代入されるが、その数値は重要ではない、原子の数だけ繰り返す
            xi = np.array([float(lines[idx].split()[i]) for i in range(3)])
            ui = np.array([float(lines[idx].split()[i + 3]) for i in range(3)])
            idx = idx + 1

            rx = xi - xc
            r = np.linalg.norm(rx)

            green = np.zeros((3, 3, 3))
            kelvin_solution(g, v, rx, green)

            for i in range(9):
                for j in range(9):
                    for k in range(3):
                        a[i][j] += (
                            green[k][iidx[i]][jidx[i]] * green[k][iidx[j]][jidx[j]]#*:ただ単に要素ごとの積、a[i][j]でうまくテンソルの形に直している
                        )

            for i in range(9):
                for j in range(3):
                    b[i] -= ui[j] * green[j][iidx[i]][jidx[i]]

    a_inv = np.linalg.inv(a)#正方行列aの逆行列inverse
    x = np.dot(a_inv, b)#行列積の計算

       


    with open(f"{sys.argv[1]}/data_eV_9.out", "w") as f:#sys.argvで指定されたフォルダ内のdata.outファイルに書き込む
        for i in range(9):
            f.write(f"{iidx[i]} {jidx[i]} {x[i]*conv_J_to_eV}\n")

       


if __name__ == "__main__":#『このファイルがコマンドラインからスクリプトとして直接実行された場合のみ以降の処理を実行する
    main()

    