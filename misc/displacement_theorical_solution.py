import numpy as np # type: ignore
import sys


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
    return green


def main():
    a = np.zeros((9, 9))
    b = np.zeros(9)

    iidx = np.array([0, 0, 0, 1, 1, 1, 2, 2, 2])
    jidx = np.array([0, 1, 2, 0, 1, 2, 0, 1, 2])

    green_tensor = []
    xj = []

    with open(f"{sys.argv[1]}/data.inp", "r") as f:#data.inpというファイルをsys.argv[1]で指定されたフォルダから読み込む、python script.py foldernameの場合、sys.argv[1]にfoldernameが入る 
        lines = f.readlines()#ファイルの内容を行ごとに読み込んで、各行を要素とするリストを返す

        g = float((lines[0].split())[0])#わかりやすいように、（）をつけた、ただのインデック操作
        v = float(lines[0].split()[1])
        xc = np.array([float(lines[1].split()[i]) for i in range(3)])#xcは欠陥の座標
        n = int(lines[2].split()[0])
        idx = 3


        for m in range(n):
            xi = np.array([float(lines[idx].split()[i]) for i in range(3)])
            xj.append(xi)
            ui = np.array([float(lines[idx].split()[i + 3]) for i in range(3)])
            idx = idx + 1

            rx = xi - xc
            green = np.zeros((3, 3, 3))

            green = kelvin_solution(g, v, rx, green)#pythonでは参照渡し
            green_tensor.append(green)

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

    with open(f"/Users/kyou/Library/CloudStorage/Box-Box/lammps_related/displacement_of_theoricalsolution.xyz","w") as f1:
        with open(f"{sys.argv[1]}/displacement_of_theoricalsolution.out","w") as f2:
            f1.write(f"{n}\n")
        # 拡張XYZのヘッダ行（各列の意味を定義）
        # species: 元素記号, pos: 座標, disp: 変位ベクトル（これがOVITOでベクトル表示される）
            f1.write("Properties=species:S:1:pos:R:3:disp:R:3\n")

            for atom_id in range(n):
                ux=uy=uz=0
                for j in range(9):
                    ux -= green_tensor[atom_id][0][iidx[j]][jidx[j]]*x[j]
                    uy -= green_tensor[atom_id][1][iidx[j]][jidx[j]]*x[j]
                    uz -= green_tensor[atom_id][2][iidx[j]][jidx[j]]*x[j]

                X,Y,Z = xj[atom_id]
                scale_factor=100
                f1.write(f"Cu {X:.5f} {Y:.5f} {Z:.5f} {ux*scale_factor:.5f} {uy*scale_factor:.5f} {uz*scale_factor:.5f}\n")#f文字列の桁数指
                scale_factor=1
                f2.write(f"{X:.5f} {Y:.5f} {Z:.5f} {ux*scale_factor:.5f} {uy*scale_factor:.5f} {uz*scale_factor:.5f}\n")

        
    


    with open(f"{sys.argv[1]}/data_p.out", "w") as f:#sys.argvで指定されたフォルダ内のdata.outファイルに書き込む
        for i in range(9):
            f.write(f"{iidx[i]} {jidx[i]} {x[i]}\n")

    with open(f"{sys.argv[1]}/green.out", "w") as f:#sys.argvで指定されたフォルダ内のdata.outファイルに書き込む
        for i in range(n):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        f.write(f"{i} {green_tensor[i][j][k][l]}")
                    f.write("\n")
                f.write("\n")



if __name__ == "__main__":#『このファイルがコマンドラインからスクリプトとして直接実行された場合のみ以降の処理を実行する
    main()



    