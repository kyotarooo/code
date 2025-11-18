#xiでソートし、原子番号をふり、そこで繰り返し処理をかける

import numpy as np # type: ignore
import sys
from tqdm import tqdm#type: ignore


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
    r_range_limited = 10*(1.0e-10)#←←←←←←←適用領域←←←←←←←←←←

    iidx = np.array([0, 0, 0, 1, 1, 1, 2, 2, 2])
    jidx = np.array([0, 1, 2, 0, 1, 2, 0, 1, 2])
    P_redidual = np.array([3451101.8, 2.9636976E-07, 1.5559328E-07, 2.9636976E-07, 3451101.8, 2.0846301E-07, 1.5559328E-07, 2.0846301E-07, 3451101.8])
    P_redidual = P_redidual * 1.0e-25

    atom_distance_list = []
    green_tensor = []
    xi_list = []
    ui_list = []

    with open(f"{sys.argv[1]}/data.inp", "r") as f:#data.inpというファイルをsys.argv[1]で指定されたフォルダから読み込む、python script.py foldernameの場合、sys.argv[1]にfoldernameが入る 
        lines = f.readlines()#ファイルの内容を行ごとに読み込んで、各行を要素とするリストを返す

        g = float((lines[0].split())[0])#わかりやすいように、（）をつけた、ただのインデック操作
        v = float(lines[0].split()[1])
        xc = np.array([float(lines[1].split()[i]) for i in range(3)])#xcは欠陥の座標
        n = int(lines[2].split()[0])
        idx = 3

        for m in range(n):
            xi = np.array([float(lines[idx].split()[i]) for i in range(3)])
            ui = np.array([float(lines[idx].split()[i + 3]) for i in range(3)])
            idx = idx + 1

            rx = xi - xc
            atom_distance_list.append((xi,ui,rx))
           
    counter = 0

    sorted_atom_distance_list = sorted(atom_distance_list, key=lambda x: np.linalg.norm(x[2]))
    for xi,ui,rx in tqdm(sorted_atom_distance_list):

        xi_list.append(xi)
        ui_list.append(ui)
        if np.linalg.norm(rx)<r_range_limited: ##適用領域
            green = np.zeros((3, 3, 3))
            counter += 1
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

    with open(f"/Users/kyou/Library/CloudStorage/Box-Box/lammps_related/disp_of_green.xyz","w") as f1:
        with open(f"{sys.argv[1]}/disp_of_green_sort.out","w") as f2:
            with open(f"{sys.argv[1]}/disp_of_MD_sort.out","w") as f3: 
                with open(f"{sys.argv[1]}/disp_of_green_r_sort.out","w") as f4:
                    f1.write(f"{n}\n")
                # 拡張XYZのヘッダ行（各列の意味を定義）
                # species: 元素記号, pos: 座標, disp: 変位ベクトル（これがOVITOでベクトル表示される）
                    f1.write("Properties=species:S:1:pos:R:3:disp:R:3\n")
                    f2.write(f"{r_range_limited}\n")
                    f3.write(f"{r_range_limited}\n")
                    f4.write(f"{r_range_limited}\n")

                    for atom_id in range(n):
                        X,Y,Z = xi_list[atom_id]
                        if atom_id < counter:
                            ux=uy=uz=0
                            ux_r = uy_r = uz_r = 0
                            for j in range(9):
                                ux -= green_tensor[atom_id][0][iidx[j]][jidx[j]]*x[j]
                                uy -= green_tensor[atom_id][1][iidx[j]][jidx[j]]*x[j]
                                uz -= green_tensor[atom_id][2][iidx[j]][jidx[j]]*x[j]

                                ux_r -= green_tensor[atom_id][0][iidx[j]][jidx[j]]*P_redidual[j]
                                uy_r -= green_tensor[atom_id][1][iidx[j]][jidx[j]]*P_redidual[j]
                                uz_r -= green_tensor[atom_id][2][iidx[j]][jidx[j]]*P_redidual[j]

                        else:
                            ux,uy,uz = ui_list[atom_id]

                        scale_factor=100
                        f1.write(f"Cu {X:.12e} {Y:.12e} {Z:.12e} {ux*scale_factor:.12e} {uy*scale_factor:.12e} {uz*scale_factor:.12e}\n")#f文字列の桁数指
                        scale_factor=1
                        f2.write(f"{X:.12e} {Y:.12e} {Z:.12e} {ux*scale_factor:.12e} {uy*scale_factor:.12e} {uz*scale_factor:.12e}\n")
                        f4.write(f"{X:.12e} {Y:.12e} {Z:.12e} {ux_r*scale_factor:.12e} {uy_r*scale_factor:.12e} {uz_r*scale_factor:.12e}\n")


                    for atom_id in range(n):
                        X,Y,Z = xi_list[atom_id]
                        ux,uy,uz = ui_list[atom_id]
                        scale_factor=1
                        f3.write(f"{X:.12e} {Y:.12e} {Z:.12e} {ux*scale_factor:.12e} {uy*scale_factor:.12e} {uz*scale_factor:.12e}\n")
# ...既存のコード...
    with open(f"{sys.argv[1]}/force_dipole_green.out", "w") as f:
        for i in range(9):
            f.write(f"{iidx[i]} {jidx[i]} {x[i]}\n")

    with open(f"{sys.argv[1]}/force_dipole_redidual.out", "w") as f:
        for i in range(9):
            f.write(f"{iidx[i]} {jidx[i]} {P_redidual[i]}\n")


    with open(f"{sys.argv[1]}/green.out", "w") as f:#sys.argvで指定されたフォルダ内のdata.outファイルに書き込む
        for i in range(counter):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        f.write(f"{i} {green_tensor[i][j][k][l]}")
                    f.write("\n")
                f.write("\n")



if __name__ == "__main__":#『このファイルがコマンドラインからスクリプトとして直接実行された場合のみ以降の処理を実行する
    main()



