import numpy as np # type: ignore
import sys
from tqdm import tqdm#type: ignore
import os

######## 環境変数の取得 ########
output_dir = os.environ.get("OUTPUT_PATH")
if output_dir is None:
    raise ValueError("環境変数 OUTPUT_PATH error")

######## 材料定数の取得 ############
lattice_const = os.environ.get("LATTICE_CONST") # ang
lattice_const = float(lattice_const)

######### 単位変換係数 ########
conv_eV_to_J = 1.602176634E-19
conv_bars_to_Pa = 1.0E+5
conv_Vang_to_Vm = 1.0E-30
conv_barVang_to_NVm = conv_bars_to_Pa * conv_Vang_to_Vm
conv_J_to_eV = 1/conv_eV_to_J
conv_ang_to_m = 1.0E-10

######## ファイル数の読み込み ########
with open(f"{output_dir}/infilename.txt", 'r') as f_num:
    lines = f_num.readlines()
    file = lines[0].split()
    file_num = len(file)
    print(f"{file_num}")

######## kelvin解 ########
def kelvin_solution(g, v, rx, green):
    r = np.linalg.norm(rx)
    rx /= r #←←微分
    c = 1.0 / (16.0 * np.pi * g * (1.0 - v) * r * r)
    d = np.zeros((3, 3))
    for i in range(3):
        d[i][i] = 1.0

    for i in range(3):
        for j in range(3):
            for k in range(3):
                green[i][j][k] = (-1) * c * (
                    -1.0 * rx[i] * d[j][k]
                    - rx[j] * d[k][i]
                    + (3.0 - 4.0 * v) * rx[k] * d[i][j]
                    + 3.0 * rx[i] * rx[j] * rx[k]
                )
    return green,r

def main():
    ######### cut offの定義 ########
    x_list = []
    limit_range = np.array([0, 0.5, 1, 1.5, 2]) * conv_ang_to_m

    for l in range(file_num):
        with open(f"{output_dir}/disp_data/MD_disp_data{l}.txt", "w") as f_MD_disp_data:
            for m in limit_range:
                r_range_limit = lattice_const * m
                a = np.zeros((9, 9))
                b = np.zeros(9)
                iidx = np.array([0, 0, 0, 1, 1, 1, 2, 2, 2])
                jidx = np.array([0, 1, 2, 0, 1, 2, 0, 1, 2])
                atom_distance_list = []

                ######### MDのデータ読み込み ########
                with open(f"{output_dir}/dump/displacement/displacement{l}.inp", "r") as f_disp: 
                    lines = f_disp.readlines()
                    g = float((lines[0].split())[0])
                    v = float(lines[0].split()[1])
                    xc = np.array([float(lines[1].split()[i]) for i in range(3)])#xcは欠陥の座標
                    n = int(lines[2].split()[0])
                    idx = 3

                    for _ in range(n):
                        xi = np.array([float(lines[idx].split()[i]) for i in range(3)])
                        ui = np.array([float(lines[idx].split()[i + 3]) for i in range(3)])
                        idx = idx + 1
                        rx = xi - xc
                        atom_distance_list.append((xi,ui,rx))

                ######### 最小二乗法 ######## 
                counter = 0
                sorted_atom_distance_list = sorted(atom_distance_list, key=lambda x: np.linalg.norm(x[2]))
                for xi,ui,rx in tqdm(sorted_atom_distance_list):
                    f_MD_disp_data.write(f"{np.linalg.norm(rx)} {np.linalg.norm(ui)} {ui[0]} {ui[1]} {ui[2]}\n")
                    if np.linalg.norm(rx)>r_range_limit:
                        green = np.zeros((3, 3, 3))
                        counter += 1
                        green,r = kelvin_solution(g, v, rx, green)

                        for i in range(9):
                            for j in range(9):
                                for k in range(3):
                                    a[i][j] += (
                                        green[k][iidx[i]][jidx[i]] * green[k][iidx[j]][jidx[j]]
                                    )

                        for i in range(9):
                            for j in range(3):
                                b[i] -= ui[j] * green[j][iidx[i]][jidx[i]]

                    else:
                        continue
                
                ######### 逆行列の計算 (特異行列の場合はnan) ########
                try:       
                    a_inv = np.linalg.inv(a)#正方行列aの逆行列inverse
                    x = np.dot(a_inv, b)#行列積の計算
                    x_list.append(x)
                except np.linalg.LinAlgError:
                    print(f"[Warning] matrix 'a' is singular at file index [{l}][{m}], skip this case.")
                    x_dammy = np.full(9, np.nan)
                    x_list.append(x_dammy)
                    continue 

    ######## ファイルに出力 "p_j", "p_eV", "green_tensor"
    with open(f"{output_dir}/force_dipole/data_p_J.inp", "w") as f_J,\
         open(f"{output_dir}/force_dipole/data_p_eV.inp", "w") as f_eV,\
         open(f"{output_dir}/force_dipole/green.txt","w") as f_g:
        for l in range(file_num*len(limit_range)):
            # xx yy zz xy xz yz
            f_J.write(f"{x_list[l][0]:.12e} {x_list[l][4]:.12e} {x_list[l][8]:.12e} {x_list[l][1]:.12e} {x_list[l][2]:.12e} {x_list[l][5]:.12e}")
            f_J.write("\n")

            f_eV.write(f"{x_list[l][0]*conv_J_to_eV:.12e} {x_list[l][4]*conv_J_to_eV:.12e} {x_list[l][8]*conv_J_to_eV:.12e} {x_list[l][1]*conv_J_to_eV:.12e} {x_list[l][2]*conv_J_to_eV:.12e} {x_list[l][5]*conv_J_to_eV:.12e}")
            f_eV.write("\n")

        for i in range(3):
            for j in range(3):
                for k in range(3):
                    f_g.write(f"G[{i}][{j}][{k}] = {green[i][j][k]}")
                    f_g.write("\n")


if __name__ == "__main__":
    main()



