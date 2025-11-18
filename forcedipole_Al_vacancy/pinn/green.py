import numpy as np # type: ignore
import sys
from tqdm import tqdm#type: ignore

######### 単位変換係数 ########
conv_eV_to_J = 1.602176634E-19
conv_bars_to_Pa = 1.0E+5
conv_Vang_to_Vm = 1.0E-30
conv_barVang_to_NVm = conv_bars_to_Pa * conv_Vang_to_Vm
conv_J_to_eV = 1/conv_eV_to_J

######## ファイル数の読み込み ########
with open(f"/Users/kyou/Library/CloudStorage/Box-Box/lammps/pointdefects_Al/infilename.txt", 'r') as f_num:
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
                green[i][j][k] = c * (
                    -1.0 * rx[i] * d[j][k]
                    - rx[j] * d[i][k]
                    + (3.0 - 4.0 * v) * rx[k] * d[i][j]
                    + 3.0 * rx[i] * rx[j] * rx[k]
                )
    return green


def main():
    for l in range(file_num):
        a = np.zeros((9, 9))
        b = np.zeros(9)
        iidx = np.array([0, 0, 0, 1, 1, 1, 2, 2, 2])
        jidx = np.array([0, 1, 2, 0, 1, 2, 0, 1, 2])
        atom_distance_list = []

        with open(f"/Users/kyou/Library/CloudStorage/Box-Box/dump/vacancy_Al/displacement/displacement{l}.inp", "r") as f_disp: 
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
            
        counter = 0
        sorted_atom_distance_list = sorted(atom_distance_list, key=lambda x: np.linalg.norm(x[2]))
        for xi,ui,rx in tqdm(sorted_atom_distance_list):

            green = np.zeros((3, 3, 3))
            counter += 1
            green = kelvin_solution(g, v, rx, green)

            for i in range(9):
                for j in range(9):
                    for k in range(3):
                        a[i][j] += (
                            green[k][iidx[i]][jidx[i]] * green[k][iidx[j]][jidx[j]]
                        )

            for i in range(9):
                for j in range(3):
                    b[i] -= ui[j] * green[j][iidx[i]][jidx[i]]
            
        a_inv = np.linalg.inv(a)#正方行列aの逆行列inverse
        x = np.dot(a_inv, b)#行列積の計算
            
        with open(f"/Users/kyou/Library/CloudStorage/Box-Box/force_dipole/Al/vacancy/displacement/data_p_J.inp", "a") as f_J,\
             open(f"/Users/kyou/Library/CloudStorage/Box-Box/force_dipole/Al/vacancy/displacement/data_p_eV.inp", "a") as f_eV:
            
            # xx yy zz xy xz yz
            f_J.write(f"{x[0]:.12e} {x[4]:.12e} {x[8]:.12e} {x[1]:.12e} {x[2]:.12e} {x[5]:.12e}")
            f_J.write("\n")

            f_eV.write(f"{x[0]*conv_J_to_eV:.12e} {x[4]*conv_J_to_eV:.12e} {x[8]*conv_J_to_eV:.12e} {x[1]*conv_J_to_eV:.12e} {x[2]*conv_J_to_eV:.12e} {x[5]*conv_J_to_eV:.12e}")
            f_eV.write("\n")

if __name__ == "__main__":
    main()



