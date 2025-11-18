import numpy as np #type: ignore

import matplotlib.pyplot as plt #type: ignore

def main():
    with open(f"/Users/kyou/Library/CloudStorage/Box-Box/force_dipole/Al/vacancy/displacement/data_p_eV.inp", "r") as f_eV,\
         open(f"/Users/kyou/Library/CloudStorage/Box-Box/force_dipole/Al/vacancy/displacement/data_p_J.inp", "r") as f_J,\
         open(f"/Users/kyou/Library/CloudStorage/Box-Box/lammps/pointdefects_Al/N_atom.txt", "r") as f_N:
        
        lines_eV = f_eV.readlines()
        lines_J = f_J.readlines()
        lines_N = f_N.readlines()
        n = len(lines_eV)
        n_N = len(lines_N)
        
        pxx_eV = []
        pyy_eV = []
        pzz_eV = []
        pxx_J = []
        pyy_J = []
        pzz_J = []
        atoms = []

        for i in range(n):
            p1_eV = float(lines_eV[i].split()[0])
            p2_eV = float(lines_eV[i].split()[1])
            p3_eV = float(lines_eV[i].split()[2])
            p1_J = float(lines_J[i].split()[0])
            p2_J = float(lines_J[i].split()[1])
            p3_J = float(lines_J[i].split()[2])
            pxx_eV.append(p1_eV)
            pyy_eV.append(p2_eV)
            pzz_eV.append(p3_eV)
            pxx_J.append(p1_J)
            pyy_J.append(p2_J)
            pzz_J.append(p3_J)

        for i in range(n_N-n, n_N): #←計算が完了してる場合
            a = int(lines_N[i].split()[0])
            atoms.append(a)
            print(a)

        # for i in range(n_N-(n+1), n_N-1):  #←計算が完了して無い場合
        #     a = int(lines_N[i].split()[0])
        #     atoms.append(a)
        #     print(a)

    plt.plot(atoms, pxx_eV, marker="^", mfc="#0075c2", markersize = 8, mec="black", label = "$P_{11}$", color="black",ls="--",linewidth =1.0)
    plt.plot(atoms, pxx_eV, marker="o", mfc="#33ff33", markersize = 8, mec="black", label = "$P_{22}$", color="black",ls=":",linewidth =1.0)
    plt.plot(atoms, pzz_eV, marker="v", mfc="#ea5549", markersize = 8, mec="black", label = "$P_{33}$", color="black",ls="-",linewidth =1.0)
    # plt.axhline(-0.65, color='black', linestyle='--', linewidth=1.0, label = "※")
    # plt.axhline(-0.79, color='black', linestyle='-', linewidth=1.0,label = "※")
    plt.legend()
    plt.xlabel("Number of atoms[-]", fontsize=14)
    plt.ylabel("Force dipole [eV]", fontsize=14)
    plt.legend(fontsize=14)
    plt.tick_params(labelsize=14)
    plt.grid(color='lightgray', linestyle=':', linewidth=0.8)
    plt.show()

if __name__ == "__main__":
    main()          
     

