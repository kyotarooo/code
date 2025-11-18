import numpy as np #type: ignore

import matplotlib.pyplot as plt #type: ignore

def main():
    with open(f"/Users/kyou/Library/CloudStorage/Box-Box/lammps/pointdefects_Zr/force_dipole_eV/force_dipole_eV.txt", "r") as f_eV,\
         open(f"/Users/kyou/Library/CloudStorage/Box-Box/lammps/pointdefects_Zr/force_dipole_J/force_dipole_J.txt", "r") as f_J,\
         open(f"/Users/kyou/Library/CloudStorage/Box-Box/lammps/pointdefects_Zr/N_atom.txt", "r") as f_N:
        
        lines_eV = f_eV.readlines()
        lines_J = f_J.readlines()
        lines_N = f_N.readlines()
        n = len(lines_eV)
        n_N = len(lines_N)
        pxx_eV = []
        pzz_eV = []
        pxx_J = []
        pzz_J = []
        atoms = []

        for i in range(n):
            p1_eV = float(lines_eV[i].split()[0])
            p3_eV = float(lines_eV[i].split()[2])
            p1_J = float(lines_J[i].split()[0])
            p3_J = float(lines_J[i].split()[2])
            pxx_eV.append(p1_eV)
            pzz_eV.append(p3_eV)
            pxx_J.append(p1_J)
            pzz_J.append(p3_J)

        for i in range(n_N-n, n_N):
            a = int(lines_N[i].split()[0])
            atoms.append(a)
            print(a)





    plt.plot(atoms, pxx_eV, marker="^", markersize = 7, label = "p11", color='#0075c2',ls="--",linewidth = 0.5)
    plt.plot(atoms, pzz_eV, marker="v", markersize = 7, label = "p33", color='#ea5549',ls="--",linewidth = 0.5)
    #EAM#3
    # plt.axhline(-5.43, color='#0068b7', linestyle='-', linewidth=0.7, label = "※")
    # plt.axhline(-5.51, color='#ea5550', linestyle='-', linewidth=0.7,label = "※")
    #EAM#2
    plt.axhline(-0.65, color='#0068b7', linestyle='-', linewidth=0.7, label = "※")
    plt.axhline(-0.79, color='#ea5550', linestyle='-', linewidth=0.7,label = "※")
    plt.legend()
    plt.xlabel("atoms")
    plt.ylabel("force dipole[eV]")
    plt.show()

if __name__ == "__main__":
    main()          
     
