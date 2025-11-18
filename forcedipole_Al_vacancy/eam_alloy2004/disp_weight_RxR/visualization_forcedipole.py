import numpy as np #type: ignore
import matplotlib.pyplot as plt #type: ignore
import os

######## 環境変数の取得 ########
output_dir = os.environ.get("OUTPUT_PATH")
if output_dir is None:
    raise ValueError("環境変数 OUTPUT_PATH error")

def main():
    # cutoff半径の数
    num_cutoff = 4

    with open(f"{output_dir}/force_dipole/data_p_eV.inp", "r") as f_eV,\
         open(f"{output_dir}/force_dipole/data_p_J.inp", "r") as f_J,\
         open(f"{output_dir}/N_atom_perfect.txt", "r") as f_N:
        
        lines_eV = f_eV.readlines()
        lines_J = f_J.readlines()
        lines_N = f_N.readlines()
        n = len(lines_eV)
        n_N = len(lines_N)
        
        pxx_eV = {i: [] for i in range(4)}
        pyy_eV = {i: [] for i in range(4)}
        pzz_eV = {i: [] for i in range(4)}
        pxx_J = {i: [] for i in range(4)}
        pyy_J = {i: [] for i in range(4)}
        pzz_J = {i: [] for i in range(4)}
        atoms = []

        for i in range(n):
            mod = i % num_cutoff
            values_eV = lines_eV[i].split()
            values_J = lines_J[i].split()

            # nanを含む場合に対応
            def safe_float(x):
                return np.nan if x.lower() == 'nan' else float(x)

            p1_eV = safe_float(values_eV[0])
            p2_eV = safe_float(values_eV[4])
            p3_eV = safe_float(values_eV[8])

            p1_J = safe_float(values_J[0])
            p2_J = safe_float(values_J[4])
            p3_J = safe_float(values_J[8])

            pxx_eV[mod].append(p1_eV)
            pyy_eV[mod].append(p2_eV)
            pzz_eV[mod].append(p3_eV)
            pxx_J[mod].append(p1_J)
            pyy_J[mod].append(p2_J)
            pzz_J[mod].append(p3_J)

        n = int(len(lines_eV)/num_cutoff)
        for i in range(n_N-n, n_N): #←計算が完了してる場合
            a = int(lines_N[i].split()[0])
            atoms.append(a)
            print(a)

        # for i in range(n_N-(n+1), n_N-1):  #←計算が完了して無い場合
        #     a = int(lines_N[i].split()[0])
        #     atoms.append(a)
        #     print(a)

        # p_cutoff = []
        # for i in range(5):
        #     if i != 1:
        #         p_cutoff.append(pxx_eV[i][2])
        
        r_cutoff = [0.0, 1.0, 1.5, 2.0]

    # plt.plot(atoms, pxx_eV, marker="^", mfc="#0075c2", markersize = 8, mec="black", label = "$P_{11}$", color="black",ls="--",linewidth =1.0)
    # plt.plot(atoms, pxx_eV, marker="o", mfc="#33ff33", markersize = 8, mec="black", label = "$P_{22}$", color="black",ls=":",linewidth =1.0)
    # plt.plot(atoms, pzz_eV, marker="v", mfc="#ea5549", markersize = 8, mec="black", label = "$P_{33}$", color="black",ls="-",linewidth =1.0)
    # plt.axhline(-0.65, color='black', linestyle='--', linewidth=1.0, label = "※")
    # plt.axhline(-0.79, color='black', linestyle='-', linewidth=1.0,label = "※")
    
    plt.plot(atoms, pxx_eV[0], marker="o", mfc="#000000", markersize = 8, label = "$r_{excl} / a = 0.0$", color="black",ls="-",linewidth =1.0)
    plt.plot(atoms, pxx_eV[1], marker="^", mfc="#808080", markersize = 8, label = "$r_{excl} / a = 1.0$", color="black",ls="-",linewidth =1.0)
    plt.plot(atoms, pxx_eV[2], marker="v", mfc="#dcdcdc", markersize = 8, label = "$r_{excl} / a = 1.5$", color="black",ls="-",linewidth =1.0)
    plt.plot(atoms, pxx_eV[3], marker="s", mfc="#f5f5f5", markersize = 8, label = "$r_{excl} / a = 2.0$", color="black",ls="-",linewidth =1.0)
    plt.axhline(-3.197, color='black', linestyle='--', linewidth=1.5, label = "$Residual-stress-method$")

    # plt.plot(r_cutoff, p_cutoff, marker="o", mfc="#000000", markersize = 8, color="black",ls="-",linewidth =2.0)
    # plt.axhline(-3.197, color='black', linestyle='--', linewidth=1.5, label = "$Residual-stress-method$")
    
    plt.legend()
    plt.xlabel("Number of atoms[-]", fontsize=14)
    plt.ylabel("Force dipole [eV]", fontsize=14)
    #plt.ylim(top=-2.0)
    plt.legend(fontsize=14)
    plt.tick_params(labelsize=14)
    plt.grid(color='lightgray', linestyle=':', linewidth=0.8)
    plt.show()


if __name__ == "__main__":
    main()          
     
