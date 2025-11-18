import numpy as np #type: ignore
import matplotlib.pyplot as plt #type: ignore
import os

######## 環境変数の取得 ########
output_dir = os.environ.get("OUTPUT_PATH")
if output_dir is None:
    raise ValueError("環境変数 OUTPUT_PATH error")

############ Get force dipole(residual stress) directory path ############
dipole_r_path = os.environ.get("DIPOLE_R_PATH")
if dipole_r_path is None:
    raise ValueError("環境変数 DIPOLE_R_PATH error")

############ Load Force Dipole tensor (residual_stress) [J] ############
with open(f"{dipole_r_path}/force_dipole/force_dipole_eV.txt", "r") as f_r:
    lines = f_r.readlines()
    P_r = []
    for i in range(len(lines)):
        values = np.array([float(v) for v in lines[i].split()])
        P_r.append(values)
    num_P_r = (len(P_r))

def main():
    # cutoff半径の数
    num_cutoff = 4

    with open(f"{output_dir}/force_dipole/data_p_eV.inp", "r") as f_eV,\
         open(f"{output_dir}/force_dipole/data_p_J.inp", "r") as f_J,\
         open(f"{dipole_r_path}/force_dipole/force_dipole_eV.txt", "r") as f_residual,\
         open(f"{output_dir}/N_atom_perfect.txt", "r") as f_N:
        
        lines_eV = f_eV.readlines()
        lines_J = f_J.readlines()
        lines_N = f_N.readlines()
        lines_R = f_residual.readlines()
        n = len(lines_eV)
        n_N = len(lines_N)
        
        pxx_eV = {i: [] for i in range(4)}
        pyy_eV = {i: [] for i in range(4)}
        pzz_eV = {i: [] for i in range(4)}
        pxx_J = {i: [] for i in range(4)}
        pyy_J = {i: [] for i in range(4)}
        pzz_J = {i: [] for i in range(4)}
        pxx_R = []
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
            
        for i in range(num_P_r):
            values_r_eV = lines_R[i].split()
            P1_residual = safe_float(values_r_eV[0])
            pxx_R.append(P1_residual)

        n = int(len(lines_eV)/num_cutoff)
        for i in range(n_N-n, n_N): #←計算が完了してる場合
            a = int(lines_N[i].split()[0])
            atoms.append(a)
            

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
    plt.figure()
    #color_list = ["#000080", "#0000ff", "#1e90ff", "#87cefa"]
    color_list= plt.cm.Blues(np.linspace(0.2, 1.0, len(r_cutoff)))
    label_list = [r"$\it{Displacement}\,:\,\,r_{excl}/a = 0.0$", r"$\it{Displacement}\,:\,\,r_{excl}/a = 1.0$", r"$\it{Displacement}\,:\,\,r_{excl}/a = 1.5$", r"$\it{Displacement}\,:\,\,r_{excl}/a = 2.0$"]
    marker_list = ["o", "^", "v", "s"]
    plt.rcParams["mathtext.fontset"] = "stix"  # STIXフォントはTimes系
    plt.rcParams["font.family"] = "STIXGeneral"
    for la in range(len(r_cutoff)):
            plt.plot(atoms, pxx_eV[la], marker_list[la], label = label_list[la], linewidth = 1.0, color = "black", markersize = 8, mfc=color_list[la], ls = "-")
    plt.plot(atoms, pxx_R, marker="h", mfc="#ff0000", markersize = 8, label=r"$\it{Residual\ Stress\ }$", color="black",ls="-",linewidth =1.0)
    plt.grid(True, color='gray', linestyle='--', linewidth=0.5)
    plt.xlabel(r"$\it{Number\, \, of\, \, atoms} \, \, \, [-]$")
    plt.ylabel(r"$\it{Force\, \, Dipole\, } \, \,[-]$")
    plt.legend(
    frameon=True,
    facecolor='white',
    framealpha=0.9,
    fontsize=10,
    loc='upper right',
    labelspacing=0.2,     # ← 行間（デフォルト: 0.5）
    handlelength=1.5,     # ← 線の長さ（デフォルト: 2.0）
    handletextpad=0.3,    # ← 線と文字の間隔（デフォルト: 0.8）
    borderaxespad=0.3,    # ← 軸との余白（デフォルト: 0.5）
    borderpad=0.3,        # ← 凡例枠内の余白（デフォルト: 0.4）
    )
    # plt.xlim([0.0, 5.5])
    #plt.ylim(top=-2.0)
    # --- make dir & save png ---
    save_dir = f"{output_dir}/force_dipole_graph"
    os.makedirs(save_dir, exist_ok=True)
    plt.savefig(f"{save_dir}/force_dipole.png", dpi=300)
    plt.close()


    # plt.plot(r_cutoff, p_cutoff, marker="o", mfc="#000000", markersize = 8, color="black",ls="-",linewidth =2.0)
    # plt.axhline(-3.197, color='black', linestyle='--', linewidth=1.5, label = "$Residual-stress-method$")
if __name__ == "__main__":
    main()          
     
