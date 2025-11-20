import numpy as np #type: ignore
import matplotlib.pyplot as plt #type: ignore 
import os

############ Get working directory path ############
output_dir = os.environ.get("OUTPUT_PATH")
if output_dir is None:
    raise ValueError("環境変数 OUTPUT_PATH error")

############ Get force dipole(residual stress) directory path ############
dipole_r_path = os.environ.get("DIPOLE_R_PATH")
if dipole_r_path is None:
    raise ValueError("環境変数 DIPOLE_R_PATH error")

######## 削除する原子の指定 ######## 
atom_type_to_delete = int(os.environ.get("ATOM"))
if atom_type_to_delete is None:
    raise ValueError("環境変数 ATOM error")

######## コア幅 ########  ←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←
a_core_width = float(os.environ.get("Burgers"))
if a_core_width is None:
    raise ValueError("環境変数 core error")

############ Define material constant ################
lattice_const = os.environ.get("LATTICE_CONST") # ang
lattice_const = float(lattice_const) * 1.0e-10 # m

############# Unit conversion constants ############
conv_eV_to_J = 1.602176634E-19
conv_bars_to_Pa = 1.0E+5
conv_Vang_to_Vm = 1.0E-30
conv_barVang_to_NVm = conv_bars_to_Pa * conv_Vang_to_Vm
conv_J_to_eV = 1/conv_eV_to_J
conv_ang_to_m = 1.0E-10

############ Load input file count ############
with open(f"{output_dir}/4hsic_q/filename.txt", 'r') as f_num:
    lines = f_num.readlines()
    file_num = len(lines)

############ Load Force Dipole tensor (Displacement)(J) ############
with open(f"{output_dir}/force_dipole/force_dipole_{atom_type_to_delete}/data_p_J.inp", "r") as f_P:
    lines = f_P.readlines()
    cutoff_num = 4
    P = {i: [] for i in range(cutoff_num)} 

    for i in range(len(lines)):
        values = np.array([float(v) for v in lines[i].split()])
        mod = i % cutoff_num
        P[mod].append(values)

                
############ Load Force Dipole tensor (residual_stress) [J] ############
with open(f"{dipole_r_path}/force_dipole/force_dipole_{atom_type_to_delete}/force_dipole_J.txt", "r") as f_r:
    lines = f_r.readlines()
    P_r = []
    for i in range(len(lines)):
        values = np.array([float(v) for v in lines[i].split()])
        P_r.append(values)
    num_P_r = (len(P_r))

    # if num_P_r != file_num:
    #     raise ValueError("The number of force-dipole is different from residual and Displacement")

############ Load The number of atoms ############
with open(f"{output_dir}/N_atom_perfect.txt", "r") as f_atoms:
    atom = f_atoms.readlines()
    for i in range(len(atom)):
        atoms = [float(line.strip()) for line in atom]

print(atoms)


############ Kelvin Solition (wei cai) ############ ←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←
def kelvin_solution(g, v, rx, green):
    r = np.sqrt(np.linalg.norm(rx) ** 2 + a_core_width ** 2)
    rx_norm = rx / r  
    c = 1.0 / (16.0 * np.pi * g * (1.0 - v) * r * r)
    d = np.eye(3)
    for i in range(3):
        for j in range(3):
            for k in range(3):
                green[i][j][k] = (-1) * c * (
                    -1.0 * rx_norm[i] * d[j][k]
                    - rx_norm[j] * d[k][i]
                    + (3.0 - 4.0 * v) * rx_norm[k] * d[i][j]
                    + 3.0 * rx_norm[i] * rx_norm[j] * rx_norm[k]
                )
    return green

def main():
    xc_list = []
    
    # --- Initialization ---
    a = np.zeros((9, 9))
    b = np.zeros(9)
    iidx = np.array([0, 0, 0, 1, 1, 1, 2, 2, 2])
    jidx = np.array([0, 1, 2, 0, 1, 2, 0, 1, 2])
    atom_distance_list = []
    mesh_atom_distance_list = []

    ############ Generate RSME graph ############
    # ====== Load MD displacement data ======
    with open(f"{output_dir}/dump/displacement/displacement_{atom_type_to_delete}/displacement{3}.inp", "r") as f_disp:
        lines = f_disp.readlines()
        g = float((lines[0].split())[0])
        v = float(lines[0].split()[1])
        xc = np.array([float(lines[1].split()[i]) for i in range(3)])#xc: defect coordinate
        xc_list.append(xc)
        atom_number = int(lines[2].split()[0])
        idx = 3
        for _ in range(atom_number):
                xi = np.array([float(lines[idx].split()[i]) for i in range(3)])
                ui = np.array([float(lines[idx].split()[i + 3]) for i in range(3)])
                idx = idx + 1
                rx = xi - xc
                atom_distance_list.append((xi,ui,rx))
    sorted_atom_distance_list = sorted(atom_distance_list, key=lambda x: np.linalg.norm(x[2]))

    # ====== Compute RMSE from Displacement ======
    rmse_disp = {i: [] for i in range(cutoff_num)}

    for cutoff in range(cutoff_num):
        for la in range(file_num):
            sum = 0
            for xi, ui, rx in sorted_atom_distance_list:      
                abs_rx = np.linalg.norm(rx)
                if abs_rx < 1e-15:  
                    continue
                green = np.zeros((3, 3, 3))
                green = kelvin_solution(g, v, rx, green)
                ux=uy=uz=0.0
                for j in range(9):
                    ux -= green[0][iidx[j]][jidx[j]]*P[cutoff][la][j]
                    uy -= green[1][iidx[j]][jidx[j]]*P[cutoff][la][j]
                    uz -= green[2][iidx[j]][jidx[j]]*P[cutoff][la][j]
                ui_green = [ux, uy, uz]
                ui_sub = ui - ui_green
                ui_sub = np.array(ui_sub) / lattice_const
                sum = sum + np.dot(ui_sub, ui_sub)
            rmse_disp[cutoff].append(np.sqrt(sum / atom_number))
        
    # ====== Compute RMSE from Residual Stress ======
    rmse_residual = []
    

    for la in range(num_P_r):
        sum = 0
        for xi, ui, rx in sorted_atom_distance_list:
            abs_rx = np.linalg.norm(rx)
            if abs_rx < 1e-15:  
                continue
            green = np.zeros((3, 3, 3))
            green = kelvin_solution(g, v, rx, green)
            ux=uy=uz=0.0
            for j in range(9):
                ux -= green[0][iidx[j]][jidx[j]]*P_r[la][j]
                uy -= green[1][iidx[j]][jidx[j]]*P_r[la][j]
                uz -= green[2][iidx[j]][jidx[j]]*P_r[la][j]
            ui_residual = [ux, uy, uz]
            ui_sub = ui - ui_residual
            ui_sub = np.array(ui_sub) / lattice_const
            sum = sum + np.dot(ui_sub, ui_sub)
        rmse_residual.append(np.sqrt(sum / atom_number))
    
    print(rmse_residual)


    # ############# Plot RMSE residual ############
    plt.figure()
    plt.rcParams["mathtext.fontset"] = "stix"  # STIXフォントはTimes系
    plt.rcParams["font.family"] = "STIXGeneral"
    plt.plot(atoms, rmse_residual, marker = "h", label=r"$\it{Residual\ Stress\ }$", linewidth = 1.0, mfc = "#ff0000", color = "black", ls = "-", markersize = 7)
    plt.xlabel(r"$\it{Number\, \, of\, \, atoms} \, \, [-]$")
    plt.ylabel(r"$\it{RMSE} \, [-]$")
    # plt.xlim([0.0, 5.5])
    # plt.ylim(top = 0.0050)
    #plt.title(f"$\\it{{Number\\ of\\ atoms}} = {number_atoms}$")
    plt.grid(True, color='gray', linestyle='--', linewidth=0.5)
    plt.tight_layout()
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

        # --- make dir & save png ---
    save_dir = f"{output_dir}/rmse_residual_graph/rmse_residual_graph_{atom_type_to_delete}"
    os.makedirs(save_dir, exist_ok=True)
    plt.savefig(f"{save_dir}/rsme_residual_disp.png", dpi=300)
    plt.close()

    # ############# Plot RMSE comparison: Kelvin(cut_off) vs residual ############
    label_list = [r"$\it{Displacement}\,:\,\,r_{excl}/a = 0.0$", r"$\it{Displacement}\,:\,\,r_{excl}/a = 1.0$", r"$\it{Displacement}\,:\,\,r_{excl}/a = 1.5$", r"$\it{Displacement}\,:\,\,r_{excl}/a = 2.0$"]
    marker_list = ["o", "^", "v", "s"]
    color_list= plt.cm.Blues(np.linspace(0.2, 1.0, cutoff_num))
    plt.figure()
    plt.rcParams["mathtext.fontset"] = "stix"  # STIXフォントはTimes系
    plt.rcParams["font.family"] = "STIXGeneral"
    for la in range(cutoff_num):
        plt.plot(atoms[2:], rmse_disp[la][2:], marker_list[la], label = label_list[la], linewidth = 1.0, color = "black",  markersize = 7,  mfc=color_list[la], ls = "-")
    plt.plot(atoms, rmse_residual, marker = "h", label=r"$\it{Residual\ Stress\ }$", linewidth = 1.0, mfc = "#ff0000", color = "black", ls = "-", markersize = 7)
    plt.xlabel("Number of atoms in the supercell [-]")
    plt.ylabel("Root mean square error (RMSE) [-]")
    plt.title("RMSE Comparison: Displacement Method (cutoff) vs Residual-Stress Method")

    # plt.xlim([0.0, 5.5])
    # plt.ylim(top = 0.0050)
    #plt.title(f"$\\it{{Number\\ of\\ atoms}} = {number_atoms}$")
    plt.grid(True, color='gray', linestyle='--', linewidth=0.5)
    plt.tight_layout()
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

        # --- make dir & save png ---
    save_dir = f"{output_dir}/rmse_graph/rmse_graph_{atom_type_to_delete}"
    os.makedirs(save_dir, exist_ok=True)
    plt.savefig(f"{save_dir}/rsme_disp.png", dpi=300)
    plt.close()
       
            

        

if __name__ == "__main__":
     main()


        
             

             

             


             

            