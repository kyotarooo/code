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
cutoff_num = 4

############ Load input file count ############
with open(f"{output_dir}/4hsic_q/filename.txt", 'r') as f_num:
    lines = f_num.readlines()
    file_num = len(lines)
                
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




############ Kelvin Solition (wei cai) ############ ←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←
def kelvin_solution(g, v, rx, green, a_core_width):
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
    a_core_width = [0, 2e-11, 5e-11, 8e-11, 1e-10, 3e-10]
    a_num = len(a_core_width)

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
        
    # ====== Compute RMSE from Residual Stress ======
    rmse_residual = {i: [] for i in range(a_num)}
    
    for a in range(len(a_core_width)):
        for la in range(num_P_r):
            sum = 0
            for xi, ui, rx in sorted_atom_distance_list:
                abs_rx = np.linalg.norm(rx)
                if abs_rx < 1e-15:  
                    continue
                green = np.zeros((3, 3, 3))
                green = kelvin_solution(g, v, rx, green, a_core_width[a])
                ux=uy=uz=0.0
                for j in range(9):
                    ux -= green[0][iidx[j]][jidx[j]]*P_r[la][j]
                    uy -= green[1][iidx[j]][jidx[j]]*P_r[la][j]
                    uz -= green[2][iidx[j]][jidx[j]]*P_r[la][j]
                ui_residual = [ux, uy, uz]
                ui_sub = ui - ui_residual
                ui_sub = np.array(ui_sub) / lattice_const
                sum = sum + np.dot(ui_sub, ui_sub)
            rmse_residual[a].append(np.sqrt(sum / atom_number))

    # ############# Plot RMSE comparison: Kelvin(cut_off) vs residual ############
    color_list= plt.cm.Blues(np.linspace(0.2, 1.0, len(a_core_width)))
    marker_list = ["^", "v", "s", "h", "8", "o"]
    plt.figure()
    plt.rcParams["mathtext.fontset"] = "stix"  # STIXフォントはTimes系
    plt.rcParams["font.family"] = "STIXGeneral"
    for a in range(len(a_core_width)):
        plt.plot(atoms, rmse_residual[a], label = f"b = {a_core_width[a]}[m]", marker = marker_list[a], linewidth = 1.0, color = "black",  markersize = 7,  mfc=color_list[a], ls = "-")
    plt.xlabel("Number of atoms in the supercell [-]")
    plt.ylabel("Root mean square error (RMSE) [-]")
    #plt.title("RMSE Comparison: Displacement Method (cutoff) vs Residual-Stress Method")

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
    save_dir = f"{output_dir}/rmse_graph/rmse_graph_{atom_type_to_delete}/b"
    os.makedirs(save_dir, exist_ok=True)
    plt.savefig(f"{save_dir}/rsme_disp.png", dpi=300)
    plt.close()
       
            

        

if __name__ == "__main__":
     main()


        
             

             

             


             

            