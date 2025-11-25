# =========================================================================================
# Non-singular Green’s function (Kelvin solution) following:
#   W. Cai et al., "A non-singular continuum theory of dislocations"
#
# Replace singular radius:
#       r = |rx|
# with non-singular form:
#       r = sqrt(|rx|^2 + a^2)
#
# Then the displacement field is blended:
#       w̃(rx) = (1 - m) w(rx; a1) + m w(rx; a2)
#
# where:
#   a1 = 0.9038 * |b|
#   a2 = 0.5451 * |b|
#   m  = 0.6575
# =========================================================================================

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

######## 削除する原子の指定 ########  
atom_type_to_delete = os.environ.get("ATOM")
if atom_type_to_delete is None:
    raise ValueError("環境変数 ATOM error")

######## コア幅 ########  ←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←
a_core_width = float(os.environ.get("Burgers"))
if a_core_width is None:
    raise ValueError("環境変数 core error")



    
######### ファイル数の読み込み ########
with open(f"{output_dir}/4hsic_vacancy/atom_type_delete_{atom_type_to_delete}/filename.txt", 'r') as f_num:
    lines = f_num.readlines()
    file = lines[0].split()
    file_num = len(file)

############ Load Force Dipole tensor (displacement) [J] ############
with open(f"{output_dir}/force_dipole/force_dipole_{atom_type_to_delete}/a_{a_core_width}/data_p_J.inp", "r") as f_P,\
    open(f"{output_dir}/force_dipole/mode.txt", "w") as f_mode:
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

############ Load The number of atoms ############
with open(f"{output_dir}/N_atom_perfect.txt", "r") as f_atoms:
    atom = f_atoms.readlines()
    for i in range(len(atom)):
        atoms = [float(line.strip()) for line in atom]

# ############ Generate supercell mesh points ############
# with open(f"{output_dir}/supercell_size.txt", "r") as f_cell:
#      lines = f_cell.readlines()
#      super_cell = []
#      n = 30
#      for line in lines[-file_num:]:
#         values = np.array([float(v) for v in line.split()]) * conv_ang_to_m
#         x = np.linspace(0, values[0], n)
#         y = np.linspace(0, values[1], n)
#         z = np.linspace(0, values[2], n)
#         X, Y, Z = np.meshgrid(x, y, z, indexing = "ij")
#         points = np.column_stack((X.ravel(), Y.ravel(), Z.ravel()))
#         super_cell.append(points)

############ Generate supercell points along diagonal ############
with open(f"{output_dir}/supercell.txt", "r") as f_cell:
    lines = f_cell.readlines()
    super_cell = []
    n = 200  
    for line in lines[-file_num:]:
        values = np.array([float(v) for v in line.split()]) * conv_ang_to_m
        cell_size = np.min(values)  
        x = np.linspace(0, cell_size, n)
        points = np.column_stack((x, x, x))
        super_cell.append(points)


############ Kelvin Solition (wei cai) ############ ←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←
def kelvin_solution(g, v, rx, green, a_cai):
    r = np.sqrt(np.linalg.norm(rx) ** 2 + a_cai ** 2)
    rx_norm = rx / r  # ← コピーを作る
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
    for l in range(file_num):
        # --- Initialization ---
        a = np.zeros((9, 9))
        b = np.zeros(9)
        iidx = np.array([0, 0, 0, 1, 1, 1, 2, 2, 2])
        jidx = np.array([0, 1, 2, 0, 1, 2, 0, 1, 2])
        atom_distance_list = []
        mesh_atom_distance_list = []

        ############ Generate diplacement graph from MD ############
        # ====== Load MD displacement data ======
        with open(f"{output_dir}/dump/displacement/displacement_{atom_type_to_delete}/displacement{l}.inp", "r") as f_disp:
            lines = f_disp.readlines()
            g = float((lines[0].split())[0])
            v = float(lines[0].split()[1])
            xc = np.array([float(lines[1].split()[i]) for i in range(3)])#xc: defect coordinate
            xc_list.append(xc)
            n = int(lines[2].split()[0])
            idx = 3
            for _ in range(n):
                    xi = np.array([float(lines[idx].split()[i]) for i in range(3)])
                    ui = np.array([float(lines[idx].split()[i + 3]) for i in range(3)])
                    idx = idx + 1
                    rx = xi - xc
                    atom_distance_list.append((xi,ui,rx))
        sorted_atom_distance_list=sorted(atom_distance_list, key=lambda x: np.linalg.norm(x[2]))

       # ====== Compute radial displacement from MD ======
        ur_list = []        
        abs_rx_list = []    
        for xi, ui, rx in sorted_atom_distance_list:
             abs_rx = np.linalg.norm(rx)
             nr = rx / abs_rx # radial unit vector
             ur_list.append((-1)*np.dot(nr,ui))
             abs_rx_list.append(abs_rx)
        ur_list = np.array(ur_list) / lattice_const
        abs_rx_list = np.array(abs_rx_list) / lattice_const

        ############ Generate supercell points along diagonal ############
        with open(f"{output_dir}/supercell.txt", "r") as f_cell:
            lines = f_cell.readlines()
            super_cell = []
            n = 200  
            values = np.array([float(v) for v in lines[l].split()]) * conv_ang_to_m
            x = np.linspace(xc[0], values[0], n)
            y = np.linspace(xc[1], values[1], n)
            z = np.linspace(xc[2], values[2], n)
            super_cell = np.column_stack((x, y, z))
            

        ############ Generate diplacement graph from Kelvin Solution ############
        # ====== Generate mesh distance list ======
        for xi in super_cell:
            mesh_atom_distance_list.append(xi - xc)
        sorted_mesh_atom_distance_list = sorted(mesh_atom_distance_list,key=lambda x: np.linalg.norm(x))
        

        # ====== Compute displacement from Kelvin solution ======
        abs_rx_list_green = {i: [] for i in range(file_num)}
        ur_list_green = {i: [] for i in range(file_num)}
        a_cai = {1: 0.9038 * a_core_width, 2: 0.5451 * a_core_width}
        m_weight = 0.6575

        for la in range(file_num):
            for rx in sorted_mesh_atom_distance_list:
                abs_rx = np.linalg.norm(rx)
                if abs_rx < 1e-15:  
                    continue
                abs_rx_list_green[la].append(abs_rx)
                nr = rx / abs_rx
                
                # w^~の適用   --------------------------------------
                green_1_d = np.zeros((3, 3, 3))
                green_1_d = kelvin_solution(g, v, rx, green_1_d, a_cai[1])
                ux_1_d=uy_1_d=uz_1_d=0.0
                for j in range(9):
                    ux_1_d -= green_1_d[0][iidx[j]][jidx[j]]*P[0][la][j]
                    uy_1_d -= green_1_d[1][iidx[j]][jidx[j]]*P[0][la][j]
                    uz_1_d -= green_1_d[2][iidx[j]][jidx[j]]*P[0][la][j]
                ui_1 = np.array([ux_1_d, uy_1_d, uz_1_d])
                
                green_2_d = np.zeros((3, 3, 3))
                green_2_d = kelvin_solution(g, v, rx, green_2_d, a_cai[2])
                ux_2_d=uy_2_d=uz_2_d=0.0
                for j in range(9):
                    ux_2_d -= green_2_d[0][iidx[j]][jidx[j]]*P[0][la][j]
                    uy_2_d -= green_2_d[1][iidx[j]][jidx[j]]*P[0][la][j]
                    uz_2_d -= green_2_d[2][iidx[j]][jidx[j]]*P[0][la][j]
                ui_2 = np.array([ux_2_d, uy_2_d, uz_2_d])
                ui = (1-m_weight) * ui_1 + m_weight * ui_2
                # -------------------------------------------------
                
                ur = (-1)*np.dot(nr, ui)
                ur_list_green[la].append(ur)
            abs_rx_list_green[la] = np.array(abs_rx_list_green[la]) / lattice_const
            ur_list_green[la] = np.array(ur_list_green[la]) / lattice_const

        # ====== Compute displacement from Residual stress ======
        abs_rx_list_residual = {i: [] for i in range(num_P_r)}
        ur_list_residual = {i: [] for i in range(num_P_r)}

        for la in range(num_P_r):
            for rx in sorted_mesh_atom_distance_list:
                abs_rx = np.linalg.norm(rx)
                if abs_rx < 1e-15:  
                    continue
                abs_rx_list_residual[la].append(abs_rx)
                nr = rx / abs_rx
                
                # w^~の適用   --------------------------------------
                green_1_r = np.zeros((3, 3, 3))
                green_1_r = kelvin_solution(g, v, rx, green_1_r, a_cai[1])
                ux_1_r=uy_1_r=uz_1_r=0.0
                for j in range(9):
                    ux_1_r -= green_1_r[0][iidx[j]][jidx[j]]*P_r[la][j]
                    uy_1_r -= green_1_r[1][iidx[j]][jidx[j]]*P_r[la][j]
                    uz_1_r -= green_1_r[2][iidx[j]][jidx[j]]*P_r[la][j]
                ui_1 = np.array([ux_1_r, uy_1_r, uz_1_r])
                
                green_2_r = np.zeros((3, 3, 3))
                green_2_r = kelvin_solution(g, v, rx, green_2_r, a_cai[2])
                ux_2_r=uy_2_r=uz_2_r=0.0
                for j in range(9):
                    ux_2_r -= green_2_r[0][iidx[j]][jidx[j]]*P_r[la][j]
                    uy_2_r -= green_2_r[1][iidx[j]][jidx[j]]*P_r[la][j]
                    uz_2_r -= green_2_r[2][iidx[j]][jidx[j]]*P_r[la][j]
                ui_2 = np.array([ux_2_r, uy_2_r, uz_2_r])
                ui = (1-m_weight) * ui_1 + m_weight * ui_2
                # -------------------------------------------------
                
                ur = (-1)*np.dot(nr, ui)
                ur_list_residual[la].append(ur)
            abs_rx_list_residual[la] = np.array(abs_rx_list_residual[la]) / lattice_const
            ur_list_residual[la] = np.array(ur_list_residual[la]) / lattice_const


        # ############# Plot comparison: MD vs Kelvin(all) ############
        label_list = []
        num_plot = 4  # プロットしたい本数
        plot_idx = np.linspace(2, file_num - 1, num_plot, dtype=int)

        for i in range(file_num):
            label_list.append(f"$\\it{{{{Number\\ of\\ atoms}}}} = {atoms[i]:.0f}$")
        plt.figure()
        plt.rcParams["mathtext.fontset"] = "stix"  # STIXフォントはTimes系
        plt.rcParams["font.family"] = "STIXGeneral"
        plt.plot(abs_rx_list, ur_list, "o",mfc="#ff0000", markersize = 2.5, mec = "black", label = "$Ur$",  color="black", mew=0.5)
        for la in plot_idx:
            plt.plot(abs_rx_list_green[la][2:], ur_list_green[la][2:], "-",label = label_list[la], linewidth=1.0)
        plt.axhline(y=0.0, color='#808080', linestyle='--', linewidth=0.5)
        plt.xlabel(r"$\it{Distance} \,/\, a \, [-]$")
        plt.ylabel(r"$\it{Displacement} \,/\, a \, [-]$")
        plt.title(f"$\\it{{Number\\ of\\ atoms}}$ = {atoms[l]}")
        plt.ylim([-0.005,0.005])
        plt.grid(True, color='gray', linestyle='--', linewidth=0.5)
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
        save_dir = f"{output_dir}/MD_disp/core_width_{a_core_width}"
        os.makedirs(save_dir, exist_ok=True)
        plt.savefig(f"{save_dir}/green_disp{atoms[l]}.png", dpi=300)
        plt.close()

        # ############# Plot comparison: MD vs Kelvin(all) vs residual ############
        label_list = []
        for i in range(file_num):
            label_list.append(f"$\\it{{{{Displacement\\ :\\ N}}}} = {atoms[i]:.0f}$")
        plt.figure()
        plt.rcParams["mathtext.fontset"] = "stix"  # STIXフォントはTimes系
        plt.rcParams["font.family"] = "STIXGeneral"
        colors = plt.cm.Greens(np.linspace(0.2, 1.0, len(atom)))
        plt.plot(abs_rx_list, ur_list, "o",mfc="#ff0000", markersize = 2.5, mec = "black", label = "$Ur\\ (MD)$",  color="black", mew=0.5)
        for la in range(file_num):
            plt.plot(abs_rx_list_green[la], ur_list_green[la], "-",label = label_list[la], linewidth=1.0, color=colors[la])
        plt.plot(abs_rx_list_residual[6], ur_list_residual[6], "-",label = "$\\it{Residual\\ Stress}$", linewidth=1.0, color = "#0000ff")
        plt.axhline(y=0.0, color='#808080', linestyle='--', linewidth=0.5)
        plt.xlabel("Distance / a [-]")
        plt.ylabel("Displacement / a [-]")
        #plt.title(fr"$\it{{Number\ of\ atoms}} = {atoms[l]}\,\,(r_{{\mathrm{{excl}}}}\,=\,0)$")
        plt.ylim([-1.5, 0.02])
        plt.xlim(left = 0.0)
        plt.grid(True, color='gray', linestyle='--', linewidth=0.5)
        plt.legend(
        frameon=True,
        facecolor='white',
        framealpha=0.9,
        fontsize=10,
        #loc='upper right',
        loc='lower right', 
        labelspacing=0.2,     # ← 行間（デフォルト: 0.5）
        handlelength=1.5,     # ← 線の長さ（デフォルト: 2.0）
        handletextpad=0.3,    # ← 線と文字の間隔（デフォルト: 0.8）
        borderaxespad=0.3,    # ← 軸との余白（デフォルト: 0.5）
        borderpad=0.3,        # ← 凡例枠内の余白（デフォルト: 0.4）
        )

        # --- make dir & save png ---
        save_dir = f"{output_dir}/MD_disp_residual/MD_disp_residual_{atom_type_to_delete}/core_width_{a_core_width}"
        os.makedirs(save_dir, exist_ok=True)
        plt.savefig(f"{save_dir}/green_disp{atoms[l]}.png", dpi=300)
        plt.close()

    
        

if __name__ == "__main__":
     main()


        
             

             

             


             

            