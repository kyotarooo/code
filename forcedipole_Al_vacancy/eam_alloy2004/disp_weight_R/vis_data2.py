import numpy as np #type: ignore
import matplotlib.pyplot as plt #type: ignore 
import os

############ Get working directory path ############
output_dir = os.environ.get("OUTPUT_PATH")
if output_dir is None:
    raise ValueError("環境変数 OUTPUT_PATH error")

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
with open(f"{output_dir}/infilename.txt", 'r') as f_num:
    lines = f_num.readlines()
    file_num = len(lines[0].split())

############ Load Force Dipole tensor(J) ############
with open(f"{output_dir}/force_dipole/data_p_J.inp", "r") as f_P,\
    open(f"{output_dir}/force_dipole/mode.txt", "w") as f_mode:
    lines = f_P.readlines()
    cutoff_num = 4
    P = {i: [] for i in range(cutoff_num)} 

    for i in range(len(lines)):
        values = np.array([float(v) for v in lines[i].split()])
        mod = i % cutoff_num
        P[mod].append(values)

############ Load The number of atoms ############
with open(f"{output_dir}/N_atom_perfect.txt", "r") as f_atoms:
    atom = f_atoms.readlines()
    for i in range(len(atom)):
        atoms = np.array([float(v) for v in atom[i].strip()])

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
with open(f"{output_dir}/supercell_size.txt", "r") as f_cell:
    lines = f_cell.readlines()
    super_cell = []
    n = 200  
    for line in lines[-file_num:]:
        values = np.array([float(v) for v in line.split()]) * conv_ang_to_m
        cell_size = np.min(values)  
        x = np.linspace(0, cell_size, n)
        points = np.column_stack((x, x, x))
        super_cell.append(points)

       
############ Kelvin Solition ############
def kelvin_solution(g, v, rx, green):
    r = np.linalg.norm(rx)
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
        with open(f"{output_dir}/dump/displacement/displacement{l}.inp", "r") as f_disp:
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


        ############ Generate diplacement graph from Kelvin Solution ############
        # ====== Generate mesh distance list ======
        for xi in super_cell[l]:
            mesh_atom_distance_list.append(xi - xc)
        sorted_mesh_atom_distance_list = sorted(mesh_atom_distance_list,key=lambda x: np.linalg.norm(x))

        # ====== Compute displacement from Kelvin solution ======
        abs_rx_list_green = {i: [] for i in range(file_num)}
        ur_list_green = {i: [] for i in range(file_num)}

        for la in range(file_num):
            for rx in sorted_mesh_atom_distance_list:
                abs_rx = np.linalg.norm(rx)
                if abs_rx < 1e-15:  
                    continue
                abs_rx_list_green[la].append(abs_rx)
                nr = rx / abs_rx
                green = np.zeros((3, 3, 3))
                green = kelvin_solution(g, v, rx, green)
                ux=uy=uz=0.0
                for j in range(9):
                    ux -= green[0][iidx[j]][jidx[j]]*P[0][la][j]
                    uy -= green[1][iidx[j]][jidx[j]]*P[0][la][j]
                    uz -= green[2][iidx[j]][jidx[j]]*P[0][la][j]
                ui = [ux, uy, uz]
                ur = (-1)*np.dot(nr, ui)
                ur_list_green[la].append(ur)
            abs_rx_list_green[la] = np.array(abs_rx_list_green[la]) / lattice_const
            ur_list_green[la] = np.array(ur_list_green[la]) / lattice_const

        # ############# Plot comparison: MD vs Kelvin() ############
        label_list = []
        for i in range(file_num):
            label_list.append(f"$\\it{{Number\\ of\\ atoms}}$ = {atom[i]}")
        plt.figure()
        plt.plot(abs_rx_list, ur_list, "o",mfc="#ff0000", markersize = 2.5, mec = "black", label = "$Ur$",  color="black", mew=0.5)
        for la in range(file_num):
            plt.plot(abs_rx_list_green[la], ur_list_green[la], "-",label = label_list[la], linewidth=1.0)
        plt.axhline(y=0.0, color='#808080', linestyle='--', linewidth=0.5)
        plt.xlabel(r"$\it{Distance} \,/\, a \, [-]$")
        plt.ylabel(r"$\it{Displacement} \,/\, a \, [-]$")
        plt.title(f"$\\it{{Number\\ of\\ atoms}}$ = {n}")
        plt.ylim([-0.001,0.014])
        plt.legend()
        # --- make dir & save png ---
        save_dir = f"{output_dir}/MD_vs_green(all)_graph"
        os.makedirs(save_dir, exist_ok=True)
        plt.savefig(f"{save_dir}/green_disp{n}.png", dpi=300)
        plt.close()

    
        

if __name__ == "__main__":
     main()


        
             

             

             


             

            