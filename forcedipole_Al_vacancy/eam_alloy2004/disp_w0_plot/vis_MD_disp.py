import numpy as np #type: ignore
import matplotlib.pyplot as plt #type: ignore
import os

######## 環境変数の取得 ########
output_dir = os.environ.get("OUTPUT_PATH")
if output_dir is None:
    raise ValueError("環境変数 OUTPUT_PATH error")

######## MD dislacement data(U)-ファイル名の取得 ########
items_U = os.listdir(f"{output_dir}/disp_data")
num = len(items_U)
Udata_name_list = []
for i in range(len(items_U)):
    Udata_name_list.append(f"MD_disp_data{i + 1}.txt")

######## theorical solution data(u)-ファイル名の取得 ########
items_u = os.listdir(f"{output_dir}/ux_data")
udata_name_list = []
for i in range(len(items_u)):
    udata_name_list.append(f"ux_data{int(i / 4)+1}_{i % 4}")
    
######## 材料定数の取得 ############
lattice_const = os.environ.get("LATTICE_CONST") # ang
lattice_const = float(lattice_const) * 1.0e-10 # m

######## main処理 ########
def main():
    # ==== list_listの定義 ====
    abs_ur_list_list = []
    abs_Ur_list_list = []
    ux_list_list = []
    uy_list_list = []
    uz_list_list = []
    Ux_list_list = []
    Uy_list_list = []
    Uz_list_list = []
    # ==== theorical dataでループ ====
    for filename in udata_name_list:
        with open(f"{output_dir}/ux_data/{filename}", "r") as f_u:

            # ==== listの定義 ====
            abs_ur_list = []
            ux_list = []
            uy_list = []
            uz_list = []
            lines = f_u.readlines()

            num_lines = len(lines)
            for i in range(num_lines):
                abs_r, rx, ry, rz, ux, uy, uz = [float(v) for v in lines[i].split()[:7]]
                abs_ur_list.append((abs_r / lattice_const))
                if rx > 0:
                    ux_list.append((-1*ux) / lattice_const)
                else:
                    ux_list.append(ux / lattice_const)
                if ry > 0:
                    uy_list.append((-1*uy) / lattice_const)
                else:
                    uy_list.append(uy / lattice_const)
                if rz > 0:
                    uz_list.append((-1*uz) / lattice_const)
                else:
                    uz_list.append(uz / lattice_const)

        abs_ur_list_list.append(abs_ur_list)
        ux_list_list.append(ux_list)
        uy_list_list.append(uy_list)
        uz_list_list.append(uz_list)

    # ==== MD dataでループ ====
    for filename in Udata_name_list:
        with open(f"{output_dir}/disp_data/{filename}", "r") as f_U:

            # ==== listの定義 ====
            abs_Ur_list = []
            Ux_list = []
            Uy_list = []
            Uz_list = []
            lines = f_U.readlines()

            num_lines = len(lines)
            for i in range(num_lines):
                abs_r, rx, ry, rz, ux, uy, uz = [float(v) for v in lines[i].split()[:7]]
                abs_Ur_list.append((abs_r / lattice_const))
                if rx > 0:
                    Ux_list.append((-1*ux) / lattice_const)
                else:
                    Ux_list.append(ux / lattice_const)
                if ry > 0:
                    Uy_list.append((-1*uy) / lattice_const)
                else:
                    Uy_list.append(uy / lattice_const)
                if rz > 0:
                    Uz_list.append((-1*uz) / lattice_const)
                else:
                    Uz_list.append(uz / lattice_const)
        
        abs_Ur_list_list.append(abs_Ur_list)
        Ux_list_list.append(Ux_list)
        Uy_list_list.append(Uy_list)
        Uz_list_list.append(Uz_list)

    for i in range(num):

        plt.plot(abs_Ur_list_list[i] ,Ux_list_list[i], "^",mfc="#0075c2", markersize = 6, mec = "black", label = "$U_{1}$",  color="black")
        # plt.plot(abs_Ur_list_list[i] ,Uy_list_list[i], "o",mfc="#33ff33", markersize = 4, mec = "black", label = "$u_{22}$",  color="black")
        # plt.plot(abs_Ur_list_list[i] ,Uz_list_list[i], "v",mfc="#ea5549", markersize = 4, mec = "black", label = "$u_{33}$",  color="black")
        plt.plot(abs_ur_list_list[i * 4 + 0], ux_list_list[i * 4 + 0], "s",mfc="#808080", markersize = 4, mec = "black", label = "$u_{1}(r_{excl}/a = 0)$",  color="black")
        plt.plot(abs_ur_list_list[i * 4 + 1], ux_list_list[i * 4 + 1], "d",mfc="#c0c0c0", markersize = 4, mec = "black", label = "$u_{1}(r_{excl}/a = 1)$",  color="black")
        plt.plot(abs_ur_list_list[i * 4 + 2], ux_list_list[i * 4 + 2], "H",mfc="#dcdcdc", markersize = 4, mec = "black", label = "$u_{1}(r_{excl}/a = 1.5)$",  color="black")
        plt.plot(abs_ur_list_list[i * 4 + 3], ux_list_list[i * 4 + 3], "o",mfc="#ffffff", markersize = 4, mec = "black", label = "$u_{1}(r_{excl}/a = 2.0)$",  color="black")
        plt.axvline(x=1.0, color='#808080', linestyle='--', linewidth=0.5)
        plt.axvline(x=1.5, color='#808080', linestyle='--', linewidth=0.5)
        plt.axvline(x=2.0, color='#808080', linestyle='--', linewidth=0.5)
        plt.xlabel("$Distance / a [-]$")
        plt.ylabel("$Displacement / a[-]$")
        plt.legend()
        # === png 保存 ===
        save_dir = f"{output_dir}/MD_data_plot"
        os.makedirs(save_dir, exist_ok=True)
        plt.savefig(f"{save_dir}/disp_comparison_{i}.png", dpi=300)
        plt.close()

if __name__ == "__main__":
    main()          
     
