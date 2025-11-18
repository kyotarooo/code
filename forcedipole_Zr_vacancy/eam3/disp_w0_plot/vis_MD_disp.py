import numpy as np #type: ignore
import matplotlib.pyplot as plt #type: ignore
import os

######## 環境変数の取得 ########
output_dir = os.environ.get("OUTPUT_PATH")
if output_dir is None:
    raise ValueError("環境変数 OUTPUT_PATH error")

######## MD dislacement data数の取得 ########
items = os.listdir(f"{output_dir}/disp_data")
data_num = len(items)

######## main処理 ########
def main():
    # ==== listの定義 ====
    abs_r_list = []
    abs_u_list = []
    ux_list = []
    uy_list = []
    uz_list = []
    # ==== MDdataでループ ====
    for filename in items:
        with open(f"{output_dir}/disp_data/{filename}", "r") as f:
            lines = f.readlines()
            num_lines = len(lines)
            for i in range(num_lines):
                abs_r, abs_u, ux, uy, uz = [float(v) for v in lines[i].split()[:5]]
                abs_r_list.append(abs_r)
                abs_u_list.append(abs_u)
                ux_list.append(ux)
                uy_list.append(uy)
                uz_list.append(uz)
    
    plt.plot(abs_r_list, ux_list, "^",mfc="#0075c2", markersize = 4, mec = "black", label = "$u_{11}$",  color="black")
    plt.plot(abs_r_list, uy_list, "o",mfc="#33ff33", markersize = 4, mec = "black", label = "$u_{22}$",  color="black")
    plt.plot(abs_r_list, uz_list, "v",mfc="#ea5549", markersize = 4, mec = "black", label = "$u_{33}$",  color="black")
    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()          
     
