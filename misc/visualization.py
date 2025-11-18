import numpy as np#type: ignore
import sys
import matplotlib.pyplot as plt#type: ignore

def main():
    r_limited = [10.0, 15.0,20.0,30.0] #ここは自分で変更していく
    r_limited_array = np.array(r_limited)
    r_limited_array = r_limited_array * 1e-10#単位はメートル
    x_MD_list = []
    y_MD_list = []
    x_green_list = []
    y_green_list = []
    x_green_r_list = []
    y_green_r_list = []
    

    a = 3.615

    for n in range(len(r_limited)):
        x_MD = []
        y_MD = []
        x_green = []
        y_green = []
        x_green_r = []
        y_green_r = []
        with open(f"{sys.argv[1]}/ur_green_r_{r_limited_array[n]}.out","r") as f_green_r,\
             open(f"{sys.argv[1]}/ur_MD_{r_limited_array[n]}.out","r") as f_MD,\
             open(f"{sys.argv[1]}/ur_green_{r_limited_array[n]}.out","r") as f_green:
                    
            lines_green_r = f_green_r.readlines()
            lines_MD = f_MD.readlines()
            lines_green = f_green.readlines()

        for line in lines_green_r:
            x_green_r.append(float(line.split()[0]))
            y_green_r.append(float(line.split()[1]))
        y_green_r = [val / a for val in y_green_r]  # ← 無次元化

        for line in lines_MD:
            x_MD.append(float(line.split()[0]))
            y_MD.append(float(line.split()[1]))
        y_MD = [val / a for val in y_MD]  # ← 無次元化

        for line in lines_green:
            x_green.append(float(line.split()[0]))
            y_green.append(float(line.split()[1]))
        y_green= [val / a for val in y_green]  # ← 無次元化

        x_green_r_list.append(x_green_r)
        y_green_r_list.append(y_green_r)
        x_MD_list.append(x_MD)
        y_MD_list.append(y_MD)
        x_green_list.append(x_green)
        y_green_list.append(y_green) 


    plt.plot(x_green_list[0], y_green_list[0], label=f'green func r={r_limited_array[0]:.2e} [m]', color='#387238')
    plt.plot(x_green_list[1], y_green_list[1], label=f'green func r={r_limited_array[1]:.2e} [m]', color='#55aa55')
    plt.plot(x_green_list[2], y_green_list[2], label=f'green func r={r_limited_array[2]:.2e} [m]', color='#8dc7aa')
    plt.plot(x_green_list[3], y_green_list[3], label=f'green func r={r_limited_array[3]:.2e} [m]', color='#c7e2d5')
    plt.plot(x_green_r_list[0], y_green_r_list[0], label=f'residual stress r={r_limited_array[0]:.2e} [m]', color='#b8d5f1')
    plt.plot(x_green_r_list[1], y_green_r_list[1], label=f'residual stress r={r_limited_array[1]:.2e} [m]', color='#95bfea')
    plt.plot(x_green_r_list[2], y_green_r_list[2], label=f'residual stress r={r_limited_array[2]:.2e} [m]', color='#2b80d5')
    plt.plot(x_green_r_list[3], y_green_r_list[3], label=f'residual stress r={r_limited_array[3]:.2e} [m]', color='#1d558d')
    plt.axhline(0, color='gray', linestyle='-', linewidth=1)  # y=0の線を追加
    plt.scatter(x_MD_list[3], y_MD_list[3], label='MD', marker=".", color='#ff2b2b')
    plt.xlabel("Distance (x-x')[m]")
    plt.ylabel('dislacement/a[-]')          
    plt.legend()
    plt.tight_layout()#グラフのレイアウトを自動的に調整、ラベルやタイトルなどがはみ出さないようにする
    plt.show()

   


if __name__ == "__main__":
    main()


    




#difはいらない気がする！！！！
# 
#       
# import numpy as np#type: ignore
# import sys
# import matplotlib.pyplot as plt#type: ignore

# def main():
#     r_limited = [10.0, 15.0,20.0,30.0] #ここは自分で変更していく
#     r_limited_array = np.array(r_limited)
#     r_limited_array = r_limited_array * 1e-10#単位はメートル
#     x_dif_list = []
#     y_dif_list = []
#     x_MD_list = []
#     y_MD_list = []
#     x_green_list = []
#     y_green_list = []
#     a = 3.615

#     for n in range(len(r_limited)):
#         x_dif = []
#         y_dif = []
#         x_MD = []
#         y_MD = []
#         x_green = []
#         y_green = []
#         with open(f"{sys.argv[1]}/dif_MD_green_{r_limited[n]}.out","r") as f_dif,\
#              open(f"{sys.argv[1]}/ur_MD_{r_limited[n]}.out","r") as f_MD,\
#              open(f"{sys.argv[1]}/ur_green_{r_limited[n]}.out","r") as f_green:
                    
#             lines_dif = f_dif.readlines()
#             lines_MD = f_MD.readlines()
#             lines_green = f_green.readlines()

#         for line in lines_dif:
#             x_dif.append(float(line.split()[0]))
#             y_dif.append(float(line.split()[1]))
#         y_dif = [val / a for val in y_dif]  # ← 無次元化

#         for line in lines_MD:
#             x_MD.append(float(line.split()[0]))
#             y_MD.append(float(line.split()[1]))
#         y_MD = [val / a for val in y_MD]  # ← 無次元化

#         for line in lines_green:
#             x_green.append(float(line.split()[0]))
#             y_green.append(float(line.split()[1]))
#         y_green= [val / a for val in y_green]  # ← 無次元化

#         x_dif_list.append(x_dif)
#         y_dif_list.append(y_dif)
#         x_MD_list.append(x_MD)
#         y_MD_list.append(y_MD)
#         x_green_list.append(x_green)
#         y_green_list.append(y_green) 

#     fig, axes = plt.subplots(1,3,sharey ="all",tight_layout = True)
#     for n in range(len(r_limited)):
#         axes[0].plot(x_dif_list[n],y_dif_list[n],'.', markersize=5.0,label=f"r_apply = {r_limited[n]}")
#         axes[1].plot(x_green_list[n],y_green_list[n],'.',markersize=5.0,label=f"r_apply = {r_limited[n]}")
#         axes[2].plot(x_MD_list[n],y_MD_list[n],'.',markersize=5.0,label=f"r_apply = {r_limited[n]}")

#     for ax in axes:
#         ax.grid(True)
#         ax.set_xlabel("Distance r (|x - x'|)")
#         ax.set_ylabel("Displacement / a")

#     axes[0].legend()
#     axes[1].legend()

#     axes[0].set_title('Difference (=|Green - MD|)')
#     axes[1].set_title("Green's Function(radial direction)")
#     axes[2].set_title('MD Simulation(radial direction)')
#     plt.show()

# if __name__ == "__main__":
#     main()




