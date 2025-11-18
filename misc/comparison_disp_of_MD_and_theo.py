import numpy as np#type: ignore
import sys
import matplotlib.pyplot as plt#type: ignore



def main():
    with open(f"{sys.argv[1]}/data.inp","r") as f:
        lines = f.readlines()
        xc = np.array([float(lines[1].split()[i]) for i in range(3)])#xcは欠陥の座標
        n = int(lines[2].split()[0])

    with open(f"{sys.argv[1]}/disp_of_green_sort.out","r") as f:
        lines = f.readlines()
        idx = 1
        ui_list_green =[]
        xi_list_green = []
        abs_ri_list = []
        abs_ui_list_green = []
        r_limited = float(lines[0].split()[0])
        ur_list_green = []

        for _ in range(n):
            xi = np.array([float(lines[idx].split()[i]) for i in range(3)])
            ui = np.array([float(lines[idx].split()[i+3]) for i in range(3)])
            idx += 1
            ri = xi - xc
            abs_ri = np.linalg.norm(xi - xc)
            nr = ri / abs_ri #半径方向単位ベクトル
            ur = (-1)*np.dot(nr, ui)

            if abs_ri>=r_limited:
                break      

            ur_list_green.append(ur)
            abs_ri_list.append(abs_ri)
            ui_list_green.append(ui)
            xi_list_green.append(xi)   
            abs_ui_list_green.append(np.linalg.norm(ui))
        
        n = len(abs_ri_list)

    with open(f"{sys.argv[1]}/disp_of_MD_sort.out","r") as f:
        lines = f.readlines()
        idx = 1
        ui_list_MD = []
        xi_list_MD = []
        abs_ui_list_MD = []
        ur_list_MD = []

        for _ in range(n):
            xi = np.array([float(lines[idx].split()[i]) for i in range(3)])
            ui = np.array([float(lines[idx].split()[i+3]) for i in range(3)])
            idx += 1
            ri = xi - xc
            abs_ri = np.linalg.norm(xi - xc)
            nr = ri / abs_ri #半径方向単位ベクトル
            ur = (-1)*np.dot(nr, ui) #中心方向を正とする
            ur_list_MD.append(ur)
            ui_list_MD.append(ui)
            xi_list_MD.append(xi)
            abs_ui_list_MD.append(np.linalg.norm(ui))

    with open(f"{sys.argv[1]}/disp_of_green_r_sort.out","r") as f:
        lines = f.readlines()
        idx = 1
        ui_r_list_green = []
        xi_r_list_green = []
        abs_ui_r_list_green = []
        ur_r_list_green = []

        for _ in range(n):
            xi = np.array([float(lines[idx].split()[i]) for i in range(3)])
            ui = np.array([float(lines[idx].split()[i+3]) for i in range(3)])
            idx += 1
            ri = xi - xc
            abs_ri = np.linalg.norm(xi - xc)
            nr = ri / abs_ri #半径方向単位ベクトル
            ur = (-1)*np.dot(nr, ui) #中心方向を正とする
            ur_r_list_green.append(ur)
            ui_r_list_green.append(ui)
            xi_r_list_green.append(xi)
            abs_ui_r_list_green.append(np.linalg.norm(ui))
        
        

    dif_u_list = []
    dif_ur_list = []
    for m in range(n):
        dif_u = ui_list_green[m]-ui_list_MD[m]
        dif_ur = ur_list_green[m] - ur_list_MD[m]
        abs_dif_u = np.linalg.norm(dif_u)
        abs_dif_ur = np.linalg.norm(dif_ur)
        dif_u_list.append(abs_dif_u)
        dif_ur_list.append(abs_dif_ur)

    with open(f"{sys.argv[1]}/dif_MD_green_{r_limited}.out","w") as f_dif,\
         open(f"{sys.argv[1]}/ur_MD_{r_limited}.out","w") as f_MD,\
         open(f"{sys.argv[1]}/ur_green_{r_limited}.out","w") as f_green,\
         open(f"{sys.argv[1]}/ur_green_r_{r_limited}.out","w") as f_green_r:
            for i in range(n):
                f_dif.write(f"{abs_ri_list[i]:.12e} {dif_ur_list[i]:.12e}\n")
                f_MD.write(f"{abs_ri_list[i]:.12e} {ur_list_MD[i]:.12e}\n")
                f_green.write(f"{abs_ri_list[i]:.12e} {ur_list_green[i]:.12e}\n")
                f_green_r.write(f"{abs_ri_list[i]:.12e} {ur_r_list_green[i]:.12e}\n")


    scall_factor = 1.0e+10
    dif_ur_array = np.array(dif_ur_list)
    ur_list_green_array = np.array(ur_list_green)
    ur_list_MD_array = np.array(ur_list_MD)
    ur_r_list_green_array = np.array(ur_r_list_green)





    plt.plot(abs_ri_list,ur_list_green_array*scall_factor, label='green func', color='#008000')
    plt.scatter(abs_ri_list,ur_list_MD_array*scall_factor, label='MD',marker=".", color='#ff2b2b')
    plt.plot(abs_ri_list,ur_r_list_green_array*scall_factor, label='residual stress', color='#002b55')
    plt.xlabel("Distance r(=|x-x'|")
    plt.ylabel('dislacement[Å]')          
    plt.legend()
    plt.grid(True)
    plt.tight_layout()#グラフのレイアウトを自動的に調整、ラベルやタイトルなどがはみ出さないようにする
    plt.show()

    #16進カラーコードはhttps://www.tagindex.com/color/color_gradation.htmlを使用






    # fig, axes = plt.subplots(1,4,sharey ="all",tight_layout = True)
    # axes[0].plot(abs_ri_list, dif_ur_array * scall_factor, '.', markersize=4.0)
    # axes[1].plot(abs_ri_list, ur_list_green_array * scall_factor, '.', markersize=4.0)
    # axes[2].plot(abs_ri_list, ur_list_MD_array * scall_factor, '.', markersize=4.0)
    # axes[3].plot(abs_ri_list, ur_r_list_green_array * scall_factor, '.', markersize=4.0)
    # axes[0].grid(True)
    # axes[1].grid(True)
    # axes[2].grid(True)
    # axes[3].grid(True)
    # axes[0].set_title('Difference (=|Green - MD|)')
    # axes[1].set_title("Green's Function(radial direction)")
    # axes[2].set_title('MD Simulation(radial direction)')
    # axes[0].set_xlabel("Distance r(=|x-x'|)")
    # axes[1].set_xlabel("Distance r(=|x-x'|)")
    # axes[2].set_xlabel("Distance r(=|x-x'|)")
    # fig.suptitle(f'r = {r_limited}')
    # plt.show()


if __name__ == "__main__":#『このファイルがコマンドラインからスクリプトとして直接実行された場合のみ以降の処理を実行する
    main()









