from lammps import lammps

# シンプルに呼び出し
lmp = lammps()

print("LAMMPS loaded successfully!")
    


# from lammps import lammps #type: ignore
# lmp = lammps(libfile="/usr/local/Cellar/lammps/20250722-update1/lib/liblammps_serial.0.dylib")


# with open(f"/Users/kyou/Library/CloudStorage/Box-Box/lammps/pointdefects_Zr/infilename.txt", 'r') as f_num:
#     f = f_num.readline()
#     filename = f.split()
    
    
# for i in range(len(filename)):
#     lmp.file(f"/Users/kyou/Library/CloudStorage/Box-Box/lammps/pointdefects_Zr/in_file/{filename[i]}")
 