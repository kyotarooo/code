from lammps import lammps #type: ignore
lmp = lammps()


with open(f"/Users/kyou/Library/CloudStorage/Box-Box/lammps/pointdefects_Zr/infilename.txt", 'r') as f_num:
    f = f_num.readline()
    filename = f.split()
    
    
for i in range(len(filename)):
    lmp.file(f"/Users/kyou/Library/CloudStorage/Box-Box/lammps/pointdefects_Zr/in_file/{filename[i]}")
    


