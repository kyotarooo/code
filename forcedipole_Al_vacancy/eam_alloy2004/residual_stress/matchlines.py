import numpy as np #type: ignore

def main():
     with open(f"/Users/kyou/Library/CloudStorage/Box-Box/lammps/pointdefects_Al/force_dipole_J/perfect_J_sumstress.txt", "r") as f_pe_J,\
         open(f"/Users/kyou/Library/CloudStorage/Box-Box/lammps/pointdefects_Al/force_dipole_eV/perfect_eV_sumstress.txt", "r") as f_pe_eV,\
         open(f"/Users/kyou/Library/CloudStorage/Box-Box/lammps/pointdefects_Al/force_dipole_J/defect_J_sumstress.txt", "r") as f_de_J,\
         open(f"/Users/kyou/Library/CloudStorage/Box-Box/lammps/pointdefects_Al/force_dipole_eV/defect_eV_sumstress.txt", "r") as f_de_eV:
        lines_perfect_J = f_pe_J.readlines()
        lines_perfect_eV = f_pe_eV.readlines()
        lines_defect_J = f_de_J.readlines()
        lines_defect_eV = f_de_eV.readlines()
        
        n_J = len(lines_perfect_J)
        n_eV = len(lines_defect_eV)

        minnum = min(n_eV, n_J)

        print (len(lines_perfect_J))
        print (len(lines_defect_J))
        print (len(lines_perfect_eV))
        print (len(lines_defect_eV))
        print (f"minnumber = {minnum}")

     with open(f"/Users/kyou/Library/CloudStorage/Box-Box/lammps/pointdefects_Al/force_dipole_J/perfect_J_sumstress_sub.txt", "w") as f_pe_J,\
         open(f"/Users/kyou/Library/CloudStorage/Box-Box/lammps/pointdefects_Al/force_dipole_eV/perfect_eV_sumstress_sub.txt", "w") as f_pe_eV,\
         open(f"/Users/kyou/Library/CloudStorage/Box-Box/lammps/pointdefects_Al/force_dipole_J/defect_J_sumstress_sub.txt", "w") as f_de_J,\
         open(f"/Users/kyou/Library/CloudStorage/Box-Box/lammps/pointdefects_Al/force_dipole_eV/defect_eV_sumstress_sub.txt", "w") as f_de_eV:

         f_pe_J.writelines(lines_perfect_J[-minnum:])
         f_pe_eV.writelines(lines_perfect_eV[-minnum:])
         f_de_J.writelines(lines_defect_J[-minnum:])
         f_de_eV.writelines(lines_defect_eV[-minnum:])
        
if __name__ == "__main__":
    main()
