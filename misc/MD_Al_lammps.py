import numpy as np #type: ignore
import subprocess

#***************材料定数の定義*******************
E = 70.6e+09 #Pa 裳華房　材料力学　P.15
nu = 0.3
g = E / (2 * (1 + nu)) # 剛性率
lattice_const = 4.04934

#************lammpsのスクリプト******************
cub_side_length = np.array(list(range(3, 31)))
lammpsfilename = "lammps_Al.in"
force_dipole_factor = 1.0e-25 #bars times Å^3　をN times mへ
vector_list = ["pxx", "pyy", "pzz", "pxy","pxz", "pyz"]

with open(f"/Users/kyou/Library/CloudStorage/Box-Box/lammps_related/{lammpsfilename}","w") as f:
    f.write('log /Users/kyou/Library/CloudStorage/Box-Box/lammps_related/log.Al\n')
    f.write('clear\n')
    f.write('units metal\n')
    f.write('dimension 3\n')
    f.write('boundary p p p\n')#←←←←  境界条件
    f.write('atom_style atomic\n')
    f.write('lattice fcc 4.04934\n')
    f.write('region Al_region block -20 20 -20 20 -20 20\n')
    f.write('create_box 1 Al_region\n')
    f.write('create_atoms 1 region Al_region\n')
    f.write('mass 1 26.98\n')
    f.write('pair_style eam/alloy\n')
    f.write('pair_coeff * * /Users/kyou/Library/CloudStorage/Box-Box/interatomic_potentials/Al-LEA.eam.alloy Al\n')
    #******原子空孔領域******    
    f.write('region vacancyspot sphere 0 0.1 0 0.5\n')
    f.write('group vacancyatom region vacancyspot\n')
    #******force dipole の適用領域　立方体で定義********
    for n in cub_side_length:
        m = n/2   
        f.write(f'region residual_applied_area_{n} block -{m} {m} -{m} {m} -{m} {m}\n')
        f.write(f'group residual_applied_atom_{n} region residual_applied_area_{n}\n')
#******原子空孔の情報を出力*******←dumpの番号を分ければ、undumpする必要ない？
    f.write('dump 1 vacancyatom custom 1 /Users/kyou/Library/CloudStorage/Box-Box/lammps_related/Al_deleted_atoms.dump id type x y z\n')
    f.write('dump_modify 1 format line "%d %d %.10f %.10f %.10f"\n')
    f.write('run 0\n')
    f.write('undump 1\n')
    #******原子を削除******
    f.write('delete_atoms group vacancyatom\n')
    f.write('dump 1 all custom 1 /Users/kyou/Library/CloudStorage/Box-Box/lammps_related/MD_Al_perfect.dump id type x y z\n')
    f.write('run 0\n')
    f.write('undump 1\n')
    #******緩和******
    f.write('minimize 1e-8 1e-10 1000 10000\n')
    f.write('run 1000\n')
    #******応力の計算******
    f.write('compute peratom all stress/atom NULL\n')
    for k in range(3,31):
        f.write(f'compute sumstress{k} residual_applied_atom_{k} reduce sum c_peratom[1] c_peratom[2] c_peratom[3] c_peratom[4] c_peratom[5] c_peratom[6]\n')
        for l in range(1,7):
            f.write(f'variable {vector_list[l-1]}{k} equal c_sumstress{k}[{l}]*{force_dipole_factor}\n')
        f.write(f'fix printstress{k} all print 1 "${{pxx{k}}} ${{pyy{k}}} ${{pzz{k}}} ${{pxy{k}}} ${{pxz{k}}} ${{pyz{k}}}" file /Users/kyou/Library/CloudStorage/Box-Box/lammps_related/sumstress_Al/stress{k}.txt screen no\n')        
    f.write('run 1\n')

    subprocess.run(["/Users/kyou/venv/bin/lmp", "-in", f"/Users/kyou/Library/CloudStorage/Box-Box/lammps_related/{lammpsfilename}"])
    
    




