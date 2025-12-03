import os
import glob

if os.path.isfile("restart.dat"):
    os.remove("restart.dat")
    print("restart.dat is removed.")

if os.path.isfile("volume.vtk"):
    os.remove("volume.vtk")
    print("volume.vtk is removed.")

if os.path.isfile("time.out"):
    os.remove("time.out")
    print("time.out is removed.")

if os.path.isfile("elem.pvd"):
    os.remove("elem.pvd")
    print("elem.pvd is removed.")

if os.path.isfile("elem.plt"):
    os.remove("elem.plt")
    print("elem.plt is removed.")
    
vtu_files = glob.glob("elem_*.vtu")
for vtu_file in vtu_files:
    if os.path.isfile(vtu_file):
        os.remove(vtu_file)
        print(vtu_file + " is removed.")

vtk_files = glob.glob("elem_*.vtk")
for vtk_file in vtk_files:
    if os.path.isfile(vtk_file):
        os.remove(vtk_file)
        print(vtk_file + " is removed.")
