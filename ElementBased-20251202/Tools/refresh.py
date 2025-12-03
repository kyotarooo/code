# ===================================================================================
# refresh.py: Delete all simulation histroy
# -----------------------------------------------------------------------------------
# This script is useful for deleting all simulation history to make a new simulation.
# -----------------------------------------------------------------------------------
# Developer: A. Takahashi (Tokyo University of Science)
# ===================================================================================
import os
import glob


def delete_files(files):
    for file in files:
        if os.path.isfile(file):
            os.remove(file)
            print(file + " was removed.")


files = (
    "restart.dat",
    "density.out",
    "mechanical_behavior.out",
    "time.out",
    "elem.pvd",
    "elem.plt",
    "direct_interaction.log"
)

delete_files(files)
delete_files(glob.glob("*.vtu"))
delete_files(glob.glob("*.vtk"))
delete_files(glob.glob("*.txt"))
