from ase.build import bulk
from ase.build import make_supercell

# === 1. Create 4H-SiC unit cell ===
sic = bulk('SiC', crystalstructure='wurtzite', a=3.073, c=10.053, orthorhombic=False)

# === 2. Define supercell repeat (nx, ny, nz) ===
# For example: 2x2x2 supercell = 8 times the atoms and volume
repeat = (2, 2, 2)  # adjust this to get desired number of atoms or volume

# Method 1: Easy way using ASE's built-in repeat method
supercell = sic.repeat(repeat)

# === 3. Check new size ===
print(f"Number of atoms: {len(supercell)}")
print(f"New volume: {supercell.get_volume():.2f} Å^3")

# === 4. Save to file (optional) ===
supercell.write("4H_SiC_supercell.xyz")
