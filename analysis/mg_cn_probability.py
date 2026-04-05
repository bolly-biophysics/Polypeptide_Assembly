import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import distances
from tqdm import tqdm
import matplotlib.pyplot as plt

# Load trajectory
u = mda.Universe('pep-pep.prmtop', ['pbc_water/fixed.xtc', 'pbc_water/fixed1.xtc', 'xtc_tmp/fixed2.xtc', 'xtc_tmp/fixed3.xtc', 'xtc_tmp/fixed4.xtc', 'xtc_tmp/fixed5.xtc', 'xtc_tmp/fixed6.xtc'])

# Define atom selections
mg_atoms = u.select_atoms("resname MG")  # Magnesium ions (compatible with AMBER/GROMACS naming)
water_oxygens = u.select_atoms("resname WAT and name O")  # Water oxygen atoms (TIP3P/TIP4P)

print(f"System contains: {len(mg_atoms)} Mg²⁺, {len(water_oxygens)} water molecules")

# Parameter settings
CUTOFF = 2.8  # Magnesium ion hydration shell distance threshold (Å)
STRIDE = 10   # Analyze every 10 frames
output_file = "cn_mg_10us.dat"  # Output filename

# Initialize statistics array
hydration_counts = []

# Analyze frame by frame (every STRIDE frames)
for ts in tqdm(u.trajectory[::STRIDE], total=len(u.trajectory)//STRIDE, desc="Analysis progress"):
    # Compute distance matrix between Mg²⁺ and water oxygen atoms
    dist_matrix = distances.distance_array(
        mg_atoms.positions,
        water_oxygens.positions,
        box=ts.dimensions  # Consider periodic boundary conditions
    )
    
    # Count water molecules around each Mg²⁺
    for mg_idx in range(len(mg_atoms)):
        n_water = np.sum(dist_matrix[mg_idx, :] < CUTOFF)
        hydration_counts.append(n_water)

# Convert to numpy array for statistics
hydration_counts = np.array(hydration_counts)

# Calculate distribution probability
max_count = np.max(hydration_counts)
bins = np.arange(-0.5, max_count + 1.5, 1)  # Bin edges (centered on integers)
hist, bin_edges = np.histogram(hydration_counts, bins=bins, density=True)

# Output results
with open(output_file, 'w') as f:
    f.write("# Hydration_Number\tProbability\n")
    for num, prob in zip(bin_edges[:-1]+0.5, hist):
        f.write(f"{int(num)}\t{prob:.4f}\n")

print(f"Analysis completed! Results saved to {output_file}")