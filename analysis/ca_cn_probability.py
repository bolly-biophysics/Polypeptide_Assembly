import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import distances
from tqdm import tqdm
import matplotlib.pyplot as plt

# Load trajectory
u = mda.Universe('pep-pep.prmtop', ['pbc_water/fixed.xtc', 'pbc_water/fixed1.xtc'])

# Define atom selections
ca_atoms = u.select_atoms("resname CA")  # Calcium ions (compatible with AMBER/GROMACS naming)
water_oxygens = u.select_atoms("resname WAT and name O")  # Water oxygen atoms (TIP3P/TIP4P)

print(f"System contains: {len(ca_atoms)} Ca²⁺, {len(water_oxygens)} water molecules")

# Parameter settings
CUTOFF = 3.2  # Hydration shell distance threshold (Å)
STRIDE = 10   # Analyze every 10 frames
output_file = "cn_ca.dat"  # Output filename

# Initialize statistics array
hydration_counts = []

# Analyze frame by frame (every STRIDE frames)
for ts in tqdm(u.trajectory[::STRIDE], total=len(u.trajectory)//STRIDE, desc="Analysis progress"):
    # Compute distance matrix between Ca²⁺ and water oxygen atoms
    dist_matrix = distances.distance_array(
        ca_atoms.positions,
        water_oxygens.positions,
        box=ts.dimensions  # Consider periodic boundary conditions
    )
    
    # Count water molecules around each Ca²⁺
    for ca_idx in range(len(ca_atoms)):
        n_water = np.sum(dist_matrix[ca_idx, :] < CUTOFF)
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