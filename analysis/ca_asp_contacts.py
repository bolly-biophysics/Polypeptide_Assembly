import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import distances  # Correct import method
from tqdm import tqdm

# Load trajectory
u = mda.Universe('nonwater.prmtop', 'nonwater.nc')
ca_atoms = u.select_atoms("resname CA")
asp_oxygens = u.select_atoms("resname ASP and name OD1 OD2")
print(f"System contains: {len(ca_atoms)} Ca²⁺ ions, {len(asp_oxygens)} ASP oxygen atoms")

# Parameter settings
CUTOFF = 3.2  # Contact distance threshold (Å)
output_file = "ca_asp_contacts.dat"  # Output filename
STRIDE = 10   # Analyze every 10 frames

# Open output file
with open(output_file, 'w') as f:
    f.write("# Frame\tTotal_Contacts\n")  # Write header
    
    # Analyze frame by frame (every STRIDE frames)
    for ts in tqdm(u.trajectory[::STRIDE], total=len(u.trajectory)//STRIDE, desc="Analysis progress"):
        # Compute all Ca-oxygen distances
        dist_matrix = distances.distance_array(
            ca_atoms.positions,
            asp_oxygens.positions,
            box=ts.dimensions
        )
        
        # Count number of contacts (distance < CUTOFF)
        total_contacts = np.sum(dist_matrix < CUTOFF)
        
        # Write current frame result
        f.write(f"{ts.frame}\t{total_contacts}\n")

print(f"Analysis completed! Results saved to {output_file}")