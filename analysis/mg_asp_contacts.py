import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import distances  # Correct import method
from tqdm import tqdm

# Load trajectory
u = mda.Universe('nonwater.prmtop', 'nonwater_10us.nc')
mg_atoms = u.select_atoms("name MG")
asp_oxygens = u.select_atoms("resname ASP and name OD1 OD2")
print(f"System contains: {len(mg_atoms)} Mg²⁺ ions, {len(asp_oxygens)} ASP oxygen atoms")

# Parameter settings
CUTOFF = 2.8  # Contact distance threshold (Å)
output_file = "mg_asp_contacts_10us.dat"  # Output filename
STRIDE = 10   # Analyze every 10 frames

# Open output file
with open(output_file, 'w') as f:
    f.write("# Frame\tTotal_Contacts\n")  # Write header
    
    # Analyze frame by frame (every STRIDE frames)
    for ts in tqdm(u.trajectory[::STRIDE], total=len(u.trajectory)//STRIDE, desc="Analysis progress"):
        # Compute all Mg-oxygen distances
        dist_matrix = distances.distance_array(
            mg_atoms.positions,
            asp_oxygens.positions,
            box=ts.dimensions
        )
        
        # Count number of contacts (distance < CUTOFF)
        total_contacts = np.sum(dist_matrix < CUTOFF)
        
        # Write current frame result
        f.write(f"{ts.frame}\t{total_contacts}\n")

print(f"Analysis completed! Results saved to {output_file}")