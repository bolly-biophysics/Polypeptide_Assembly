import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import distances
from tqdm import tqdm

# Load trajectory
u = mda.Universe('nonwater.prmtop', 'nonwater_10us.nc')

# Define atom selections
arg_nitrogens = u.select_atoms("resname ARG and name NH1 NH2")  # Positively charged nitrogen atoms of ARG
asp_oxygens = u.select_atoms("resname ASP and name OD1 OD2")   # Negatively charged oxygen atoms of ASP

print(f"System contains: {len(arg_nitrogens)} ARG nitrogen atoms, {len(asp_oxygens)} ASP oxygen atoms")

# Parameter settings
CUTOFF = 3.5  # Salt bridge contact distance threshold (Å), recommended 3.5-4.0 Å
output_file = "arg_asp_contacts_10us.dat"  # Output filename
STRIDE = 10  # Analyze every 10 frames

# Open output file
with open(output_file, 'w') as f:
    f.write("# Frame\tTotal_Contacts\n")  # Write header
    
    # Analyze frame by frame (every STRIDE frames)
    for ts in tqdm(u.trajectory[::STRIDE], total=len(u.trajectory)//STRIDE, desc="Analysis progress"):
        # Compute distance matrix between all ARG(NH1/NH2) and ASP(OD1/OD2)
        dist_matrix = distances.distance_array(
            arg_nitrogens.positions,
            asp_oxygens.positions,
            box=ts.dimensions  # Consider periodic boundary conditions
        )
        
        # Count number of contacts (distance < CUTOFF)
        total_contacts = np.sum(dist_matrix < CUTOFF)
        
        # Write current frame result
        f.write(f"{ts.frame}\t{total_contacts}\n")

print(f"Analysis completed! Results saved to {output_file}")