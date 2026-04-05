#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to calculate time-dependent Ca-Cl contact count
Function: Count Ca-Cl contacts (distance < 3.5Å) for each frame in the trajectory
Output: Time series data of frame number and corresponding contact count
"""
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import distances
from tqdm import tqdm  # Progress bar tool

# ====================== Parameter settings ======================
prmtop = 'nonwater.prmtop'      # Topology file
traj = 'nonwater.nc'            # Trajectory file
output_file = 'ca_cl_contacts.dat'  # Output file

ca_selection = "resname CA"  # Calcium ion selection statement
cl_selection = "name Cl- or element Cl"  # Chloride ion selection statement
CUTOFF = 3.5                     # Contact distance threshold (Å)
STRIDE = 10                      # Frame skipping interval

# ====================== Main program ======================
def main():
    # 1. Load trajectory
    u = mda.Universe(prmtop, traj)
    ca_atoms = u.select_atoms(ca_selection)
    cl_atoms = u.select_atoms(cl_selection)
    
    print(f"System contains: {len(ca_atoms)} Ca²⁺ ions, {len(cl_atoms)} Cl⁻ ions")
    print(f"Analysis parameters: cutoff distance={CUTOFF}Å, frame skip interval={STRIDE} frames")

    # 2. Open output file
    with open(output_file, 'w') as f:
        f.write("# Frame\tTotal_Contacts\n")  # Write header
        
        # 3. Analyze frame by frame (with progress bar)
        for ts in tqdm(u.trajectory[::STRIDE], 
                      total=len(u.trajectory)//STRIDE,
                      desc="Analysis progress"):
            # Calculate distance matrix (automatically handles periodicity)
            dist_matrix = distances.distance_array(
                ca_atoms.positions,
                cl_atoms.positions,
                box=ts.dimensions  # Key: consider periodic boundary conditions
            )
            
            # Count number of contacts (distance < CUTOFF)
            total_contacts = np.sum(dist_matrix < CUTOFF)
            
            # Write result
            f.write(f"{ts.frame}\t{total_contacts}\n")

    print(f"\nAnalysis completed! Results saved to {output_file}")

if __name__ == "__main__":
    main()