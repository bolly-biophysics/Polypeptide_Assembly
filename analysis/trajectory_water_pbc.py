import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import distances
from MDAnalysis.lib.distances import transform_RtoS, transform_StoR

# ====================== Parameter settings ======================
prmtop = "../pep-pep.prmtop"
nc_traj = "../md.nc" 
output_traj = "fixed.xtc"

mg_selection = "resname MG"            # Mg²⁺ selection statement
water_selection = "resname WAT"        # Select complete water molecules (including O and H)
r_cutoff = 2.8                         # Strict hydration shell distance threshold (Å), no tolerance

# ====================== Main program ======================
def main():
    u = mda.Universe(prmtop, nc_traj)
    print(f"Original box type: {u.trajectory[0].dimensions}")

    mg_ions = u.select_atoms(mg_selection)
    waters = u.select_atoms(water_selection)
    water_oxygens = waters.select_atoms("name O")
    water_dict = {res.resindex: res.atoms.indices for res in waters.residues}
    
    n_mg = len(mg_ions)
    n_waters = len(water_oxygens)
    mg_positions = np.zeros((n_mg, 3))
    water_o_positions = water_oxygens.positions.copy()
    
    with mda.Writer(output_traj, u.atoms.n_atoms) as w:
        for ts in u.trajectory:
            print(f"Processing frame {ts.frame}...")
            box = ts.dimensions
            mg_positions[:] = mg_ions.positions
            water_o_positions[:] = water_oxygens.positions
            
            # Record corrected water molecules (avoid duplicate processing)
            corrected = np.zeros(n_waters, dtype=bool)
            
            for i in range(n_mg):
                dist = distances.distance_array(
                    mg_positions[i].reshape(1, 3),
                    water_o_positions,
                    box=box
                )

                in_hydration = (dist[0] <= r_cutoff) & ~corrected

                # Apply PBC correction to water molecules entering the hydration shell
                for idx in np.where(in_hydration)[0]:
                    res_idx = water_oxygens.resindices[idx]
                    water_atoms = u.atoms[water_dict[res_idx]]
                    r = water_atoms.positions - mg_positions[i]
                    s = transform_RtoS(r, box)
                    s_wrapped = s - np.round(s)
                    water_atoms.positions = mg_positions[i] + transform_StoR(s_wrapped, box)
                    corrected[idx] = True  # Mark as corrected

            w.write(u.atoms)
    
    print(f"Trajectory saved to {output_traj}")

if __name__ == "__main__":
    main()