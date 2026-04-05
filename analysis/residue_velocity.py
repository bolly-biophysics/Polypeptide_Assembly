#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Complete script for peptide side chain velocity analysis
Features:
1. Trajectory centering with PBC correction
2. Calculate ARG CZ and ASP CG side chain velocities
3. Generate per-chain and overall velocity distribution plots
4. Automatically filter out abnormally high velocity values
"""

import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import align
from tqdm import tqdm

# ======================
# 1. Parameter settings
# ======================
top_file = '../nonwater.prmtop'      # Topology file path
traj_file = 'nonwater.nc'           # Trajectory file path
dt = 40.0                           # Time step (ps)
n_chains = 10                        # Number of peptide chains
residues_per_chain = 22              # Number of residues per chain
output_prefix = 'peptide_analysis'   # Output file prefix
max_plausible_speed = 0.2           # Upper limit of reasonable speed (Å/ps)

# ======================
# 2. Trajectory processing functions
# ======================
def enhanced_trajectory_processing(universe, n_chains, residues_per_chain):
    """
    Enhanced trajectory centering processing
    Includes: PBC correction, global translation removal, per-chain rotation removal
    """
    # Create AtomGroup for each chain
    chains = [universe.select_atoms(f"resid {i*residues_per_chain+1}:{(i+1)*residues_per_chain}") 
              for i in range(n_chains)]
    
    # Store reference structure (first frame)
    ref_positions = [chain.positions.copy() for chain in chains]
    ref_com = sum(chain.center_of_mass() for chain in chains) / n_chains
    
    for ts in tqdm(universe.trajectory, desc="Processing trajectory"):
        # 1. Ensure molecules are fully within the primary box
        universe.atoms.wrap(compound='fragments')
        
        # 2. Remove global translation of the system
        current_com = sum(chain.center_of_mass() for chain in chains) / n_chains
        universe.atoms.positions -= (current_com - ref_com)
        
        # 3. Remove rotation for each chain independently
        for i, chain in enumerate(chains):
            if ts.frame > 0:  # First frame as reference
                R, _ = align.rotation_matrix(chain.positions, ref_positions[i])
                chain.positions = np.dot(chain.positions, R.T)
        
        # 4. Wrap again to ensure correct positions
        universe.atoms.wrap(compound='fragments')

# ======================
# 3. Velocity calculation functions
# ======================
def pbc_corrected_displacement(current_pos, prev_pos, box):
    """Displacement correction considering periodic boundary conditions"""
    delta = current_pos - prev_pos
    return delta - np.round(delta / box) * box

def calculate_velocities_with_pbc(universe, n_chains, residues_per_chain, dt):
    """Calculate side chain velocities for each chain (with PBC correction)"""
    # Initialize data structures
    arg_vel = [[] for _ in range(n_chains)]
    asp_vel = [[] for _ in range(n_chains)]
    
    # Create atom groups
    chains = []
    for i in range(n_chains):
        chain = universe.select_atoms(f"resid {i*residues_per_chain+1}:{(i+1)*residues_per_chain}")
        chains.append({
            'arg': chain.select_atoms("resname ARG and name CZ"),
            'asp': chain.select_atoms("resname ASP and name CG"),
            'backbone': chain.select_atoms("backbone")  # For reference
        })
    
    # Calculate velocities
    box = universe.dimensions[:3]
    for i in tqdm(range(len(universe.trajectory)-1), desc="Calculating velocities"):
        universe.trajectory[i]
        prev_pos = {
            'arg': [c['arg'].positions.copy() for c in chains],
            'asp': [c['asp'].positions.copy() for c in chains],
            'backbone': [c['backbone'].positions.copy() for c in chains]
        }
        
        universe.trajectory[i+1]
        
        for ch in range(n_chains):
            # Calculate backbone center of mass drift (after PBC correction)
            com_drift = pbc_corrected_displacement(
                chains[ch]['backbone'].center_of_mass(),
                prev_pos['backbone'][ch].mean(axis=0),
                box
            )
            
            # ARG CZ velocity (relative to backbone motion)
            if len(chains[ch]['arg']) > 0:
                dr = pbc_corrected_displacement(
                    chains[ch]['arg'].positions,
                    prev_pos['arg'][ch],
                    box
                ) - com_drift
                v = np.linalg.norm(dr, axis=1) / dt
                arg_vel[ch].extend(v)
            
            # ASP CG velocity (relative to backbone motion)
            if len(chains[ch]['asp']) > 0:
                dr = pbc_corrected_displacement(
                    chains[ch]['asp'].positions,
                    prev_pos['asp'][ch],
                    box
                ) - com_drift
                v = np.linalg.norm(dr, axis=1) / dt
                asp_vel[ch].extend(v)
    
    return arg_vel, asp_vel

# ======================
# 4. Data post-processing
# ======================
def filter_velocities(velocities, max_speed):
    """Filter out unreasonably high speed values"""
    return [np.array(v)[np.array(v) < max_speed] for v in velocities]

def save_velocity_data(arg_vel, asp_vel, prefix):
    """Save velocity data"""
    # Save as NumPy binary format
    np.savez(f"{prefix}_velocity_data.npz",
             arg_velocities=arg_vel,
             asp_velocities=asp_vel)
    
    # Save statistical summary
    with open(f"{prefix}_velocity_stats.txt", "w") as f:
        f.write("Chain\tARG_median\tASP_median\tARG_mean\tASP_mean\n")
        for ch in range(len(arg_vel)):
            arg_med = np.median(arg_vel[ch]) if len(arg_vel[ch]) > 0 else np.nan
            asp_med = np.median(asp_vel[ch]) if len(asp_vel[ch]) > 0 else np.nan
            arg_mean = np.mean(arg_vel[ch]) if len(arg_vel[ch]) > 0 else np.nan
            asp_mean = np.mean(asp_vel[ch]) if len(asp_vel[ch]) > 0 else np.nan
            f.write(f"{ch+1}\t{arg_med:.4f}\t{asp_med:.4f}\t{arg_mean:.4f}\t{asp_mean:.4f}\n")
        
        # Add global statistics
        all_arg = np.concatenate(arg_vel)
        all_asp = np.concatenate(asp_vel)
        f.write("\nOverall Statistics:\n")
        f.write(f"ARG global median: {np.median(all_arg):.4f}\n")
        f.write(f"ASP global median: {np.median(all_asp):.4f}\n")
        f.write(f"ARG global mean: {np.mean(all_arg):.4f}\n")
        f.write(f"ASP global mean: {np.mean(all_asp):.4f}\n")

# ======================
# 5. Visualization functions
# ======================
def plot_velocity_distributions(arg_vel, asp_vel, output_prefix):
    """Plot two types of velocity distribution plots"""
    # Combine all data
    all_arg = np.concatenate([v for v in arg_vel if len(v) > 0])
    all_asp = np.concatenate([v for v in asp_vel if len(v) > 0])
    
    # Set global style
    plt.style.use('seaborn')
    plt.rcParams.update({'font.size': 12, 'figure.dpi': 300})
    
    # 1. Per-chain velocity distribution plot
    plt.figure(figsize=(14, 10))
    
    # ARG per-chain distribution
    plt.subplot(2, 1, 1)
    vp = plt.violinplot(arg_vel, showmedians=True, showextrema=False)
    # Set colors
    for pc in vp['bodies']:
        pc.set_facecolor('royalblue')
        pc.set_alpha(0.7)
    vp['cmedians'].set_color('red')
    
    plt.title('ARG CZ Velocity Distribution by Chain (Å/ps)', pad=20)
    plt.ylabel('Velocity (Å/ps)')
    plt.xticks(range(1, n_chains+1), [f'Chain {i+1}' for i in range(n_chains)])
    plt.ylim(0, max_plausible_speed)
    plt.grid(True, linestyle='--', alpha=0.6)
    
    # Annotate medians
    for i, data in enumerate(arg_vel):
        if len(data) > 0:
            median = np.median(data)
            plt.text(i+1, median+0.005, f'{median:.3f}', 
                     ha='center', va='bottom', fontsize=10, color='red')
    
    # ASP per-chain distribution
    plt.subplot(2, 1, 2)
    vp = plt.violinplot(asp_vel, showmedians=True, showextrema=False)
    # Set colors
    for pc in vp['bodies']:
        pc.set_facecolor('forestgreen')
        pc.set_alpha(0.7)
    vp['cmedians'].set_color('red')
    
    plt.title('ASP CG Velocity Distribution by Chain (Å/ps)', pad=20)
    plt.ylabel('Velocity (Å/ps)')
    plt.xticks(range(1, n_chains+1), [f'Chain {i+1}' for i in range(n_chains)])
    plt.ylim(0, max_plausible_speed)
    plt.grid(True, linestyle='--', alpha=0.6)
    
    # Annotate medians
    for i, data in enumerate(asp_vel):
        if len(data) > 0:
            median = np.median(data)
            plt.text(i+1, median+0.005, f'{median:.3f}', 
                     ha='center', va='bottom', fontsize=10, color='red')
    
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_by_chain.eps", format='eps', dpi=300)
    plt.savefig(f"{output_prefix}_by_chain.png", dpi=300)
    plt.close()
    
    # 2. Overall velocity distribution plot
    plt.figure(figsize=(12, 6))
    
    # ARG overall distribution
    plt.subplot(1, 2, 1)
    vp = plt.violinplot([all_arg], showmedians=True, showextrema=False)
    vp['bodies'][0].set_facecolor('royalblue')
    vp['bodies'][0].set_alpha(0.7)
    vp['cmedians'].set_color('red')
    
    plt.title('Overall ARG CZ Velocity (Å/ps)')
    plt.ylabel('Velocity (Å/ps)')
    plt.xticks([1], ['All Chains'])
    plt.ylim(0, max_plausible_speed)
    plt.grid(True, linestyle='--', alpha=0.6)
    
    # Add statistical annotation
    stats_text = f"Median: {np.median(all_arg):.3f}\nMean: {np.mean(all_arg):.3f}"
    plt.text(1, max_plausible_speed*0.8, stats_text, 
             ha='center', va='top', bbox=dict(facecolor='white', alpha=0.8))
    
    # ASP overall distribution
    plt.subplot(1, 2, 2)
    vp = plt.violinplot([all_asp], showmedians=True, showextrema=False)
    vp['bodies'][0].set_facecolor('forestgreen')
    vp['bodies'][0].set_alpha(0.7)
    vp['cmedians'].set_color('red')
    
    plt.title('Overall ASP CG Velocity (Å/ps)')
    plt.ylabel('Velocity (Å/ps)')
    plt.xticks([1], ['All Chains'])
    plt.ylim(0, max_plausible_speed)
    plt.grid(True, linestyle='--', alpha=0.6)
    
    # Add statistical annotation
    stats_text = f"Median: {np.median(all_asp):.3f}\nMean: {np.mean(all_asp):.3f}"
    plt.text(1, max_plausible_speed*0.8, stats_text, 
             ha='center', va='top', bbox=dict(facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_overall.eps", format='eps', dpi=300)
    plt.savefig(f"{output_prefix}_overall.png", dpi=300)
    plt.close()

# ======================
# 6. Main execution flow
# ======================
if __name__ == "__main__":
    print("=== Peptide Side Chain Velocity Analysis ===")
    print(f"Parameter settings:\n- Number of chains: {n_chains}\n- Residues per chain: {residues_per_chain}")
    print(f"- Time step: {dt} ps\n- Velocity upper limit: {max_plausible_speed} Å/ps")
    
    # 1. Load trajectory
    print("\n1. Loading trajectory file...")
    u = mda.Universe(top_file, traj_file)
    print(f"System information:\n- Total atoms: {len(u.atoms)}\n- Total frames: {len(u.trajectory)}")
    print(f"- Time range: {u.trajectory[0].time/1000:.1f}-{u.trajectory[-1].time/1000:.1f} ns")
    
    # 2. Trajectory preprocessing
    print("\n2. Trajectory centering processing...")
    enhanced_trajectory_processing(u, n_chains, residues_per_chain)
    
    # 3. Calculate velocities
    print("\n3. Calculating side chain velocities...")
    arg_vel, asp_vel = calculate_velocities_with_pbc(u, n_chains, residues_per_chain, dt)
    
    # 4. Data post-processing
    print("\n4. Data processing...")
    arg_vel = filter_velocities(arg_vel, max_plausible_speed)
    asp_vel = filter_velocities(asp_vel, max_plausible_speed)
    save_velocity_data(arg_vel, asp_vel, output_prefix)
    
    # 5. Visualization
    print("\n5. Generating plots...")
    plot_velocity_distributions(arg_vel, asp_vel, output_prefix)
    
    # Results summary
    all_arg = np.concatenate([v for v in arg_vel if len(v) > 0])
    all_asp = np.concatenate([v for v in asp_vel if len(v) > 0])
    
    print("\n=== Analysis Results Summary ===")
    print(f"ARG CZ median velocity: {np.median(all_arg):.4f} ± {np.std(all_arg):.4f} Å/ps")
    print(f"ASP CG median velocity: {np.median(all_asp):.4f} ± {np.std(all_asp):.4f} Å/ps")
    
    print(f"\nAnalysis completed! Output files:")
    print(f"1. Plot files:")
    print(f"- {output_prefix}_by_chain.eps/png (Per-chain velocity distribution)")
    print(f"- {output_prefix}_overall.eps/png (Overall velocity distribution)")
    print(f"2. Data files:")
    print(f"- {output_prefix}_velocity_data.npz (Raw velocity data)")
    print(f"- {output_prefix}_velocity_stats.txt (Detailed statistics)")