import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import distances
from tqdm import tqdm
import warnings

# Set matplotlib to use English labels
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.unicode_minus'] = False

# Load trajectory
u = mda.Universe("../nonwater.prmtop", "nonwater_10us.nc")

# Parameters
n_chains = 10                  # Number of peptide chains
residues_per_chain = 22        # Residues per chain
contact_cutoff = 4.5           # Contact distance threshold (Å)
min_contacts = 5               # Minimum contact pairs
box = u.dimensions if hasattr(u, 'dimensions') else None

def get_chain_atom_groups():
    """Get heavy atom groups for each chain"""
    chains = []
    protein = u.select_atoms("protein")
    
    # Method 1: Use segid if available
    if hasattr(protein, 'segids'):
        segids = np.unique(protein.segids)
        if len(segids) == n_chains:
            return [u.select_atoms(f"segid {segid} and not name H*") for segid in segids]
    
    # Method 2: Automatic segmentation by resid
    resids = np.unique(protein.resids)
    chain_resids = np.array_split(np.sort(resids), n_chains)
    return [u.select_atoms(f"resid {' '.join(map(str, res))} and not name H*") 
            for res in chain_resids]

chains = get_chain_atom_groups()

def analyze_frame():
    # Calculate contact matrix
    contact_matrix = np.zeros((n_chains, n_chains))
    
    for i in range(n_chains):
        for j in range(i+1, n_chains):
            dist = distances.distance_array(chains[i].positions, 
                                         chains[j].positions,
                                         box=box)
            contact_matrix[i,j] = contact_matrix[j,i] = np.sum(dist < contact_cutoff)
    
    # Cluster identification
    parent = list(range(n_chains))
    def find(u):
        while parent[u] != u:
            parent[u] = parent[parent[u]]
            u = parent[u]
        return u
    
    pairs = np.argwhere(contact_matrix >= min_contacts)
    for i,j in pairs:
        if i < j:
            parent[find(j)] = find(i)
    
    return len(set(find(i) for i in range(n_chains))), contact_matrix

# Main analysis
cluster_counts = []
contact_history = []

# Create output file
with open("cluster_evolution_10us.dat", "w") as f:
    f.write("# Frame  ClusterCount\n")  # Header
    
    for ts in tqdm(u.trajectory, desc="Analyzing"):
        frame = ts.frame
        n_clusters, contacts = analyze_frame()
        cluster_counts.append(n_clusters)
        contact_history.append(contacts)
        
        # Write to file each frame
        f.write(f"{frame:6d}  {n_clusters:4d}\n")

# Analysis summary
print("\n=== Results Summary ===")
print(f"Frames analyzed: {len(u.trajectory)}")
print(f"Initial clusters: {cluster_counts[0]}")
print(f"Final clusters: {cluster_counts[-1]}")
print(f"Results saved to 'cluster_evolution.dat'")

# Plotting
plt.figure(figsize=(12,6))
plt.plot(cluster_counts)
plt.xlabel("Frame")
plt.ylabel("Number of Clusters")
plt.title("Cluster Evolution")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("cluster_evolution.png", dpi=300)
plt.show()
