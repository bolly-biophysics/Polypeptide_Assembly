import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import distances
from tqdm import tqdm

# Load trajectory
u = mda.Universe('../nonwater.prmtop', '../nonwater_10us.nc')
mg_atoms = u.select_atoms("name MG")
asp_oxygens = u.select_atoms("resname ASP and name OD1 OD2")

# State definition: 0asp ~ 6asp (number of bound ASP residues)
STATE_NAMES = {i: f"{i}asp" for i in range(7)}

# Initialization
transitions = np.zeros((7, 7), dtype=int)  # Transition matrix
state_counts = np.zeros(7, dtype=int)      # State occurrence counts
prev_states = np.zeros(len(mg_atoms), dtype=int)  # Initial state is 0asp
rare_transitions = []                      # Rare transition events (count < 50)
transition_counts = {}                     # Transition type counter

print("Analyzing trajectory...")
for ts in tqdm(u.trajectory[1:], total=len(u.trajectory)-1):
    # Calculate distance matrix between Mg²⁺ and all ASP carboxyl oxygens
    dist_matrix = distances.distance_array(
        mg_atoms.positions,
        asp_oxygens.positions,
        box=ts.dimensions
    )
    
    for mg_idx in range(len(mg_atoms)):
        prev_state = prev_states[mg_idx]
        
        # Step 1: Find all ASP carboxyl oxygens with distance < 2.8Å
        bound_oxy_indices = np.where(dist_matrix[mg_idx] < 2.8)[0]
        
        # Step 2: Count involved ASP residues (deduplicate)
        bound_asp_residues = set()
        for oxy_idx in bound_oxy_indices:
            asp_residue = asp_oxygens[oxy_idx].residue.resid
            bound_asp_residues.add(asp_residue)
        
        # Current state = number of bound ASP residues (max 6)
        curr_state = min(len(bound_asp_residues), 6)
        
        # Record transition
        transitions[prev_state, curr_state] += 1
        trans_key = (prev_state, curr_state)
        transition_counts[trans_key] = transition_counts.get(trans_key, 0) + 1
        
        # Record rare transition events
        if prev_state != curr_state and transition_counts[trans_key] < 50:
            rare_transitions.append({
                "frame": ts.frame,
                "mg_index": mg_idx,
                "from_state": prev_state,
                "to_state": curr_state,
                "asp_residues": list(bound_asp_residues),
                "transition_count": transition_counts[trans_key]
            })
        
        prev_states[mg_idx] = curr_state

# Count state occurrences
state_counts = np.sum(transitions, axis=0)

# Output results (same as before)
print("\n=== Absolute Transition Count Matrix ===")
abs_transitions = pd.DataFrame(
    transitions,
    index=[f"From {s}" for s in STATE_NAMES.values()],
    columns=[f"To {s}" for s in STATE_NAMES.values()]
)
pd.set_option('display.max_columns', None)
print(abs_transitions)

print("\n=== State Statistics ===")
print(f"Total frames × Mg²⁺ count: {len(u.trajectory) * len(mg_atoms)}")
for s in STATE_NAMES:
    print(f"{STATE_NAMES[s]}: {state_counts[s]}")

print("\n=== Detailed Transition Event Statistics ===")
for from_state in range(7):
    for to_state in range(7):
        if from_state != to_state and transitions[from_state, to_state] > 0:
            print(f"{STATE_NAMES[from_state]} → {STATE_NAMES[to_state]}: {transitions[from_state, to_state]}")

# Save rare events
if rare_transitions:
    with open("rare_transitions.log", "w") as f:
        f.write("Frame\tMg²⁺ Index\tFrom State\tTo State\tASP Residue List\tCumulative Count\n")
        for event in rare_transitions:
            f.write(f"{event['frame']}\t{event['mg_index']}\t"
                   f"{STATE_NAMES[event['from_state']]}\t"
                   f"{STATE_NAMES[event['to_state']]}\t"
                   f"{event['asp_residues']}\t"
                   f"{event['transition_count']}\n")
    print(f"\nRare transition events saved to rare_transitions.log (total {len(rare_transitions)} entries)")
else:
    print("\nNo events with transition count less than 50 found.")