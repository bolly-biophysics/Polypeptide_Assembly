import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import distances
from tqdm import tqdm

# Load trajectory
u = mda.Universe('../nonwater.prmtop', '../nonwater.nc')
ca_atoms = u.select_atoms("resname CA")  # Modified to select calcium ions
asp_oxygens = u.select_atoms("resname ASP and name OD1 OD2")

# State definition (extended to S8)
STATE_NAMES = {
    0: "S0 (8H₂O)",
    1: "S1 (7H₂O+1COO⁻)",
    2: "S2 (6H₂O+2COO⁻)",
    3: "S3 (5H₂O+3COO⁻)",
    4: "S4 (4H₂O+4COO⁻)",
    5: "S5 (3H₂O+5COO⁻)",
    6: "S6 (2H₂O+6COO⁻)",
    7: "S7 (1H₂O+7COO⁻)",
    8: "S8 (0H₂O+8COO⁻)"
}

# Initialization
n_states = 9  # Number of states increased from 7 to 9
transitions = np.zeros((n_states, n_states), dtype=int)  # 9x9 transition count matrix
state_counts = np.zeros(n_states, dtype=int)             # Total occurrences of each state
rare_transitions = []                                    # Record rare transition events (count < 10)
prev_states = np.zeros(len(ca_atoms), dtype=int)         # Track previous state for each Ca²⁺ (initial S0)
transition_counts = {}                                   # Track cumulative count for each transition type

print("Analyzing trajectory...")
for ts in tqdm(u.trajectory[1:], total=len(u.trajectory)-1):  # Start from second frame
    dist_matrix = distances.distance_array(
        ca_atoms.positions,
        asp_oxygens.positions,
        box=ts.dimensions
    )
    
    for ca_idx in range(len(ca_atoms)):
        prev_state = prev_states[ca_idx]  # Get previous state
        n_bound = np.sum(dist_matrix[ca_idx] < 3.2)  # Modified distance threshold to 3.2Å
        curr_state = min(n_bound, 8)                # Current state (max is 8)
        
        # Record state transition (prev_state -> curr_state)
        transitions[prev_state, curr_state] += 1
        
        # Maintain a counter for each transition type
        trans_key = (prev_state, curr_state)
        transition_counts[trans_key] = transition_counts.get(trans_key, 0) + 1
        
        # Only record state changes with cumulative count < 10
        if prev_state != curr_state and transition_counts[trans_key] < 10:
            # Get involved ASP residues (residues of carboxyl oxygens with distance < 3.2Å)
            bound_asp_indices = []
            for oxy_idx in np.where(dist_matrix[ca_idx] < 3.2)[0]:
                asp_residue = asp_oxygens[oxy_idx].residue.resid
                bound_asp_indices.append(asp_residue)
            
            rare_transitions.append({
                "frame": ts.frame,
                "ca_index": ca_idx,  # Modified to ca_index
                "from_state": prev_state,
                "to_state": curr_state,
                "asp_residues": bound_asp_indices,
                "transition_count": transition_counts[trans_key]  # Current cumulative count
            })
        
        # Update current state as previous state for next frame
        prev_states[ca_idx] = curr_state

# Calculate total occurrences of each state (derived from transition matrix)
state_counts = np.sum(transitions, axis=0)

# Output absolute transition count matrix
print("\n=== Absolute Transition Count Matrix ===")
abs_transitions = pd.DataFrame(
    transitions,
    index=[f"From {s}" for s in STATE_NAMES.values()],
    columns=[f"To {s}" for s in STATE_NAMES.values()]
)
pd.set_option('display.max_columns', None)
pd.set_option('display.expand_frame_repr', False)
print(abs_transitions)

# Output state statistics
print("\n=== State Statistics ===")
print(f"Total frames × Ca²⁺ count: {len(u.trajectory) * len(ca_atoms)}")
print("Occurrences per state:")
for s in STATE_NAMES:
    print(f"{STATE_NAMES[s]}: {state_counts[s]}")

# Output detailed transition event statistics (all possible transition combinations)
print("\n=== Detailed Transition Event Statistics ===")
for from_state in range(n_states):
    for to_state in range(n_states):
        if from_state != to_state and transitions[from_state, to_state] > 0:
            print(f"{STATE_NAMES[from_state]} → {STATE_NAMES[to_state]}: {transitions[from_state, to_state]}")

# Calculate and output net binding/dissociation trends
print("\n=== Net Change Trends ===")
for state in range(1, n_states):
    incoming = np.sum(transitions[:, state]) - transitions[state, state]
    outgoing = np.sum(transitions[state, :]) - transitions[state, state]
    net_change = incoming - outgoing
    print(f"{STATE_NAMES[state]}: Net change = {net_change} (incoming {incoming} times, outgoing {outgoing} times)")

# Output rare transition events to file
if rare_transitions:
    with open("rare_transitions.log", "w") as f:
        f.write("Frame\tCa²⁺ Index\tFrom State\tTo State\tInvolved ASP Residues\tCumulative Count\n")
        for event in rare_transitions:
            f.write(f"{event['frame']}\t{event['ca_index']}\t" 
                   f"{STATE_NAMES[event['from_state']]}\t"
                   f"{STATE_NAMES[event['to_state']]}\t"
                   f"{event['asp_residues']}\t"
                   f"{event['transition_count']}\n")
    print(f"\nRare transition events saved to rare_transitions.log (total {len(rare_transitions)} entries)")
else:
    print("\nNo events with transition count less than 10 found.")