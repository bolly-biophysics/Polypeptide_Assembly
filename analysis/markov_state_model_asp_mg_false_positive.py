import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import distances
from tqdm import tqdm

# Load trajectory
u = mda.Universe('../nonwater.prmtop', '../nonwater_10us.nc')
mg_atoms = u.select_atoms("name MG")
asp_oxygens = u.select_atoms("resname ASP and name OD1 OD2")

# State definition
STATE_NAMES = {
    0: "S0 (6H₂O)",
    1: "S1 (5H₂O+1COO⁻)",
    2: "S2 (4H₂O+2COO⁻)",
    3: "S3 (3H₂O+3COO⁻)",
    4: "S4 (2H₂O+4COO⁻)",
    5: "S5 (1H₂O+5COO⁻)",
    6: "S6 (0H₂O+6COO⁻)"
}

# Initialization
transitions = np.zeros((7, 7), dtype=int)  # 7x7 transition count matrix
state_counts = np.zeros(7, dtype=int)      # Total occurrences of each state
rare_transitions = []                      # Record rare transition events (count < 10)
prev_states = np.zeros(len(mg_atoms), dtype=int)  # Track previous state for each Mg²⁺ (initial S0)
transition_counts = {}                     # Track cumulative count for each transition type

print("Analyzing trajectory...")
for ts in tqdm(u.trajectory[1:], total=len(u.trajectory)-1):  # Start from second frame
    dist_matrix = distances.distance_array(
        mg_atoms.positions,
        asp_oxygens.positions,
        box=ts.dimensions
    )
    
    for mg_idx in range(len(mg_atoms)):
        prev_state = prev_states[mg_idx]  # Get previous state
        n_bound = np.sum(dist_matrix[mg_idx] < 2.8)  # Number of bound carboxyl oxygens
        curr_state = min(n_bound, 6)                # Current state (max is 6)
        
        # Record state transition (prev_state -> curr_state)
        transitions[prev_state, curr_state] += 1
        
        # Maintain a counter for each transition type
        trans_key = (prev_state, curr_state)
        transition_counts[trans_key] = transition_counts.get(trans_key, 0) + 1
        
        # Only record state changes with cumulative count < 10
        if prev_state != curr_state and transition_counts[trans_key] < 10:
            # Get involved ASP residues (residues of carboxyl oxygens with distance < 2.8Å)
            bound_asp_indices = []
            for oxy_idx in np.where(dist_matrix[mg_idx] < 2.8)[0]:
                asp_residue = asp_oxygens[oxy_idx].residue.resid
                bound_asp_indices.append(asp_residue)
            
            rare_transitions.append({
                "frame": ts.frame,
                "mg_index": mg_idx,
                "from_state": prev_state,
                "to_state": curr_state,
                "asp_residues": bound_asp_indices,
                "transition_count": transition_counts[trans_key]  # Current cumulative count
            })
        
        # Update current state as previous state for next frame
        prev_states[mg_idx] = curr_state

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
print(f"Total frames × Mg²⁺ count: {len(u.trajectory) * len(mg_atoms)}")
print("Occurrences per state:")
for s in STATE_NAMES:
    print(f"{STATE_NAMES[s]}: {state_counts[s]}")

# Output detailed transition event statistics (all possible transition combinations)
print("\n=== Detailed Transition Event Statistics ===")
for from_state in range(7):
    for to_state in range(7):
        if from_state != to_state and transitions[from_state, to_state] > 0:
            print(f"{STATE_NAMES[from_state]} → {STATE_NAMES[to_state]}: {transitions[from_state, to_state]}")

# Calculate and output net binding/dissociation trends
print("\n=== Net Change Trends ===")
for state in range(1, 7):
    incoming = np.sum(transitions[:, state]) - transitions[state, state]
    outgoing = np.sum(transitions[state, :]) - transitions[state, state]
    net_change = incoming - outgoing
    print(f"{STATE_NAMES[state]}: Net change = {net_change} (incoming {incoming} times, outgoing {outgoing} times)")

# Output rare transition events to file
if rare_transitions:
    with open("rare_transitions.log", "w") as f:
        f.write("Frame\tMg²⁺ Index\tFrom State\tTo State\tInvolved ASP Residues\tCumulative Count\n")
        for event in rare_transitions:
            f.write(f"{event['frame']}\t{event['mg_index']}\t" 
                   f"{STATE_NAMES[event['from_state']]}\t"
                   f"{STATE_NAMES[event['to_state']]}\t"
                   f"{event['asp_residues']}\t"
                   f"{event['transition_count']}\n")
    print(f"\nRare transition events saved to rare_transitions.log (total {len(rare_transitions)} entries)")
else:
    print("\nNo events with transition count less than 10 found.")