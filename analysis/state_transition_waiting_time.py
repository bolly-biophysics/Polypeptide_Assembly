import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import distances
from tqdm import tqdm

# Load the trajectory
u = mda.Universe('../nonwater.prmtop', '../nonwater_10us.nc')
mg_atoms = u.select_atoms("name MG")
asp_oxygens = u.select_atoms("resname ASP and name OD1 OD2")

# Time between frames in ps (40ps as mentioned)
dt = 40  # ps

# Data structures to track binding events
binding_records = {i: [] for i in range(1, 4)}  # For 1st, 2nd, 3rd binding events
mg_history = {mg_idx: {'binding_steps': [], 'current_asp': set()} for mg_idx in range(len(mg_atoms))}

print("Analyzing trajectory...")
for ts in tqdm(u.trajectory, total=len(u.trajectory)):
    # Calculate distance matrix between Mg²⁺ and ASP oxygens
    dist_matrix = distances.distance_array(
        mg_atoms.positions,
        asp_oxygens.positions,
        box=ts.dimensions
    )
    
    for mg_idx in range(len(mg_atoms)):
        # Get current Mg²⁺ info
        mg_info = mg_history[mg_idx]
        prev_asp = mg_info['current_asp'].copy()
        
        # Find ASP oxygens within 2.8Å
        bound_oxy_indices = np.where(dist_matrix[mg_idx] < 2.8)[0]
        current_asp = set()
        
        # Get ASP residues involved in current binding
        for oxy_idx in bound_oxy_indices:
            asp_residue = asp_oxygens[oxy_idx].residue.resid
            current_asp.add(asp_residue)
        
        # Update current ASP set
        mg_info['current_asp'] = current_asp
        
        # Check for new binding events
        new_asp = current_asp - prev_asp
        
        if new_asp:
            binding_step = len(mg_info['binding_steps']) + 1
            if binding_step <= 3:  # We're interested in up to 3rd binding
                # Record the binding event
                event = {
                    'frame': ts.frame,
                    'time_ps': ts.frame * dt,
                    'mg_index': mg_idx,
                    'asp_residues': list(new_asp),
                    'total_bound_asp': list(current_asp)
                }
                
                # Calculate time since previous binding (if applicable)
                if mg_info['binding_steps']:
                    last_event = mg_info['binding_steps'][-1]
                    time_since_last = event['time_ps'] - last_event['time_ps']
                    event['time_since_last_ps'] = time_since_last
                else:
                    event['time_since_last_ps'] = 0
                
                mg_info['binding_steps'].append(event)
                
                # Add to our records if it's 1st, 2nd or 3rd binding
                if binding_step in binding_records:
                    binding_records[binding_step].append(event)

# Function to write binding events to file
def write_binding_events(events, filename):
    with open(filename, 'w') as f:
        f.write("Frame\tTime(ps)\tMg_Index\tNew_ASP_Residues\tAll_Bound_ASP\tTime_Since_Last(ps)\n")
        for event in events:
            f.write(f"{event['frame']}\t{event['time_ps']}\t{event['mg_index']}\t"
                   f"{event['asp_residues']}\t{event['total_bound_asp']}\t"
                   f"{event.get('time_since_last_ps', 0)}\n")

# Write the binding events to separate files
for step, events in binding_records.items():
    if events:
        filename = f"binding_step_{step}.log"
        write_binding_events(events, filename)
        print(f"Written {len(events)} {step}st/nd/rd binding events to {filename}")

# Additional analysis: Time intervals between binding steps
time_intervals = {1: [], 2: []}  # 1->2 and 2->3 intervals

for mg_idx, mg_info in mg_history.items():
    steps = mg_info['binding_steps']
    if len(steps) >= 2:
        time_intervals[1].append(steps[1]['time_since_last_ps'])
    if len(steps) >= 3:
        time_intervals[2].append(steps[2]['time_since_last_ps'])

# Write time interval statistics to file
with open("binding_time_intervals.log", 'w') as f:
    f.write("Binding_Step_Interval\tAverage_Time(ps)\tMedian_Time(ps)\tMin_Time(ps)\tMax_Time(ps)\tCount\n")
    for interval, times in time_intervals.items():
        if times:
            f.write(f"{interval}->{interval+1}\t{np.mean(times):.1f}\t{np.median(times):.1f}\t"
                   f"{min(times)}\t{max(times)}\t{len(times)}\n")

print("Analysis complete.")
