#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Enhanced Mg-Cl contact analysis script
New features: Track contact IDs and statistical information for specific ion pairs
"""
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import distances
from MDAnalysis.lib.distances import transform_RtoS, transform_StoR
from collections import defaultdict, Counter
import time

# ====================== Parameter settings ======================
prmtop = "../nonwater.prmtop"
traj = "../nonwater_10us.nc"
output_file = "mg_cl_contact_results.dat"
pair_stats_file = "ion_pairs_frequency.csv"  # New: ion pair statistics file

mg_selection = "name MG or element Mg"
cl_selection = "name Cl- or element Cl"
contact_cutoff = 3.0  # Å
skip_frames = 1
min_duration = 40  # ps

# ====================== New module: Ion pair tracking ======================
class IonPairTracker:
    def __init__(self):
        self.pair_durations = defaultdict(list)  # Record duration for each ion pair
        self.pair_occurrences = Counter()       # Record occurrence count for each ion pair
    
    def add_contact(self, mg_id, cl_id, duration):
        """Record ion pair contact"""
        pair = (int(mg_id), int(cl_id))  # Convert to atom ID tuple
        self.pair_durations[pair].append(duration)
        self.pair_occurrences[pair] += 1
    
    def save_stats(self, filename):
        """Save ion pair statistics"""
        with open(filename, 'w') as f:
            f.write("Mg_ID,Cl_ID,Occurrences,Avg_Duration(ps),Total_Duration(ps)\n")
            for (mg_id, cl_id), counts in self.pair_occurrences.most_common():
                avg_dur = np.mean(self.pair_durations[(mg_id, cl_id)])
                total_dur = np.sum(self.pair_durations[(mg_id, cl_id)])
                f.write(f"{mg_id},{cl_id},{counts},{avg_dur:.1f},{total_dur:.1f}\n")

# ====================== Main analysis function ======================
def analyze_ion_contacts():
    print("=== Mg-Cl Contact Analysis (with ion pair tracking) ===")
    start_time = time.time()
    u = mda.Universe(prmtop, traj)
    
    # Initialize ion pair tracker
    pair_tracker = IonPairTracker()  # New
    
    # Select ions
    mg = u.select_atoms(mg_selection)
    cl = u.select_atoms(cl_selection)
    print(f"System contains: {len(mg)} Mg²⁺ (IDs: {mg.ids.tolist()})")
    print(f"          {len(cl)} Cl⁻ (IDs: {cl.ids.tolist()})")

    # Periodic boundary handling (omitted, same as previous script)
    # ...
    
    # Contact analysis
    contact_pairs = defaultdict(list)
    global_durations = []
    frame_time = u.trajectory.dt * skip_frames
    
    for i, ts in enumerate(u.trajectory):
        if i % skip_frames != 0:
            continue
            
        # Distance calculation
        dist = distances.distance_array(mg.positions, cl.positions, box=ts.dimensions)
        current_contacts = set(zip(*np.where(dist < contact_cutoff)))
        
        # Update contact status
        for pair in list(contact_pairs.keys()):
            mg_idx, cl_idx = pair
            if pair in current_contacts:
                contact_pairs[pair].append(i)
            else:
                duration = len(contact_pairs.pop(pair)) * frame_time
                if duration >= min_duration:
                    global_durations.append(duration)
                    # New: record specific ion pair
                    mg_id = mg[pair[0]].id
                    cl_id = cl[pair[1]].id
                    pair_tracker.add_contact(mg_id, cl_id, duration)
        
        # Add new contacts
        for pair in current_contacts:
            if pair not in contact_pairs:
                contact_pairs[pair] = [i]

    # Process final contacts
    for pair in contact_pairs:
        duration = len(contact_pairs[pair]) * frame_time
        if duration >= min_duration:
            global_durations.append(duration)
            mg_id = mg[pair[0]].id
            cl_id = cl[pair[1]].id
            pair_tracker.add_contact(mg_id, cl_id, duration)  # New

    # Output results
    if global_durations:
        global_durations = np.array(global_durations)
        print(f"\nTotal contact events: {len(global_durations)}")
        print(f"Average duration: {np.mean(global_durations):.1f} ps")
        
        # Save ion pair statistics
        pair_tracker.save_stats(pair_stats_file)  # New
        print(f"\nIon pair statistics saved to {pair_stats_file}")
        
        # Example: Display top 3 most common ion pairs
        top_pairs = pair_tracker.pair_occurrences.most_common(3)
        print("\nMost common ion pairs:")
        for (mg_id, cl_id), count in top_pairs:
            avg = np.mean(pair_tracker.pair_durations[(mg_id, cl_id)])
            print(f"Mg{mg_id}-Cl{cl_id}: {count} occurrences, average duration {avg:.1f}ps")

if __name__ == "__main__":
    analyze_ion_contacts()