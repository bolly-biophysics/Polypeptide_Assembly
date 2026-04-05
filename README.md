# Timescale Modulation of Charged Polypeptide Assembly through Divalent Cation Hydration Dynamics

This repository contains the complete analysis code for the research:

> **"Timescale Modulation of Charged Polypeptide Assembly through Divalent Cation Hydration Dynamics"**

The code supports the key finding that **Mg²⁺** acts as a slow structural organizer via stepwise dehydration, while **Ca²⁺** functions as a fast phase separation trigger through rapid coordination switching.

---

## 📁 Repository Structure
├── README.md
│
├── cmd_amber_input_files/ # Conventional MD input files (AMBER)
│ ├── README.md
│ ├── min1.in # Initial minimization
│ ├── min2.in # Second minimization
│ ├── heat.in # Heating step
│ ├── equil_0.in # Equilibration 0
│ ├── equil_1.in # Equilibration 1
│ ├── equil_10.in # Equilibration 10
│ ├── equil_100.in # Equilibration 100
│ └── md.in # Production MD
│
├── force_field_construction/ # Force field parameter generation
│ ├── README.md
│ ├── ADD_CMAP.py # CMAP correction addition
│ ├── leap_126_tip3p.in # LEAP input (12-6 LJ potential, TIP3P water)
│ ├── leap_1264_spce.in # LEAP input (12-6-4 potential, SPC/E water)
│ ├── leap_1264_spce_cl_merz.in # LEAP input (12-6-4 potential, Cl⁻ Merz charges)
│ └── leap_1264_tip3p.in # LEAP input (12-6-4 potential, TIP3P water)
│
├── metad_gromacs_input_files/ # Metadynamics input files (GROMACS)
│ ├── README.md
│ ├── 10_polypeptides.gro # Coordinate file (10 chains)
│ ├── 10_polypeptides.top # Topology file
│ ├── atom_groups.dat # PLUMED atom indices (Mg²⁺, Asp, Arg)
│ ├── plumed.dat # PLUMED metadynamics input
│ ├── em.mdp # Energy minimization parameters
│ ├── nvt.mdp # NVT equilibration parameters
│ ├── npt.mdp # NPT equilibration parameters
│ ├── md.mdp # Production MD parameters
│ └── md.tpr # Run input file (portable)
│
├── polypeptides_system_construction/ # System construction and ion addition
│ ├── README.md
│ ├── single_polypeptides.pdb # Single polypeptide chain
│ ├── 10_polypeptides.pdb # 10 polypeptide chains
│ ├── 10_polypeptides_MgCl2.pdb # 10 chains with Mg²⁺/Cl⁻ ions
│ ├── 10_polypeptides_construction.cpp # C++ code for system construction
│ └── add_salt_scattered.cpp # Random salt addition
│
├── Analysis Scripts
│ │
│ ├── Trajectory Preprocessing
│ │ └── trajectory_water_pbc.py # PBC correction: wrap water around ions
│ │
│ ├── Contact Analysis
│ │ ├── mg_asp_contacts.py # Mg²⁺-Asp contact count (2.8Å)
│ │ ├── ca_asp_contacts.py # Ca²⁺-Asp contact count (3.2Å)
│ │ ├── mg_cl_contacts.py # Mg²⁺-Cl⁻ contact count (3.0Å)
│ │ ├── ca_cl_contacts.py # Ca²⁺-Cl⁻ contact count (3.5Å)
│ │ ├── arg_asp_contacts.py # ARG-Asp salt bridge contacts
│ │ ├── cluster_count.py # Cluster counting analysis
│ │ └── bridged_cation.cpp # Bridged cation analysis
│ │
│ ├── Coordination Number & Hydration
│ │ ├── mg_cn_probability.py # Mg²⁺ water coordination number distribution
│ │ ├── ca_cn_probability.py # Ca²⁺ water coordination number distribution
│ │ └── hydration_shell_exchange_rate.cpp # Mg²⁺ hydration shell water exchange (OpenMP)
│ │
│ ├── Markov State Model (State Transition)
│ │ ├── markov_state_model_asp_mg.py # Mg²⁺ state transitions (S0–S6)
│ │ ├── markov_state_model_asp_ca.py # Ca²⁺ state transitions (S0–S8)
│ │ ├── markov_state_model_asp_mg_false_positive.py # False positive analysis for Mg²⁺
│ │ ├── markov_state_model_asp_ca_false_positive.py # False positive analysis for Ca²⁺
│ │ ├── state_transition_waiting_time.py # Waiting time distribution
│ │ ├── plot_markov_state_transition_mg.sh # Plot script for Mg²⁺
│ │ └── plot_markov_state_transition_ca.sh # Plot script for Ca²⁺
│ │
│ ├── Velocity Analysis
│ │ ├── residue_velocity.py # ARG CZ / ASP CG side-chain velocities
│ │ ├── cation_velocity_cn.cpp # 2D: ion velocity vs coordination number
│ │ └── mg_cl_contact_detail.py # Detailed Mg-Cl contact with pair tracking
│ │
│ ├── Free Energy Surface
│ │ └── free_energy_metad.py # WT-MTD → FES (kJ/mol ↔ kcal/mol)
│ │
│ └── Utilities
│ └── (additional scripts as needed)
text


---

## ⚛️ Force Field Naming Convention

| Filename | Potential Form | Water Model | Description |
|----------|----------------|-------------|-------------|
| `leap_126_tip3p.in` | 12-6 Lennard-Jones | TIP3P | Standard non-bonded interactions |
| `leap_1264_spce.in` | 12-6-4 (LJ + ion-induced dipole) | SPC/E | Improved ion description (e.g., Ca²⁺) |
| `leap_1264_spce_cl_merz.in` | 12-6-4 with Merz charges | SPC/E | Cl⁻ parameter optimization |
| `leap_1264_tip3p.in` | 12-6-4 (LJ + ion-induced dipole) | TIP3P | Improved ion description |

**Note:** The 12-6-4 potential includes an additional \( r^{-4} \) term to account for ion-induced dipole interactions, which is essential for accurately describing divalent cations like Ca²⁺.

---

## 🧪 Dependencies

### Python (3.6+)
```bash
pip install MDAnalysis numpy pandas tqdm matplotlib

C++ (C++11 or later)
bash

# OpenMP required for hydration_shell_exchange_rate.cpp
g++ -fopenmp -O2 hydration_shell_exchange_rate.cpp -o hydration_exchange

AMBER (for conventional MD)

    AMBER 18 or later recommended

GROMACS (for metadynamics)

    GROMACS 2019 or later

    PLUMED 2.8+ with METAD module

🚀 Usage Examples
1. Force field construction
bash

cd force_field_construction/
# Add CMAP corrections
python ADD_CMAP.py
# Run LEAP to generate topology
tleap -f leap_1264_tip3p.in

2. System construction
bash

cd polypeptides_system_construction/
# Build 10 polypeptide chains
g++ 10_polypeptides_construction.cpp -o build_system
./build_system
# Add Mg²⁺ and Cl⁻ ions
g++ add_salt_scattered.cpp -o add_salt
./add_salt

3. Conventional MD (AMBER)
bash

cd cmd_amber_input_files/
# Minimization
mpirun -np 4 pmemd.MPI -i min1.in -o min1.out -p system.prmtop -c system.inpcrd
mpirun -np 4 pmemd.MPI -i min2.in -o min2.out -p system.prmtop -c min1.rst
# Heating
mpirun -np 4 pmemd.MPI -i heat.in -o heat.out -p system.prmtop -c min2.rst
# Equilibration
mpirun -np 4 pmemd.MPI -i equil_0.in -o equil_0.out -p system.prmtop -c heat.rst
# Production MD
mpirun -np 4 pmemd.MPI -i md.in -o md.out -p system.prmtop -c equil_100.rst

4. Metadynamics (GROMACS + PLUMED)
bash

cd metad_gromacs_input_files/
# Energy minimization
gmx mdrun -deffnm em -v -c 10_polypeptides.gro -s em.tpr -o em.trr -e em.edr -g em.log
# NVT equilibration
gmx mdrun -deffnm nvt -v -c em.gro -s nvt.tpr -o nvt.trr -e nvt.edr -g nvt.log
# NPT equilibration
gmx mdrun -deffnm npt -v -c nvt.gro -s npt.tpr -o npt.trr -e npt.edr -g npt.log
# Production metadynamics with PLUMED
gmx mdrun -deffnm md -v -c npt.gro -s md.tpr -o md.trr -e md.edr -g md.log -plumed plumed.dat

5. Trajectory preprocessing (PBC correction)
bash

python trajectory_water_pbc.py

6. Contact analysis
bash

python mg_asp_contacts.py          # Mg²⁺-Asp contacts
python ca_cl_contacts.py           # Ca²⁺-Cl⁻ contacts
python arg_asp_contacts.py         # ARG-Asp salt bridges

7. Coordination number distribution
bash

python mg_cn_probability.py        # Mg²⁺ water CN
python ca_cn_probability.py        # Ca²⁺ water CN

8. Markov State Model (state transition)
bash

# Mg²⁺ state transitions (S0–S6)
python markov_state_model_asp_mg.py
# Ca²⁺ state transitions (S0–S8)
python markov_state_model_asp_ca.py
# Waiting time analysis
python state_transition_waiting_time.py

9. Water exchange rate (Mg²⁺ hydration shell)
bash

g++ -fopenmp hydration_shell_exchange_rate.cpp -o hydration_exchange
./hydration_exchange

10. Free energy reconstruction (WT-MTD)
bash

python free_energy_metad.py \
    --negbias HILLS \
    --gamma 15.0 \
    --input-unit kJ/mol \
    --output-unit kcal/mol \
    --plot --contour --format png,eps

11. Velocity analysis
bash

python residue_velocity.py                # Side-chain velocities
python mg_cl_contact_detail.py            # Mg-Cl contact details

12. Plot transition matrices
bash

bash plot_markov_state_transition_mg.sh
bash plot_markov_state_transition_ca.sh

📊 Key Parameters (as used in the manuscript)
Analysis	Cutoff (Å)	Notes
Mg²⁺–Asp contact	2.8	First coordination shell
Ca²⁺–Asp contact	3.2	First coordination shell
Mg²⁺–water (hydration)	2.8	First hydration shell
Ca²⁺–water (hydration)	3.2	First hydration shell
Mg²⁺–Cl⁻ contact	3.0	Ion pair distance
Ca²⁺–Cl⁻ contact	3.5	Ion pair distance
ARG–Asp salt bridge	3.5	N–O distance

Simulation parameters:

    Time step per frame: 40 ps

    Total simulation: 10 µs

    Analysis stride: 10 frames (400 ps)

📊 Output Files
Script	Output file	Description
mg_asp_contacts.py	mg_asp_contacts.dat	Frame, total contacts
ca_asp_contacts.py	ca_asp_contacts.dat	Frame, total contacts
mg_cl_contacts.py	mg_cl_contacts.dat	Frame, total contacts
ca_cl_contacts.py	ca_cl_contacts.dat	Frame, total contacts
arg_asp_contacts.py	arg_asp_contacts.dat	Frame, total contacts
mg_cn_probability.py	cn_mg.dat	Hydration number, probability
ca_cn_probability.py	cn_ca.dat	Hydration number, probability
markov_state_model_asp_mg.py	mg_transition_matrix.npy	Transition matrix
markov_state_model_asp_ca.py	ca_transition_matrix.npy	Transition matrix
state_transition_waiting_time.py	waiting_times.log	Waiting time distribution
hydration_shell_exchange_rate.cpp	water_exchanges.log	Exchange events
free_energy_metad.py	*_fes_kcal/mol.dat	Free energy surface
mg_cl_contact_detail.py	ion_pairs_frequency.csv	Pair statistics
residue_velocity.py	velocity_stats.txt	Velocity statistics
cluster_count.py	cluster_counts.dat	Cluster statistics
bridged_cation.cpp	bridge_MG.dat, bridge_CL.dat	Bridging cation counts
📄 Notes

    All Python scripts are compatible with Python 3.6+

    All Python scripts include periodic boundary condition (PBC) handling via box=ts.dimensions

    The C++ XYZ reader assumes format: ID x y z per line

    WT-MTD conversion factor: γ/(γ-1) with default γ = 15

    Unit conversion: 1 kJ/mol = 0.239006 kcal/mol

📝 License & Citation

This code is released for reproducibility of the associated manuscript.
If you use any part of it, please cite the original publication.

For questions, please open an issue or contact the corresponding author.
