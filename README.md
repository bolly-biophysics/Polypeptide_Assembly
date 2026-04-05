# Timescale Modulation of Charged Polypeptide Assembly through Divalent Cation Hydration Dynamics

This repository contains the complete analysis code for the manuscript:

> **"Timescale Modulation of Charged Polypeptide Assembly through Divalent Cation Hydration Dynamics"**  
> *Submitted to Nature Communications*

The code supports the key finding that **Mg²⁺** acts as a slow structural organizer via stepwise dehydration, while **Ca²⁺** functions as a fast phase separation trigger through rapid coordination switching.

---

## 📁 Repository Structure
├── README.md
├── atom_groups.dat # PLUMED atom indices (Mg²⁺, Asp, Arg)
├── plumed_input.dat # PLUMED metadynamics input
│
├── add_ligands.cpp # Add polyproline ligands to aptamer
├── add_ions.cpp # Add Mg²⁺/Cl⁻ ions to system
├── pbc_correction.py # Wrap water molecules around Mg²⁺
│
├── mg_asp_contacts.py # Mg²⁺-Asp contact count (2.8Å)
├── ca_asp_contacts.py # Ca²⁺-Asp contact count (3.2Å)
├── mg_cl_contacts.py # Mg²⁺-Cl⁻ contact count (3.0Å)
├── ca_cl_contacts.py # Ca²⁺-Cl⁻ contact count (3.5Å)
├── arg_asp_contacts.py # ARG-Asp salt bridge contacts
├── arg_asp_contacts_10us.py # ARG-Asp contacts over 10µs
│
├── mg_hydration.py # Mg²⁺ water coordination number
├── ca_hydration.py # Ca²⁺ water coordination number
├── cn_distribution.cpp # Probability distribution of CN
│
├── mg_state_transition.py # Mg²⁺ (0-6 Asp-bound) transitions
├── ca_state_transition.py # Ca²⁺ (0-8 Asp-bound) transitions
├── bridge_analysis.cpp # Mg²⁺/Cl⁻ bridging between charge centers
│
├── peptide_velocity_analysis.py # ARG CZ / ASP CG side-chain velocities
├── velocity_cn_distribution.cpp # 2D: velocity vs coordination number
├── velocity_data_io.cpp # Velocity & CN data I/O
│
├── exchange_rate_mg.cpp # Mg²⁺ hydration shell water exchange (OpenMP)
│
├── fes_reconstructor.py # WT-MTD → FES (kJ/mol ↔ kcal/mol)
│
├── ion_pair_tracking.py # Mg-Cl pair duration & frequency
│
└── contact_count.cpp # Mg²⁺/Cl⁻ contact counts
text


---

## 🧪 Dependencies

### Python (3.8+)
```bash
pip install MDAnalysis numpy pandas tqdm matplotlib

C++ (C++11 or later)
bash

# OpenMP required for exchange_rate_mg.cpp
g++ -fopenmp -O2 exchange_rate_mg.cpp -o exchange_rate_mg

PLUMED (optional, for metadynamics)

    Version 2.8+ with METAD module

🚀 Usage Examples
1. Trajectory preprocessing
bash

# Add ligands to aptamer PDB
g++ add_ligands.cpp -o add_ligands
./add_ligands

# Add Mg²⁺/Cl⁻ ions
g++ add_ions.cpp -o add_ions
./add_ions

# Apply PBC correction (wrap water around Mg²⁺)
python pbc_correction.py

2. Contact analysis
bash

python mg_asp_contacts.py      # Mg²⁺-Asp contacts
python ca_cl_contacts.py       # Ca²⁺-Cl⁻ contacts
python arg_asp_contacts.py     # ARG-Asp salt bridges

3. Hydration number distribution
bash

python mg_hydration.py      # Mg²⁺ water CN
python ca_hydration.py      # Ca²⁺ water CN

4. State transition analysis
bash

python mg_state_transition.py   # Mg²⁺: S0–S6
python ca_state_transition.py   # Ca²⁺: S0–S8

5. Water exchange rate (Mg²⁺)
bash

g++ -fopenmp exchange_rate_mg.cpp -o exchange_rate
./exchange_rate

6. Free energy reconstruction
bash

python fes_reconstructor.py \
    --negbias HILLS \
    --gamma 15.0 \
    --input-unit kJ/mol \
    --output-unit kcal/mol \
    --plot --contour --format png,eps

7. Ion pair tracking
bash

python ion_pair_tracking.py

8. Peptide side-chain velocity
bash

python peptide_velocity_analysis.py

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
mg_hydration.py	cn_mg.dat	Hydration number, probability
mg_state_transition.py	rare_transitions.log	Rare transition events
exchange_rate_mg.cpp	water_exchanges_mg.log	Exchange events (frame, Mg, water, type)
fes_reconstructor.py	*_fes_kcal/mol.dat	Free energy surface
ion_pair_tracking.py	ion_pairs_frequency.csv	Pair statistics
peptide_velocity_analysis.py	peptide_analysis_velocity_stats.txt	Velocity stats
📄 Notes

    All Python scripts include periodic boundary condition (PBC) handling via box=ts.dimensions

    The C++ XYZ reader assumes format: ID x y z per line

    WT-MTD conversion factor: γ/(γ-1) with default γ = 15

    Unit conversion: 1 kJ/mol = 0.239006 kcal/mol

📝 License & Citation

This code is released for reproducibility of the above manuscript.
If you use any part of it, please cite the original publication (to be added upon acceptance).

For questions, please open an issue or contact the corresponding author.
