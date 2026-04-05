/**
 * @file exchange_rate_stride_mg_final.cpp
 * @brief Statistical analysis of water exchange in magnesium ion hydration shell (full output version)
 * 
 * Features:
 * 1. Count complete exchange cycles of water molecules in the first hydration shell of Mg ions (entry + exit = 1 exchange)
 * 2. Explicit initialization of frame 0 state to avoid false counting
 * 3. Support multi-threaded parallel computing
 * 4. Output system configuration and complete statistical results
 * 5. Support reading trajectory frames with stride
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <chrono>
#include <omp.h>

using namespace std;

// ================= System parameter configuration =============================
const int PROTEIN_ATOMS  = 3720;     // Number of protein atoms (0 if no protein)
const int MG_IONS        = 180;      // Number of magnesium ions
const int CL_IONS        = 360;      // Number of chloride ions
const int WATER_OXYGENS  = 64553;    // Number of water molecules (counted by oxygen atoms)
const int TOTAL_ATOMS    = PROTEIN_ATOMS + MG_IONS + CL_IONS + WATER_OXYGENS;

// Hydration shell parameters (adjusted according to literature)
const double rb = 2.8;        // Inner boundary of hydration shell (Å), distances less than this are definitely inside
const double rub = 2.8;       // Outer boundary of hydration shell (Å), distances greater than this are definitely outside
const double lambda = 1.00;   // Center position parameter for smoothing function
const double beta = 5.0;      // Steepness parameter for smoothing function

// Analysis parameters
const double dt = 40.0;       // Time interval per frame (ps)
const int n_frames = 50000;   // Total number of frames to analyze
const int stride = 10;        // Frame reading stride (default is 1, meaning read all frames)
const string history_file = "water_exchanges_mg_final_400ps.log";  // Exchange event record file

// ================= Data structures ================================
struct Atom {
    int id;      // Atom ID
    double x, y, z;  // 3D coordinates (Å)
};

/**
 * @brief Water molecule state tracking structure
 * Records whether currently inside hydration shell and associated Mg ion ID
 */
struct WaterState {
    bool current_state = false;  // true indicates inside hydration shell
    int current_mg_id = -1;      // Currently associated Mg ion ID (-1 indicates no association)
};

// ================= Core functions ================================
/**
 * @brief Calculate distance between two atoms
 * @param a Atom 1
 * @param b Atom 2
 * @return Euclidean distance (Å)
 */
double distance(const Atom& a, const Atom& b) {
    return sqrt((a.x-b.x)*(a.x-b.x) + 
           (a.y-b.y)*(a.y-b.y) + 
           (a.z-b.z)*(a.z-b.z));
}

/**
 * @brief Determine whether a water molecule is in the hydration shell of a magnesium ion
 * @param water Water oxygen atom
 * @param mg Magnesium ion
 * @return true if inside hydration shell
 * 
 * Uses sigmoid function for smooth transition region:
 * - Distance < rb: definitely inside
 * - Distance > rub: definitely outside
 * - Middle region: smooth transition
 */
bool inFirstHydrationShell(const Atom& water, const Atom& mg) {
    double r = distance(water, mg);
    if (r < rb) return true;
    if (r > rub) return false;
    return (1.0 + exp(beta * (r - lambda * rb))) < 2.0;
}

/**
 * @brief Read one frame of data from trajectory file
 * @param fin Input file stream
 * @return Array of atom coordinates
 * 
 * File format requirements:
 * ID x y z
 * ...
 */
vector<Atom> read_frame(ifstream& fin) {
    vector<Atom> frame(TOTAL_ATOMS);
    for (int i = 0; i < TOTAL_ATOMS; ++i) {
        fin >> frame[i].id >> frame[i].x >> frame[i].y >> frame[i].z;
    }
    return frame;
}

// ================= Main processing logic ================================
/**
 * @brief Process a single frame and count exchange events
 * @param frame_id Current frame number
 * @param frame Atom coordinate data
 * @param waters Water molecule state array
 * @param total_exchanges Exchange count accumulator
 * @param logfile Log file output stream
 * 
 * Core logic:
 * 1. Detect whether each water molecule is in the hydration shell of any Mg ion
 * 2. Record exchange events when state changes (0.5 exchange per crossing)
 * 3. Use OpenMP for parallel acceleration
 */
void process_frame(int frame_id, 
                  const vector<Atom>& frame,
                  vector<WaterState>& waters,
                  double& total_exchanges,
                  ofstream& logfile) 
{
    vector<Atom> mg_ions(MG_IONS);
    vector<Atom> water_oxygens(WATER_OXYGENS);

    // Extract magnesium ion coordinates (skip protein atoms)
    for (int i = 0; i < MG_IONS; ++i) {
        mg_ions[i] = frame[PROTEIN_ATOMS + i];
    }

    // Extract water oxygen atom coordinates (skip protein and ions)
    for (int i = 0; i < WATER_OXYGENS; ++i) {
        water_oxygens[i] = frame[PROTEIN_ATOMS + MG_IONS + CL_IONS + i];
    }

    double frame_exchanges = 0.0;

    // Parallel processing of all water molecules
    #pragma omp parallel for reduction(+:frame_exchanges)
    for (int w = 0; w < WATER_OXYGENS; ++w) {
        bool new_state = false;
        int new_mg_id = -1;

        // Check if within hydration shell of any magnesium ion
        for (int m = 0; m < MG_IONS; ++m) {
            if (inFirstHydrationShell(water_oxygens[w], mg_ions[m])) {
                new_state = true;
                new_mg_id = m;
                break;
            }
        }

        WaterState& ws = waters[w];
        
        // Record exchange event when state changes
        if (new_state != ws.current_state) {
            frame_exchanges += 0.5;  // 0.5 exchange per crossing
            #pragma omp critical
            {
                if (new_state) {
                    // Enter event: record associated Mg ion ID
                    logfile << frame_id << "," << new_mg_id << "," 
                           << w << ",enter\n";
                } else {
                    // Exit event: record Mg ion ID from previous frame (from ws.current_mg_id)
                    logfile << frame_id << "," << ws.current_mg_id << "," 
                           << w << ",exit\n";
                }
            }
            // Update state
            ws.current_state = new_state;
            ws.current_mg_id = new_mg_id;
        }
    }

    total_exchanges += frame_exchanges;
}

// ================= Main program ==================================
int main() {
    // Open trajectory file
    ifstream fin("system_all.xyz");
    if (!fin) {
        cerr << "Error: Unable to open trajectory file!" << endl;
        return 1;
    }

    // Output system configuration information
    cout << "============================================\n"
         << "      Mg Ion Hydration Shell Exchange Statistics\n"
         << "============================================\n"
         << "System configuration:\n"
         << "  - Number of Mg ions:       " << MG_IONS << "\n"
         << "  - Number of water molecules: " << WATER_OXYGENS << "\n"
         << "  - Total number of atoms:   " << TOTAL_ATOMS << "\n"
         << "Hydration shell parameters:\n"
         << "  - Inner boundary (rb):     " << rb << " Å\n"
         << "  - Outer boundary (rub):    " << rub << " Å\n"
         << "  - Transition zone center:  " << lambda*rb << " Å\n"
         << "Analysis parameters:\n"
         << "  - Total frames:            " << n_frames << "\n"
         << "  - Frame interval:          " << dt << " ps\n"
         << "  - Stride:                  " << stride << "\n"
         << "  - Effective frames:        " << n_frames/stride << "\n"
         << "  - Total duration:          " << n_frames*dt/1000 << " ns\n"
         << "  - Effective duration:      " << (n_frames/stride)*dt/1000 << " ns\n"
         << "  - Parallel threads:        " << omp_get_max_threads() << "\n"
         << "============================================\n";

    // Initialize output file and state array
    ofstream logfile(history_file);
    logfile << "frame,mg_id,water_id,event_type\n";
    vector<WaterState> waters(WATER_OXYGENS);
    double total_exchanges = 0.0;
    auto t_start = chrono::high_resolution_clock::now();

    // ===== Initialize frame 0 state =====
    auto init_frame = read_frame(fin); 
    int initial_hydrated = 0;
    for (int w = 0; w < WATER_OXYGENS; ++w) {
        for (int m = 0; m < MG_IONS; ++m) {
            if (inFirstHydrationShell(init_frame[PROTEIN_ATOMS + MG_IONS + CL_IONS + w],
                                    init_frame[PROTEIN_ATOMS + m])) {
                waters[w].current_state = true;
                waters[w].current_mg_id = m;  // Initialize associated Mg ion ID
                initial_hydrated++;
                break;
            }
        }
    }
    cout << "Initialization complete: Frame 0 has " << initial_hydrated 
         << " water molecules in hydration shell (theoretical coordination number: " 
         << MG_IONS*6 << ")\n" << endl;

    // ===== Main processing loop =====
    cout << "Starting trajectory data processing..." << endl;
    for (int f = 1; f < n_frames; ++f) {
        auto frame = read_frame(fin);
        if (fin.eof()) break;

        // Only process frames that meet the stride condition
        if (f % stride == 0) {
            process_frame(f, frame, waters, total_exchanges, logfile);
        }

        // Output progress every 1000 frames
        if (f % 1000 == 0) {
            auto duration = chrono::duration_cast<chrono::minutes>(
                chrono::high_resolution_clock::now() - t_start);
            cout << "Processed " << f << "/" << n_frames << " frames | "
                 << "Effective frames: " << f/stride << " | "
                 << "Elapsed time: " << duration.count() << " minutes | "
                 << "Exchange count: " << total_exchanges << endl;
        }
    }

    // ===== Final results output =====
    double total_ns = (n_frames / stride) * dt / 1000;  // Effective duration considering stride
    auto t_end = chrono::high_resolution_clock::now();
    auto total_time = chrono::duration_cast<chrono::minutes>(t_end - t_start);

    cout << "\n============================================\n"
         << "           Analysis Results Summary\n"
         << "============================================\n"
         << "Total exchanges:        " << total_exchanges << "\n"
         << "Average exchange rate:  " << total_exchanges/total_ns << " /ns\n"
         << "Exchange rate per ion:  " << total_exchanges/total_ns/MG_IONS << " /ns/ion\n"
         << "Total computation time: " << total_time.count() << " minutes\n"
         << "Event log saved to:     " << history_file << "\n"
         << "============================================\n";

    return 0;
}