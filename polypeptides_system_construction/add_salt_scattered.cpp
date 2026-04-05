#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#define N_BioSystem_atom 3720
#define N_Cation_added 180
#define N_Anion_added 360
#define d_cut 10
using namespace std;

struct Atom_info
{
    string ATOM, atomName, resName, chainID;
    int serial, resSeq;
    double x, y, z, ocpy, bfac;
};

Atom_info BioSystem[4000];

int main()
{
    ifstream infile;
    infile.open("10_polypeptides_MgCl2.pdb", ios::in);

    if (!infile)
    {
        cerr << "Open infile failure!" << endl;
        return -1;
    }

    // Format biological system information into structure array
    int i = 1, j;
    while (!infile.eof())
    {
        infile >> BioSystem[i].ATOM >> BioSystem[i].serial >> BioSystem[i].atomName >> BioSystem[i].resName >> BioSystem[i].chainID >> BioSystem[i].resSeq >> BioSystem[i].x >> BioSystem[i].y >> BioSystem[i].z >> BioSystem[i].ocpy >> BioSystem[i].bfac;
        i++;
    }
    infile.close();

    // Calculate geometric center coordinates of biological system
    double BioSystem_center_x, BioSystem_center_y, BioSystem_center_z;
    for (i = 1; i <= N_BioSystem_atom; i++)
    {
        BioSystem_center_x += BioSystem[i].x;
        BioSystem_center_y += BioSystem[i].y;
        BioSystem_center_z += BioSystem[i].z;
    }
    BioSystem_center_x /= N_BioSystem_atom;
    BioSystem_center_y /= N_BioSystem_atom;
    BioSystem_center_z /= N_BioSystem_atom;

    // Calculate maximum distance from biological system center to boundary
    double d_max_BioSystem = 10;
    for (i = 1; i <= N_BioSystem_atom; i++)
    {
        double d = sqrt((BioSystem[i].x - BioSystem_center_x) * (BioSystem[i].x - BioSystem_center_x) + (BioSystem[i].y - BioSystem_center_y) * (BioSystem[i].y - BioSystem_center_y) + (BioSystem[i].z - BioSystem_center_z) * (BioSystem[i].z - BioSystem_center_z));
        if (d > d_max_BioSystem)
            d_max_BioSystem = d;
    }

    // Add TER at the end of the original biological system PDB file
    ofstream outfile;
    outfile.open("10_polypeptides_MgCl2.pdb", ios::out | ios::app);
    outfile << "TER" << endl;
    outfile.close();

    // Randomly generate new cation geometric centers under given constraints, and append their full information to the biological system PDB file
    double Cation_center_tmp_x, Cation_center_tmp_y, Cation_center_tmp_z;
    double d_Cation_BioSystem_center, d_Cation_BioSystem_atom, d_Cation_Cation;
    double Cation_center_loc_x[1000], Cation_center_loc_y[1000], Cation_center_loc_z[1000];
    srand((unsigned int)time(NULL));
    for (i = 1; i <= N_Cation_added; i++)
    {
        bool gen_Cation_center_failed = true;
        while (gen_Cation_center_failed)
        {
            bool Cation_Cation_is_too_close = false;
            bool Cation_BioSystem_is_too_close = false;
            // Randomly generate cation coordinates
            Cation_center_tmp_x = rand() % 2000 - 1000;
            Cation_center_tmp_y = rand() % 2000 - 1000;
            Cation_center_tmp_z = rand() % 2000 - 1000;
            // Calculate distance between this coordinate and the biological system geometric center
            d_Cation_BioSystem_center = sqrt((Cation_center_tmp_x - BioSystem_center_x) * (Cation_center_tmp_x - BioSystem_center_x) + (Cation_center_tmp_y - BioSystem_center_y) * (Cation_center_tmp_y - BioSystem_center_y) + (Cation_center_tmp_z - BioSystem_center_z) * (Cation_center_tmp_z - BioSystem_center_z));
            // If the coordinate falls within the preset range, assign its index and proceed to the next step
            if (d_Cation_BioSystem_center > 0 && d_Cation_BioSystem_center < d_max_BioSystem)
            {
                Cation_center_loc_x[i] = Cation_center_tmp_x;
                Cation_center_loc_y[i] = Cation_center_tmp_y;
                Cation_center_loc_z[i] = Cation_center_tmp_z;
            }
            // Otherwise, re-run the loop to generate cation coordinates
            else
            {
                continue;
            }
            // Check whether the newly added cation is too close to any atom of the biological system
            for (j = 1; j <= N_BioSystem_atom; j++)
            {
                d_Cation_BioSystem_atom = sqrt((Cation_center_loc_x[i] - BioSystem[j].x) * (Cation_center_loc_x[i] - BioSystem[j].x) + (Cation_center_loc_y[i] - BioSystem[j].y) * (Cation_center_loc_y[i] - BioSystem[j].y) + (Cation_center_loc_z[i] - BioSystem[j].z) * (Cation_center_loc_z[i] - BioSystem[j].z));
                if (d_Cation_BioSystem_atom < d_cut)
                {
                    Cation_BioSystem_is_too_close = true;
                    break;
                }
            }
            // If the cation is too close to the biological system, re-run the loop to generate cation coordinates
            if (Cation_BioSystem_is_too_close)
            {
                continue;
            }
            // Check whether the newly added cation is too close to any previously added cation
            for (j = 1; j < i; j++)
            {
                d_Cation_Cation = sqrt((Cation_center_loc_x[i] - Cation_center_loc_x[j]) * (Cation_center_loc_x[i] - Cation_center_loc_x[j]) + (Cation_center_loc_y[i] - Cation_center_loc_y[j]) * (Cation_center_loc_y[i] - Cation_center_loc_y[j]) + (Cation_center_loc_z[i] - Cation_center_loc_z[j]) * (Cation_center_loc_z[i] - Cation_center_loc_z[j]));
                if (d_Cation_Cation < 5)
                {
                    Cation_Cation_is_too_close = true;
                    break;
                }
            }
            // If cations are too close to each other, re-run the loop to generate cation coordinates
            if (Cation_Cation_is_too_close)
            {
                continue;
            }
            // If all the above checks pass successfully, append the cation coordinates to the biological system PDB file
            outfile.open("10_polypeptides_MgCl2.pdb", ios::out | ios::app);
            outfile << setw(6) << "HETATM" << setw(5) << i << "  " << setw(3) << setiosflags(ios::left) << "MG" << " " << setw(3) << resetiosflags(ios::left) << "MG" << " " << "X" << setw(4) << resetiosflags(ios::right) << i << "    " << setw(8) << fixed << setprecision(3) << Cation_center_loc_x[i] << setw(8) << fixed << setprecision(3) << Cation_center_loc_y[i] << setw(8) << fixed << setprecision(3) << Cation_center_loc_z[i] << setw(6) << fixed << setprecision(2) << "1.00" << setw(6) << fixed << setprecision(2) << endl;
            outfile.close();
            gen_Cation_center_failed = false;
        }
    }

    // Randomly generate new anion geometric centers under given constraints, and append their full information to the biological system PDB file
    double Anion_center_tmp_x, Anion_center_tmp_y, Anion_center_tmp_z;
    double d_Anion_BioSystem_center, d_Anion_BioSystem_atom, d_Anion_Cation, d_Anion_Anion;
    double Anion_center_loc_x[1000], Anion_center_loc_y[1000], Anion_center_loc_z[1000];
    srand((unsigned int)time(NULL));
    for (i = 1; i <= N_Anion_added; i++)
    {
        bool gen_Anion_center_failed = true;
        while (gen_Anion_center_failed)
        {
            bool Anion_Anion_is_too_close = false;
            bool Anion_Cation_is_too_close = false;
            bool Anion_BioSystem_is_too_close = false;
            // Randomly generate anion coordinates
            Anion_center_tmp_x = rand() % 2000 - 1000;
            Anion_center_tmp_y = rand() % 2000 - 1000;
            Anion_center_tmp_z = rand() % 2000 - 1000;
            // Calculate distance between this coordinate and the biological system geometric center
            d_Anion_BioSystem_center = sqrt((Anion_center_tmp_x - BioSystem_center_x) * (Anion_center_tmp_x - BioSystem_center_x) + (Anion_center_tmp_y - BioSystem_center_y) * (Anion_center_tmp_y - BioSystem_center_y) + (Anion_center_tmp_z - BioSystem_center_z) * (Anion_center_tmp_z - BioSystem_center_z));
            // If the coordinate falls within the preset range, assign its index and proceed to the next step
            if (d_Anion_BioSystem_center > 0 && d_Anion_BioSystem_center < d_max_BioSystem)
            {
                Anion_center_loc_x[i] = Anion_center_tmp_x;
                Anion_center_loc_y[i] = Anion_center_tmp_y;
                Anion_center_loc_z[i] = Anion_center_tmp_z;
            }
            // Otherwise, re-run the loop to generate anion coordinates
            else
            {
                continue;
            }
            // Check whether the newly added anion is too close to any atom of the biological system
            for (j = 1; j <= N_BioSystem_atom; j++)
            {
                d_Anion_BioSystem_atom = sqrt((Anion_center_loc_x[i] - BioSystem[j].x) * (Anion_center_loc_x[i] - BioSystem[j].x) + (Anion_center_loc_y[i] - BioSystem[j].y) * (Anion_center_loc_y[i] - BioSystem[j].y) + (Anion_center_loc_z[i] - BioSystem[j].z) * (Anion_center_loc_z[i] - BioSystem[j].z));
                if (d_Anion_BioSystem_atom < d_cut)
                {
                    Anion_BioSystem_is_too_close = true;
                    break;
                }
            }
            // If the anion is too close to the biological system, re-run the loop to generate anion coordinates
            if (Anion_BioSystem_is_too_close)
            {
                continue;
            }
            // Check whether the newly added anion is too close to any previously added cation
            for (j = 1; j <= N_Cation_added; j++)
            {
                d_Anion_Cation = sqrt((Anion_center_loc_x[i] - Cation_center_loc_x[j]) * (Anion_center_loc_x[i] - Cation_center_loc_x[j]) + (Anion_center_loc_y[i] - Cation_center_loc_y[j]) * (Anion_center_loc_y[i] - Cation_center_loc_y[j]) + (Anion_center_loc_z[i] - Cation_center_loc_z[j]) * (Anion_center_loc_z[i] - Cation_center_loc_z[j]));
                if (d_Anion_Cation < 15)
                {
                    Anion_Cation_is_too_close = true;
                    break;
                }
            }
            // If the anion is too close to a cation, re-run the loop to generate anion coordinates
            if (Anion_Cation_is_too_close)
            {
                continue;
            }
            // Check whether the newly added anion is too close to any previously added anion
            for (j = 1; j < i; j++)
            {
                d_Anion_Anion = sqrt((Anion_center_loc_x[i] - Anion_center_loc_x[j]) * (Anion_center_loc_x[i] - Anion_center_loc_x[j]) + (Anion_center_loc_y[i] - Anion_center_loc_y[j]) * (Anion_center_loc_y[i] - Anion_center_loc_y[j]) + (Anion_center_loc_z[i] - Anion_center_loc_z[j]) * (Anion_center_loc_z[i] - Anion_center_loc_z[j]));
                if (d_Anion_Anion < 5)
                {
                    Anion_Anion_is_too_close = true;
                    break;
                }
            }
            // If anions are too close to each other, re-run the loop to generate anion coordinates
            if (Anion_Anion_is_too_close)
            {
                continue;
            }
            // If all the above checks pass successfully, append the anion coordinates to the biological system PDB file
            outfile.open("10_polypeptides_MgCl2.pdb", ios::out | ios::app);
            outfile << setw(6) << "HETATM" << setw(5) << i << "  " << setw(3) << setiosflags(ios::left) << "Cl-" << " " << setw(3) << resetiosflags(ios::left) << "Cl-" << " " << "X" << setw(4) << resetiosflags(ios::right) << i << "    " << setw(8) << fixed << setprecision(3) << Anion_center_loc_x[i] << setw(8) << fixed << setprecision(3) << Anion_center_loc_y[i] << setw(8) << fixed << setprecision(3) << Anion_center_loc_z[i] << setw(6) << fixed << setprecision(2) << "1.00" << setw(6) << fixed << setprecision(2) << endl;
            outfile.close();
            gen_Anion_center_failed = false;
        }
    }

    // Add END at the end of the total simulation system PDB file
    outfile.open("10_polypeptides_MgCl2.pdb", ios::out | ios::app);
    outfile << std::left << "END" << endl;
    outfile.close();

    return 0;
}