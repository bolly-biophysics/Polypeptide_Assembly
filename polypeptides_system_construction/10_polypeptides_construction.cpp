#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#define N_Apt_atom 372
#define N_Lig_atom 372
#define N_Lig_add 9
using namespace std;

struct Molecule_info
{
    string ATOM, atomName, resName, chainID, element;
    int serial, resSeq;
    double x, y, z, ocpy, bfac;
};

Molecule_info Apt[2000];
Molecule_info Lig[2000];

int main()
{
    ifstream infile1, infile2;
    infile1.open("single_polypeptides.pdb", ios::in);
    infile2.open("single_polypeptides.pdb", ios::in);

    if (!infile1)
    {
        cerr << "Open infile1 failure!" << endl;
        return -1;
    }

    // Format aptamer information into structure array
    int i = 1, j, k, t;
    while (!infile1.eof())
    {
        infile1 >> Apt[i].ATOM >> Apt[i].serial >> Apt[i].atomName >> Apt[i].resName >> Apt[i].resSeq >> Apt[i].x >> Apt[i].y >> Apt[i].z >> Apt[i].ocpy >> Apt[i].bfac;
        i++;
    }
    t = i - 2;
    cout << t << endl;
    infile1.close();

    if (!infile2)
    {
        cerr << "Open infile2 failure!" << endl;
        return -1;
    }

    // Format ligand information into structure array
    i = 1;
    while (!infile2.eof())
    {
        infile2 >> Lig[i].ATOM >> Lig[i].serial >> Lig[i].atomName >> Lig[i].resName >> Lig[i].resSeq >> Lig[i].x >> Lig[i].y >> Lig[i].z >> Lig[i].ocpy >> Lig[i].bfac;
        i++;
    }
    t = i - 2;
    cout << t << endl;
    infile2.close();

    // Calculate geometric center coordinates of aptamer
    double Apt_center_x, Apt_center_y, Apt_center_z;
    for (i = 1; i <= N_Apt_atom; i++)
    {
        Apt_center_x += Apt[i].x;
        Apt_center_y += Apt[i].y;
        Apt_center_z += Apt[i].z;
    }
    Apt_center_x /= N_Apt_atom;
    Apt_center_y /= N_Apt_atom;
    Apt_center_z /= N_Apt_atom;

    // Calculate aptamer direction vector
    double vec_Apt_x, vec_Apt_y, vec_Apt_z;
    vec_Apt_x = Apt[345].x - Apt[9].x;
    vec_Apt_y = Apt[345].y - Apt[9].y;
    vec_Apt_z = Apt[345].z - Apt[9].z;

    // Calculate maximum distance from aptamer center to boundary
    double d_max = 10;
    for (i = 1; i <= N_Apt_atom; i++)
    {
        double d = sqrt((Apt[i].x - Apt_center_x) * (Apt[i].x - Apt_center_x) + (Apt[i].y - Apt_center_y) * (Apt[i].y - Apt_center_y) + (Apt[i].z - Apt_center_z) * (Apt[i].z - Apt_center_z));
        if (d > d_max)
            d_max = d;
    }
    cout << d_max << endl;
    
    // Calculate geometric center coordinates of ligand
    double Lig_center_x, Lig_center_y, Lig_center_z;
    for (i = 1; i <= N_Lig_atom; i++)
    {
        Lig_center_x += Lig[i].x;
        Lig_center_y += Lig[i].y;
        Lig_center_z += Lig[i].z;
    }
    Lig_center_x /= N_Lig_atom;
    Lig_center_y /= N_Lig_atom;
    Lig_center_z /= N_Lig_atom;

    // Calculate offset vector between each ligand atom and its geometric center
    double offset_vec_x[2000], offset_vec_y[2000], offset_vec_z[2000];
    for (i = 1; i <= N_Lig_atom; i++)
    {
        offset_vec_x[i] = Lig[i].x - Lig_center_x;
        offset_vec_y[i] = Lig[i].y - Lig_center_y;
        offset_vec_z[i] = Lig[i].z - Lig_center_z;
    }

    // Add TER at the end of the original aptamer PDB file
    ofstream outfile;
    outfile.open("10_polypeptides.pdb", ios::out | ios::app);
    outfile << "TER" << endl;
    outfile.close();

    // Randomly generate new ligand geometric centers under given constraints, and append their full information to the aptamer PDB file
    double Lig_center_tmp_x, Lig_center_tmp_y, Lig_center_tmp_z;
    double vec_Lig_Apt_x, vec_Lig_Apt_y, vec_Lig_Apt_z;
    double vec_product, d_Lig_Apt, d_ij;
    double Lig_center_loc_x[2000], Lig_center_loc_y[2000], Lig_center_loc_z[2000];
    double Lig_atom_loc_x[2000], Lig_atom_loc_y[2000], Lig_atom_loc_z[2000];
    srand((unsigned int)time(NULL));
    for (i = 1; i <= N_Lig_add; i++)
    {
        bool gen_center_failed = true;
        while (gen_center_failed)
        {
            bool is_too_close = false;
            Lig_center_tmp_x = rand() % 200 - 100;
            Lig_center_tmp_y = rand() % 200 - 100;
            Lig_center_tmp_z = rand() % 200 - 100;
            vec_Lig_Apt_x = Lig_center_tmp_x - Apt_center_x;
            vec_Lig_Apt_y = Lig_center_tmp_y - Apt_center_y;
            vec_Lig_Apt_z = Lig_center_tmp_z - Apt_center_z;
            vec_product = vec_Apt_x * vec_Lig_Apt_x + vec_Apt_y * vec_Lig_Apt_y + vec_Apt_z * vec_Lig_Apt_z;
            d_Lig_Apt = sqrt((Lig_center_tmp_x - Apt_center_x) * (Lig_center_tmp_x - Apt_center_x) + (Lig_center_tmp_y - Apt_center_y) * (Lig_center_tmp_y - Apt_center_y) + (Lig_center_tmp_z - Apt_center_z) * (Lig_center_tmp_z - Apt_center_z));
            if (d_Lig_Apt > d_max && d_Lig_Apt < d_max + 20 && vec_product < 0.0001)
            {
                Lig_center_loc_x[i] = Lig_center_tmp_x;
                Lig_center_loc_y[i] = Lig_center_tmp_y;
                Lig_center_loc_z[i] = Lig_center_tmp_z;
            }
            else
            {
                continue;
            }
            for (j = 1; j < i; j++)
            {
                d_ij = sqrt((Lig_center_loc_x[i] - Lig_center_loc_x[j]) * (Lig_center_loc_x[i] - Lig_center_loc_x[j]) + (Lig_center_loc_y[i] - Lig_center_loc_y[j]) * (Lig_center_loc_y[i] - Lig_center_loc_y[j]) + (Lig_center_loc_z[i] - Lig_center_loc_z[j]) * (Lig_center_loc_z[i] - Lig_center_loc_z[j]));
                if (d_ij < d_max)
                {
                    is_too_close = true;
                    break;
                }
            }
            if (is_too_close)
            {
                continue;
            }
            else
            {
                for (k = 1; k <= N_Lig_atom; k++)
                {
                    Lig_atom_loc_x[k] = Lig_center_loc_x[i] + offset_vec_x[k];
                    Lig_atom_loc_y[k] = Lig_center_loc_y[i] + offset_vec_y[k];
                    Lig_atom_loc_z[k] = Lig_center_loc_z[i] + offset_vec_z[k];
                    ofstream outfile;
                    outfile.open("10_polypeptides.pdb", ios::out | ios::app);
                    if (Lig[k].atomName.length() <= 3)
                        outfile << setw(6) << "ATOM  " << setw(5) << Lig[k].serial + i * 372 << "  " << setw(3) << setiosflags(ios::left) << Lig[k].atomName << " " << setw(3) << Lig[k].resName << " " << " " << setw(4) << setiosflags(ios::right) << Lig[k].resSeq << "    " << setw(8) << fixed << setprecision(3) << Lig_atom_loc_x[k] << setw(8) << fixed << setprecision(3) << Lig_atom_loc_y[k] << setw(8) << fixed << setprecision(3) << Lig_atom_loc_z[k] << setw(6) << fixed << setprecision(2) << Lig[k].ocpy << setw(6) << fixed << setprecision(2) << Lig[k].bfac << endl;
                    else
                        outfile << setw(6) << "ATOM  " << setw(5) << Lig[k].serial + i * 372 << " " << setw(4) << setiosflags(ios::left) << Lig[k].atomName << " " << setw(3) << Lig[k].resName << " " << " " << setw(4) << setiosflags(ios::right) << Lig[k].resSeq << "    " << setw(8) << fixed << setprecision(3) << Lig_atom_loc_x[k] << setw(8) << fixed << setprecision(3) << Lig_atom_loc_y[k] << setw(8) << fixed << setprecision(3) << Lig_atom_loc_z[k] << setw(6) << fixed << setprecision(2) << Lig[k].ocpy << setw(6) << fixed << setprecision(2) << Lig[k].bfac << endl;
                    outfile.close();
                }
                ofstream outfile;
                outfile.open("10_polypeptides.pdb", ios::out | ios::app);
                outfile << "TER" << endl;
                outfile.close();
                gen_center_failed = false;
            }
        }
    }

    // Add END at the end of the aptamer PDB file
    outfile.open("10_polypeptides.pdb", ios::out | ios::app);
    outfile << "END" << endl;
    outfile.close();

    return 0;
}