#include<iostream>
#include<fstream>
#include<cmath>
#define N_tot 940
#define N_PEP 400
#define N_MG 180
#define N_CL 360
#define R_MG 2.8
#define R_CL 4
using namespace std;

float sc1_x[26000][500], sc1_y[26000][500], sc1_z[26000][500], sc2_x[26000][500], sc2_y[26000][500], sc2_z[26000][500];
float density_MG[500][500][26000], density_MG_avg[500][500], density_MG_std[500][500];
float density_CL[500][500][26000], density_CL_avg[500][500], density_CL_std[500][500];

struct example
{
    int serial;
    float x, y, z;
};

example atom_info[26000][1000];

int main()
{
    ifstream infile;
    infile.open("ion_env.xyz");

    if (!infile)
    {
        cerr << "Open infile failure!" << endl;
        return -1;
    }

    // Format atom coordinates into structure array
    int i = 1, j, k, l, t;
    float aa;
    while (!infile.eof())
    {
        aa = (i + 0.0) / (N_tot + 0.0); 
        k = fmod(i, N_tot);
        if (k == 0) 
        {
            k = N_tot;
            j = floor(aa);
        }
        else 
            j = floor(aa) + 1;
        infile >> atom_info[j][k].serial >> atom_info[j][k].x >> atom_info[j][k].y >> atom_info[j][k].z;
        i++;
    }
    t = j - 1;
    cout << t << endl;
    infile.close();

    ofstream outfile_1, outfile_2;
    outfile_1.open("bridge_MG.dat", ios::out | ios::ate);
    outfile_2.open("bridge_CL.dat", ios::out | ios::ate);

    // Calculate side chain charge center 1 over time
    for (i = 1; i <= t; i++)
    {
        for (j = 1; j <= N_PEP; j = j + 2)
        {
            sc1_x[i][(j + 1) / 2] = atom_info[i][j].x;
            sc1_y[i][(j + 1) / 2] = atom_info[i][j].y;
            sc1_z[i][(j + 1) / 2] = atom_info[i][j].z;
        }
    }

    // Calculate side chain charge center 2 over time
    for (i = 1; i <= t; i++)
    {
        for (j = 2; j <= N_PEP; j = j + 2)
        {
            sc2_x[i][(j + 0) / 2] = atom_info[i][j].x;
            sc2_y[i][(j + 0) / 2] = atom_info[i][j].y;
            sc2_z[i][(j + 0) / 2] = atom_info[i][j].z;
        }
    }
    
    // Calculate number of bridging MG ions between each pair of charge centers over time
    float d_i1, d_i2, d_j1, d_j2;
    for (i = 1; i <= 200; i++)
    {
        for (j = 1; j <= 200; j++)
        {
            for (k = 1; k <= t; k++)
            {
                density_MG[i][j][k] = 0;
                for (l = N_PEP + 1; l <= N_PEP + N_MG; l++)
                {
                    d_i1 = sqrt((atom_info[k][l].x - sc1_x[k][i]) * (atom_info[k][l].x - sc1_x[k][i]) + (atom_info[k][l].y - sc1_y[k][i]) * (atom_info[k][l].y - sc1_y[k][i]) + (atom_info[k][l].z - sc1_z[k][i]) * (atom_info[k][l].z - sc1_z[k][i]));
                    d_i2 = sqrt((atom_info[k][l].x - sc2_x[k][i]) * (atom_info[k][l].x - sc2_x[k][i]) + (atom_info[k][l].y - sc2_y[k][i]) * (atom_info[k][l].y - sc2_y[k][i]) + (atom_info[k][l].z - sc2_z[k][i]) * (atom_info[k][l].z - sc2_z[k][i]));
                    d_j1 = sqrt((atom_info[k][l].x - sc1_x[k][j]) * (atom_info[k][l].x - sc1_x[k][j]) + (atom_info[k][l].y - sc1_y[k][j]) * (atom_info[k][l].y - sc1_y[k][j]) + (atom_info[k][l].z - sc1_z[k][j]) * (atom_info[k][l].z - sc1_z[k][j]));
                    d_j2 = sqrt((atom_info[k][l].x - sc2_x[k][j]) * (atom_info[k][l].x - sc2_x[k][j]) + (atom_info[k][l].y - sc2_y[k][j]) * (atom_info[k][l].y - sc2_y[k][j]) + (atom_info[k][l].z - sc2_z[k][j]) * (atom_info[k][l].z - sc2_z[k][j]));
                    if (d_i1 < R_MG || d_i2 < R_MG)
                    {
                        if (d_j1 < R_MG || d_j2 < R_MG)
                            density_MG[i][j][k]++;
                    }
                }
            }
        }
    }

    // Calculate average number of bridging MG ions between each pair of charge centers
    for (i = 1; i <= 200; i++)
    {
        for (j = 1; j <= 200; j++)
        {
            density_MG_avg[i][j] = 0;
            for (k = 1; k <= t; k++)
            {
                density_MG_avg[i][j] += density_MG[i][j][k];
            }
            density_MG_avg[i][j] = density_MG_avg[i][j] / t;
        }
    }

    // Calculate standard deviation of bridging MG ions between each pair of charge centers
    for (i = 1; i <= 200; i++)
    {
        for (j = 1; j <= 200; j++)
        {        
            density_MG_std[i][j] = 0;
            for (k = 1; k <= t; k++)
            {
                density_MG_std[i][j] += (density_MG[i][j][k] - density_MG_avg[i][j]) * (density_MG[i][j][k] - density_MG_avg[i][j]);
            }
            density_MG_std[i][j] = sqrt(density_MG_std[i][j] / t);
            outfile_1 << i << " " << j << " " << density_MG_avg[i][j] << " " << density_MG_std[i][j] << endl;
        }
        outfile_1 << endl;
    }
    outfile_1.close();

    // Calculate number of bridging CL ions between each pair of charge centers over time
    for (i = 1; i <= 200; i++)
    {
        for (j = 1; j <= 200; j++)
        {
            for (k = 1; k <= t; k++)
            {
                density_CL[i][j][k] = 0;
                for (l = N_PEP + N_MG + 1; l <= N_PEP + N_MG + N_CL; l++)
                {
                    d_i1 = sqrt((atom_info[k][l].x - sc1_x[k][i]) * (atom_info[k][l].x - sc1_x[k][i]) + (atom_info[k][l].y - sc1_y[k][i]) * (atom_info[k][l].y - sc1_y[k][i]) + (atom_info[k][l].z - sc1_z[k][i]) * (atom_info[k][l].z - sc1_z[k][i]));
                    d_i2 = sqrt((atom_info[k][l].x - sc2_x[k][i]) * (atom_info[k][l].x - sc2_x[k][i]) + (atom_info[k][l].y - sc2_y[k][i]) * (atom_info[k][l].y - sc2_y[k][i]) + (atom_info[k][l].z - sc2_z[k][i]) * (atom_info[k][l].z - sc2_z[k][i]));
                    d_j1 = sqrt((atom_info[k][l].x - sc1_x[k][j]) * (atom_info[k][l].x - sc1_x[k][j]) + (atom_info[k][l].y - sc1_y[k][j]) * (atom_info[k][l].y - sc1_y[k][j]) + (atom_info[k][l].z - sc1_z[k][j]) * (atom_info[k][l].z - sc1_z[k][j]));
                    d_j2 = sqrt((atom_info[k][l].x - sc2_x[k][j]) * (atom_info[k][l].x - sc2_x[k][j]) + (atom_info[k][l].y - sc2_y[k][j]) * (atom_info[k][l].y - sc2_y[k][j]) + (atom_info[k][l].z - sc2_z[k][j]) * (atom_info[k][l].z - sc2_z[k][j]));
                    if (d_i1 < R_CL || d_i2 < R_CL)
                    {
                        if (d_j1 < R_CL || d_j2 < R_CL)
                            density_CL[i][j][k]++;
                    }
                }
            }
        }
    }

    // Calculate average number of bridging CL ions between each pair of charge centers
    for (i = 1; i <= 200; i++)
    {
        for (j = 1; j <= 200; j++)
        {
            density_CL_avg[i][j] = 0;
            for (k = 1; k <= t; k++)
            {
                density_CL_avg[i][j] += density_CL[i][j][k];
            }
            density_CL_avg[i][j] = density_CL_avg[i][j] / t;
        }
    }

    // Calculate standard deviation of bridging CL ions between each pair of charge centers
    for (i = 1; i <= 200; i++)
    {
        for (j = 1; j <= 200; j++)
        {        
            density_CL_std[i][j] = 0;
            for (k = 1; k <= t; k++)
            {
                density_CL_std[i][j] += (density_CL[i][j][k] - density_CL_avg[i][j]) * (density_CL[i][j][k] - density_CL_avg[i][j]);
            }
            density_CL_std[i][j] = sqrt(density_CL_std[i][j] / t);
            outfile_2 << i << " " << j << " " << density_CL_avg[i][j] << " " << density_CL_std[i][j] << endl;
        }
        outfile_2 << endl;
    }
    outfile_2.close();

    return 0;
}