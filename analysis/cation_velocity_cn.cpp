#include<iostream>
#include<fstream>
#include<cmath>
#define N_tot 68813
#define N_PEP 3720
#define N_MG 180
#define N_CL 360
#define r_MG 2.8
#define r_CL 4
using namespace std;

float v_MG[26000][200], v_CL[26000][400];
float cn_MG[26000][200], cn_CL[26000][400];
float P_v_cn[100][100], P_v_cn_CL[100][100];

struct example
{
    int serial;
    float x, y, z;
};

example atom_info[26000][70000];

int main()
{
    ifstream infile;
    infile.open("system_all.xyz");

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
    outfile_1.open("p_v_cn.dat", ios::out | ios::ate);
    outfile_2.open("p_v_cn_CL.dat", ios::out | ios::ate);

    // Calculate Mg ion velocity over time
    for (i = 1; i <= t - 1; i++)
    {
        for (j = N_PEP + 1; j <= N_PEP + N_MG; j++)
        {
            v_MG[i][j - N_PEP] = sqrt((atom_info[i + 1][j].x - atom_info[i][j].x) * (atom_info[i + 1][j].x - atom_info[i][j].x) + (atom_info[i + 1][j].y - atom_info[i][j].y) * (atom_info[i + 1][j].y - atom_info[i][j].y) + (atom_info[i + 1][j].z - atom_info[i][j].z) * (atom_info[i + 1][j].z - atom_info[i][j].z));
        }
    }

    // Calculate water coordination number of Mg ions
    int count;
    float dist_O_WAT;
    for (i = 1; i <= t; i++)
    {
        for (j = N_PEP + 1; j <= N_PEP + N_MG; j++)
        {
            count = 0;
            for (k = N_PEP + N_MG + N_CL + 1; k <= N_tot; k++)
            {
                dist_O_WAT = sqrt((atom_info[i][k].x - atom_info[i][j].x) * (atom_info[i][k].x - atom_info[i][j].x) + (atom_info[i][k].y - atom_info[i][j].y) * (atom_info[i][k].y - atom_info[i][j].y) + (atom_info[i][k].z - atom_info[i][j].z) * (atom_info[i][k].z - atom_info[i][j].z));
                if (dist_O_WAT < r_MG)
                    count++;
            }
            cn_MG[i][j - N_PEP] = count;
        }
    }
    
    // Calculate 2D distribution of Mg ion velocity vs water coordination number
    for (i = 0; i <= 40; i++)
    {
        for (j = 0; j <= 7; j++)
        {
            count = 0;
            for (k = 1; k <= t - 1; k++)
            {
                for (l = 1; l <= N_MG; l++)
                {
                    if (v_MG[k][l] >= i && v_MG[k][l] < i + 1)
                    {
                        if (cn_MG[k][l] >= (float)j && cn_MG[k][l] < (float)(j + 1))
                            count++;
                    }
                }
            }
            P_v_cn[i][j] = (float)count / (t - 1) / N_MG;
            outfile_1 << i << " " << j << " " << P_v_cn[i][j] << endl;
        }
        outfile_1 << endl;
    }
    outfile_1.close();

    // Calculate Cl ion velocity over time
    for (i = 1; i <= t - 1; i++)
    {
        for (j = N_PEP + N_MG + 1; j <= N_PEP + N_MG + N_CL; j++)
        {
            v_CL[i][j - N_PEP - N_MG] = sqrt((atom_info[i + 1][j].x - atom_info[i][j].x) * (atom_info[i + 1][j].x - atom_info[i][j].x) + (atom_info[i + 1][j].y - atom_info[i][j].y) * (atom_info[i + 1][j].y - atom_info[i][j].y) + (atom_info[i + 1][j].z - atom_info[i][j].z) * (atom_info[i + 1][j].z - atom_info[i][j].z));
        }
    }

    // Calculate water coordination number of Cl ions
    for (i = 1; i <= t; i++)
    {
        for (j = N_PEP + N_MG + 1; j <= N_PEP + N_MG + N_CL; j++)
        {
            count = 0;
            for (k = N_PEP + N_MG + N_CL + 1; k <= N_tot; k++)
            {
                dist_O_WAT = sqrt((atom_info[i][k].x - atom_info[i][j].x) * (atom_info[i][k].x - atom_info[i][j].x) + (atom_info[i][k].y - atom_info[i][j].y) * (atom_info[i][k].y - atom_info[i][j].y) + (atom_info[i][k].z - atom_info[i][j].z) * (atom_info[i][k].z - atom_info[i][j].z));
                if (dist_O_WAT < r_CL)
                    count++;
            }
            cn_CL[i][j - N_PEP - N_MG] = count;
        }
    }
    
    // Calculate 2D distribution of Cl ion velocity vs water coordination number
    for (i = 0; i <= 50; i++)
    {
        for (j = 0; j <= 10; j++)
        {
            count = 0;
            for (k = 1; k <= t - 1; k++)
            {
                for (l = 1; l <= N_CL; l++)
                {
                    if (v_CL[k][l] >= i && v_CL[k][l] < i + 1)
                    {
                        if (cn_CL[k][l] >= (float)j && cn_CL[k][l] < (float)(j + 1))
                            count++;
                    }
                }
            }
            P_v_cn_CL[i][j] = (float)count / (t - 1) / N_CL;
            outfile_2 << i << " " << j << " " << P_v_cn_CL[i][j] << endl;
        }
        outfile_2 << endl;
    }
    outfile_2.close();

    return 0;
}