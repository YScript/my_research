#include <iostream>
#include <fstream>

using namespace std;

struct coor
{
	int x;
	int y;
	int z;
	int type;
};

int ndx = 4;

int len_free_a = 6;
int len_free_b = 15;
int num_free_chains = 600;
int len_grafted_hp = 11;

int main(int argc, char const *argv[])
{
	int lx = 60,ly = 60,lz = 60;
	int lxyz = lx*ly*lz;
	int num_NearPoint = 18;
	int num_grafted_chains = (lx*ly)/(ndx^2);
	int len_free_chains = len_free_a + len_free_b;
	int num_total_chains = num_grafted_chains + num_free_chains;
	
	int nchain[num_total_chains][len_free_chains];
	for (int i = 0; i < num_total_chains; ++i)
	{
		for (int j = 0; j < len_free_chains; ++j)
		{
			nchain[i][j] = 0;
		}
	}
	coor *atom = new coor[lxyz];
	double *densA_z = new double[lz];
	double *densB_z = new double[lz];
	for (int i = 0; i < lz; ++i)
	{
		densA_z[i] = 0.0;
		densB_z[i] = 0.0;
	}
	cout <<num_grafted_chains<<endl;
	int nct = 0;
	for (int i = 0; i < lx; ++i)
	{
		for (int j = 0; j < ly; ++j)
		{
			for (int k = 0; k < lz; ++k)
			{
				atom[nct].x = i;
				atom[nct].y = j;
				atom[nct].z = k;
				nct++;
			}
		}
	}
	for (int i = 0; i < lxyz; ++i)
	{
		atom[i].type = 3;
	}
	int natom = nct;

	ifstream fin;
	fin.open("chains.txt");
/*	cout <<"fin open"<<endl;*/
	nct = 0;
	int ichain,id;
	for (int i = 0; i < num_total_chains; ++i)
	{
		fin >>ichain;
		if (ichain <num_grafted_chains)
		{
			for (int j = 0; j < len_grafted_hp; ++j)
			{
				fin >>id;
				id--;
				atom[id].type = 0;
			}	
		}
		else
		{
			for (int k = 0; k < len_free_chains; ++k)
			{	
				fin >>id;
				id--;
				if (k <len_free_a)
				{	
					atom[id].type = 1;
				}
				else
				{
					atom[id].type = 2;
				}
			}
		}
	}
	fin.close();
/*	cout <<"fin close"<<endl;*/
	ofstream fout;
	fout.open("vertical_density.txt");
	fout <<"lz"<<"\t"<<"density_A_zdir"<<"\t"<<"density_B_zdir"<<endl;
	nct = 0;
	int zc;
	int num_A_monomer,num_B_monomer;

	for (int j = 0; j < lz; ++j)
	{
		for (int i = 0; i < natom; ++i)
		{

			if (atom[i].z == j )
			{
				if (atom[i].type == 1)
				{
					densA_z[j]++;
				}
				else if (atom[i].type ==2)
				{
					densB_z[j]++;
				}
			}
			
			
		}
		/*densA_z[j] = (double)densA_z[j]/(double)(lx*ly);
		densB_z[j] = (double)densB_z[j]/(double)(lx*ly);*/
	}
	for (int i = 0; i < lz; ++i)
	{
		fout <<i+1<<"\t"<<densA_z[i]<<"\t"<<densB_z[i]<<endl;
	}
	fout.close();

	delete []densA_z;
	delete []densB_z;
	delete []atom;


	return 0;

}


/*
(-40*4^(-(0.02*(-90+0.1*(i-1)))^2)+(-20+0.1*(j-1)))*(10*4^(-(0.01*(-90+0.1*(i-1)))^2))
*/