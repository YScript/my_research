
#include <iostream>
#include <fstream>

using namespace std;

/*
	this is  a little produce to handle the problem about the contant numbers
between A /B monomer in BCP and solvent;
*/
/*	
	by___fengyuan /2017/3/20
*/


struct coor
{
	int x;
	int y;
	int z;
	int type;
};
int main(int argc, char const *argv[])
{
	int lx = 60,ly = 60,lz = 60;
	int lxyz = lx* ly *lz ;
	double contact_Number = 0.0;
	double cntnumb_of_AandB,cntnumb_of_AandS,cntnumb_of_BandS;

/*	init the coordinates system;*/
	coor *atoms = new coor [lxyz];
	for (int i = 0; i < lxyz; ++i)
	{
		atoms[i].type = 3;
	}
	int count = 0;
	for (int i = 0; i < lx; ++i)
	{
		for (int j = 0; j < ly; ++j)
		{
			for (int k = 0; k < lz; ++k)
			{
				atoms[count].x = i;
				atoms[count].y = j;
				atoms[count].z = k;
				count++;
			}
		}
	}
	for (int i = 0; i < lxyz; ++i)
	{
		if (atoms[i].z == 0)
		{
			atoms[i].type = 5;  /* the bottom flat type is 5*/
		}
		else if (atoms[i].z == lx-1)
		{
			atoms[i].type == 6; /*the top flat type is 6*/
		}
	}
	cout <<"after init coordinates"<<endl;
/*	build the near points array;*/
	int const number_of_nearpoints = 18;
	double const len_of_nearpoints = 2.0;
	coor *near_coor = new coor[number_of_nearpoints];
	count = 0;
	int xi,yi,zi,rr;
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			for (int k = 0; k < 3; ++k)
			{	
				xi = near_coor[count].x;
				yi = near_coor[count].y;
				zi = near_coor[count].z;
				rr = xi *xi +yi *yi + zi *zi ;
			//	rr = near_coor[count].x * near_coor[count].x + near_coor[count].y*near_coor[count].y+near_coor[count].z*near_coor[count].z;
				if (rr <= len_of_nearpoints)			
				{
					near_coor[count].x = i-1;
					near_coor[count].y = j-1;
					near_coor[count].z = k-1;					
					count ++;
				}
				
			}
		}
	}
	int number_np = count;
	cout <<"number_of_nearpoints"<<"\t"<<number_np<<endl;
	int id;
	int **id_near_point;
	id_near_point = new int *[lxyz];
	for (int i = 0; i < lxyz; ++i)
	{
		id_near_point[i] = new int [number_np];
	}
	for (int i = 0; i < lxyz; ++i)
	{
		for (int j = 0; j < number_np; ++j)
		{
			xi = atoms[i].x + near_coor[j].x;
			yi = atoms[i].y + near_coor[j].y;
			zi = atoms[i].z + near_coor[j].z;
			id = xi*ly*lz + yi*lz + zi;
			id_near_point[i][j] = id;		
		}
	}

/*	read the data of particle's coor and type*/
	ifstream fin;
	fin.open("chains.txt");
		/* doing  with all atoms is not optimum, so we do it with chains infomation*/
	fin.close();

	cntnumb_of_AandB = 0.0;
	cntnumb_of_AandS = 0.0;
	cntnumb_of_BandS = 0.0;
	ofstream fout;
	fout.open("contact_Number.txt");
	fout <<"particles_type"<<"\t"<<contact_Number<<endl;

	fout << "cntnumb_of_AandB"<<"\t"<<cntnumb_of_AandB<<endl;
	fout << "cntnumb_of_AandS"<<"\t"<<cntnumb_of_AandS<<endl;
	fout << "cntnumb_of_BandS"<<"\t"<<cntnumb_of_BandS<<endl;
	fout.close();
	delete  []atoms;
	atoms = NULL;
	delete  []near_coor;
	near_coor = NULL;
	for (int i = 0; i < lxyz; ++i)
	{
		delete []id_near_point[i];
	}
	delete	[]id_near_point;
/*
	system("pause");*/
	return 0;
}