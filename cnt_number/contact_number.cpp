
/*
	this is  a little produce to handle the problem about the contant numbers
between A /B monomer in BCP and solvent;
*/
/*	
	by___fengyuan /2017/3/20
*/

#include <iostream>
#include <fstream>

using namespace std;

int ndx = 4;
int number_bcp = 600;
int len_fa = 11;
int len_fb = 9;

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
	int len_ga = 11;
	int number_ghp = lx*ly/(ndx*ndx);
	int number_total_chains = number_ghp + number_bcp;
	double contact_Number = 0.0;
	double cntnumb_of_AandB,cntnumb_of_AandS,cntnumb_of_BandS;

/*	init the coordinates system;*/
	coor *atoms = new coor [lxyz];
	for (int i = 0; i < lxyz; ++i)
	{
		atoms[i].type = 3;
	}
	int **id_near_point= NULL;
        id_near_point = new int*[lxyz];
 	int np = 18;
	for (int i = 0; i < lxyz; ++i)
        {
               id_near_point[i] = new int [np];
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
			atoms[i].type = 6; /*the top flat type is 6*/
		}
		else
		{
			atoms[i].type =3;
		}
	}
/*	test init coor and type infomation */
	/*ofstream fo1;
	fo1.open("fo1");
	for (int i = 0; i < lxyz; ++i)
	{
		if (atoms[i].type ==5)		
		{
			fo1 <<i<<"\t"<<atoms[i].x<<"\t"<<atoms[i].y<<"\t"<<atoms[i].type<<endl;
		}
	}
	fo1.close();*/


/*	build the near points coor array;*/
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
				xi = i-1;
				yi = j-1;
				zi = k-1;
				rr = xi *xi +yi *yi + zi *zi ;
			//	rr = near_coor[count].x * near_coor[count].x + near_coor[count].y*near_coor[count].y+near_coor[count].z*near_coor[count].z;
				if (rr <= len_of_nearpoints && rr >0)			
				{
					near_coor[count].x = xi;
					near_coor[count].y = yi;
					near_coor[count].z = zi;
					near_coor[count].type = 0;			
					count ++;
				}
				
			}
		}
	}
	/*for (int i = 0; i < count; ++i)
	{
		cout <<i<<"\t"<<near_coor[i].x<<"\t"<<near_coor[i].y<<"\t"<<near_coor[i].z<<endl;
	}*/
	int number_np = count;
//	cout <<"number_of_nearpoints"<<"\t"<<np<<"\t"<<number_np<<endl;
	int id;

/*	build the id of near points*/	
	for (int i = 0; i < lxyz; ++i)
	{
		for (int j = 0; j < number_np; ++j)
		{
			xi = atoms[i].x + near_coor[j].x;
			yi = atoms[i].y + near_coor[j].y;
			zi = atoms[i].z + near_coor[j].z;
			if(xi > (lx-1))
			{	
				xi = xi -lx;
			}
			else if(xi < 0)
			{
				xi = xi + lx;
			}
			if(yi >( ly-1))
			{
				yi = yi - ly;
			}
			else if(yi < 0)
			{
				yi = yi +ly;		
			}
			if(zi > (lz-1))
			{
				zi = zi -lz;
			}
			else if(zi < 0)
			{
				zi = zi +lz;
			}

			id = xi*ly*lz + yi*lz + zi;
			id_near_point[i][j] = id;		
		}
	}

/*	read the data of particle's coor and type*/
	ifstream fin;
	fin.open("chains.txt");
		/* doing  with all atoms is not optimum, so we do it with chains infomation*/
	int id_of_chain,len;
	int id_of_atom;
	int nb = 0;
	for(int i= 0; i<number_total_chains; ++i)
	{
		fin >> id_of_chain;
		if(i <number_ghp)
		{
			len = len_ga;
			for(int j = 0;j <len;++j)
			{
				fin >>id_of_atom;
				id_of_atom--;
				atoms[id_of_atom].type = 1;
			}	
		}
		else
		{
			len = len_fa + len_fb;
			for(int k = 0;k<len;++k)
			{
				fin >>id_of_atom;
				id_of_atom--;
				if(k <len_fa)
				{
					atoms[id_of_atom].type = 1;
//					cout <<id_of_atom<<"type A"<<endl;
				}
				else
				{
					atoms[id_of_atom].type = 2;
//					cout <<id_of_atom<<"type A"<<endl;
					nb++;
				}		
			}
		}
	}
	fin.close();
	int natom_total_number = number_ghp*len_ga+number_bcp*(len_fa+len_fb);
/*	cout <<natom_total_number<<endl;*/

/*	ouput section*/
	cntnumb_of_AandB = 0.0;
	cntnumb_of_AandS = 0.0;
	cntnumb_of_BandS = 0.0;
	double numaa,numbb,numas,numaw,numab,numbs,numbw;
	numas = 0.0;
	numaa = 0.0;
	numab = 0.0;
	numaw = 0.0;
	numbb = 0.0;
	numbs = 0.0;
	numbw = 0.0;
	ofstream fout;
	fout.open("contact_Number.txt");
	fout << "ctnber_AA\tctnber_AB\tctnber_AS\tctnber_AW\tctnber_BB\tctnber_BS\tctnber_BW"<<endl;
	
	int nct = 0;
	for(int i= 0;i <lxyz; i++)
	{	
		if(atoms[i].type == 1)
		{
		//	fout <<nct<<"\t"<<atoms[i].x<<"\t"<<atoms[i].y<<"\t"<<atoms[i].z<<"\t"<<atoms[i].type<<endl;
			for (int j = 0; j <number_np ; ++j)
			{	
				if (atoms[id_near_point[i][j]].type == 1)
				{
					numaa += 1.0;
				}
				else if (atoms[id_near_point[i][j]].type == 2)
				{
					numab += 1.0;
				}
				else if (atoms[id_near_point[i][j]].type == 3)
				{
					numas += 1.0;
				}
				else if (atoms[id_near_point[i][j]].type == 5)
				{
					numaw += 1.0;
				}
			}		
		}
		else if (atoms[i].type == 2)
		{
			for (int j = 0; j <number_np ; ++j)
			{
				if (atoms[id_near_point[i][j]].type == 2)
				{
					numbb += 1.0;
				}
				if (atoms[id_near_point[i][j]].type == 3)
				{
					numbs += 1.0;
				}
				else if (atoms[id_near_point[i][j]].type == 5)
				{
					numbw += 1.0;
				}
			}
		}
	}
	double numab1;
	int number_a = number_bcp*len_fa+number_ghp*len_ga;
	int number_b = number_bcp*len_fb;

	
	numaa /= (double)(number_a*number_np);	
	numab /= (double)(number_a*number_np);
	numas /= (double)(number_a*number_np);
	numaw /= (double)(number_a*number_np);
	numbb /= (double)(number_b*number_np);
	numbs /= (double)(number_b*number_np);
	numbw /= (double)(number_b*number_np);
	fout <<numaa<<"\t"<<numab<<"\t"<<numas<<"\t"<<numaw<<"\t"<<numbb<<"\t"<<numbs<<"\t"<<numbw<<endl;
	fout.close();
/*	relese the memory*/	
	delete  []atoms;
	atoms = NULL;
	delete  []near_coor;
	near_coor = NULL;
	for (int i = 0; i < lxyz; ++i)
	{
		delete []id_near_point[i];
	}
	delete	[]id_near_point;
	cout <<"finished"<<endl;
	return 0;
}
