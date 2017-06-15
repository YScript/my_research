
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

int ndx = 3;
int number_bcp = 400;
int len_fa = 11;
int len_fb = 9;

int distance_r1_r2(int ,int ,int ,int ,int ,int ,int);

struct coor
{
	int x;
	int y;
	int z;
	int type;
};

int main(int argc, char const *argv[])
{	
	int zc = 11;
	int lx = 60,ly = 60,lz = 60;
	int lxyz = lx* ly *lz ;
	int len_ga = 11;
	int len_ft = len_fa +len_fb;
	int number_ghp = lx*ly/(ndx*ndx);	
	int number_total_chains = number_bcp + number_ghp;
	int nchain[number_total_chains][len_ft];
	for (int i = 0; i < number_total_chains; ++i)
    	{	
    		for (int j = 0; j < len_ft; ++j)
    		{
    			nchain[i][j] = 0;
    		}
    	}

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

	ifstream fin;
	fin.open("chains.txt");
		/* doing  with all atoms is not optimum, so we do it with chains infomation*/
	int id_of_chain,len;
	int id_of_atom;
	int nb = 0;
	
	int id;

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
				nchain[i][j] = id_of_atom;
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
				}
				else
				{
					atoms[id_of_atom].type = 2;
				}
				nchain[i][k] = id_of_atom;	
			}
		}

	}
	fin.close();
//	cout <<"read the data of particle's coor and type"<<endl;
	int natom_total_number = number_ghp*len_ga+number_bcp*(len_fa+len_fb);
/*	deal with data*/

	int xk,yk,zk,xj,yj,zj;
	double array_rg2_ghp[number_ghp];
	for (int i = 0; i < number_ghp; ++i)
	{
		array_rg2_ghp[i] = 0.0;
	}
	double array_rg2_fa[number_bcp];
	for (int i = 0; i < number_bcp; ++i)
	{
		array_rg2_fa[i] = 0.0;
	}
	double rg2_ghp = 0.0;		//	mean square radius of gyraitom of grafted homopolymer;
	double rg2_fa = 0.0;		// mean square radius of gyraitom of A segment of free copolymer;
	double sum_rr2_ghp,sum_rr2_fa ;
	double sum_ree2_ghp,sum_ree2_fa;
	int rr2;
	int id_chains = 0;
	int nct = 0;
	double array_ree2_ghp[number_ghp];
	double array_ree2_bcp[number_bcp];
	for (int i = 0; i < number_ghp; ++i)
	{
		array_ree2_ghp[i] = 0.0;
	}
	for (int i = 0; i < number_bcp; ++i)
	{
		array_ree2_bcp[i] = 0.0;
	}
	double ree2_ghp = 0.0;
	double ree2_bcp = 0.0;
/*	caculate the mean square radius of gyration of A monomer belong BCP*/
	for (int i = 0; i < number_ghp; ++i)
	{	
		sum_rr2_ghp = 0.0;
		for (int j = 0; j < len_ga; ++j)
		{
			for (int k = 0; k < len_ga; ++k)
			{	
				xj =atoms[nchain[i][j]].x;
				yj =atoms[nchain[i][j]].y;
				zj =atoms[nchain[i][j]].z;

				xk =atoms[nchain[i][k]].x;
				yk =atoms[nchain[i][k]].y;
				zk =atoms[nchain[i][k]].z;

				sum_rr2_ghp += distance_r1_r2(xj,yj,zj,xk,yk,zk,lx);
			}
		}

		array_rg2_ghp[i] =sum_rr2_ghp/(double)(len_ga*len_ga);
	}

	cout <<"caculate the mean square radius of gyration of A monomer belong BCP"<<endl;

/*	caculate the mean square radius of gyration of A monomer belong BCP*/
	id_of_chain = 0;
	for (int i = 0; i < number_bcp; ++i)
	{	
		sum_rr2_fa = 0.0;
		id_of_chain =i+ number_ghp;

		for (int j = 0; j < len_fa; ++j)
		{	

			xj =atoms[nchain[id_of_chain][j]].x;
			yj =atoms[nchain[id_of_chain][j]].y;
			zj =atoms[nchain[id_of_chain][j]].z;

			for (int k = 0; k < len_fa; ++k)
			{	
				xk =atoms[nchain[id_of_chain][k]].x;
				yk =atoms[nchain[id_of_chain][k]].y;
				zk =atoms[nchain[id_of_chain][k]].z;

				sum_rr2_fa += distance_r1_r2(xj,yj,zj,xk,yk,zk,lx);
			}
		}
		array_rg2_fa[i] =sum_rr2_fa /(double)(len_fa*len_fa);
	}

	cout <<"caculate the mean square radius of gyration of A monomer belong BCP"<<endl;
/*	ouput section*/
	ofstream fout;
	fout.open("rg2_ghp.txt");
	for (int i = 0; i < number_ghp; ++i)
	{
		rg2_ghp += array_rg2_ghp[i]; 
	}
	rg2_ghp /=(double)number_ghp; 
	fout <<"rg2_ghp"<<"\t"<<rg2_ghp<<endl;
	for (int i = 0; i < number_ghp; ++i)
	{
		fout <<i<<"\t"<<array_rg2_ghp[i]<<endl;
	}
	fout.close();
	
	ofstream fp;
	fp.open("rg2_bcp.txt");
	for (int i = 0; i < number_bcp; ++i)
	{
		rg2_fa += array_rg2_fa[i]; 
	}
	rg2_fa /=(double)number_bcp; 
	fp <<"rg2_bcp"<<"\t"<<rg2_fa<<endl;
	for (int i = 0; i < number_bcp; ++i)
	{
		fp <<i<<"\t"<<array_rg2_fa[i]<<endl;
	}
	fp.close();

cout <<"rg2_ghp\t"<<rg2_ghp<<endl;
cout <<"rg2_bcp\t"<<rg2_fa <<endl;
	/*	relese the memory*/	
	delete  []atoms;
	atoms = NULL;
	for (int i = 0; i < lxyz; ++i)
	{
		delete []id_near_point[i];
	}
	delete	[]id_near_point;

	cout <<"finished"<<endl;
	system("pause");
	return 0;
}

int distance_r1_r2(int x1,int y1,int z1,int x2,int y2,int z2,int l){
			int x = x1- x2;
			int y = y1- y2;
			int z = z1- z2;
			int rr2;
			if(x > l/2)
			{	
				x = x -l;
			}
			else if(x < -l/2)
			{
				x = x + l;
			}
			if(y >l/2)
			{
				y = y - l;
			}
			else if(y < -l/2)
			{
				y = y +l;		
			}
			if(z > l/2)
			{
				z = z -l;
			}
			else if(z < -l/2)
			{
				z = z +l;
			}
			rr2 = x*x + y*y+ z*z;
	return rr2;
}

