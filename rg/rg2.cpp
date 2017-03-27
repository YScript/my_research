/*
	this produce aim to caculate the rg2
*/
#include <iostream>
#include <fstream>

using namespace std;

int ndx = 3;
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
	int len_ghp = 11;
	int number_bcp;
	int number_ghp = lx*ly/(ndx*ndx);


	ifstream fin;
	fin.open("chains.txt");

/*	output section*/
	ofstream fout;
	fout.open("rg2.txt");
	
	return 0;
}

