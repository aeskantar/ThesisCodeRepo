#pragma once

#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

struct gridInfo2D_Structured {
	int Imax = 0, Jmax = 0, bountaryStart = 0, bountaryEnd = 0;
	int totalNodes;
	int** idMatrix = 0;
	vector<double> x0, y0;
	vector<double> xf, yf;
};

struct gridInfo2D_Unstructured {
	vector<int> nu, logfr;
	double** coor, **coorp;
	int np, nq, nall, nentries, nprev, ns, ns2;
};


class Structures2D
{
public:
	Structures2D();
	Structures2D(int ns, int np, int nq, vector <int> & nu);
	void Create();
	
	int get_ns(){ return ns; }
	int get_nseg() { return nseg; }
	int get_np() { return np; }
	int get_nq() { return nq; }
	int get_nbseg() { return nbseg; }
	vector <int> get_ndeg() { return ndeg; }
	vector <int> get_jaret() { return jaret; }
	~Structures2D();

private:
	int **nubo;
	int ns, nseg, np, nq, nbseg, nvseg;
	vector <int> ndeg, jaret, nu, nusg, ibsg2tr;

	int maxnod = 50000;
	int maxpyr = 4 * maxnod, nvmaxall = 700000, maxnu = 300000,
		maxseg = 4 * maxnod, maxbseg = 5000;

	void numsegs();
	void fjaret();
};
