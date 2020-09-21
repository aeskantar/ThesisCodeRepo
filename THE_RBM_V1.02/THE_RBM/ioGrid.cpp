#include "ioGrid.h"

/*=====================================================
  =====================================================
  Created by Alexandros Eskantar and
			 Ahmet-Triantafilos Naidi
  Master students of Computational Mechanics,
  National Technichal University of Athens. All 
  rights reserved to Professor Kyriakos C. Giannakoglou
  and the Lab. of  THERMAL TURBOMACHINES PARALLEL CFD 
						& OPTIMIZATION UNIT
  ======================================================
  ======================================================
*/

ioGrids::ioGrids() {}

void ioGrids::readDeformed2D_Unstructured(string defFile)
{
	int dum;

	std::ifstream def;
	//Check if .ele file exists
	if (def.fail()) {

		cout << "Non-existent def file" << endl;
		exit;
	}

	gInfo_un.coorp = new double*[2];
	for (unsigned int i = 0; i < 2; ++i)
		gInfo_un.coorp[i] = new double[gInfo_un.ns];
	//Read.def file
	def.open(defFile);
	def >> gInfo_un.ns2;
	if (gInfo_un.ns != gInfo_un.ns2) 
	{ cout << "Incompatible derfomation file" << endl; }
	for (int i = 0; i < gInfo_un.ns; i++)
		def >> dum;
	for (int i = 0; i < gInfo_un.ns; i++){ def >> gInfo_un.coorp[0][i];}
		
	for (int i = 0; i < gInfo_un.ns; i++){ def >> gInfo_un.coorp[1][i];}
		
	def.close();

}

void ioGrids::readInitial2D_Unstructured(string eleFile, string nodFile)
{

	std::ifstream ele;
	std::ifstream nod;
	//Check if .ele file exists
	if (ele.fail()) {
		cout << "Non-existent ele file" << endl;
		exit;
	}
	//Read .ele file
	ele.open(eleFile);
	ele >> gInfo_un.np;
	//ele >> gInfo_un.nq;
	gInfo_un.nq = 0;
	gInfo_un.nall = gInfo_un.np + gInfo_un.nq; // total number of elements
	gInfo_un.nentries = gInfo_un.np * 3 + gInfo_un.nq * 4; // entries to NU
	gInfo_un.nu.resize(gInfo_un.nentries);
	//Read entries
	if (gInfo_un.np > 0)
	{
		for (unsigned int i = 0; i < gInfo_un.np * 3; i++) ele >> gInfo_un.nu[i];
	}
	gInfo_un.nprev = gInfo_un.np * 3;
	if (gInfo_un.nq > 0)
	{
		for (unsigned int i = 0; i < gInfo_un.nq * 4; i++) ele >> gInfo_un.nu[gInfo_un.nprev + i];
	}
	ele.close();
	
	// Check if.nod file exists
	if (nod.fail())
	{
		cout << "Non-existent nod file" << endl;
		exit;
	}
	//Read .nod file
	nod.open(nodFile);
	nod >> gInfo_un.ns;
	gInfo_un.logfr.resize(gInfo_un.ns);
	for (unsigned int i = 0; i < gInfo_un.ns; i++) nod >> gInfo_un.logfr[i];

	gInfo_un.coor = new double*[2];
	for (unsigned int i = 0; i < 2; ++i)
		gInfo_un.coor[i] = new double[gInfo_un.ns];
	
	for (unsigned int i = 0; i < gInfo_un.ns; ++i) nod >> gInfo_un.coor[0][i];
	for (unsigned int i = 0; i < gInfo_un.ns; ++i) nod >> gInfo_un.coor[1][i];
}

void ioGrids::readInitial3D_Unstructured(string eleFile,string hybFile, string nodFile)
{
	std::ifstream ele;
	std::ifstream nod;
	std::ifstream hyb;

	hyb.open(hybFile);
	ele.open(eleFile);

	//Check if .hyb file exists

	if (hyb.is_open())
	{
		hyb >> gInfo3D_un.ntet>> gInfo3D_un.npyr>> gInfo3D_un.npri>> gInfo3D_un.nhex;
	}
	else //Check if .ele file exists
	{
		if (ele.fail()) { std::cout << "Non-existent ele, hyb files" << endl; exit; }
		ele >> gInfo3D_un.ntet;
		gInfo3D_un.npyr = 0;
		gInfo3D_un.npri = 0;
		gInfo3D_un.nhex = 0;
	}
	//Read .ele file

	gInfo3D_un.nall = gInfo3D_un.ntet + gInfo3D_un.npyr + gInfo3D_un.npri + gInfo3D_un.nhex; 
	gInfo3D_un.nentries = gInfo3D_un.ntet * 4 + gInfo3D_un.npyr * 5 + gInfo3D_un.npri * 6 + gInfo3D_un.nhex * 8; // entries to NU array
	gInfo3D_un.nentrsg = gInfo3D_un.ntet * 6 + gInfo3D_un.npyr * 8 + gInfo3D_un.npri * 9 + gInfo3D_un.nhex * 12; // ! entries of nusg
	gInfo3D_un.nu.resize(gInfo3D_un.nentries);

	//Read entries
	if (gInfo3D_un.ntet > 0)
	{
		for (unsigned int i = 0; i < gInfo3D_un.ntet * 4; i++) ele >> gInfo3D_un.nu[i];
	}
	int nprev = gInfo3D_un.ntet * 4;
	if (gInfo3D_un.npyr > 0)
	{
		for (unsigned int i = 0; i < gInfo3D_un.npyr * 5; i++) ele >> gInfo3D_un.nu[nprev + i];
	}
	nprev += gInfo3D_un.npyr * 5;
	if (gInfo3D_un.npri > 0)
	{
		for (unsigned int i = 0; i < gInfo3D_un.npri * 6; i++) ele >> gInfo3D_un.nu[nprev + i];
	}
	nprev += gInfo3D_un.npri * 6;
	if (gInfo3D_un.nhex > 0)
	{
		for (unsigned int i = 0; i < gInfo3D_un.nhex * 8; i++) ele >> gInfo3D_un.nu[nprev + i];
	}
	ele.close();

	// Check if.nod file exists
	if (nod.fail())
	{
		std::cout << "Non-existent nod file" << endl;
		exit;
	}
	//Read .nod file
	nod.open(nodFile);
	nod >> gInfo3D_un.ns;
	gInfo3D_un.logfr.resize(gInfo3D_un.ns);
	for (unsigned int i = 0; i < gInfo3D_un.ns; i++) nod >> gInfo3D_un.logfr[i];

	gInfo3D_un.coor = new double* [3];
	for (unsigned int i = 0; i < 3; ++i)
		gInfo3D_un.coor[i] = new double[gInfo3D_un.ns];

	for (unsigned int i = 0; i < gInfo3D_un.ns; ++i) nod >> gInfo3D_un.coor[0][i];
	for (unsigned int i = 0; i < gInfo3D_un.ns; ++i) nod >> gInfo3D_un.coor[1][i];
	for (unsigned int i = 0; i < gInfo3D_un.ns; ++i) nod >> gInfo3D_un.coor[2][i];
	nod.close();
}	

void ioGrids::write2D_Unstructured(double **coorp, double **coor, vector<int> &logfr, int ns, string gridName)
{
	std::ofstream nod;
	nod.open(gridName);
	nod << ns << endl;
	for (int i = 0; i < ns; i++){ nod << logfr[i] << " "; }
	nod << endl;
	for (int i = 0; i < ns; i++){ nod << coorp[0][i] << " "; }
	nod << endl;
	for (int i = 0; i < ns; i++){ nod << coorp[1][i] << " "; }
    
	//for (int i = 0; i < ns; i++)
	//{ 
	//	nod << coorp[0][i] << "	" << coorp[1][i] <<endl;
	//}
	//nod.close();
	
}

void ioGrids::meshQuality(double** coor, int np, vector<int> &nu)
{
	double qMin, qMax, qMean, qLoc, qSTD;
	double riza3 = sqrt(3.0);
	double wr11 = 1.;
	double wr12 = -1. / riza3;
	double wr21 = 0.;
	double wr22 = 2. / riza3;
	double a11, a12, a21, a22, s11, s12, s21, s22;
	double aja, trace;
	int is1, is2, is3;
	vector<int> nod; nod.resize(3);
	vector <double> qT, qual; qT.resize(3); qual.resize(np);
	int ind = 0;
	int kountInv = 0; // number of inverted triangles

	for (unsigned int ip = 0; ip < np; ip++)
	{
		nod[0] = nu[ind];
		nod[1] = nu[ind];
		nod[2] = nu[ind];
		is1 = nod[0]; is2 = nod[1]; is3 = nod[2];
		ind += 3;

		qT[0] = 10.; qT[1] = 10.; qT[2] = 10.; // set a high (bad) value
//
//		the three nodes of the triangle
		for (unsigned int j = 0; j < 3; j++)
		{
			if (j == 2){ is1 = nod[0]; is1 = nod[1]; is1 = nod[2]; }
			if (j == 3){ is1 = nod[0]; is1 = nod[1]; is1 = nod[2]; }

			a11 = coor[0][is2] - coor[0][is1];
			a12 = coor[0][is3] - coor[0][is1];
			a21 = coor[1][is2] - coor[1][is1];
			a22 = coor[1][is3] - coor[1][is1];

			s11 = a11*wr11 + a12*wr21;
			s12 = a11*wr12 + a12*wr22;
			s21 = a21*wr11 + a22*wr21;
			s22 = a21*wr12 + a22*wr22;

			aja = s11*s22 - s12*s21;
			trace = s11*s11 + s12*s12 + s21*s21 + s22*s22;

			if (aja < .0)
			{
				kountInv++;
				exit;
			}
			qT[j] = aja / trace;
		}

		qual[ip] = 2.*(qT[0] + qT[1] + qT[2]) / 3.;
	}

	qMin = 100.;
	qMax = -100.;
	qMean = 0.;
	for (unsigned int ip = 0; ip < np; ip++)
	{
		qLoc = qual[ip];
		qMean = qMean + qLoc;
		if (qMin > qLoc) qMin = qLoc;
		if (qMax < qLoc) qMax = qLoc;
	}
	qMean = qMean / np;

	qSTD = 0.;
	for (unsigned int ip = 0; ip < np; ip++)
	{
		qSTD = qSTD + (qual[ip] - qMean)*(qual[ip] - qMean);
	}
	qSTD = sqrt(qSTD/np);

	cout << " Quality Metrics  " << endl;
	cout << " ---------------- " << endl;
	cout << " Min: " << qMin << "	Max: " << qMax << endl;
	cout << " Mean: " << qMean << "	std: " << qSTD << endl;
	cout << "  ------  ------  ------  ------  ------ " << endl;
}

void ioGrids::getGeometry() {

	ofstream geometry;
	geometry.open("geom2D.dat");

	for (int i = 0; i < gInfo_un.ns; i++)
	{
		if (gInfo_un.logfr[i]==3)
		{
			geometry << gInfo_un.coor[0][i] << "	" << gInfo_un.coor[1][i] << endl;
		}
	}
}

ioGrids::~ioGrids() {}

