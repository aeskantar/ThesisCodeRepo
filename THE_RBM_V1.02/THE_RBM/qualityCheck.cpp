#include "qualityCheck.h"

qualityCheck::qualityCheck() {};
qualityCheck::qualityCheck(double** coorp_,vector<int> nu_,int np_) 
	: coorp(coorp_), nu(nu_), np(np_){}

bool qualityCheck::invertedElements(double** coorp) {
	bool invElnts = false;
	int is1, is2, is3, ind = 0;
	double dx12, dx13, dy12, dy13, criterion;
	double dxa, dya, dza = 0.;
	int kountInv = 0;  // number of inverted triangles
	ind = 0;
	for (int ip = 0; ip < np; ip++)
	{
		is1 = nu[ind + 1] - 1; // exw valei -1
		is2 = nu[ind + 2] - 1;// exw valei -1
		is3 = nu[ind + 3] - 1;// exw valei -1
		dx12 = coorp[0][is2] - coorp[0][is1];
		dx13 = coorp[0][is3] - coorp[0][is1];
		dy12 = coorp[1][is2] - coorp[1][is1];
		dy13 = coorp[1][is3] - coorp[1][is1];
		vecProd(dx12, dy12, 0., dx13, dy13, 0., dxa, dya, dza);
		criterion = dza * .5;
		if (criterion < 0.){ kountInv++; }
		ind = ind + 3;
	}
	if (kountInv > 0) {
		cout << "Inverted triangles : "<< kountInv << endl;
		invElnts = true;
	}
	return invElnts;
}

void qualityCheck::vecProd(double v1x, double v1y, double v1z, double v2x,
	double v2y, double v2z, double& xre, double& yre, double& zre) {

	xre = v1y * v2z - v1z * v2y;
	yre = v1z * v2x - v1x * v2z;
	zre = v1x * v2y - v1y * v2x;

}

void qualityCheck::meshQuality(double** coor, vector<int> &nu)
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
		nod[1] = nu[ind + 1];
		nod[2] = nu[ind + 2];
		is1 = nod[0]; is2 = nod[1]; is3 = nod[2];
		ind += 3;

		qT[0] = 100.; qT[1] = 100.; qT[2] = 100.; // set a high (bad) value
//
//		the three nodes of the triangle
		for (unsigned int j = 0; j < 3; j++)
		{
			if (j == 1) { is1 = nod[0]; is2 = nod[1]; is3 = nod[2]; }
			if (j == 2) { is1 = nod[0]; is2 = nod[1]; is3 = nod[2]; }

			a11 = coor[0][is2-1] - coor[0][is1-1];
			a12 = coor[0][is3-1] - coor[0][is1-1];
			a21 = coor[1][is2-1] - coor[1][is1-1];
			a22 = coor[1][is3-1] - coor[1][is1-1];

			s11 = a11 * wr11 + a12 * wr21;
			s12 = a11 * wr12 + a12 * wr22;
			s21 = a21 * wr11 + a22 * wr21;
			s22 = a21 * wr12 + a22 * wr22;

			aja = s11 * s22 - s12 * s21;
			trace = s11 * s11 + s12 * s12 + s21 * s21 + s22 * s22;

			if (aja < 0)
			{
				kountInv++;
				exit;
			}
			qT[j] = aja / trace;
		}

		qual[ip] = 2. * (qT[0] + qT[1] + qT[2]) / 3.;
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
	double a, b;
	a = qual[0];
	b = qual[np - 1];
	qMean = qMean / np;

	qSTD = 0.;
	for (unsigned int ip = 0; ip < np; ip++)
	{
		qSTD = qSTD + (qual[ip] - qMean) * (qual[ip] - qMean);
	}
	qSTD = sqrt(qSTD / np);

	cout << " Quality Metrics  " << endl;
	cout << " ---------------- " << endl;
	cout << " Min: " << qMin << "	Max: " << qMax << endl;
	cout << " Mean: " << qMean << "	std: " << qSTD << endl;
	cout << "  ------  ------  ------  ------  ------ " << endl;
}


qualityCheck::	~qualityCheck() {};