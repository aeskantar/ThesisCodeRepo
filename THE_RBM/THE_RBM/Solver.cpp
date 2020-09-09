#include "Solver.h"
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

Numerics::Numerics(){}

Numerics::Numerics(double **coorp_, double **coor_, vector<int> &logfr_, vector <int> &ndeg_, vector <int> &jaret_, 
				   vector <int> &nu_, int ns_, int np_) :
coorp(coorp_), coor(coor_), logfr(logfr_), ndeg(ndeg_), jaret(jaret_), nu(nu_), ns(ns_), np(np_){}

void Numerics::Solver2D(int Iterations)
{
	std::vector<int>::iterator it;
	it = nu.begin();
	nu.insert(it, 0);
	qualityCheck elements(coorp, nu, np);
	bool testInt;
	ndim = 2;
	dxi = new double*[ndim];
	for (unsigned int i = 0; i < ndim; ++i)
		dxi[i] = new double[ns];
//
//	- Start outer iterations loop
//-------------------------------
	int kIterOut = 0;
//*********************
	while (kIterOut <= Iterations)
	{
		kIterOut++;
//	write(*, *)' External Iteration ', kIterOut
		for (int i = 0; i < ndim; i++)
		{
			for (int j = 0; j < ns; j++)
			{
				dxi[i][j] = coorp[i][j];
			}
		}
		
//
//	 Sweep all internal grid nodes and compute deformations & rotations
		for (int is = 1; is <= ns; is++)
		{
			if (logfr[is - 1] != 0){ continue; } // only internal nodes are of interest
			mID = is; // current mID node
//
//			Count neighbouring nodes
			neiTot = kountNeis(mID);
//			write(*,*)' Node, Number of neis =',is,neiTot
//		
//			Global to Local w.r.t. point M 
			globalToLocal(ndim, cM, mID);
//
//			Compute deformations & rotations
			theta = 0.;	// initialization
			iter = 0;
			dxOld = 20000.;
			dyOld = 30000.;
			thOld = 0;
			while (iter < Iterations)
			{
				iter++;
			//-> Update X, Y
				sx = 0.; cc = cos(theta);
				sy = 0.; ss = sin(theta);
				for (int k = ndeg[mID - 1] + 1; k <= ndeg[mID]; k++)
				{
					inei = jaret[k];
					sx = sx + coorp[0][inei - 1] - cc*coor[0][inei - 1] - ss*coor[1][inei - 1]; //exw valei -1 sto inei
					sy = sy + coorp[1][inei - 1] + ss*coor[0][inei - 1] - cc*coor[1][inei - 1]; //exw valei -1 sto inei
				}
				dx = sx / neiTot;
				dy = sy / neiTot;
			//-> Update theta
				aTerm = 0.;
				bTerm = 0.;
				for (int k = ndeg[mID - 1] + 1; k <= ndeg[mID]; k++)
				{
					inei = jaret[k];
					xnei = coor[0][inei - 1]; ynei = coor[1][inei - 1]; //exw valei -1 sto inei
					xneiP = coorp[0][inei - 1]; yneiP = coorp[1][inei - 1]; //exw valei -1 sto inei
					aTerm = aTerm - xnei*dx + xnei*xneiP - ynei*dy + ynei*yneiP;
					bTerm = bTerm + ynei*dx - ynei*xneiP - xnei*dy + xnei*yneiP;
				}
				fun = aTerm*ss + bTerm*cc;
				dfun = aTerm*cc - bTerm*ss;
				if (abs(dfun) < tiny){ break; }
				theta = theta - fun / dfun;
//
			//-> Check convergence
				if (converged(dxOld, dyOld, thOld, dx, dy, theta))
				{ 
					break; 
				}
				dxOld = dx;
				dyOld = dy;
				thOld = theta;	
//
				//if (mod(iter,10)==0 .or. iter==1) write(*,*) iter,dx,dy,theta
				//write(*,*) iter,dx,dy,theta
				if (iter > Iterations){ break; }
			}
//
			//Update coordinates of point M
			coorp[0][mID - 1] = cM[0] + dx;  //exw valei -1 sto mID
			coorp[1][mID - 1] = cM[1] + dy; //exw valei -1 sto mID
//
//			Local to Global w.r.t. point M
			localToGlobal(ndim, cM, mID);
		}
//		
		xrms = 0.; yrms = 0.;
		xErrMax = -1e10;
		yErrMax = -1e10;
		for (int is = 0; is < ns; is++)
		{
			xi = dxi[0][is] - coorp[0][is];
			yi = dxi[1][is] - coorp[1][is];
			xi2r = sqrt(xi*xi);
			yi2r = sqrt(yi*yi);
			if (xi2r > xErrMax) { xErrMax = xi2r; }
			if (yi2r > yErrMax) { yErrMax = yi2r; }
			xrms = xrms + xi*xi;
			yrms = yrms + yi*yi;
		}
		xrms = sqrt(xrms / ns);
		yrms = sqrt(yrms / ns);
		xrms = log10(xrms);
		yrms = log10(yrms);
		xErrMax = log10(xErrMax);
		yErrMax = log10(yErrMax);
		if (kIterOut <= 10. || kIterOut % 100 == 0.){ cout << kIterOut << "	" << xrms << "	" << yrms << "	" << xErrMax << "	" << yErrMax << endl; }
//		
		if (kIterOut % 100 == 0)
		{
			testInt = elements.invertedElements(coorp);
		}
	}
}

int Numerics::kountNeis(int nodeId) {

	int neiTot = 0;

	for (int i = ndeg[nodeId - 1] + 1; i <= ndeg[nodeId]; i++)
	{
		neiTot++;
	}


	return neiTot;
}

void Numerics::globalToLocal(int ndim, vector<double> &cM, int nodeID) {

	int inei;
	cM.resize(ndim);

	for (int id = 0; id < ndim; id++)
	{
		cM[id] = coor[id][nodeID-1]; // exw valei -1 sto coor[id][nodeID]
	}

	for (int k = ndeg[nodeID - 1] + 1; k <= ndeg[nodeID]; k++)
	{
		inei = jaret[k];
		for (int id = 0; id < ndim; id++)
		{
			coor[id][inei - 1] = coor[id][inei - 1] - cM[id];     // exw valei -1 sto inei
			coorp[id][inei - 1] = coorp[id][inei - 1] - cM[id];	// exw valei -1 sto inei
		}
	}
}

void Numerics::localToGlobal(int ndim, vector<double> &cM, int nodeID) {

	int inei;
	for (int k = ndeg[nodeID - 1] + 1; k <= ndeg[nodeID]; k++)
	{
		inei = jaret[k];
		for (int id = 0; id < ndim; id++)
		{
			coor[id][inei - 1] = coor[id][inei - 1] + cM[id]; // exw valei -1 sto inei
			coorp[id][inei - 1] = coorp[id][inei - 1] + cM[id];// exw valei -1 sto inei
		}
	}
}

bool Numerics::converged(double dxOld, double dyOld, double thOld, double dx, double dy, double theta) {

	double kconvT = 0;
	bool convergence = false;

	if (abs(dxOld - dx) < 1.e-14) { kconvT++; }
	if (abs(dyOld - dy) < 1.e-14) { kconvT++; }
	if (abs(thOld - theta)< 1.e-14){ kconvT++; }

	if (kconvT == 3) { convergence = true; }

	return convergence;
}

Numerics::~Numerics(){}