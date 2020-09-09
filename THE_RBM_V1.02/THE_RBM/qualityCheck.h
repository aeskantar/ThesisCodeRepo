#pragma once
#include"DataStructures.h"

using namespace std;

class qualityCheck
{
public:
	qualityCheck();
	qualityCheck(double** coorp_,vector<int> nu_,int np_);
	bool invertedElements(double** coorp);
	void meshQuality(double** coor, vector<int>& nu);
	~qualityCheck();

private:
	int np;
	vector<int> nu;
	double** coorp;
	void vecProd(double v1x, double v1y, double v1z, double v2x,
		double v2y, double v2z, double& xre, double& yre, double& zre);
};

