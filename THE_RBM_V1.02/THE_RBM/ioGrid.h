#pragma once

#include <vector>
#include <iostream>
#include <fstream>
#include "DataStructures.h"

using namespace std;


class ioGrids
{
public:
	ioGrids();

	void readInitial2D_Structured(string fileName);
	void write2D_Structured(gridInfo2D_Structured gI, string gridName);

	void readInitial2D_Unstructured(string eleFile, string nodFile);
	void readDeformed2D_Unstructured(string defFile);
	void write2D_Unstructured(double ** coorp, double **coor, vector<int> &logfr, int ns, string gridName);
	
	gridInfo2D_Structured getInfo2D_Structured() { return gInfo_st; }
	gridInfo2D_Unstructured getInfo2D_Unstructured() { return gInfo_un; }
	
	~ioGrids();
private:
	void stringIndexing2D();
	void meshQuality(double** coor, int np, vector<int> &nu);
	
	string fileName;
	gridInfo2D_Structured gInfo_st;
	gridInfo2D_Unstructured gInfo_un;	
	
};
