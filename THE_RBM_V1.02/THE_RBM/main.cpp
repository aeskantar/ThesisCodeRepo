#include <iostream>
#include <math.h>
#include "DataStructures.h"
#include "ioGrid.h"
#include "qualityCheck.h"
#include "Solver.h"

using namespace std;


int main()
{

	gridInfo2D_Unstructured ginfo_un;
	ioGrids ionGrid;
	ionGrid.readInitial2D_Unstructured("hybrid3.ele", "hybrid3.nod");
	ginfo_un = ionGrid.getInfo2D_Unstructured();
	qualityCheck gQuality(ginfo_un.coor, ginfo_un.nu, ginfo_un.np);
	gQuality.meshQuality(ginfo_un.coor, ginfo_un.nu);
	ionGrid.readDeformed2D_Unstructured("hybrid3.def");
	ginfo_un = ionGrid.getInfo2D_Unstructured();
	Structures2D dataStructure(ginfo_un.ns, ginfo_un.np, ginfo_un.nq, ginfo_un.nu);
	dataStructure.Create();
	Numerics newtonRaphson2D(ginfo_un.coorp, ginfo_un.coor, ginfo_un.logfr, dataStructure.get_ndeg(), 
							 dataStructure.get_jaret(), ginfo_un.nu, ginfo_un.ns, ginfo_un.np);
	ginfo_un = ionGrid.getInfo2D_Unstructured();
	newtonRaphson2D.Solver2D(200000);
	ionGrid.write2D_Unstructured(ginfo_un.coorp, ginfo_un.coor,ginfo_un.logfr, ginfo_un.ns, "new.nod");
	

	cin.get();
	return 0;
	
}