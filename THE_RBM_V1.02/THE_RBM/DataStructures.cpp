#include "DataStructures.h"

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

Structures2D::Structures2D(){}

Structures2D::Structures2D(int ns_, int np_, int nq_, vector <int> & nu_) : ns(ns_), np(np_), nq(nq_), nu(nu_){}

void Structures2D::Create()
{
	if (np + nq > maxpyr)
	{
		cout << "geom2d: Increase maxpyr at DataStructures.h" << endl;
		exit(0);
	}
	if (ns > maxnod)
	{
		cout << "geom2d: Increase maxnod at DataStructures.h" << endl;
	}

//	--------------------------------------------------
//	Numerate the segments(nubo, nusg, ibsg2tr)
//	 -------------------------------------------------
	numsegs(); // provisionally uses jaret, ndeg
	cout << "numsegs completed...." << endl;
	fjaret();
	cout << "fjaret completed...." << endl;

//	  ------
//	  Ending
//	  ------
	cout << " Number of nodes                : " << ns << endl;
	cout << " Number of segments			 : " << nseg << endl;
	cout << " Number of triangles			 : " << np << endl;
	cout << " Number of quadrangles			 : " << nq << endl;
	cout << " Number of boundary segments	 : " << nbseg << endl;

	cout << " Geom2d completed ..." << endl;
}

void Structures2D::numsegs()
{
	int is, kpoi_1, is1, is2, iprmid, kount, k1, k2,
		is3, is4, ivseg1, ivseg2;
	int ind = 0, ind2 = 0;
	int **iquadseg;
	iquadseg = new int*[4];//side enumeration for quadrangles
	for (unsigned int i = 0; i < 4; ++i)
		iquadseg[i] = new int[2];
	iquadseg[0][0] = 2; iquadseg[1][0] = 3; iquadseg[2][0] = 4; iquadseg[3][0] = 1;
	iquadseg[0][1] = 3; iquadseg[1][1] = 4; iquadseg[2][1] = 1; iquadseg[3][1] = 2;
	vector <int> nod;
	nod.resize(3);
	ndeg.resize(ns + 1);
	nusg.resize(maxnu);
	jaret.resize(nvmaxall);
	ibsg2tr.resize(maxbseg);
	int **iauxn; //opposite nodes of segments for triangles
	iauxn = new int*[3];
	for (unsigned int i = 0; i < 3; ++i)
		iauxn[i] = new int[2];
	iauxn[0][0] = 2; iauxn[1][0] = 1; iauxn[2][0] = 1;
	iauxn[0][1] = 3; iauxn[1][1] = 3; iauxn[2][1] = 2;
	nubo = new int*[2];
	for (unsigned int i = 0; i < 2; ++i)
		nubo[i] = new int[maxseg];
	nseg = 0;
	nbseg = 0;
//
//	Only in this subroutine nu turns -ve and
//	jaret is provisionally the triangles around a node
	for (int i = 0; i < nvmaxall; i++) jaret[i] = 0;
//
//  ----------------------------------------------------------------
//	Find the elements around a node (ATTENTION: in bounds,#elements 
//									 around a node is different from
//									 #segments around a node!!)
//	----------------------------------------------------------------
	for (int i = 0; i < np; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			is = nu[ind];
			ndeg[is] = ndeg[is] + 1;
			ind += 1;
		}
	}
	for (int i = np; i < np + nq; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			is = nu[ind];
			ndeg[is] = ndeg[is] + 1;
			ind += 1;
		}
	}
//
	for (int k = 1; k <= ns; k++)
	{
		if (ndeg[k] == 0)
		{
			cout << endl;
			cout << "ERROR: hanging node, index: " << k << endl;
			exit(0);
		}
		ndeg[k] = ndeg[k] + ndeg[k - 1]; // build index
		jaret[ndeg[k]] = ndeg[k - 1];	 // provisory
	}
	if (2 * ndeg[ns] > nvmaxall)
	{
		cout << "numsegs:Increase nvmaxall at DataStructures.h-->" << 2 * ndeg[ns] << "++" << endl;
		exit(0);
	}
//
	ind = 0;
	for (int i = 1; i <= np; i++)
	{
		for (int j = 1; j <= 3; j++)
		{
			is = nu[ind];
			jaret[ndeg[is]] = jaret[ndeg[is]] + 1;
			kpoi_1 = jaret[ndeg[is]];
			jaret[kpoi_1] = i;	// Element
			ind++;
		}
	}
	for (int i = np + 1; i <= np + nq; i++)
	{
		for (int j = 1; j <= 4; j++)
		{
			is = nu[ind];
			jaret[ndeg[is]] = jaret[ndeg[is]] + 1;
			kpoi_1 = jaret[ndeg[is]];
			jaret[kpoi_1] = i;	// Element
			ind++;
		}
	}

//
//	-----------------------------------
//	Numerate the segments and find NUBO
//	-----------------------------------
	ind = 0;
	std::vector<int>::iterator it;
	it = nu.begin();
	nu.insert(it,0);
	for (int ip = 1; ip <= np; ip++)
	{
		for (int m = 0; m < 3; m++)
		{
//
//			Each segment is identidief by its opposite node, in the triangle:	
			ind++;
			if (nu[ind] < 0){ continue; } // this segment has been treated before
//
//			Else, this is a NEW segment, then:
			ind2 = (ip - 1) * 3;
			is1 = abs(nu[ind2 + iauxn[m][0]]);
			is2 = abs(nu[ind2 + iauxn[m][1]]);
			
			if (nseg + 2 * nq > maxseg)
			{
				cout << "numsegs: Increase maxseg at DataStructures.h" << endl;
				exit(0);
			}
			if (nbseg > maxbseg)
			{
				cout << "numsegs: Increase maxbseg at DataStructures.h" << endl;
				exit(0);
			}
			ibsg2tr[nbseg] = nseg;
			nubo[0][nseg] = is1;
			nubo[1][nseg] = is2;

			nusg[ind] = nseg;
			nu[ind] = -nu[ind]; // attention, nu turn -ve
//
//          Find adjacent segments to avoid duplicating.
//			Also mark the segments that lie over the boundary
			for (int ik = ndeg[is1 - 1] + 1; ik <= ndeg[is1]; ik++)
			{
				iprmid = jaret[ik];
//				
				if (iprmid == ip) continue;
				if (iprmid <= np)
				{
					kount = 0;
					ind2 = (iprmid - 1) * 3;
					for (int m1 = 0; m1 < 3; m1++)
					{						
						k1 = abs(nu[ind2 + m1 + 1]);
						k2 = 1;
						if (k1 == is1 || k1 == is2) { k2 = 0; }
						if (k2 != 0)
						{
							kount++;
							nod[0] = m1 + 1;
						}
					}
					if (kount == 2) { continue; }
					if (kount == 0) { cout << "numsegs Error, with data !" << endl; exit(0); }

					nu[ind2 + nod[0]] = -nu[ind2 + nod[0]]; // for this also,nu turns-ve
					int check;
					check = nu[ind2 + nod[0]];
					nusg[ind2 + nod[1]] = nseg;
					nbseg = nbseg - 1; // This was a boundary segment
					break;
				}
				else
				{
					kount = 0;
					ind2 = np * 3 + (iprmid - np - 1) * 4;
					for (int m1 = 1; m1 <= 4; m1++)
					{
						k1 = abs(nu[ind2 + m1]);
						k2 = 1;
						if (k1 == is1 || k1 == is2){ k2 = 0; }
						if (k2 != 0)
						{
							nod[kount] = m1;
							kount++;
						}
					}
					if (kount == 3){ continue; }
					if (kount <= 1)
					{
						cout << "numsegs: Error, with data !" << endl;
						exit(0);
					}
					if (nod[0] == 1 && nod[1] == 4){ nod[1] = 1; }
					nu[ind2 + nod[1]] = -nu[ind2 + nod[1]];
					nusg[ind2 + nod[1]] = nseg;
					nbseg = nbseg - 1; // This was not a boundary segment
					break;
				}
			}
			nseg++;
			nbseg++;
		}
	}
//
	for (int ip = np + 1; ip <= np + nq; ip++)	// Loop on quandrangle
	{
		for (int m = 0; m < 4; m++)	//Loop on the four nodes of each quadrangle
		{
			ind++;
//			Each segment is identified according to the enumeration in iquadseg:
			if (nu[ind] < 0){ continue; }
//			Else, this is a NEW segment, then:
			ind2 = np * 3 + (ip - np - 1) * 4;
			is1 = abs(nu[ind2 + iquadseg[m][0]]);
			is2 = abs(nu[ind2 + iquadseg[m][1]]);
			if (nseg + 2 * nq > maxseg)
			{
				cout << "numsegs: Increase maxseg at DataStructures.h " << endl; // EULER
				exit(0);
			}
			if (nbseg > maxbseg)
			{
				cout << "numsegs: Increase maxbseg at DataStructures.h " << endl;
			}
			ibsg2tr[nbseg] = nseg;
			nubo[0][nseg] = is1;
			nubo[1][nseg] = is2;

			nusg[ind] = nseg;
			nu[ind] = -nu[ind]; // attention, nu turns -ve
//
//			Find adjacent segment to avoid duplicating.
//			Also mark the segments that lie over the boundary.
			for (int ik = ndeg[is1 - 1] + 1; ik <= ndeg[is1]; ik++)
			{
				iprmid = jaret[ik];
//
				if (iprmid == ip){ continue; } // skip the quadrangle at hand
				if (iprmid <= np)
				{
					kount = 0;
					ind2 = (iprmid - 1) * 3;
					for (int m1 = 1; m1 <= 3; m1++) // Loop on the 3 nodes of the triangle
					{
						k1 = abs(nu[ind2 + m1]);
						k2 = 1;
						if (k1 == is1 || k1 == is2){ k2 = 0; }
						if (k2 != 0)
						{
							kount++;
							nod[0] = m1;
						}
					}
					if (kount == 2){ continue; }
					if (kount == 0)
					{
						cout << "numsegs: Error, with data !" << endl;
						exit(0);
					}
					nu[ind2 + nod[1]] = -nu[ind2 + nod[0]]; // for this also, nu turns-ve
					nusg[ind2 + nod[0]] = nseg;
					nbseg = nbseg - 1; // This was not a boundary segment
					break;
				}
				else
				{
					kount = 0;
					ind2 = np * 3 + (iprmid - np - 1) * 4;
					for (int m1 = 1; m1 <= 4; m1++) // Loop on the 4 nodes of the quadrangle
					{
						k1 = abs(nu[ind2 + m1]);
						k2 = 1;
						if (k1 == is1 || k1 == is2) { k2 = 0; }
						if (k2 != 0)
						{
							nod[kount] = m1;
							kount++;
						}
					}
					if (kount == 3) { continue; }
					if (kount <= 1)
					{
						cout << "numsegs: Error, with data !" << endl;
						exit(0);
					}
					if (nod[0] == 1 && nod[1] == 4){ nod[1] = 1; }
					nu[ind2 + nod[1]] = -nu[ind2 + nod[1]]; //for this also, nu turn-ve
					int check;
					check = nu[ind2 + nod[1]];
					nusg[ind2 + nod[1]] = nseg;
					nbseg = nbseg - 1; // This was a boundary segment
					break;
				}
			}
			nseg++;
			nbseg++;
		}
	}
//
//     -------------------------
//		Reconstitute positive NU
//     -------------------------
	kount = 0;
	ind = 0;
	for (int ip = 1; ip <= np; ip++)
	{
		for (int i = 1; i <= 3; i++)
		{
			nu[ind + i] = -nu[ind + i];
			if (nu[ind + i] < 0)
			{
				cout << ind + i << endl;
				cout << "tr:" << ip << " Nd:" << nu[ind + 1] << " " << nu[ind + 2] << " " << nu[ind + 3] << endl;
				kount++;
			}
		}
		ind = ind + 3;
	}
	for (int ip = 1; ip <= nq; ip++)
	{
		for (int i = 1; i <= 4; i++)
		{
			nu[ind + i] = -nu[ind + i];
			if (nu[ind + i] < 0)
			{
				cout << "quad: " << ip << " Nd:" << nu[ind + 1] << nu[ind + 2] << nu[ind + 3] << nu[ind + 4] << endl;
				kount++;
			}
		}
		ind = ind + 4;
	}

	if (kount > 0)
	{
		cout << "numsegs: Error with data !" << endl;
		exit(0);
	}
//
//     ------------------------------
//	   Build NUBO of VIRTUAL segments
//	   ------------------------------
	ind = np * 3;
	nvseg = 0;
	for (int ip = np + 1; ip <= np + nq; ip++)
	{
		is1 = nu[ind + 1];
		is2 = nu[ind + 2];
		is3 = nu[ind + 3];
		is4 = nu[ind + 4];
		ivseg1 = nseg + nvseg;
		ivseg2 = nseg + nvseg + 1;
		nubo[0][ivseg1] = is3;
		nubo[1][ivseg1] = is1;
		nubo[0][ivseg2] = is4;
		nubo[1][ivseg2] = is2;
		ind = ind + 4;
		nvseg = nvseg + 2;
	}
	int dum2;
	dum2 = 1;
}

void Structures2D::fjaret()
{
	int inod1, inod2, kpoi_1, kpoi_2, max_str;
//
//	Jaret: Serial Storage of Jaret
	for (int i = 0; i <= ns; i++){ ndeg[i] = 0; }
	for (int i = 0; i < nvmaxall; i++){ jaret[i] = 0; }
//	
//	Find # of nodes - segments around each node
//	ATT: Jaret also includes neighbours due to virtual segments
	for (int iseg = 0; iseg < nseg + nvseg; iseg++)
	{
		inod1 = nubo[0][iseg];
		inod2 = nubo[1][iseg];
		ndeg[inod1] = ndeg[inod1] + 1;
		ndeg[inod2] = ndeg[inod2] + 1;
	}
//
	for (int k = 1; k <= ns; k++)
	{
		ndeg[k] = ndeg[k - 1] + ndeg[k]; // build index
		jaret[ndeg[k]] = ndeg[k - 1];	 // provisory
	}
//
// Calculate Max Length of Jaret(String), Only Node Storing
	for (int iseg = 0; iseg < nseg + nvseg; iseg++)
	{
		inod1 = nubo[0][iseg];
		inod2 = nubo[1][iseg];
		jaret[ndeg[inod1]] = jaret[ndeg[inod1]] + 1;
		kpoi_1 = jaret[ndeg[inod1]];
		jaret[kpoi_1] = inod2;
		jaret[ndeg[inod2]] = jaret[ndeg[inod2]] + 1;
		kpoi_2 = jaret[ndeg[inod2]];
		jaret[kpoi_2] = inod1;
	}
//
//	Create Extra Segments
	max_str = ndeg[ns];	// Length of String when all
						// JARET(Nodes) are Stored
//
	if (2 * ndeg[ns] > nvmaxall)
	{
		cout << "fjaret:Increase nvmaxall at DataStructures.h-->" << 2 * ndeg[ns] << "++" << endl;
	}
	for (int k = 1; k <= ns; k++) { jaret[ndeg[k] + max_str] = max_str + ndeg[k - 1]; }
//
	for (int iseg = 0; iseg < nseg + nvseg; iseg++)
	{
		inod1 = nubo[0][iseg];
		inod2 = nubo[1][iseg];
		jaret[ndeg[inod1] + max_str] = jaret[ndeg[inod1] + max_str] + 1;
		kpoi_1 = jaret[ndeg[inod1] + max_str];
		jaret[kpoi_1] = iseg;
		jaret[ndeg[inod2] + max_str] = jaret[ndeg[inod2] + max_str] + 1;
		kpoi_2 = jaret[ndeg[inod2] + max_str];
		jaret[kpoi_2] = iseg + 1;
	}
}

Structures2D::~Structures2D(){}