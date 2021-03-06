/*
 ** This is a simple program to test the Tillotson EOS library.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <assert.h>
#include "tillotson.h"

#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) > (B) ? (B) : (A))

void main(int argc, char **argv) {
	/*
	** Calculate the cold curve.
	*/
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
	double rhomax = 25.0;
	int i = 0;

	TILLMATERIAL *granite;

	granite = tillInitMaterial(GRANITE, dKpcUnit, dMsolUnit, rhomax);
	
	fprintf(stderr, "Init cold curve...\n");
	tillInitColdCurve(granite);
	fprintf(stderr, "Done.\n");

	fprintf(stderr,"nTable: %i, nTableMax: %i\n",granite->nTable, granite->nTableMax);

	for (i=0; i < granite->nTable; i++)
	{
		/*fprintf(stderr,"nTable: %i, i: %i\n",granite->nTable,i);*/
		printf("%.30f %.30f %.30f\n", granite->cold[i].rho, granite->cold[i].u, granite->cold[i].dudrho);
	}

	tillFinalizeMaterial(granite);
}
