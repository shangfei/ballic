/*
 * Copyright (c) 2016 Christian Reinhardt.
 *
 * This program solves the structure equations to obtain an equilibrium model
 * for different values rhos=rho(r=R) and us=u(r=R).
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <assert.h>
#include "tillotson/tillotson.h"
#include "solve_model_array.h"

/*
 * Initialize the variables and allocate memory.
 */
MODEL *modelInit(int iMat)
{
    /* Initialize the model */
	MODEL *model;
    
    model = (MODEL *) malloc(sizeof(MODEL));
    assert(model != NULL);

    model->dKpcUnit = 2.06701e-13;
    model->dMsolUnit = 4.80438e-08;

	/*
	 * Initialize one material.
	 * i=0: Granite
	 * i=1: Iron
	 * i=2: Basalt
	 * i=3: Ice
	 */
	model->tillMat = tillInitMaterial(iMat, model->dKpcUnit, model->dMsolUnit, 1000, 1000, 100.0, 1200.0, 1);

	// Debug information
	fprintf(stderr,"\n");	
	fprintf(stderr,"Material: %i\n",iMat);	
	fprintf(stderr,"a: %g\n", model->tillMat->a);
	fprintf(stderr,"b: %g\n", model->tillMat->b);
	fprintf(stderr,"A: %g\n", model->tillMat->A);
	fprintf(stderr,"B: %g\n", model->tillMat->B);
	fprintf(stderr,"rho0: %g\n", model->tillMat->rho0);
	fprintf(stderr,"u0: %g\n", model->tillMat->u0);
	fprintf(stderr,"us: %g\n", model->tillMat->us);
	fprintf(stderr,"us2: %g\n", model->tillMat->us2);
	fprintf(stderr,"alpha: %g\n", model->tillMat->alpha);
	fprintf(stderr,"beta: %g\n", model->tillMat->beta);
	fprintf(stderr,"cv: %g\n", model->tillMat->cv);
	fprintf(stderr,"\n");

    model->nTableMax = 10000; 
    model->M = (double *) malloc(model->nTableMax*sizeof(double));
    assert(model->M != NULL);
    model->rho = (double *) malloc(model->nTableMax*sizeof(double));
    assert(model->rho != NULL);
    model->u = (double *) malloc(model->nTableMax*sizeof(double));
    assert(model->u != NULL);
    model->r = (double *) malloc(model->nTableMax*sizeof(double));
    assert(model->r != NULL);
    model->dr =  0.0;
    model->nTable = 0;
    
    return(model);
}

/*
 * Here we assume an isentropic internal energy profile
 *
 * du/drho = P(rho,u)/rho^2
 *
 * from the first law of thermodynamics.
 */
double dudrho(MODEL *model,double rho,double u)
{
	return(tillPressure(model->tillMat,rho,u)/(rho*rho));
}

double dudr(MODEL *model,double r,double rho,double M,double u)
{
	return(dudrho(model,rho,u)*drhodr(model,r,rho,M,u));
}

/*
 * The derivative  of the density with respect to r is obtained
 * from the equation for hydrostatic equilibrium.
 */
double drhodr(MODEL *model,double r,double rho,double M,double u)
{
	double dPdrho,dPdu;
	
	/*
	 * Derivatives of P with respect to rho and u.
	 */
	dPdrho=tilldPdrho(model->tillMat, rho, u); // dP/drho at u=const.
	dPdu = tilldPdu(model->tillMat, rho, u);; // dP/du at rho=const.

	/*
	 * drho/dr = -G*M*rho/(dPdrho+dPdu*dudrho)
	 */
	assert(r >= 0.0);
	if (r > 0.0)
	{
		// We assume G=1
		return(-M*rho/(r*r*(dPdrho + dPdu*dudrho(model,rho,u))));
	} else {
		// Avoid problems for r=0.
		return(0.0);
	}
}

/*
 * This derivative is independent of the model and only involves geometry.
 */
double dMdr(double r,double rho)
{
	assert(r >= 0.0);
	return(4.0*M_PI*r*r*rho);
}

/*
 * This function solves the model as an initial value problem with
 *
 * rho_initial=rho, M_initial=0 and u_initial=u at r=0
 *
 * and returns the mass when rho == model->tillMat[i]->rho0. The
 * final radius of the model is returned in *pR. For bSetModel=1
 * the results are stored in the look up table. If an error occurs
 * the function returns M < 0.
 */
double midPtRK(MODEL *model,int bSetModel,double rho,double u,double h,double *pR, double *us)
{
    FILE *fp;
    double M = 0.0;
    double r = 0.0;
    double k1rho,k1M,k1u,k2rho,k2M,k2u,x;
    int i;

	/*
	 * Check, if the initial values are below the cold curve.
	 */
	if (tillIsBelowColdCurve(model->tillMat, rho, u))
	{
		fprintf(stderr,"rhoc= %g uc= %g below the cold curve.\n", rho, u);
		M = -1.0;
		return(M);
	}

    if (bSetModel)
	{
		i = 0;
		model->rho[i] = rho;
		model->M[i] = M;
		model->u[i] = u;
		model->r[i] = r;
		fp = fopen("ballic.model","w");
		assert(fp != NULL);
		fprintf(fp,"%g %g %g %g\n",r,rho,M,u);
		++i;
	}

    while (rho > model->tillMat->rho0)
	{
		/*
		* Midpoint Runga-Kutta (2nd order).
		*/
		k1rho = h*drhodr(model,r,rho,M,u);
		k1M = h*dMdr(r,rho);
		k1u = h*dudr(model,r,rho,M,u);

		k2rho = h*drhodr(model,r+0.5*h,rho+0.5*k1rho,M+0.5*k1M,u+0.5*k1u);
		k2M = h*dMdr(r+0.5*h,rho+0.5*k1rho);
		k2u = h*dudr(model,r+0.5*h,rho+0.5*k1rho,M+0.5*k1M,u+0.5*k1u);

		rho += k2rho;
		M += k2M;
		u += k2u;
		r += h;

		if (u < 0.0)
		{
			M = -1.0;
			return(M);
		}

//		printf("%g %g %g %g\n",r,rho,M,u);
		if (bSetModel)
		{
	    	model->rho[i] = rho;
			model->M[i] = M;
			model->u[i] = u;
			model->r[i] = r;
			fprintf(fp,"%g %g %g %g\n",r,rho,M,u);
			++i;
	    }
	}

    /*
	 * Now do a linear interpolation to rho == rho0.
	 */
	x = (model->tillMat->rho0 - rho)/k2rho;
	//	assert(x <= 0.0);
	if (x > 0.0)
	{
		M = -1.0;
		return(M);
	}

	r += h*x;
	M += k2M*x;
	rho += k2rho*x;
	u += k2u*x;

	if (bSetModel)
	{
		--i;
		model->M[i] = M;
		model->r[i] = r;
		model->rho[i] = rho;
		model->u[i] = u;

		fprintf(fp,"%g %g %g %g\n",r,rho,M,u);
		fclose(fp);
		++i;

		model->nTable = i;
		model->dr = h;
	}

    if (us != NULL) *us = u;
    *pR = r;
    return(M);
}


double modelSolve(MODEL *model,double M, double uc) {
    const int nStepsMax = 10000;
    int bSetModel;
    double rmax;
    double dr,R;
    double a,Ma,b,Mb,c,Mc;

	/*
	 * First estimate the maximum possible radius.
	 */
	R = cbrt(3.0*M/(4.0*M_PI*model->tillMat->rho0));
	dr = R/nStepsMax;
	a = 1.01*model->tillMat->rho0; /* starts with 1% larger central density */
	Ma = midPtRK(model,bSetModel=0,a,uc,dr,&R,NULL);
	fprintf(stderr,"first Ma:%g R:%g\n",Ma,R);
	b = a;
	Mb = 0.5*M;

	while (Ma > M)
	{
		b = a;
		Mb = Ma;
		a = 0.5*(model->tillMat->rho0 + a);
		Ma = midPtRK(model,bSetModel=0,a,uc,dr,&R,NULL);
	}
    while (Mb < M)
	{
		b = 2.0*b;
	   	Mb = midPtRK(model,bSetModel=0,b,uc,dr,&R,NULL);	
		fprintf(stderr,"first Mb:%g R:%g\n",Mb,R);
	}

	// (CR) Debug
	fprintf(stderr,"Root bracketed.\n");

    /*
	 * Root bracketed by (a,b).
	 */
    while (Mb-Ma > 1e-10*Mc)
	{
		c = 0.5*(a + b);
		Mc = midPtRK(model,bSetModel=0,c,uc,dr,&R,NULL);	
		if (Mc < M)
		{
			a = c;
			Ma = Mc;
	    } else {
			b = c;
			Mb = Mc;
		}
//	fprintf(stderr,"c:%.10g Mc:%.10g R:%.10g\n",c/model->tillMat[0]->rho0,Mc,R);
	}

    /*
	 * Solve it once more setting up the lookup table.
	 */
	fprintf(stderr,"rho_core: %g cv: %g uc: %g (in system units)\n",c,model->tillMat->cv,uc);
	Mc = midPtRK(model,bSetModel=1,c,uc,dr,&R,NULL);
	model->R = R;
	return c;
}

void main(int argc, char **argv)
{
	const int nStepsMax = 10000;
//    double rhoCenter, uCenter, mTot;
	// Model
    MODEL *model;
	double rho, rhomin, rhomax;
	double u, umin, umax;
	int nRho, nU;
	double dr,R, M, mTot, us, u_desired;
	int iMat;
    FILE *fp;
    int i,j;

#if 0
	/*
	 * Read command line arguments.
	 */
    if (argc != 3)
	{
		fprintf(stderr,"Usage: modelsolve <TotalMass> <iMat>\n");
		exit(1);
	}
    mTot = atof(argv[1]);
    iMat = atoi(argv[2]);

	/* Assure that the parameters make sense. */

    assert(nDesired > 0);
    assert(mTot > 0.0);
    assert(ucore > 0.0);
    assert(iMat >= 0);
#endif

	// Hard code all values
	mTot = 62.366;			// About 1 Earth mass
	iMat = 0;				// Granite
    u_desired = 0.4;        // Desired internal energy at the surface
	rhomin = 1.05*7.33;		// 1.05*rho0
	rhomax = 3.0*7.33;		// 3.0*rho0
	umin = 1e-3;
	umax = 22.0;			// umax > us2
	
	// Do a nRho x nU grid
	nRho = 100;
	nU = 100;

	// Initialize model
	model = modelInit(iMat);
	tillInitLookup(model->tillMat);

	// Open output file
    fp = fopen("modelsolve.txt","w");
    assert(fp != NULL);

	M = mTot;
	/*
	 * First estimate the maximum possible radius.
	 */
	R = cbrt(3.0*M/(4.0*M_PI*model->tillMat->rho0));
	dr = R/nStepsMax;
	fprintf(stderr,"First guess: R= %g dr= %g\n", R, dr);

#if 0
	/*
	 * Some special case to debug the code.
	 */
	rho = 19.6054;
	u = 10.2;
	M = midPtRK(model,1,rho,u,dr,&R);
#endif
	for (i=0; i<=nRho; i++)
	{
		rho = rhomin + (rhomax-rhomin)/nRho*i;

		for (j=0; j<=nU; j++)
		{
			u = umin + (umax-umin)/nU*j;
			// Solve the model for rho, M and u
			M = midPtRK(model,0,rho,u,dr,&R,&us);

//			printf("i=%i j=%i rho=%15.7E u=%15.7E M=%15.7E\n",i,j,rho,u,M);
//			fprintf(fp,"%15.7E",M);
#if 0
            /*
             * Mark models that agree up to 1% with the desired mass.
             */
			if (M > 0.0 && fabs(M-mTot) < 0.01*mTot)
			{
				fprintf(fp,"%3i",1);
			} else if (M < 0.0) {
				fprintf(fp,"%3i",-1);
			} else {
				fprintf(fp,"%3i",0);
			}
#endif
			if (M > 0.0 && fabs(M-mTot) < 0.01*mTot && fabs(us-u_desired) < 5e-2)
			{
				fprintf(fp,"%3i",1);
			} else if (M < 0.0) {
				fprintf(fp,"%3i",-1);
			} else {
				fprintf(fp,"%3i",0);
			}
#if 0
            /*
             * Store (M-Mtot)/Mtot.
             */
			fprintf(fp,"%15.7E", (M-mTot)/mTot);
#endif
		}

		fprintf(fp,"\n");
	}
	fclose(fp);
}
