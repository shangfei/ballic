/*
 * Copyright (c) 2016 Christian Reinhardt.
 *
 * This program solves the structure equations to obtain an equilibrium model
 * for given mass M, surface density rho0 and internal energy us using the
 * matching point method.
 *
 * The free variable is the enclosed mass m = m(r), so the boundary conditions
 * are r(m=M)=0 at the core and rho(m=M)=rho0, u(m=M)=us at the surface.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <assert.h>
#include "tillotson/tillotson.h"
#include "modelsolve.h"

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
double dudrho(MODEL *model, double rho, double u)
{
	return(tillPressure(model->tillMat,rho,u)/(rho*rho));
}

/*
 * We have
 *
 * du/dm = du/drho*drho/dm
 *
 * since u = u(rho) for constant entropy.
 */
double dudm(MODEL *model, double M, double r, double rho, double u)
{
	return(dudrho(model, rho, u)*drhodm(model, M, r, rho, u));
}

/*
 * The derivative  of the density with respect to m is obtained
 * from the equation for hydrostatic equilibrium.
 */
double drhodm(MODEL *model, double M, double r, double rho, double u)
{
	double dPdrho, dPdu;
	/*
	 * Derivatives of P with respect to rho and u.
	 */
	dPdrho = tilldPdrho(model->tillMat, rho, u); // dP/drho at u=const.
	dPdu = tilldPdu(model->tillMat, rho, u); // dP/du at rho=const.

	/*
	 * drho/dm = -G*M/(r^4*(dPdrho+dPdu*dudrho)
	 */
	assert(M>= 0.0);

	if (M > 0.0)
	{
		// We assume G=1
		return(-M/(r*r*r*r*(dPdrho + dPdu*dudrho(model,rho,u))));
    } else {
		// Avoid problems for m=0.
		return(0.0);
	}
}

/*
 * This derivative is independent of the model and only involves geometry.
 */
double drdm(double r, double rho)
{
	assert(r >= 0.0);
	return(1.0/(4.0*M_PI*r*r*rho));
}

/*
 * This function solves the equations of hydrostatic equlibrium from the
 * surface to the matching point using M as a free variable with b.c.
 *
 * rho(m=M)=rho0 and u(m=M)=us
 *
 * where r(m=M)=R is undetermined and given as a parameter. It returns an error
 * code (0: ok (TRUE), 1: fail (FALSE)) and the values
 *
 * r(m=M_mid)=r_L, rho(m=M_mid) = rho_L and u(m=M_mid) = u_L
 *
 * at the midpoint M_mid in *pR, *pRho and *pU. The parameter h determines
 * the step size used in the integration.
 */
int midPtRKIn(MODEL *model, double Mtot, double R, double rho0, double us, double h, double M_mid, double *pR, double *pRho, double *pU)
{
    double M = Mtot;
    double r = R;
    double rho = rho0;
    double u = us;
    double k1r, k1rho, k1u, k2r, k2rho, k2u, x;
	/*
	 * Check, if the initial values are sensible.
	 */
    assert(M > 0.0);
    assert(M_mid >= 0.0);
    assert(u > 0.0);
    assert(r > 0.0);
    assert(h > 0.0);

    while (M > M_mid)
	{
		/*
		* Midpoint Runga-Kutta (2nd order).
		*/
		k1r = h*drdm(r, rho);
        k1rho = h*drhodm(model, M, r, rho, u);
		k1u = h*dudm(model, M, r, rho, u);

		k2r = h*drdm(r-0.5*k1r, rho-0.5*k1rho);
		k2rho = h*drhodm(model, M-0.5*h, r-0.5*k1r, rho-0.5*k1rho, u-0.5*k1u);
		k2u = h*dudm(model, M-0.5*h, r-0.5*k1r, rho-0.5*k1rho, u-0.5*k1u);

		r -= k2r;
		rho -= k2rho;
		u -= k2u;
		M -= h;

        printf("%15.7E %15.7E %15.7E %15.7E\n", r, rho, M, u);
	}

    /*
	 * Now do a linear interpolation to M == M_mid0.
     *
     * y = A*y_n-1 + B*yn where A=(M_n-M_mid)/(M_n-M_n-1) and B=1-A
     *
     * so we get
     *
     * y = yn - x*k2y where k2y = y_n-y_n-1
	 */
    x = (M-M_mid)/h;
    assert(x <= 0.0);

	M -= h*x;
	r -= k2r*x;
	rho -= k2rho*x;
	u -= k2u*x;

    printf("%15.7E %15.7E %15.7E %15.7E\n", r, rho, M, u);
    printf("\n");

    // Return R, rho and u
    *pR = r;
    *pRho=rho;
    *pU=u;

    return(0);
}

/*
 * This function solves the equations of hydrostatic equlibrium from the core 
 * to the matching point using M as a free variable with b.c.
 *
 * r(m=0)=0
 *
 * where rho(m=0)=rhoc and u(m=0)=uc are undetermined and given as a parameter.
 * It returns an error code (0: ok (TRUE), 1: fail (FALSE)) and the values
 * 
 * r(m=M_mid)=r_L, rho(m=M_mid) = rho_L and u(m=M_mid) = u_L
 *
 * at the midpoint M_mid in *pR, *pRho and *pU. The parameter h determines the
 * step size used in the integration.
 */
int midPtRKOut(MODEL *model, int bSetModel, double rhoc, double uc, double h, double M_mid, double *pR, double *pRho, double *pU)
{
    FILE *fp;
    double M = 0.0;
    double r = 0.0;
    double rho = rhoc;
    double u = uc;
    double k1r, k1rho, k1u, k2r, k2rho, k2u,x;
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
		fprintf(fp,"%15.7E %15.7E %15.7E %15.7E\n", r, rho, M, u);
		++i;
	}

    while (M <= M_mid)
	{
		/*
         * Midpoint Runga-Kutta (2nd order).
		 */
		k1r = h*drdm(r, rho);
        k1rho = h*drhodm(model, M, r, rho, u);
		k1u = h*dudm(model, M, r, rho, u);

		k2r = h*drdm(r+0.5*k1r, rho+0.5*k1rho);
		k2rho = h*drhodm(model, M+0.5*h, r+0.5*k1r, rho+0.5*k1rho, u+0.5*k1u);
		k2u = h*dudm(model, M+0.5*h, r+0.5*k1r, rho+0.5*k1rho, u+0.5*k1u);

		r += k2r;
		rho += k2rho;
		u += k2u;
		M += h;

		if (bSetModel)
		{
            model->M[i] = M;
			model->r[i] = r;
            model->rho[i] = rho;
			model->u[i] = u;

		    fprintf(fp,"%15.7E %15.7E %15.7E %15.7E\n", r, rho, M, u);
			++i;
	    }
	}

    /*
	 * Now do a linear interpolation to M == M_mid as described in
     * midPtRKOut().
	 */
    x = (M-M_mid)/h;
    assert(x <= 0.0);

	M += h*x;
	r += k2r*x;
	rho += k2rho*x;
	u += k2u*x;

	if (bSetModel)
	{
		--i;
		model->M[i] = M;
		model->r[i] = r;
		model->rho[i] = rho;
		model->u[i] = u;

		fprintf(fp,"%15.7E %15.7E %15.7E %15.7E\n",r,rho,M,u);
		fclose(fp);
		++i;

		model->nTable = i;
		model->dr = h;
	}

    // Return values at the midpoint
    *pR = r;
    *pRho = rho;
    *pU = u;
    return(0);
}
#if 0
NOT IMPLEMENTED YET.
/*
 * Solve the model for a given mass M and internal enery us
 * at the surface. To obtain rho(r=0) and u(r=0) which are
 * needed to currecly solve the structure equations, we first
 * integrate from the surface in to obtain rhoc and uc and then
 * use these values to calculate the final model.
 */
double modelSolve(MODEL *model,double M, double us) {
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
	Ma = midPtRK(model,bSetModel=0,a,uc,dr,&R);
	fprintf(stderr,"first Ma:%g R:%g\n",Ma,R);
	b = a;
	Mb = 0.5*M;

	while (Ma > M)
	{
		b = a;
		Mb = Ma;
		a = 0.5*(model->tillMat->rho0 + a);
		Ma = midPtRK(model,bSetModel=0,a,uc,dr,&R);
	}
    while (Mb < M)
	{
		b = 2.0*b;
	   	Mb = midPtRK(model,bSetModel=0,b,uc,dr,&R);	
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
		Mc = midPtRK(model,bSetModel=0,c,uc,dr,&R);	
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
	Mc = midPtRK(model,bSetModel=1,c,uc,dr,&R);
	model->R = R;
	return c;
}
#endif
#if 0
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
	Ma = midPtRK(model,bSetModel=0,a,uc,dr,&R);
	fprintf(stderr,"first Ma:%g R:%g\n",Ma,R);
	b = a;
	Mb = 0.5*M;

	while (Ma > M)
	{
		b = a;
		Mb = Ma;
		a = 0.5*(model->tillMat->rho0 + a);
		Ma = midPtRK(model,bSetModel=0,a,uc,dr,&R);
	}
    while (Mb < M)
	{
		b = 2.0*b;
	   	Mb = midPtRK(model,bSetModel=0,b,uc,dr,&R);	
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
		Mc = midPtRK(model,bSetModel=0,c,uc,dr,&R);	
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
	Mc = midPtRK(model,bSetModel=1,c,uc,dr,&R);
	model->R = R;
	return c;
}
#endif
void main(int argc, char **argv)
{
	const int nStepsMax = 10000;
//    double rhoCenter, uCenter, mTot;
	// Model
    MODEL *model;
	double rho, rhoc, rhos;
	double u, uc, us;
	double dm, R, M, mTot;
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

    /*
     * Hard code the values for a 1ME model.
     */
	// Hard code all values
	mTot    = 62.366;			// About 1 Earth mass
	iMat    = 0;				// Granite
    rhoc    = 19.6054;
    uc      = 10.2;
    rhos    = 7.32745;
    us      = 0.333586;
    R       = 1.04211;

    double R_L, rho_L, u_L;

	// Initialize model
	model = modelInit(iMat);
//	tillInitLookup(model->tillMat);

    dm = mTot/(nStepsMax);
    midPtRKIn(model, mTot, R, rhos, us, dm, 0.1*mTot, &R_L, &rho_L, &u_L);
    M = 0.1*mTot;

    fprintf(stderr, "%15.7E ", mTot);
    fprintf(stderr, "%15.7E ", R);
    fprintf(stderr, "%15.7E ", rhos);
    fprintf(stderr, "%15.7E ", us);

    fprintf(stderr, "\n");

    fprintf(stderr, "%15.7E ", rhoc);
    fprintf(stderr, "%15.7E ", uc);
    fprintf(stderr, "\n");

    fprintf(stderr, "%15.7E ", M);
    fprintf(stderr, "%15.7E ", R_L);
    fprintf(stderr, "%15.7E ", rho_L);
    fprintf(stderr, "%15.7E ", u_L);

    fprintf(stderr, "\n");

}
