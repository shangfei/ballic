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
#include "nr/nrutil.h"

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

	if (r > 0.0)
	{
		// We assume G=1
		return(-M/(4.0*M_PI*r*r*r*r*(dPdrho + dPdu*dudrho(model,rho,u))));
    } else {
		// Avoid problems for m=0.
		return(0.0);
//        return(-1.0/(4.0*M_PI*(dPdrho + dPdu*dudrho(model,rho,u)))*pow((4.0*M_PI*rho/3.0),4.0/3.0)*pow(M,-1.0/3.0));
	}
}

/*
 * This derivative is independent of the model and only involves geometry.
 */
double drdm(double r, double rho)
{
	assert(r >= 0.0);
    if (r > 0.0)
    {
        return(1.0/(4.0*M_PI*r*r*rho));
    } else {
        /*
         * Use an analytic expression (assuming rho=const.) in the center to
         * start integration. If the user did not do that crash.
         */
        assert(0);
    }
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

    printf("%15.7E %15.7E %15.7E %15.7E\n", r, rho, M, u);

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
    
    /*
     * Stimmt da wirklich so, wir integrieren ja von aussen nach innen....
     */
    fprintf(stderr,"x= %15.7E\n", x);
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
/*
    if (tillIsBelowColdCurve(model->tillMat, rho, u))
	{
		fprintf(stderr,"rhoc= %g uc= %g below the cold curve.\n", rho, u);
		M = -1.0;
		return(M);
	}
*/
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
	} else {
        printf("%15.7E %15.7E %15.7E %15.7E\n", r, rho, M, u);
    }

    /*
     * Because the derivatives diverge for M=0 and r=0 we use analytic
     * expressions for r and rho to start the integration.
     */
    M = 1e-0*h;
    r = pow(3.0*M/(4.0*M_PI*rhoc),1.0/3.0);

    double dPdrho = tilldPdrho(model->tillMat, rhoc, uc); // dP/drho at u=const.
    double dPdu = tilldPdu(model->tillMat, rhoc, uc); // dP/du at rho=const.
    rho = rhoc - 3.0/(8.0*M_PI*(dPdrho + dPdu*dudrho(model, rhoc, uc)))*pow((4.0*M_PI*rhoc/3.0), 4.0/3.0)*pow(M, 2.0/3.0);
    u = uc - tillPressure(model->tillMat,rhoc,uc)/(rhoc*rhoc)*3.0/(8.0*M_PI*(dPdrho + dPdu*dudrho(model, rhoc, uc)))*pow((4.0*M_PI*rhoc/3.0), 4.0/3.0)*pow(M, 2.0/3.0);
    //u +=  1e-3*h*dudm(model, M, r, rho, u);
/*
    u +=  1e-3*h*dudm(model, M, r, rho, u);
    rho +=  1e-3*h*drhodm(model, M, r, rho, u);
    r = pow(3.0*M/(4.0*M_PI*rho),1.0/3.0);
*/

    printf("%15.7E %15.7E %15.7E %15.7E\n", r, rho, M, u);

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
	    } else {
            printf("%15.7E %15.7E %15.7E %15.7E\n", r, rho, M, u);
        }
	}

    /*
	 * Now do a linear interpolation to M == M_mid as described in
     * midPtRKOut().
	 */
    x = (M-M_mid)/h;
    //assert(x <= 0.0);

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
	} else {
        printf("%15.7E %15.7E %15.7E %15.7E\n", r, rho, M, u);
    }

    // Return values at the midpoint
    *pR = r;
    *pRho = rho;
    *pU = u;
    return(0);
}

#if 0
void LinEquationSolver(MODEL *model, double *a, double *b, int n)
{
    float **A = matrix(1, n, 1, n);
    // Setting m=1 since we only have one right hand side vector b
    float **B = matrix(1, n, 1, 1);
    int i, j;

    fprintf(stderr,"n= %i\n", n);
#if 0
    for (i=0; i<n; i++)
    {
//        B[i+1][1] = b[i];
        B[i+1][1] = 1;
        for (j=0; j<n; j++)
        {
            /*
             * Remember that NR arrays start from one while C arrays always
             * start from 0.
             */
//            A[i+1][j+1] = a[i*n + j];
            A[i+1][j+1] = 1/(i+j;
        }
    }
#endif
    for (i=1; i<=n; i++)
    {
        B[i][1] = 0.0;
        for (j=1; j<=n; j++)
        {
            A[i][j] = 1.0/(i+j);
        }
    }
/*
    n = 2;
    A[1][1] = 2.0;
    A[1][2] = 1.0;
    A[2][1] = 1.0;
    A[2][2] = 2.0;
*/
    /*
     * Print the matrix A and vector b.
     */
    for (i=1; i<=n; i++)
    {
        printf("A= ");
        for (j=1; j<=n; j++)
        {
            printf("%15.7E  ", A[i][j]);
        }

        printf(" b= %15.7E\n", B[i][1]);
    }
    printf("\n");
    printf("Calling gaussj()...\n");
    printf("\n");

    gaussj(A, n, B, 1);

    /*
     * Print the matrix A and vector b.
     */
    for (i=1; i<=n; i++)
    {
/*        printf("A= ");
        for (j=1; j<=n; j++)
        {
            printf("%15.7E  ", A[i][j]);
        }
*/
        printf(" b= %15.7E\n", B[i][1]);
    }

}
#endif

void GaussJordanTest(int n)
{
#if 0
    float **A = matrix(1, n, 1, n);
    float **A_old = matrix(1, n, 1, n);
    float **Id = matrix(1, n, 1, n);
    // Setting m=1 since we only have one right hand side vector b
    float **b = matrix(1, n, 1, 1);
#endif
    double **A = dmatrix(1, n, 1, n);
    double **A_old = dmatrix(1, n, 1, n);
    double **Id = dmatrix(1, n, 1, n);
    // Setting m=1 since we only have one right hand side vector b
    double **b = dmatrix(1, n, 1, 1);
    int i, j, k;

    fprintf(stderr,"n= %i\n", n);
    
    /*
     * Setup the matrix A and the vector b.
     */
    for (i=1; i<=n; i++)
    {
        b[i][1] = 1.0;

        for (j=1; j<=n; j++)
        {
            A[i][j] = 1.0/(i+j);
            A_old[i][j] = A[i][j];
        }
    }
/*
    n = 2;
    A[1][1] = 2.0;
    A[1][2] = 1.0;
    A[2][1] = 1.0;
    A[2][2] = 2.0;
*/
    /*
     * Print the matrix A and vector b.
     */
    printf("\n");
    printf("Original A and b:\n");
    printf("\n");
    for (i=1; i<=n; i++)
    {
        printf("A= ");
        for (j=1; j<=n; j++)
        {
            printf("%15.7E  ", A[i][j]);
        }

        printf(" b= %15.7E\n", b[i][1]);
    }

    printf("\n");
    printf("Calling gaussj()...\n");
    printf("\n");

    gaussj(A, n, b, 1);

    /*
     * Print the matrix A and vector b.
     */
    printf("\n");
    printf("Inverse of A and solution b:\n");
    printf("\n");
    for (i=1; i<=n; i++)
    {
        printf("A= ");
        for (j=1; j<=n; j++)
        {
            printf("%15.7E  ", A[i][j]);
        }

        printf(" b= %15.7E\n", b[i][1]);
    }

    /*
     * Calculate A*A_old to see if we get the identity matrix Id.
     */
    for (i=1; i<=n; i++)
    {
        for (j=1; j<=n; j++)
        {
            Id[i][j] = 0.0;
            for (k=1; k<=n; k++)
            {
                Id[i][j] += A[i][k]*A_old[k][j];
            }
        }
    }

    /*
     * Print the identity matrix Id.
     */
    printf("\n");
    printf("Id:\n");
    printf("\n");
    for (i=1; i<=n; i++)
    {
        for (j=1; j<=n; j++)
        {
            printf("%15.7E  ", Id[i][j]);
        }
        printf("\n");
    }
}

/*
 * This function uses the difference
 * ddx = (x2_R-x2_L)-(x1_R-x1_L)
 * between two lightly different initial guesses dx of initial values of for R,
 * rhoc and uc to predict for which values the difference x_R-x_L of the
 * solutions at the matching point will vanish solving a linear system of
 * equations with the Gauss-Jordan method. The estimated corrections for R,
 * rhoc and uc are returned in pDelta.
 *
 * dR = R_new - R_old
 * drhoc = rhoc_new - rhoc_old
 * duc = uc_new - uc_old
 *
 * Y_R = R_L - R_R (values of the old solution?)
 * Y_rho = rho_L - rho_R
 * Y_u=u_L - u_R
 *
 * dY_R = Y_R_new-Y_R_old           (where Y_R=R_L-R_R)
 * dY_rho = Y_rho_new-Y_rho_old     (where Y_rho=rho_L-rho_R)
 * dY_u = Y_u_new-Y_u_old           (where Y_u=u_L-u_R)
 */
void modelGetNewInitialGuess(double dR, double drhoc, double duc, double Y_R, double Y_rho, double Y_u, double dY_R, double dY_rho, double dY_u, double *pDelta)
{
   /*
    * The matrix A contains the (numerical) derivatives of dx with respect to
    * the variables R, rhoc and uc.
    */
    double **A = dmatrix(1, 3, 1, 3);
    double **b = dmatrix(1, 3, 1, 1);
    int i,j;
    // Note that NR addresses arrays from one while in C arrays start at zero.
    A[1][1] = dY_R/dR;
    A[1][2] = dY_R/drhoc;
    A[1][3] = dY_R/duc;

    A[2][1] = dY_rho/dR;
    A[2][2] = dY_rho/drhoc;
    A[2][3] = dY_rho/duc;

    A[3][1] = dY_u/dR;
    A[3][2] = dY_u/drhoc;
    A[3][3] = dY_u/duc;

    /*
     * Set b to (Y_R_old, Y_rho_old, Y_u_old).
     */
    b[1][1] = -Y_R;
    b[2][1] = -Y_rho;
    b[3][1] = -Y_u;

    /*
    fprintf(stderr, "A:\n");
    for (i=1; i<=3; i++)
    {
        for (j=1; j<=3; j++)
        {
            fprintf(stderr, "%i,%i: %15.7E  ", i, j, A[i][j]);
        }
        fprintf(stderr,"\n");
    }
    */

    fprintf(stderr, "A:\n");
    for (i=1; i<=3; i++)
    {
        for (j=1; j<=3; j++)
        {
            fprintf(stderr, "%i,%i: %15.7g  ", i, j, A[i][j]);
        }
        fprintf(stderr,"\n");
    }

    fprintf(stderr, "Calling gaussj()\n");
    /*
     * Call the Gauss-Jordan method from NR, the solution is returned in b.
     */
    gaussj(A, 3, b, 1);

    fprintf(stderr, "Called gaussj()\n");

    fprintf(stderr, "A:\n");
    for (i=1; i<=3; i++)
    {
        for (j=1; j<=3; j++)
        {
            fprintf(stderr, "%i,%i: %15.7E  ", i, j, A[i][j]);
        }
        fprintf(stderr,"\n");
    }

    fprintf(stderr, "\n");
    fprintf(stderr, "b:\n");
    for (i=1; i<=3; i++)
    {
        fprintf(stderr, "%i: %15.7E  ", i,  b[i][1]);
    }
    fprintf(stderr, "\n");

    pDelta[0] = b[1][1];
    pDelta[1] = b[2][1];
    pDelta[2] = b[3][1];
    
    // Free memory
    free_dmatrix(A, 1, 3, 1, 3);
    free_dmatrix(b, 1, 3, 1, 1);
}

void modelSolveToMatchingPoint(MODEL *model, double rhoc, double uc, double M, double R, double rhos, double us, double M_mid, double dm, double *dR, double *drho, double *du)
{
    double R_L, R_R, rho_L, rho_R, u_L, u_R;
    int bSetModel;

    /*
     * First integrate from the surface to the matching point M_mid. The values
     * of R, rho and u at M_mid are returned in R_L, rho_L and u_L.
     */
    midPtRKIn(model, M, R, rhos, us, dm, M_mid, &R_L, &rho_L, &u_L);

    /*
     * Then integrate from the core (m=0) to the matching point.
     */
    midPtRKOut(model, bSetModel=0, rhoc, uc, dm, M_mid, &R_R, &rho_R, &u_R);

    /*
     * Return the differences between R, rho and u at the matching point.
     */
    *dR = R_L-R_R;
    *drho = rho_L-rho_R;
    *du = u_L-u_R;
}

/*
 * Solve a single component model for a given mass M, density rhos and internal
 * energy us at the surface.
 *
 * Because the b.c are given both at the core and the surface and the
 * underlaying ODE's do not behave well at the core (r=0 and m=0) we use the
 * matching point method (ref?) to solve the models. The basic idea here is to
 * solve the ODE's from and surface and the core and make the solutions match
 * at a given point. Since the mass is well defined both at the core and the
 * surface we use it as an independent variable.
 */
void modelSolveSingle(MODEL *model, double Mtot, double rhos, double us)
{
    double M_mid, R_old, R_new, rhoc_old, rhoc_new, uc_old, uc_new, dm;
    double R_L, R_R, rho_L, rho_R, u_L, u_R;
    double Y_R_old, Y_R_new, Y_rho_old, Y_rho_new, Y_u_old, Y_u_new;
    double dY_R, dY_rho, dY_u;
    double *Delta;
    int bSetModel;

    // This vector stores the suggested changes in R, rhoc and uc
    Delta = (double *) malloc(3*sizeof(double));

    /*
     * As a first step the model is solved for two slighly different initial
     * guesses for the parameters R, rhoc and uc. As initial guesses for R,
     * rhoc, and uc we use:
     *
     * R = (3.0*Mtot/(4.0*M_PI*rho0))
     * rhoc = 1.1*rho_L (obtained from the first integration from the surface)
     * uc = 1.1*uc      (obtained from the first integration from the surface)
     */
    R_old = cbrt(3.0*Mtot/(4.0*M_PI*model->tillMat->rho0));

    /*
     * CR: Try to get a better first guess by using rho_mean rather than rho0.
     */
    R_old = cbrt(3.0*Mtot/(4.0*M_PI*13.16));
    R_old = 1.043;
    R_old = 1.2;

    M_mid = 0.5*Mtot;
    dm = Mtot/model->nTableMax;

    midPtRKIn(model, Mtot, R_old, rhos, us, dm, M_mid, &R_L, &rho_L, &u_L);

    /*
     * Since rho(m=0)=rhoc and u(m=0)=uc are unknown we use the left hand side
     * values at the matching point to get a first guess.
     */
    rhoc_old = 3.0*rho_L;
    uc_old = 2.5*u_L;
    rhoc_old = 19.6054;
    uc_old = 10.2;

    midPtRKOut(model, bSetModel=0, rhoc_old, uc_old, dm, M_mid, &R_R, &rho_R, &u_R);

    /*
     * Y_i = xi_L - xi_R
     */
    Y_R_old = R_L-R_R;
    Y_rho_old = rho_L-rho_R;
    Y_u_old = u_L-u_R;

    /*
     * The second step to initialize the iteration is to slightly perturb the
     * initial values R_old, rhoc_old and uc_old to obtain
     */
/*    R_new = 0.9999*R_old;
    rhoc_new = 0.9999*rhoc_old;
    uc_new = 0.9999*uc_old;
*/

    R_new = 1.0001*R_old;
    rhoc_new = 1.0001*rhoc_old;
    uc_new = 1.0001*uc_old;
    modelSolveToMatchingPoint(model, rhoc_new, uc_new, Mtot, R_new, rhos, us, M_mid, dm, &Y_R_new, &Y_rho_new, &Y_u_new);
    
    /*
     * Now that we have initial values for all the free parameters (R, rhoc,
     * uc) and their correspondinding differences at the matching point we
     * can inerate until the solutions at M_mid are matching
     * (dR_new-dR_old < eps, drho_new-drho_old < eps, du_new-du_old < eps).
     */
    dY_R = Y_R_new - Y_R_old;
    dY_rho = Y_rho_new - Y_rho_old;
    dY_u = Y_u_new - Y_u_old;

    fprintf(stderr,"Starting interation:\n");
    fprintf(stderr,"R_old=  %15.7E rhoc_old=  %15.7E uc_old=  %15.7E\n", R_old, rhoc_old, uc_old);
    fprintf(stderr,"R_new=  %15.7E rhoc_new=  %15.7E uc_new=  %15.7E\n", R_new, rhoc_new, uc_new);
    fprintf(stderr,"dR=  %15.7E drhoc=  %15.7E duc=  %15.7E\n", R_new-R_old, rhoc_new-rhoc_old, uc_new-uc_old);
    fprintf(stderr,"Y_R_old= %15.7E Y_rho_old= %15.7E Y_u_old= %15.7E\n", Y_R_old, Y_rho_old, Y_u_old);
    fprintf(stderr,"Y_R_new= %15.7E Y_rho_new= %15.7E Y_u_new= %15.7E\n", Y_R_new, Y_rho_new, Y_u_new);
    fprintf(stderr,"\n");

    fprintf(stderr,"dY_R= %15.7E dY_rho= %15.7E dY_u= %15.7E\n", dY_R, dY_rho, dY_u);
    fprintf(stderr,"\n");

    while (sqrt(dY_R*dY_R+dY_rho*dY_rho+dY_u*dY_u) > 1e-10)
#define EPS 1e-3
//    while (ddR > EPS && ddrho > EPS && ddu > EPS)
    {
        fprintf(stderr,"Interation:\n");
        fprintf(stderr,"Y_R_old= %15.7E Y_rhoc_old= %15.7E Y_uc_old= %15.7E\n", Y_R_old, Y_rho_old, Y_u_old);
        fprintf(stderr,"Y_R_new= %15.7E Y_rhoc_new= %15.7E Y_uc_new= %15.7E\n", Y_R_new, Y_rho_new, Y_u_new);

        modelGetNewInitialGuess(R_new-R_old, rhoc_new-rhoc_old, uc_new-uc_old, Y_R_old, Y_rho_old, Y_u_old, dY_R, dY_rho, dY_u, Delta);
        fprintf(stderr,"Delta= %15.7E, %15.7E, %15.7E\n", Delta[0], Delta[1], Delta[2]); 

        /*
         * Store the values from the last iteration.
         */
        R_old = R_new;
        rhoc_old = rhoc_new;
        uc_old = uc_new;
        Y_R_old = Y_R_new;
        Y_rho_old = Y_rho_new;
        Y_u_old = Y_u_new;

        /*
         * Update the guess for R, rhoc and uc and solve the model using these
         * new values.
         */
        R_new = R_old+Delta[0];
        rhoc_new = rhoc_old+Delta[1];
        uc_new = uc_old+Delta[2];
        
        modelSolveToMatchingPoint(model, rhoc_new, uc_new, Mtot, R_new, rhos, us, M_mid, dm, &Y_R_new, &Y_rho_new, &Y_u_new);

        fprintf(stderr,"Y_R_old= %15.7E Y_rhoc_old= %15.7E Y_uc_old= %15.7E\n", Y_R_old, Y_rho_old, Y_u_old);
        fprintf(stderr,"Y_R_new= %15.7E Y_rhoc_new= %15.7E Y_uc_new= %15.7E\n", Y_R_new, Y_rho_new, Y_u_new);
        fprintf(stderr,"\n");

        dY_R = Y_R_new - Y_R_old;
        dY_rho = Y_rho_new - Y_rho_old;
        dY_u = Y_u_new - Y_u_old;
    }

    fprintf(stderr,"Final values:\n");
    fprintf(stderr,"R= %15.7E rhoc= %15.7E uc= %15.7E\n", R_new, rhoc_new, uc_new);
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
    int bSetModel;
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

    double R_L, rho_L, u_L, R_R, rho_R, u_R;
    double M_mid;
	// Initialize model
	model = modelInit(iMat);
//	tillInitLookup(model->tillMat);

    dm = mTot/(nStepsMax);
    M_mid = 0.5*mTot;

    /*
     * Test the solver for a single component model.
     */
    modelSolveSingle(model, mTot, rhos, us);
    exit(1);

    /*
     * Test the Gauss-Jordan elimination method.
     */
    GaussJordanTest(3);
    exit(1);
//#if 0
    /*
     * Solve the model from the surface.
     */
    midPtRKIn(model, mTot, R, rhos, us, dm, M_mid, &R_L, &rho_L, &u_L);

    fprintf(stderr, "mTot= %15.7E ", mTot);
    fprintf(stderr, "R= %15.7E ", R);
    fprintf(stderr, "rhos= %15.7E ", rhos);
    fprintf(stderr, "us= %15.7E ", us);

    fprintf(stderr, "\n");

    fprintf(stderr, "rhoc= %15.7E ", rhoc);
    fprintf(stderr, "uc= %15.7E ", uc);
    fprintf(stderr, "\n");

    fprintf(stderr, "M_mid= %15.7E ", M_mid);
    fprintf(stderr, "R_L= %15.7E ", R_L);
    fprintf(stderr, "rho_L= %15.7E ", rho_L);
    fprintf(stderr, "u_L= %15.7E ", u_L);

    fprintf(stderr, "\n");
//#endif

//#if 0
    /*
     * Solve the model from the core.
     */
    //M_mid = 1.0*mTot;
    midPtRKOut(model, bSetModel=0, rhoc, uc, dm/100.0, M_mid, &R_R, &rho_R, &u_R);
    fprintf(stderr, "mTot= %15.7E ", mTot);
    fprintf(stderr, "R= %15.7E ", R);
    fprintf(stderr, "rhos= %15.7E ", rhos);
    fprintf(stderr, "us= %15.7E ", us);

    fprintf(stderr, "\n");

    fprintf(stderr, "rhoc= %15.7E ", rhoc);
    fprintf(stderr, "uc= %15.7E ", uc);
    fprintf(stderr, "\n");

    fprintf(stderr, "M_mid= %15.7E ", M_mid);
    fprintf(stderr, "R_R= %15.7E ", R_R);
    fprintf(stderr, "rho_R= %15.7E ", rho_R);
    fprintf(stderr, "u_R= %15.7E ", u_R);

    fprintf(stderr, "\n");
//#endif
}
