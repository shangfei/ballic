/*
 * This file provides all functions needed to solve the Lane-Emden equations
 * to generate a polytrope.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <assert.h>
#include "lane-emden.h"

MODEL *modelInit(double rhoc, double ucore, double gamma)
{
    /* Initialize the model */
	MODEL *model;
    
    model = malloc(sizeof(MODEL));
    assert(model != NULL);

    model->dKpcUnit = 2.06701e-13;
    model->dMsolUnit = 4.80438e-08;

    model->rhoc = rhoc;
    model->uc = ucore;
    model->gamma = gamma;

    /*
     * gamma = (n+1)/n
     */
    model->n = 1.0/(model->gamma-1.0);

    /*
     * K = (gamma-1)*uc/rhoc^(1/n)
     */
    model->K = (model->gamma-1.0)*model->uc/pow(model->rhoc, 1.0/model->n);

    /*
     * alpha = sqrt((n+1)/(4*M_PI*G)*K*rhoc^((1-n)/n))
     *
     * Remember that G==1 in code units.
     */
    model->alpha = sqrt((model->n+1.0)/(4.0*M_PI)*model->K*pow(model->rhoc, (1.0-model->n)/model->n));

    model->nTableMax = 10000; 

    model->w1 = malloc(model->nTableMax*sizeof(double));
    assert(model->w1 != NULL);

    model->z = malloc(model->nTableMax*sizeof(double));
    assert(model->z != NULL);

    model->M = malloc(model->nTableMax*sizeof(double));
    assert(model->M != NULL);

    model->rho = malloc(model->nTableMax*sizeof(double));
    assert(model->rho != NULL);

    model->u = malloc(model->nTableMax*sizeof(double));
    assert(model->u != NULL);

    model->r = malloc(model->nTableMax*sizeof(double));
    assert(model->r != NULL);

    model->dr =  0.0;
    model->nTable = 0;
    
    return(model);
}

/*
 * Calculate the pressure for an ideal gas,
 */
double Pressure(MODEL *model, double rho, double u)
{
    return ((model->gamma-1.0)*rho*u);
}

/*
 * We use dimensionless variables defined by: r = alpha*xi, rho = rhoc*theta^2.
 *
 * In order to obtain two first order ODEs the variables
 *
 * z = xi, w1 = theta, w2 = dtheta/dz
 *
 * are used.
 */
double dw1dz(MODEL *model, double z, double w1, double w2)
{
    return(w2);
}

double dw2dz(MODEL *model, double z, double w1, double w2)
{
    if (z > 0)
    {
//        return(-pow(w1, model->n) - 2.0/z*w2);
        return(-w1 - 2.0/z*w2);
    } else {
        return(0.0);
    }
}

double dmudz(MODEL *model, double z, double w1, double w2)
{
    return(w1*z*z);
}

/*
 * Calculate the gravitational accelleration due to the enclosed mass M(R).
 */
double CalcGrav(double r, double M)
{
	assert(r>=0.0);
	
	if (r > 0.0)
	{
		// Again assuming G=1
		return(-M/(r*r));
	} else {
		return(0.0);
	}
}

/*
 * rho = rhoc*theta^n
 */
double RhoTheta(MODEL *model, double w1)
{
    return(model->rhoc*pow(w1,model->n));
}

/*
 * u = K/(gamma-1)*rho^((gamma-1)
 */
double UTheta(MODEL *model, double w1)
{
    return(model->K/(model->gamma-1.0)*pow(model->rhoc, 1.0/model->n)*w1);
}

/*
 * r = alpha * z
 */
double RZ(MODEL *model, double z)
{
    return(model->alpha*z);
}

/*
 * M = 4.0*M_PI*rhoc*alpha^3*mu(z)
 */
double MMu(MODEL *model, double mu)
{
    return(4.0*M_PI*model->rhoc*pow(model->alpha, 3.0)*mu);
}

int WriteModeltoFile(MODEL *model, const char *Name)
{
    FILE *fp;
    int i;

	fp = fopen(Name, "w");

    for (i=0; i<model->nTable; i++)
    {
//        fprintf(fp,"%15.7E %15.7E\n", model->z[i], model->w1[i]);
        fprintf(fp,"%15.7E %15.7E %15.7E %15.7E %15.7E %15.7E %15.7E\n", model->z[i], model->w1[i], model->r[i], model->rho[i], model->M[i], model->u[i], CalcGrav(model->r[i], model->M[i]));
    }

    fclose(fp);
    return(0);
}

/*
 * This function solves the Lane-Emden equations using a 2nd order Runge-Kutta
 * with initial values
 *
 * z = 0, w1 = 1, w2 = 0
 *
 * using a step size of h. If bSetModel=1 the model is stored (r, rho and u)
 * and the total mass and radius are returned.
 */
double LESolver(MODEL *model, int bSetModel, double h, double *pR)
{
    FILE *fp;
    double w1 = 1.0;
    double w2 = 0.0;
    double z = 0.0;
    // Mass in dimensionless variables
    double mu = 0.0;
    double k1w1, k1w2, k1mu, k2w1, k2w2, k2mu, x;
    int i;

    if (bSetModel) {
		i = 0;
		model->w1[i] = w1;
		model->z[i] = z;

		model->rho[i] = RhoTheta(model, w1);
		model->M[i] = MMu(model, mu);
		model->u[i] = UTheta(model, w1);
		model->r[i] = RZ(model, z);

		fp = fopen("lane-emden.model", "w");
		assert(fp != NULL);
		fprintf(fp,"%15.7E %15.7E %15.7E %15.7E %15.7E %15.7E %15.7E\n", z, w1, model->r[i], model->rho[i], model->M[i], model->u[i], CalcGrav(model->r[i], model->M[i]));
		++i;
	}

    while (w1 > 0.0)
    {
		/*
         * Midpoint Runga-Kutta (2nd order).
         */
		k1w1 = h*dw1dz(model, z, w1, w2);
		k1w2 = h*dw2dz(model, z, w1, w2);
		k1mu = h*dmudz(model, z, w1, w2);

		k2w1 = h*dw1dz(model, z+0.5*h, w1+0.5*k1w1, w2+0.5*k1w2);
		k2w2 = h*dw2dz(model, z+0.5*h, w1+0.5*k1w1, w2+0.5*k1w2);
		k2mu = h*dmudz(model, z+0.5*h, w1+0.5*k1w1, w2+0.5*k1w2);

		w1 += k2w1;
		w2 += k2w2;
        mu += k2mu;
		z += h;
		
        printf("%g %g %g %g %g %g\n", z, w1, w2, sin(z)/z, 1.0/(z*z)*(cos(z)*z-sin(z)), h);

        if (bSetModel)
        {
	        model->w1[i] = w1;
            model->z[i] = z;

            model->rho[i] = RhoTheta(model, w1);
    		model->M[i] = MMu(model, mu);
    		model->u[i] = UTheta(model, w1);
    		model->r[i] = RZ(model, z);

            fprintf(fp,"%15.7E %15.7E %15.7E %15.7E %15.7E %15.7E %15.7E\n", z, w1, model->r[i], model->rho[i], model->M[i], model->u[i], CalcGrav(model->r[i], model->M[i]));
			++i;
		}
	}

    /*
     * Now do a linear interpolation to rho == 0.0.
     */
    x = -w1/k2w1;

    /*
     * Remember: x = - (y-y_i+1)/(y_i+1-y_1) = -A.
     */
    assert(x <= 0.0);
    z += h*x;
    w1 += k2w1*x;
    w2 += k2w2*x;
	mu += k2mu*x;

	if (bSetModel) {
		--i;

		model->w1[i] = w1;
		model->z[i] = z;

        model->rho[i] = RhoTheta(model, w1);
    	model->M[i] = MMu(model, mu);
    	model->u[i] = UTheta(model, w1);
    	model->r[i] = RZ(model, z);

        fprintf(fp,"%15.7E %15.7E %15.7E %15.7E %15.7E %15.7E %15.7E\n", z, w1, model->r[i], model->rho[i], model->M[i], model->u[i], CalcGrav(model->r[i], model->M[i]));

		fclose(fp);
		++i;
		model->nTable = i;
		model->dr = h;
	}

    *pR = RZ(model, z);
    return(MMu(model, mu));
}

void main(int argc, char **argv)
{
	const int nStepsMax = 10000;
    int bSetModel;
	// Model
    MODEL *model;
	double rhoc, uc, gamma;
	double dz, R, M;

    /*
     * Initialise density and internal energy in the core.
     */
    rhoc = 1.0;
    uc = 1.0;

    /*
     * Solve a model with n=1 (gamma=2).
     */
    gamma = 2.0;
    dz = M_PI/(nStepsMax-1);

    model = modelInit(rhoc, uc, gamma);

    fprintf(stderr,"Model initialised.\n");
    fprintf(stderr,"gamma= %g n= %g\n", model->gamma, model->n);
    fprintf(stderr,"\n");

    fprintf(stderr,"Solving LE equation: dz= %g\n", dz);
    M = LESolver(model, bSetModel=1, dz, &R);

//    M = LESolver(model, bSetModel=0, dz, &R);
    WriteModeltoFile(model, "model.txt");
    
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

}


