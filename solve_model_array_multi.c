#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <assert.h>
#include "omp.h"
#include "tipsy.h"
#include "tillotson/tillotson.h"
#include "solve_model_array_multi.h"

MODEL *modelInit(double M,double ucore)
{
    const double KBOLTZ = 1.38e-16;      /* bolzman constant in cgs */
    const double MHYDR = 1.67e-24;       /* mass of hydrogen atom in grams */
    const double MSOLG = 1.99e33;        /* solar mass in grams */
    const double GCGS = 6.67e-8;         /* G in cgs */
    const double KPCCM = 3.085678e21;    /* kiloparsec in centimeters */
	int i;
    /* Initialize the model */
	MODEL *model;
    
    model = (MODEL *) malloc(sizeof(MODEL));
    assert(model != NULL);

    model->dKpcUnit = 2.06701e-13;
    model->dMsolUnit = 4.80438e-08;
	
	/* Hard coded */
	model->nMat = 6;
	assert(model->nMat == TILL_N_MATERIAL_MAX);

	model->tillMat = (TILLMATERIAL **) malloc(model->nMat*sizeof(TILLMATERIAL *));
	assert(model->tillMat != NULL);
	
	model->iMatOrder = (int *) malloc(model->nMat*sizeof(int));
	assert(model->iMatOrder != NULL);
	
	model->fM = (double *) malloc(model->nMat*sizeof(double));
	assert(model->fM != NULL);

	/* Hard coded too */
	// Earth like planet
#if 0
	model->iMatOrder[0] = IRON;
	model->iMatOrder[1] = GRANITE;
//	model->iMatOrder[1] = BASALT;
	// It might be better so save M in model and use only mass fractions in fM
	model->fM[0] = 0.3*M;
	model->fM[1] = 0.7*M;
#endif

#if 0
	// Try 10% rocky core and 90% ice mantle
	model->iMatOrder[0] = GRANITE;
	model->iMatOrder[1] = ICE;
	model->fM[0] = 0.1*M;
	model->fM[1] = 0.9*M;
#endif
#if 0
	// Try 50% rocky core and 50% ice mantle
	model->iMatOrder[0] = GRANITE;
	model->iMatOrder[1] = ICE;
	model->fM[0] = 0.5*M;
	model->fM[1] = 0.5*M;
#endif
#if 0
	// Try 99% rocky core and 1% water ocean
	model->iMatOrder[0] = GRANITE;
	model->iMatOrder[1] = WATER;
	model->fM[0] = 0.99*M;
	model->fM[1] = 0.01*M;
#endif

	// Earth like planet with an atmosphere
#if 0
	model->iMatOrder[0] = IRON;
	model->iMatOrder[1] = GRANITE;
	model->iMatOrder[2] = IDEALGAS;
	// It might be better so save M in model and use only mass fractions in fM
	model->fM[0] = 0.3*M;
	model->fM[1] = 0.6*M;
	model->fM[2] = 0.1*M;
#endif

	// Ice giant (10% rock core, 75.19% ice mantle and 14.81% atmosphere)
//#if 0
	model->iMatOrder[0] = GRANITE;
	model->iMatOrder[1] = ICE;
	model->iMatOrder[2] = IDEALGAS;
	// It might be better so save M in model and use only mass fractions in fM
	model->fM[0] = 0.1*M;
	model->fM[1] = 0.7519*M;
	model->fM[2] = 0.1481*M;
//#endif

	fprintf(stderr,"Initializing model:\n");
	fprintf(stderr,"Mtot=%g ucore=%g\n",M,ucore);
	fprintf(stderr,"iMatOrder[%i, %i]\n",model->iMatOrder[0],model->iMatOrder[1]);
	fprintf(stderr,"fM[%g, %g]\n",model->fM[0],model->fM[1]);

	for (i=0; i<model->nMat; i++)
	{
		/*
         * Initialize one material:
         * 
         * i=0: Ideal gas
         * i=1: Granite
         * i=2: Iron
         * i=3: Basalt
         * ...
         */
		model->tillMat[i] = tillInitMaterial(i, model->dKpcUnit, model->dMsolUnit, 100, 100, 100.0, 1200.0, 1);

		// Debug information
        if (model->tillMat[i]->iMaterial == IDEALGAS)
        {
            fprintf(stderr,"\n");	
    		fprintf(stderr,"Material: %i\n",i);	
    		fprintf(stderr,"dConstGamma: %g\n", model->tillMat[i]->dConstGamma);
    		fprintf(stderr,"dMeanMolMass: %g\n", model->tillMat[i]->dMeanMolMass);
    		fprintf(stderr,"rho0: %g\n", model->tillMat[i]->rho0);
    		fprintf(stderr,"cv: %g\n", model->tillMat[i]->cv);
    		fprintf(stderr,"\n");
        } else {
            fprintf(stderr,"\n");	
    		fprintf(stderr,"Material: %i\n",i);	
    		fprintf(stderr,"a: %g\n", model->tillMat[i]->a);
    		fprintf(stderr,"b: %g\n", model->tillMat[i]->b);
    		fprintf(stderr,"A: %g\n", model->tillMat[i]->A);
    		fprintf(stderr,"B: %g\n", model->tillMat[i]->B);
    		fprintf(stderr,"rho0: %g\n", model->tillMat[i]->rho0);
    		fprintf(stderr,"u0: %g\n", model->tillMat[i]->u0);
    		fprintf(stderr,"us: %g\n", model->tillMat[i]->us);
    		fprintf(stderr,"us2: %g\n", model->tillMat[i]->us2);
    		fprintf(stderr,"alpha: %g\n", model->tillMat[i]->alpha);
    		fprintf(stderr,"beta: %g\n", model->tillMat[i]->beta);
    		fprintf(stderr,"cv: %g\n", model->tillMat[i]->cv);
    		fprintf(stderr,"\n");
        }

		/* Generate the look up table needed for tillColdULookup(). */
        if (model->tillMat[i]->iMaterial != IDEALGAS) tillInitLookup(model->tillMat[i]);
	}

    /*
     * For an ideal gas tillInitLookup() does not work and screws up cv!.
     */
    model->uc = ucore;

    model->nTableMax = 100000; 
    model->M = (double *) malloc(model->nTableMax*sizeof(double));
    assert(model->M != NULL);
    model->rho =(double *)  malloc(model->nTableMax*sizeof(double));
    assert(model->rho != NULL);
    model->u =(double *)  malloc(model->nTableMax*sizeof(double));
    assert(model->u != NULL);
    model->r =(double *) malloc(model->nTableMax*sizeof(double));
    assert(model->r != NULL);
    model->mat = (int *) malloc(model->nTableMax*sizeof(int));
    assert(model->mat != NULL);
    model->dr =  0.0;
    model->nTable = 0;
    
	return(model);
}

/*
** dudrho depends on the internal energy profile that we choose!
*/
double dudrho(MODEL *model,int iMat,double rho,double u)
{
	return(eosPressure(model->tillMat[iMat],rho,u)/(rho*rho));
}

/*
** Calculate dudrho to solve for the equilibrium model.
*/
double drhodr(MODEL *model, int iMat, double r,double rho,double M,double u)
{
    double dPdrho,dPdu;

	dPdrho = eosdPdrho(model->tillMat[iMat], rho, u); // dP/drho at u=const.
	dPdu = eosdPdu(model->tillMat[iMat], rho, u);; // dP/du at rho=const.

	/*
	** drho/dr = -G*M*rho/(dPdrho+dPdu*dudrho)
	*/
	assert(r >= 0.0);
	if (r > 0.0) {
		// We assume G=1
		return(-M*rho/(r*r*(dPdrho + dPdu*dudrho(model,iMat,rho,u))));
	}
	else {
		return(0.0);
	}
}

double dudr(MODEL *model,int iMat,double r,double rho,double M,double u)
{
  return(dudrho(model,iMat,rho,u)*drhodr(model,iMat,r,rho,M,u));
}

/*
** This derivative is independent of the model and only involves geometry.
*/
double dMdr(double r,double rho)
{
	assert(r >= 0.0);
	return(4.0*M_PI*r*r*rho);
}

/*
** Solve for rho2 and u2 using the b.c. P1=P2 and T1=T2.
*/
void modelSolveBC(MODEL *model, double *prho, double *pu, int iMat1, int iMat2)
{
	double P1, T1, P2, T2;
    double rho1,u1,rho2,u2;
    double rhoa,ua,rhob,ub;

	assert(prho != NULL);
	assert(pu != NULL);

	rho1 = *prho;
	u1 = *pu;

	P1 = tillPressure(model->tillMat[iMat1], rho1, u1);
	T1 = tillTempRhoU(model->tillMat[iMat1], rho1, u1);

	/*
	** We use rho1 as an upper limit for rho2 assuming that the denser component is in the inner shell.
	*/
	rhoa = rho1;
	ua = tillURhoTemp(model->tillMat[iMat2],rhoa,T1);

	rhob = 0.0;
	ub = tillURhoTemp(model->tillMat[iMat2],rhob,T1);

	fprintf(stderr,"\n");
	fprintf(stderr,"*****************************************************************\n");
	fprintf(stderr,"modelSolveBC:\n");
	fprintf(stderr,"iMat1=%i iMat2=%i\n",iMat1,iMat2);
	fprintf(stderr,"rho1=%g u1=%g P1=%g T1=%g\n",rho1,u1,P1,T1);
	fprintf(stderr,"rhoa=%g ua=%g Pa=%g Ta=%g\n",rhoa,ua,tillPressure(model->tillMat[iMat2], rhoa, ua),tillTempRhoU(model->tillMat[iMat2], rhoa, ua));
   	fprintf(stderr,"rhob=%g ub=%g Pb=%g Tb=%g\n",rhob,ub,tillPressure(model->tillMat[iMat2], rhob, ub),tillTempRhoU(model->tillMat[iMat2], rhob, ub));
	
	tillSolveBC(model->tillMat[iMat1],model->tillMat[iMat2],rho1,u1,&rho2,&u2);
	fprintf(stderr,"modelSolveBC: rho1: %g, u1: %g, rho2:%g, u2:%g\n",rho1,u1,rho2,u2);
	fprintf(stderr,"*****************************************************************\n");
	fprintf(stderr,"\n");

	/*
	** Return values.
	*/
	*prho = rho2;
	*pu = u2; 
//	assert(0);
}

/*
 * This function integrates the structure equations with b.c.:
 *
 * rho_initial=rho1, u_initial=u1 and M_initial=M1
 *
 * until either M=M2 or (if bLastLayer=1) rho=rho0. The final values for rho
 * and u are returned (so the original value will be overwritten!!). The
 * parameter h sets the stepsize for the RK2 algorithm and for bSetModel=1 the
 * results are saved in the lookup table (starting from index pIndex). If
 * bLastLayer=1 then the algorithm enforces rho(r=R)=rho0.
 */
int modelSolveComponent(MODEL *model, int iMat, int bSetModel, int bLastLayer, int *pIndex, double h, double *prho, double *pu,double *pM, double M_component, double *pR)
{
	double r, rho, M, u;
    double k1rho,k1M,k1u,k2rho,k2M,k2u,x;
    int i, iRet;

    /*
     * Initial values.
     */
    rho = *prho;
    u = *pu;
    M = *pM;
    r = *pR;

    if (bSetModel)
    {
        i = *pIndex;
    }

    iRet = 0;

	/*
     * Output some diagnostic information.
     */
	fprintf(stderr,"\n");
	fprintf(stderr,"******************************************************************\n");
	fprintf(stderr,"modelSolveComponent (inital values):\n");
	fprintf(stderr,"iMat: %i, Index: %i, h: %g\n",iMat, i, h);
	fprintf(stderr,"rho: %g, u: %g, M:%g, r:%g\n",rho, u, M, r);
	fprintf(stderr,"bSetModel: %i, bLastLayer: %i\n", bSetModel, bLastLayer);
	fprintf(stderr,"******************************************************************\n");
	fprintf(stderr,"\n");


	if (bSetModel)
    {
		model->rho[i] = rho;
		model->M[i] = M;
		model->u[i] = u;
		model->r[i] = r;
		model->mat[i] = iMat;
		++i;
	}

    /*
     * Check, if the given initial conditions are physical.
     */
    if (tillIsBelowColdCurve(model->tillMat[iMat], rho, u))
    {
        iRet = 1;
        return(iRet);
    }

    if (rho <= model->tillMat[iMat]->rho0)
    {
        iRet = 1;
        return(iRet);
    }

	if (bLastLayer != 1)
	{
        /*
         * For the inner layers integrate until M == M_component and assert that
         * rho >= rho0.
         */
		while (M < M_component)
        {
			/*
             * Midpoint Runga-Kutta (2nd order).
             */
			k1rho = h*drhodr(model, iMat, r, rho, M, u);
			k1M = h*dMdr(r, rho);
			k1u = h*dudr(model, iMat, r, rho, M, u);

			k2rho = h*drhodr(model, iMat, r+0.5*h, rho+0.5*k1rho, M+0.5*k1M, u+0.5*k1u);
			k2M = h*dMdr(r+0.5*h, rho+0.5*k1rho);
			k2u = h*dudr(model, iMat, r+0.5*h, rho+0.5*k1rho, M+0.5*k1M, u+0.5*k1u);

			rho += k2rho;
			M += k2M;
			u += k2u;
			r += h;

            /*
             * Assert that rho >= rho0.
             */
            if (rho <= model->tillMat[iMat]->rho0 && M < M_component)
            {
                iRet = 1;
                return(iRet);
            }

            printf("%15.7E %15.7E %15.7E %15.7E\n", r, rho, M, u);

			if (bSetModel)
            {
				model->rho[i] = rho;
				model->M[i] = M;
				model->u[i] = u;
				model->r[i] = r;
				model->mat[i] = iMat;
				++i;
		    }
		}
        
        /*
         * Now do a linear interpolation to M == M_component.
         */
        x = (M_component - M)/k2M;
//        printf(stderr,"M2=%g, M=%g, x=%g\n",Mc,M,x);
        assert(x <= 0.0);
        r += h*x;
        M += k2M*x;
        rho += k2rho*x;
        u += k2u*x;
//        fprintf(stderr,"After correction: M2=%g, M=%g, x=%g\n",Mc,M,x);

        if (bSetModel) {
            --i;
            model->M[i] = M;
            model->r[i] = r;
            model->rho[i] = rho;
            model->u[i] = u;
            model->mat[i] = iMat;
            ++i;
        }
	} else {
        /*
         * In case of the most outer layer integrate until rho == rho0.
         */
        while (rho > model->tillMat[iMat]->rho0)
        {
            /*
             * Midpoint Runga-Kutta (2nd order).
             */
			k1rho = h*drhodr(model, iMat, r, rho, M, u);
			k1M = h*dMdr(r, rho);
			k1u = h*dudr(model, iMat, r, rho, M, u);

			k2rho = h*drhodr(model, iMat, r+0.5*h, rho+0.5*k1rho, M+0.5*k1M, u+0.5*k1u);
			k2M = h*dMdr(r+0.5*h, rho+0.5*k1rho);
			k2u = h*dudr(model, iMat, r+0.5*h, rho+0.5*k1rho, M+0.5*k1M, u+0.5*k1u);

			rho += k2rho;
			M += k2M;
			u += k2u;
			r += h;

            printf("%15.7E %15.7E %15.7E %15.7E\n", r, rho, M, u);

			if (bSetModel)
            {
				model->rho[i] = rho;
				model->M[i] = M;
				model->u[i] = u;
				model->r[i] = r;
				model->mat[i] = iMat;
				++i;
		    }
		}

        /*
         * Now do a linear interpolation to rho == rho0.
         */
        x = (model->tillMat[iMat]->rho0 - rho)/k2rho;
        assert(x <= 0.0);
        
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
            model->mat[i] = iMat;
            ++i;
        }
    } // if (bSetModel != 1)


    /*
     * Return values.
     */
    *pR = r;
    *prho = rho;
    *pM = M;
    *pu = u;

    if (bSetModel)
    {
        *pIndex = i;
    }

    return (iRet);
}

/*
 * The same as modelSolveComponent() but it uses a 4th order Runge-Kutta to
 * integrate the structure equations.
 */
int modelSolveComponentRK4(MODEL *model, int iMat, int bSetModel, int bLastLayer, int *pIndex, double h, double *prho, double *pu,double *pM, double M_component, double *pR)
{
	double r, rho, M, u;
    double k1rho, k1M, k1u, k2rho, k2M, k2u, x;
    double k3rho, k3M, k3u, k4rho, k4M, k4u;
    int i, iRet;

    /*
     * Initial values.
     */
    rho = *prho;
    u = *pu;
    M = *pM;
    r = *pR;

    if (bSetModel)
    {
        i = *pIndex;
    }

    iRet = 0;

	/*
     * Output some diagnostic information.
     */
	fprintf(stderr,"\n");
	fprintf(stderr,"******************************************************************\n");
	fprintf(stderr,"modelSolveComponent (inital values):\n");
	fprintf(stderr,"iMat: %i, Index: %i, h: %g\n",iMat, i, h);
	fprintf(stderr,"rho: %g, u: %g, M:%g, r:%g\n",rho, u, M, r);
	fprintf(stderr,"bSetModel: %i, bLastLayer: %i\n", bSetModel, bLastLayer);
	fprintf(stderr,"******************************************************************\n");
	fprintf(stderr,"\n");


	if (bSetModel)
    {
		model->rho[i] = rho;
		model->M[i] = M;
		model->u[i] = u;
		model->r[i] = r;
		model->mat[i] = iMat;
		++i;
	}

    /*
     * Check, if the given initial conditions are physical.
     */
    if (tillIsBelowColdCurve(model->tillMat[iMat], rho, u))
    {
        iRet = 1;
        return(iRet);
    }

    if (rho <= model->tillMat[iMat]->rho0)
    {
        iRet = 1;
        return(iRet);
    }

	if (bLastLayer != 1)
	{
        /*
         * For the inner layers integrate until M == M_component and assert that
         * rho >= rho0.
         */
		while (M < M_component)
        {
			/*
             * 4th order Runge-Kutta.
             */
			k1rho = h*drhodr(model, iMat, r, rho, M, u);
			k1M = h*dMdr(r, rho);
			k1u = h*dudr(model, iMat, r, rho, M, u);

			k2rho = h*drhodr(model, iMat, r+0.5*h, rho+0.5*k1rho, M+0.5*k1M, u+0.5*k1u);
			k2M = h*dMdr(r+0.5*h, rho+0.5*k1rho);
			k2u = h*dudr(model, iMat, r+0.5*h, rho+0.5*k1rho, M+0.5*k1M, u+0.5*k1u);

			k3rho = h*drhodr(model, iMat, r+0.5*h, rho+0.5*k2rho, M+0.5*k2M, u+0.5*k2u);
			k3M = h*dMdr(r+0.5*h, rho+0.5*k2rho);
			k3u = h*dudr(model, iMat, r+0.5*h, rho+0.5*k2rho, M+0.5*k2M, u+0.5*k2u);

			k4rho = h*drhodr(model, iMat, r+h, rho+k3rho, M+k3M, u+k3u);
			k4M = h*dMdr(r+h, rho+k3rho);
			k4u = h*dudr(model, iMat, r+h, rho+k3rho, M+k3M, u+k3u);

            rho += k1rho/6.0+k2rho/3.0+k3rho/3.0+k4rho/6.0;
			M += k1M/6.0+k2M/3.0+k3M/3.0+k4M/6.0;
			u += k1u/6.0+k2u/3.0+k3u/3.0+k4u/6.0;
			r += h;

            /*
             * Assert that rho >= rho0.
             */
            if (rho <= model->tillMat[iMat]->rho0 && M < M_component)
            {
                iRet = 1;
                return(iRet);
            }

            printf("%15.7E %15.7E %15.7E %15.7E\n", r, rho, M, u);

			if (bSetModel)
            {
				model->rho[i] = rho;
				model->M[i] = M;
				model->u[i] = u;
				model->r[i] = r;
				model->mat[i] = iMat;
				++i;
		    }
		}
        
        /*
         * Now do a linear interpolation to M == M_component.
         */
        x = (M_component - M)/k2M;
//        printf(stderr,"M2=%g, M=%g, x=%g\n",Mc,M,x);
        assert(x <= 0.0);
        r += h*x;
        M += k2M*x;
        rho += k2rho*x;
        u += k2u*x;
//        fprintf(stderr,"After correction: M2=%g, M=%g, x=%g\n",Mc,M,x);

        if (bSetModel) {
            --i;
            model->M[i] = M;
            model->r[i] = r;
            model->rho[i] = rho;
            model->u[i] = u;
            model->mat[i] = iMat;
            ++i;
        }
	} else {
        /*
         * In case of the most outer layer integrate until rho == rho0.
         */
        while (rho > model->tillMat[iMat]->rho0)
        {
            /*
             * 4th order Runge-Kutta.
             */
			k1rho = h*drhodr(model, iMat, r, rho, M, u);
			k1M = h*dMdr(r, rho);
			k1u = h*dudr(model, iMat, r, rho, M, u);

			k2rho = h*drhodr(model, iMat, r+0.5*h, rho+0.5*k1rho, M+0.5*k1M, u+0.5*k1u);
			k2M = h*dMdr(r+0.5*h, rho+0.5*k1rho);
			k2u = h*dudr(model, iMat, r+0.5*h, rho+0.5*k1rho, M+0.5*k1M, u+0.5*k1u);

			k3rho = h*drhodr(model, iMat, r+0.5*h, rho+0.5*k2rho, M+0.5*k2M, u+0.5*k2u);
			k3M = h*dMdr(r+0.5*h, rho+0.5*k2rho);
			k3u = h*dudr(model, iMat, r+0.5*h, rho+0.5*k2rho, M+0.5*k2M, u+0.5*k2u);

			k4rho = h*drhodr(model, iMat, r+h, rho+k3rho, M+k3M, u+k3u);
			k4M = h*dMdr(r+h, rho+k3rho);
			k4u = h*dudr(model, iMat, r+h, rho+k3rho, M+k3M, u+k3u);

            rho += k1rho/6.0+k2rho/3.0+k3rho/3.0+k4rho/6.0;
			M += k1M/6.0+k2M/3.0+k3M/3.0+k4M/6.0;
			u += k1u/6.0+k2u/3.0+k3u/3.0+k4u/6.0;
			r += h;

            printf("%15.7E %15.7E %15.7E %15.7E\n", r, rho, M, u);

			if (bSetModel)
            {
				model->rho[i] = rho;
				model->M[i] = M;
				model->u[i] = u;
				model->r[i] = r;
				model->mat[i] = iMat;
				++i;
		    }
		}

        /*
         * Now do a linear interpolation to rho == rho0.
         */
        x = (model->tillMat[iMat]->rho0 - rho)/k2rho;
        assert(x <= 0.0);
        
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
            model->mat[i] = iMat;
            ++i;
        }
    } // if (bSetModel != 1)


    /*
     * Return values.
     */
    *pR = r;
    *prho = rho;
    *pM = M;
    *pu = u;

    if (bSetModel)
    {
        *pIndex = i;
    }

    return (iRet);
}
#if 0
/*
** This function integrates the ODEs with b.c. rho_initial=rho1, u_initial=u1
** and M_initial=M1 until M=M2. The final values for rho and u are returned
** (so the original value will be overwritten!!). The parameter h sets the
** stepsize for the RK2 algorithm and for bSetModel=1 the results are saved
** in the lookup table (starting from index pIndex). If bLastLayer=1 then
** the algorithm enforces rho(r=R)=rho0.
*/
void modelSolveComponent(MODEL *model,int iMat,int bSetModel,int bLastLayer,int *pIndex,double h,double *prho1,double *pu1,double *pM1,double M2,double *pR)
{
    FILE *fp;
	// Set inital values rho1, u1, M1
	double rho=*prho1;
	double u = *pu1;
    double M = *pM1;
	double r = *pR;

    double k1rho,k1M,k1u,k2rho,k2M,k2u,x;

	// (CR) Debug information
	fprintf(stderr,"\n");
	fprintf(stderr,"******************************************************************\n");
	fprintf(stderr,"modelSolveComponent (inital values):\n");
	fprintf(stderr,"iMat: %i, Index: %i, h: %g\n",iMat,*pIndex,h);
	fprintf(stderr,"rho: %g, u: %g, M:%g, r:%g\n",rho,u,M,r);
	fprintf(stderr,"bSetModel: %i, bLastLayer: %i\n",bSetModel,bLastLayer);
	fprintf(stderr,"******************************************************************\n");
	fprintf(stderr,"\n");

	if (bSetModel) {
		model->rho[*pIndex] = rho;
		model->M[*pIndex] = M;
		model->u[*pIndex] = u;
		model->r[*pIndex] = r;
		model->mat[*pIndex] = iMat;
		++*pIndex;
	}

	if (bLastLayer != 1)
	{
		/* Integrate from M1 to M2 for the inner layers. */
		while (M < M2) {
			/*
			** Midpoint Runga-Kutta (2nd order).
			*/
			k1rho = h*drhodr(model,iMat,r,rho,M,u);
			k1M = h*dMdr(r,rho);
			k1u = h*dudr(model,iMat,r,rho,M,u);

			k2rho = h*drhodr(model,iMat,r+0.5*h,rho+0.5*k1rho,M+0.5*k1M,u+0.5*k1u);
			k2M = h*dMdr(r+0.5*h,rho+0.5*k1rho);
			k2u = h*dudr(model,iMat,r+0.5*h,rho+0.5*k1rho,M+0.5*k1M,u+0.5*k1u);

			rho += k2rho;
			M += k2M;
			u += k2u;
			r += h;
		
			if (bSetModel) {
				model->rho[*pIndex] = rho;
				model->M[*pIndex] = M;
				model->u[*pIndex] = u;
				model->r[*pIndex] = r;
				model->mat[*pIndex] = iMat;
	//			fprintf(fp,"%g %g %g %g\n",r,rho,M,u);
				++*pIndex;
		    }
		}

		/*
		** Now do a linear interpolation to M == M2.
		*/
		x = (M2 - M)/k2M;
		fprintf(stderr,"M2=%g, M=%g, x=%g\n",M2,M,x);
		assert(x <= 0.0);
		r += h*x;
		M += k2M*x;
		rho += k2rho*x;
		u += k2u*x;

		if (bSetModel) {
			--*pIndex;
			model->M[*pIndex] = M;
			model->r[*pIndex] = r;
			model->rho[*pIndex] = rho;
			model->u[*pIndex] = u;
			model->mat[*pIndex] = iMat;
			++*pIndex;
		}
		/* Make sure that the material is in the condensed state */
		assert(rho >= model->tillMat[iMat]->rho0);
	} else {
		// For the last layer we integrate until rho == rho0. */
		while (rho > model->tillMat[iMat]->rho0) {
			/*
			** Midpoint Runga-Kutta (2nd order).
			*/
			k1rho = h*drhodr(model,iMat,r,rho,M,u);
			k1M = h*dMdr(r,rho);
			k1u = h*dudr(model,iMat,r,rho,M,u);

			k2rho = h*drhodr(model,iMat,r+0.5*h,rho+0.5*k1rho,M+0.5*k1M,u+0.5*k1u);
			k2M = h*dMdr(r+0.5*h,rho+0.5*k1rho);
			k2u = h*dudr(model,iMat,r+0.5*h,rho+0.5*k1rho,M+0.5*k1M,u+0.5*k1u);

			rho += k2rho;
			M += k2M;
			u += k2u;
			r += h;
	
//			fprintf(stderr,"r=%g rho=%g M=%g u=%g k1rho=%g k1M=%g k1u=%g k2rho=%g k2M=%g k2u=%g\n",
//							r, rho, M, u, k1rho, k1M, k1u, k2rho, k2M, k2u);	
			if (bSetModel) {
				model->M[*pIndex] = M;
				model->r[*pIndex] = r;
				model->rho[*pIndex] = rho;
				model->u[*pIndex] = u;
				model->mat[*pIndex] = iMat;
				++*pIndex;
			}
		}
		
		/*
		** Now do a linear interpolation to rho == rho0.
		*/
		x = (model->tillMat[iMat]->rho0 - rho)/k2rho;
		fprintf(stderr,"iMat=%i, rho0=%g, rho=%g, x=%g\n",iMat,model->tillMat[iMat]->rho0,rho,x);
		fprintf(stderr,"rho0-rho=%g, k2rho0=%g\n",model->tillMat[iMat]->rho0-rho,k2rho);
		assert(x <= 0.0);
		r += h*x;
		M += k2M*x;
		rho += k2rho*x;
		u += k2u*x;
		
		if (bSetModel) {
			--*pIndex;
			model->M[*pIndex] = M;
			model->r[*pIndex] = r;
			model->rho[*pIndex] = rho;
			model->u[*pIndex] = u;
			model->mat[*pIndex] = iMat;
			++*pIndex;
		}	
	}

	// Return values
	*prho1 = rho;
	*pu1 = u;
	*pM1 = M;
    *pR = r;
}
#endif
/*
 * This function integrates the ODEs for a two component model with b.c.
 * rho_initial=rho, u_initial=u until rho(r=R)=rho0. It returns the total
 * mass of the model. The parameter h sets the stepsize for the RK2 algorithm
 * and for bSetModel=1 the results are saved.
 */
double modelSolveTwoComponent(MODEL *model,int bSetModel,double rho,double u,double h,double *pR,double *us)
{
//    FILE *fp;
	// Set inital values for M1 and R
    double M = 0.0;
	double r = 0.0;
	double Mc = model->fM[0];
	int iMat;
    double k1rho,k1M,k1u,k2rho,k2M,k2u,x;
	int i = 0;

	// (CR) Debug information
	fprintf(stderr,"\n");
	fprintf(stderr,"******************************************************************\n");
	fprintf(stderr,"modelSolveTwoComponent (inital values):\n");
	fprintf(stderr,"rho: %g, u: %g, M:%g, r:%g\n",rho,u,M,r);
	fprintf(stderr,"bSetModel: %i, h: %g\n",bSetModel,h);
	fprintf(stderr,"******************************************************************\n");
	fprintf(stderr,"\n");

    /* Set iMat to the material for the core (usually iron). */
	iMat = model->iMatOrder[0];

    /*
     * Check, if the given initial conditions are physical.
     */
    if (tillIsBelowColdCurve(model->tillMat[iMat], rho, u))
    {
        return(-1.0);
    }

    if (rho <= model->tillMat[iMat]->rho0)
    {
        return(-1.0);
    }

	if (bSetModel) {
		model->rho[i] = rho;
		model->M[i] = M;
		model->u[i] = u;
		model->r[i] = r;
		model->mat[i] = iMat;
		++i;
	}

	/*
     * First integrate the core until M == Mcore. But at the same time we have
     * to make sure, that rho >= rho0.
     */
	while (M < Mc) {
		/*
		** Midpoint Runga-Kutta (2nd order).
		*/
		k1rho = h*drhodr(model,iMat,r,rho,M,u);
		k1M = h*dMdr(r,rho);
		k1u = h*dudr(model,iMat,r,rho,M,u);

		k2rho = h*drhodr(model,iMat,r+0.5*h,rho+0.5*k1rho,M+0.5*k1M,u+0.5*k1u);
		k2M = h*dMdr(r+0.5*h,rho+0.5*k1rho);
		k2u = h*dudr(model,iMat,r+0.5*h,rho+0.5*k1rho,M+0.5*k1M,u+0.5*k1u);

		rho += k2rho;
		M += k2M;
		u += k2u;
		r += h;

        /*
         * Check, if the initial density is too low, so that the density falls
         * below rho0 before we reach the desired mass.
         */
        if (rho <= model->tillMat[iMat]->rho0 && M > Mc)
        {
            return(-1.0);
        }

		if (bSetModel) {
			model->rho[i] = rho;
			model->M[i] = M;
			model->u[i] = u;
			model->r[i] = r;
			model->mat[i] = iMat;
			++i;
		}
	}

	/*
	** Now do a linear interpolation to M == Mcore.
	*/
	x = (Mc - M)/k2M;
	fprintf(stderr,"M2=%g, M=%g, x=%g\n",Mc,M,x);
	assert(x <= 0.0);
	r += h*x;
	M += k2M*x;
	rho += k2rho*x;
	u += k2u*x;
	fprintf(stderr,"After correction: M2=%g, M=%g, x=%g\n",Mc,M,x);

	if (bSetModel) {
		--i;
		model->M[i] = M;
		model->r[i] = r;
		model->rho[i] = rho;
		model->u[i] = u;
		model->mat[i] = iMat;
		++i;
	}
	
	fprintf(stderr,"\n");
	fprintf(stderr,"******************************************************************\n");
	fprintf(stderr,"modelSolveTwoComponent (core/mantle boundary):\n");
	fprintf(stderr,"Core: iMat=%i rho=%g u=%g M=%g r=%g P=%g T=%g\n",iMat,rho,u,M,r,tillPressure(model->tillMat[iMat],rho,u),tillTempRhoU(model->tillMat[iMat],rho,u));

	/*
	** Now calculate rho and u for the mantle using b.c. T=const and P=const.
	*/
	tillSolveBC(model->tillMat[model->iMatOrder[0]],model->tillMat[model->iMatOrder[1]],rho,u,&rho,&u);

	iMat = model->iMatOrder[1];

	fprintf(stderr,"Mantle: iMat=%i rho=%g u=%g M=%g r=%g P=%g T=%g\n",iMat,rho,u,M,r,tillPressure(model->tillMat[iMat],rho,u),tillTempRhoU(model->tillMat[iMat],rho,u));
	fprintf(stderr,"******************************************************************\n");
	fprintf(stderr,"\n");
	
	/*
	** The integrate the mantle until rho(r=R)=rho0.
	*/
	while (rho > model->tillMat[iMat]->rho0) {
		/*
		** Midpoint Runga-Kutta (2nd order).
		*/
		k1rho = h*drhodr(model,iMat,r,rho,M,u);
		k1M = h*dMdr(r,rho);
		k1u = h*dudr(model,iMat,r,rho,M,u);

		k2rho = h*drhodr(model,iMat,r+0.5*h,rho+0.5*k1rho,M+0.5*k1M,u+0.5*k1u);
		k2M = h*dMdr(r+0.5*h,rho+0.5*k1rho);
		k2u = h*dudr(model,iMat,r+0.5*h,rho+0.5*k1rho,M+0.5*k1M,u+0.5*k1u);

		rho += k2rho;
		M += k2M;
		u += k2u;
		r += h;
	
//		fprintf(stderr,"r=%g rho=%g M=%g u=%g k1rho=%g k1M=%g k1u=%g k2rho=%g k2M=%g k2u=%g\n",
//							r, rho, M, u, k1rho, k1M, k1u, k2rho, k2M, k2u);	
		if (bSetModel) {
			model->M[i] = M;
			model->r[i] = r;
			model->rho[i] = rho;
			model->u[i] = u;
			model->mat[i] = iMat;
			++i;
			model->nTable = i;
			model->dr = h;
		}
	}

	/*
	** Now do a linear interpolation to rho == rho0.
	*/
	x = (model->tillMat[iMat]->rho0 - rho)/k2rho;
	fprintf(stderr,"iMat=%i, rho0=%g, rho=%g, x=%g\n",iMat,model->tillMat[iMat]->rho0,rho,x);
	fprintf(stderr,"rho0-rho=%g, k2rho0=%g\n",model->tillMat[iMat]->rho0-rho,k2rho);
	assert(x <= 0.0);
	r += h*x;
	M += k2M*x;
	rho += k2rho*x;
	u += k2u*x;
		
	if (bSetModel) {
		--i;
		model->M[i] = M;
		model->r[i] = r;
		model->rho[i] = rho;
		model->u[i] = u;
		model->mat[i] = iMat;
		++i;
	}	

	// Return values
    if (us != NULL) *us = u;

    *pR = r;
	return(M);
}

/*
 * Calculate the density and internal energy of the atmosphere assuming
 * P=const and T=const for an ideal gas.
 */
void modelSolveBCAtmosphere(MODEL *model, TILLMATERIAL *mat, double rho1, double u1, double *rho2, double *u2)
{
    double Ps, Ts, rho_atm, u_atm;

    assert(mat->iMaterial > 0);

    /*
     * Calculate pressure and temperature at the surface (not including the atmosphere).
     */
    Ps = tillPressure(mat, rho1, u1);
	Ts = tillTempRhoU(mat, rho1, u1);

    assert(model->tillMat[0]->iMaterial == IDEALGAS);

    /*
     * Use cv and gamma from the ideal gas EOS.
     */
    u_atm = Ts*model->tillMat[0]->cv;

//    fprintf(stderr,"iMat= %3i cv= %g.\n", model->tillMat[0]->iMaterial, model->tillMat[0]->cv);
//    exit(1);
    rho_atm = Ps/((model->tillMat[0]->dConstGamma-1.0)*u_atm);

#if 0
    printf("\n");
    printf("rho= %g u= %g dConstGamma= %g (dConstGamma-1.0)= %g P= %g T= %g\n", rho_atm, u_atm, model->tillMat[0]->dConstGamma, model->tillMat[0]->dConstGamma-1.0, (model->tillMat[0]->dConstGamma-1.0)*rho_atm*u_atm, u_atm/model->tillMat[0]->cv);
    printf("\n");

    printf("mantle:   rho= %g u= %g P= %g T= %g\n", rho1, u1, Ps, Ts);
    printf("envelope: rho= %g u= %g P= %g T= %g\n", rho_atm, u_atm, (model->tillMat[0]->dConstGamma-1.0)*rho_atm*u_atm, u_atm/model->tillMat[0]->cv);
//    exit(1);
#endif
    *rho2 = rho_atm;
    *u2 = u_atm;
}

/*
 * This function integrates the ODEs for a three component model with b.c.
 *
 * rho_initial=rho, u_initial=u, M_initial=0
 *
 * For the inner layers it integrates until M=M_i (M_i: desired mass of
 * layer i) and fails if rho < rho0_i. The last layer is solved until rho=rho_s
 * where rho_s is the desired density at the surface. It returns the total mass
 * and internal energy at the surface of the model.
 *
 * The parameter h sets the stepsize for the RK2 algorithm and for bSetModel=1
 * the results are saved.
 */
double modelSolveThreeComponent(MODEL *model,int bSetModel,double rho,double u,double h,double *pR,double *us)
{
//    FILE *fp;
	// Set inital values for M1 and R
    double M = 0.0;
	double r = 0.0;
    // The cumulative mass of the core
	double Mc = model->fM[0];
    // The cumulative mass of the mantle (Mcore+Mmantle)
	double Mm = Mc + model->fM[1];
	int iMat;
    double k1rho,k1M,k1u,k2rho,k2M,k2u,x;
	int i = 0;

	// (CR) Debug information
	fprintf(stderr,"\n");
	fprintf(stderr,"******************************************************************\n");
	fprintf(stderr,"modelSolveThreeComponent (inital values):\n");
	fprintf(stderr,"rho: %g, u: %g, M:%g, r:%g rhos:%g\n",rho,u,M,r,model->tillMat[0]->rho0);
	fprintf(stderr,"bSetModel: %i, h: %g\n",bSetModel,h);
	fprintf(stderr,"iMat: %i %i %i\n",model->iMatOrder[0],model->iMatOrder[1],model->iMatOrder[2]);
    fprintf(stderr,"Mc: %g Mm: %g Matm: %g\n",model->fM[0],model->fM[1],model->fM[2]);
    fprintf(stderr,"******************************************************************\n");
	fprintf(stderr,"\n");

    /* Set iMat to the material for the core (usually iron). */
	iMat = model->iMatOrder[0];

    /*
     * Check, if the given initial conditions are physical.
     */
    if (tillIsBelowColdCurve(model->tillMat[iMat], rho, u))
    {
        return(-1.0);
    }

    if (rho <= model->tillMat[iMat]->rho0)
    {
        return(-1.0);
    }

	if (bSetModel) {
		model->rho[i] = rho;
		model->M[i] = M;
		model->u[i] = u;
		model->r[i] = r;
		model->mat[i] = iMat;
		++i;
	}

	/*
     * Step 1: Integrate the core until M == Mcore. But at the same time we
     * have to make sure, that rho >= rho0.
     */
	while (M < Mc) {
		/*
		** Midpoint Runga-Kutta (2nd order).
		*/
		k1rho = h*drhodr(model,iMat,r,rho,M,u);
		k1M = h*dMdr(r,rho);
		k1u = h*dudr(model,iMat,r,rho,M,u);

		k2rho = h*drhodr(model,iMat,r+0.5*h,rho+0.5*k1rho,M+0.5*k1M,u+0.5*k1u);
		k2M = h*dMdr(r+0.5*h,rho+0.5*k1rho);
		k2u = h*dudr(model,iMat,r+0.5*h,rho+0.5*k1rho,M+0.5*k1M,u+0.5*k1u);

		rho += k2rho;
		M += k2M;
		u += k2u;
		r += h;

        /*
         * Check, if the initial density is too low, so that the density falls
         * below rho0 before we reach the desired mass.
         */
        if (rho <= model->tillMat[iMat]->rho0 && M < Mc)
        {
            return(-1.0);
        }

//      printf("%15.7E %15.7E %15.7E %15.7E %3i %15.7E\n", r, rho, M, u, iMat, model->tillMat[iMat]->rho0);

		if (bSetModel) {
			model->rho[i] = rho;
			model->M[i] = M;
			model->u[i] = u;
			model->r[i] = r;
			model->mat[i] = iMat;
			++i;
		}
	}

	/*
	** Now do a linear interpolation to M == Mcore.
	*/
	x = (Mc - M)/k2M;
	fprintf(stderr,"M2=%g, M=%g, x=%g\n",Mc,M,x);
	assert(x <= 0.0);
	r += h*x;
	M += k2M*x;
	rho += k2rho*x;
	u += k2u*x;
	fprintf(stderr,"After correction: M2=%g, M=%g, x=%g\n",Mc,M,x);

	if (bSetModel) {
		--i;
		model->M[i] = M;
		model->r[i] = r;
		model->rho[i] = rho;
		model->u[i] = u;
		model->mat[i] = iMat;
		++i;
	}
	
	fprintf(stderr,"\n");
	fprintf(stderr,"******************************************************************\n");
	fprintf(stderr,"modelSolveThreeComponent (core/mantle boundary):\n");
	fprintf(stderr,"Core: iMat=%i rho=%g u=%g M=%g r=%g P=%g T=%g\n",iMat,rho,u,M,r,tillPressure(model->tillMat[iMat],rho,u),tillTempRhoU(model->tillMat[iMat],rho,u));

	/*
     * Step 2: Calculate rho and u for the mantle using b.c. T=const and P=const.
     */
	tillSolveBC(model->tillMat[model->iMatOrder[0]], model->tillMat[model->iMatOrder[1]], rho, u, &rho, &u);

	iMat = model->iMatOrder[1];

	fprintf(stderr,"Mantle: iMat=%i rho=%g u=%g M=%g r=%g P=%g T=%g\n",iMat,rho,u,M,r,tillPressure(model->tillMat[iMat],rho,u),tillTempRhoU(model->tillMat[iMat],rho,u));
	fprintf(stderr,"******************************************************************\n");
	fprintf(stderr,"\n");

    if (tillIsBelowColdCurve(model->tillMat[iMat], rho, u))
    {
        return(-1.0);
    }

	/*
     * Step 3: Integrate the mantle until M == Mmantle. But at the same time we
     * have to make sure, that rho >= rho0.
     */
	while (M < Mm) {
		/*
		** Midpoint Runga-Kutta (2nd order).
		*/
		k1rho = h*drhodr(model,iMat,r,rho,M,u);
		k1M = h*dMdr(r,rho);
		k1u = h*dudr(model,iMat,r,rho,M,u);

		k2rho = h*drhodr(model,iMat,r+0.5*h,rho+0.5*k1rho,M+0.5*k1M,u+0.5*k1u);
		k2M = h*dMdr(r+0.5*h,rho+0.5*k1rho);
		k2u = h*dudr(model,iMat,r+0.5*h,rho+0.5*k1rho,M+0.5*k1M,u+0.5*k1u);

		rho += k2rho;
		M += k2M;
		u += k2u;
		r += h;

        /*
         * Check, if the initial density is too low, so that the density falls
         * below rho0 before we reach the desired mass.
         */
        if (rho <= model->tillMat[iMat]->rho0 && M < Mm)
        {
            return(-1.0);
        }

//      printf("%15.7E %15.7E %15.7E %15.7E %3i %15.7E\n", r, rho, M, u, iMat, model->tillMat[iMat]->rho0);

		if (bSetModel) {
			model->rho[i] = rho;
			model->M[i] = M;
			model->u[i] = u;
			model->r[i] = r;
			model->mat[i] = iMat;
			++i;
		}
	}

	/*
	** Now do a linear interpolation to M == Mmantle.
	*/
	x = (Mm - M)/k2M;
	fprintf(stderr,"M2=%g, M=%g, x=%g\n",Mm,M,x);
	assert(x <= 0.0);
	r += h*x;
	M += k2M*x;
	rho += k2rho*x;
	u += k2u*x;
	fprintf(stderr,"After correction: M2=%g, M=%g, x=%g\n",Mm,M,x);

	fprintf(stderr,"rho= %g u= %g\n",rho,u);
	
    if (bSetModel) {
		--i;
		model->M[i] = M;
		model->r[i] = r;
		model->rho[i] = rho;
		model->u[i] = u;
		model->mat[i] = iMat;
		++i;
	}

    /*
     * Step 4: Calculate rho_atm and u_atm from rho_surface and u_surface.
     */
	fprintf(stderr,"\n");
	fprintf(stderr,"******************************************************************\n");
	fprintf(stderr,"modelSolveThreeComponent (mantle/atmosphere boundary):\n");
	fprintf(stderr,"Mantle: iMat=%i rho=%g u=%g M=%g r=%g P=%g T=%g\n",iMat,rho,u,M,r,tillPressure(model->tillMat[iMat],rho,u),tillTempRhoU(model->tillMat[iMat],rho,u));

#if 0
    /*
     * Some code to check, why we have T<0 in some cases.
     */
    if (rho < model->tillMat[iMat]->rho0)
        printf("rho= %15.7E rho0= %15.7E\n", rho, model->tillMat[iMat]->rho0);

    if (tillTempRhoU(model->tillMat[iMat],rho,u) < 0.0)
        printf("T= %g rho= %g u= %g uc=%g\n", tillTempRhoU(model->tillMat[iMat],rho,u), rho, u, tillColdULookup(model->tillMat[iMat], rho));

    if (tillIsBelowColdCurve(model->tillMat[iMat], rho, u))
        printf("Value rho= %g u= %g is below the cold curve (uc= %g).\n", rho, u, tillColdULookup(model->tillMat[iMat], rho)); 
#endif

    double rho_old = rho;
    modelSolveBCAtmosphere(model, model->tillMat[iMat], rho, u, &rho, &u);

//#if 0
    if (rho < model->tillMat[model->iMatOrder[2]]->rho0)
    {
		if (bSetModel) {
			model->M[i] = M;
			model->r[i] = r;
			model->rho[i] = rho;
			model->u[i] = u;
			model->mat[i] = iMat;
			model->nTable = i;
		}
    }

    if (rho_old < rho)
        return(-2);
//#endif

    /*
     * Set the material for the atmosphere.
     */
	iMat = model->iMatOrder[2];

	fprintf(stderr,"Atmosphere: iMat=%i rho=%g u=%g M=%g r=%g P=%g T=%g\n",iMat,rho,u,M,r,eosPressure(model->tillMat[iMat],rho,u),eosTempRhoU(model->tillMat[iMat],rho,u));
	fprintf(stderr,"******************************************************************\n");
	fprintf(stderr,"\n");

	/*
     * Step 5: Integrate the atmosphere until rho(r=R)=rho_s.
     */
	while (rho > model->tillMat[model->iMatOrder[2]]->rho0) {
		/*
		** Midpoint Runga-Kutta (2nd order).
		*/
		k1rho = h*drhodr(model,iMat,r,rho,M,u);
		k1M = h*dMdr(r,rho);
		k1u = h*dudr(model,iMat,r,rho,M,u);

		k2rho = h*drhodr(model,iMat,r+0.5*h,rho+0.5*k1rho,M+0.5*k1M,u+0.5*k1u);
		k2M = h*dMdr(r+0.5*h,rho+0.5*k1rho);
		k2u = h*dudr(model,iMat,r+0.5*h,rho+0.5*k1rho,M+0.5*k1M,u+0.5*k1u);

		rho += k2rho;
		M += k2M;
		u += k2u;
		r += h;

        /*
         * Fail, if the mass is more than three times the total desired mass.
         */
        if (M > 3.0*(model->fM[0]+model->fM[1]+model->fM[2]))
        {
            return(-1.0);
        }

//      printf("%15.7E %15.7E %15.7E %15.7E %3i %15.7E\n", r, rho, M, u, iMat, model->tillMat[iMat]->rho0);
//		fprintf(stderr,"r=%g rho=%g M=%g u=%g k1rho=%g k1M=%g k1u=%g k2rho=%g k2M=%g k2u=%g\n",
//							r, rho, M, u, k1rho, k1M, k1u, k2rho, k2M, k2u);	
		if (bSetModel) {
			model->M[i] = M;
			model->r[i] = r;
			model->rho[i] = rho;
			model->u[i] = u;
			model->mat[i] = iMat;
			++i;
			model->nTable = i;
			model->dr = h;
		}
	}

	/*
	** Now do a linear interpolation to rho == rhos.
	*/
	x = (model->tillMat[model->iMatOrder[2]]->rho0 - rho)/k2rho;
	fprintf(stderr,"iMat=%i, rhos=%g, rho=%g, x=%g\n",iMat,model->tillMat[model->iMatOrder[2]]->rho0,rho,x);
	fprintf(stderr,"rhos-rho=%g, k2rho=%g\n",model->tillMat[model->iMatOrder[2]]->rho0,k2rho);
	assert(x <= 0.0);
	r += h*x;
	M += k2M*x;
	rho += k2rho*x;
	u += k2u*x;
		
	if (bSetModel) {
		--i;
		model->M[i] = M;
		model->r[i] = r;
		model->rho[i] = rho;
		model->u[i] = u;
		model->mat[i] = iMat;
		++i;
	}	

	// Return values
    if (us != NULL) *us = u;

    *pR = r;
	return(M);
}
/*
 * Write the lookup table to ballic.model.
 */ 
void modelWriteToFile(MODEL *model)
{
	FILE *fp;
	int i;

    fprintf(stderr, "nTable= %i nTableMax= %i\n", model->nTable, model->nTableMax);
	fp = fopen("ballic.model", "w");
	assert(fp != NULL);

    assert(model->nTable <= model->nTableMax);

	fprintf(fp,"#R  rho  M  u  mat tillPressure  tillTempRhoU\n");
	for (i=0; i<model->nTable; i++)
	{
//        printf("%i\n", i);
		fprintf(fp,"%g %g %g %g %i %g %g\n",model->r[i],model->rho[i],model->M[i],model->u[i],model->mat[i],
			eosPressure(model->tillMat[model->mat[i]], model->rho[i], model->u[i]),
			eosTempRhoU(model->tillMat[model->mat[i]], model->rho[i], model->u[i]));
//		fprintf(fp,"%15.7E %15.7E %15.7E %15.7E %i\n",model->r[i],model->rho[i],model->M[i],model->u[i],model->mat[i]);
	}
	fclose(fp);
}

/*
 * We allocate an array as described in the Numerical Recipes.
 */	
double **modelMatrixAlloc(int nRow, int nCol)
{
	double **matrix;
	int i;

	matrix = (double **) malloc(nRow*sizeof(double*));
	matrix[0] = (double *) malloc(nRow*nCol*sizeof(double));

	assert(matrix != NULL);

	/* Set a pointer to each row. */
	for (i=1; i<nRow; i++)
	{
		matrix[i] = matrix[i-1]+nCol;
	}

	return (matrix);
}

int modelWriteArrayToFile(FILE *fp,double **array, int nRho, int nU)
{
    int i, j;

    assert(fp != NULL);

	for (i=0; i<nRho; i++)
	{
		for (j=0; j<nU; j++)
		{
            fprintf(fp, "%15.7E", array[i][j]);
		}

		fprintf(fp,"\n");
	}

    return(0);
}

/*
 * Build an array of models for different initial rhoc and uc and mark those,
 * that have the proper mass and surface temperature.
 */
void main(int argc, char **argv) {
	const int nStepsMax = 10000;
	// Model
    MODEL *model;
	double rho, rhomin, rhomax;
	double u, umin, umax;
	int nRho, nU;
	double dr,R, M, mTot, us, u_desired;
    double **M_array, **rhoc_array, **uc_array, **us_array;
	int nLayer, nThreads;
    FILE *fp;
    int i,j;

	// Hard code all values
	mTot = 62.366;			// About 1 Earth mass
    u_desired = 1.0;        // Desired internal energy at the surface
//	rhomin = 1.05*7.33;		// 1.05*rho0
//	rhomax = 3.0*7.33;		// 3.0*rho0
//	umin = 1e-3;
//	umax = 22.0;			// umax > us2
	
	rhomin = 7.33;		// 1.05*rho0
    rhomax = 100;		// 3.0*rho0
	umin = 0.0;
	umax = 100.0;			// umax > us2

    /*
     * A model for Uranus.
     */
	mTot = 844.393;		// 13.5 Earth masses

#if 0
    /*
     * Zoom in for the two component models.
     */
	rhomin = 21.0;
    rhomax = 54.0;
	umin = 0.0;
	umax = 25.0;
#endif
#if 0
    /*
     * Zoom in for the three component models.
     */
	rhomin = 30.0;
    rhomax = 50.0;
	umin = 5.0;
	umax = 25.0;
#endif

#if 0
    /*
     * Zoom in for the Uranus (~14 ME) three component models.
     */
	rhomin = 30.0;
    rhomax = 45.0;
	umin = 60.0;
	umax = 75.0;
#endif

//#if 0
    /*
     * Zoom in for a Uranus (13.5 ME) model:
     *  
     *  951   641   4.42650000e+01   6.96150000e+01   8.44393180e+02   2.90598500e-02
     */
    
	rhomin = 42.0;
    rhomax = 47.0;
	umin = 69.0;
	umax = 74.0;
//#endif

    // Do a nRho x nU grid
	nRho = 100;
	nU = 100;

    M_array = modelMatrixAlloc(nRho, nU);
    rhoc_array = modelMatrixAlloc(nRho, nU);
    uc_array = modelMatrixAlloc(nRho, nU);
    us_array = modelMatrixAlloc(nRho, nU);

    // Initialize model (M and uc are not used!)
	model = modelInit(mTot, 0.0);

    nLayer = 3;

	// Open output file
    fp = fopen("solve_model_array_multi.txt","w");
    assert(fp != NULL);

	M = mTot;
	/*
	 * First estimate the maximum possible radius.
	 */
//	R = cbrt(3.0*M/(4.0*M_PI*model->tillMat[nLayer-1]->rho0));
	R = cbrt(3.0*M/(4.0*M_PI*model->tillMat[model->iMatOrder[nLayer-2]]->rho0));
    dr = R/nStepsMax;
	fprintf(stderr,"First guess: R= %g dr= %g iMat= %i\n", R, dr, model->iMatOrder[nLayer-2]);

#if 0
	/*
	 * Some special case to debug the code.
	 */
	rho = 19.6054;
	u = 10.2;
	M = midPtRK(model,1,rho,u,dr,&R);
#endif

#if 0
//M = modelSolveThreeComponent(model, 0, 60.0, 80.0, dr, &R, &us);
//M = modelSolveThreeComponent(model, 0, 31.23, 8.25, dr, &R, &us);
//M = modelSolveThreeComponent(model, 0, 27.27, 22.75, dr, &R, &us);
//  35   9   3.9764500E+01   9.0000000E+00   6.5255092E+01   2.4232745E-01
//    M = modelSolveThreeComponent(model, 1, 3.9764500E+01, 9.0000000E+00, dr, &R, &us);
    /*
     * Model with M=62.366 and three components:
     * Core:     Iron
     * Mantle:   Granite
     * Envelope: Ideal gas
     *
     * 589 176   4.1780000E+01   8.5200000E+00   6.2366582E+01   6.8733496E-02
     */
    M = modelSolveThreeComponent(model, 1, 4.1780000E+01, 8.5200000E+00, dr, &R, &us);
    fprintf(stderr, "Writing model to file...\n");
    modelWriteToFile(model);
    exit(1);
#endif

//#if 0
    /*
     * Test modelSolveComponent() by comparing the results to a known solution.
     */
    M = 62.366;
	R = cbrt(3.0*M/(4.0*M_PI*model->tillMat[GRANITE]->rho0));
	dr = R/10000.0;

    R = 0.0;
    rho = 19.6054;
	u = 10.2;
    M = 0.0;

    fprintf(stderr, "Running modelSolveComponent() with parameters: r= %g, rho=%g, M=%g and u=%g\n", R, rho, M, u);

    // With bLastLayer = 0
//    modelSolveComponent(model, GRANITE, 0, 0, NULL, dr, &rho, &u, &M, 62.366, &R);
//    modelSolveComponentRK4(model, GRANITE, 0, 0, NULL, dr, &rho, &u, &M, 62.366, &R);

    // With bLastLayer = 1
//    modelSolveComponent(model, GRANITE, 0, 1, NULL, dr, &rho, &u, &M, 62.366, &R);
    modelSolveComponentRK4(model, GRANITE, 0, 1, NULL, dr, &rho, &u, &M, 62.366, &R);

    fprintf(stderr,"R= %15.7E rho= %15.7E M= %15.7E u= %15.7E\n");
    exit(1);
//#endif

#if 0
    /*
     * Investigate why we sometimes get weird return values for M.
     */
    //M= -6976.14, -6976 rho=  4.4950000E+01 u=  6.9450000E+01
    
    M = modelSolveThreeComponent(model, 1, 4.4950000E+01, 6.9450000E+01, dr, &R, &us);
    printf("nTable= %i\n", model->nTable);
    modelWriteToFile(model);
    printf("M= %g, %3i rho=%15.7E u=%15.7E\n", M, (int) M, rho, u);

    // Now solve just the ice independently
    modelSolveComponent(model, ICE, 0, 0, NULL, dr, &rho, &u, &M, 
            
            double *prho, double *pu,double *pM, double M_component, double *pR)
    exit(1);
#endif

#ifdef _OPENMP
    fprintf(stderr,"Code was compiled with OpenMP.\n");
//    omp_set_num_threads(16);
#else
    fprintf(stderr,"Code was compiled for single core.\n");
//    assert(0);
#endif

    /*
     * Solve the array of models and store the results in M_array.
     */
//#pragma omp parallel private(i, j, rho, u, M, R, us, nThreads)
#pragma omp parallel private(i, j, rho, u, M, R, us, nThreads)
{
#ifdef _OPENMP
    nThreads = omp_get_num_threads();
    fprintf(stderr, "nThreads= %i.\n", nThreads);
#else
    fprintf(stderr, "Not compiled with OpenMP.\n");
#endif


#pragma omp for private(i, j, rho, u, M, R, us, nThreads)
    for (i=0; i<nRho; i++)
	{
		rho = rhomin + (rhomax-rhomin)/nRho*i;

		for (j=0; j<nU; j++)
		{
			u = umin + (umax-umin)/nU*j;
            rhoc_array[i][j] = rho;
            uc_array[i][j] = u;
            us = 0.0;

			// Solve the model for rho, M and u
//            M = modelSolveTwoComponent(model, 0, rho, u, dr, &R, &us);
            M = modelSolveThreeComponent(model, 0, rho, u, dr, &R, &us);
            us_array[i][j] = us;

//            printf("i= %i j= %i M= %g\n", i, j, M);

            if (M < 0.0) printf("M= %g, %3i rho=%15.7E u=%15.7E\n", M, (int) M, rho, u);
            M_array[i][j] = M;
        }
    }
}
    fprintf(stderr, "All models solved.\n");

    printf("#  i   j            rhoc              uc               M              us\n");
    /*
     * Print the results.
     */
	for (i=0; i<nRho; i++)
	{
		for (j=0; j<nU; j++)
		{
#if 0
            /*
             * Mark models that agree up to 1% with the desired mass.
             */
			if (M_array[i][j] > 0.0 && fabs(M_array[i][j]-mTot) < 0.01*mTot)
			{
				fprintf(fp,"%3i",1);
			} else if (M_array[i][j] < 0.0) {
				fprintf(fp,"%3i",-1);
			} else {
				fprintf(fp,"%3i",0);
			}
#endif
//#if 0
//			if (M_array[i][j] > 0.0 && fabs(M_array[i][j]-mTot) < 0.01*mTot && fabs(us_array[i][j]-u_desired) < 1e-1)

			if (M_array[i][j] > 0.0 && fabs(M_array[i][j]-mTot) < 0.001*mTot)
            {
				fprintf(fp,"%3i",1);
//                printf("i= %i j= %i: rhoc= %15.7E uc= %15.7E M= %15.7E us= %15.7E\n", i, j, rhoc_array[i][j], uc_array[i][j], M, us_array[i][j]);
                printf(" %3i %3i %15.7E %15.7E %15.7E %15.7E\n", i, j, rhoc_array[i][j], uc_array[i][j], M_array[i][j], us_array[i][j]);
			} else if (M_array[i][j] < 0.0) {
//				fprintf(fp,"%3i",-1);
				fprintf(fp,"%3i", (int) M_array[i][j]);
            } else {
				fprintf(fp,"%3i",0);
			}
//#endif
#if 0
            /*
             * Store (M-Mtot)/Mtot.
             */
			fprintf(fp,"%15.7E", (M_array[i][j]-mTot)/mTot);
#endif
		}

		fprintf(fp,"\n");
	}
	fclose(fp);

    /*
     * Write the density and internal energy to a file too.
     */
    fp = fopen("solve_model_array_multi.rhoc","w");
    modelWriteArrayToFile(fp, rhoc_array,  nRho, nU);
    fclose(fp);

    fp = fopen("solve_model_array_multi.uc","w");
    modelWriteArrayToFile(fp, uc_array,  nRho, nU);
    fclose(fp);

    fp = fopen("solve_model_array_multi.us","w");
    modelWriteArrayToFile(fp, us_array,  nRho, nU);
    fclose(fp);

    fp = fopen("solve_model_array_multi.mass","w");
    modelWriteArrayToFile(fp, M_array,  nRho, nU);
    fclose(fp);
}
