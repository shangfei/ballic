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
	int i;
    /* Initialize the model */
	MODEL *model;
    
    model = (MODEL *) malloc(sizeof(MODEL));
    assert(model != NULL);

    model->dKpcUnit = 2.06701e-13;
    model->dMsolUnit = 4.80438e-08;
	
	/* Hard coded */
	model->nMat = 5;
	assert(model->nMat == TILL_N_MATERIAL_MAX);

	model->tillMat = (TILLMATERIAL **) malloc(model->nMat*sizeof(TILLMATERIAL *));
	assert(model->tillMat != NULL);
	
	model->iMatOrder = (int *) malloc(model->nMat*sizeof(int));
	assert(model->iMatOrder != NULL);
	
	model->fM = (double *) malloc(model->nMat*sizeof(double));
	assert(model->fM != NULL);

	/* Hard coded too */
	// Earth like planet
//#if 0
	model->iMatOrder[0] = IRON;
	model->iMatOrder[1] = GRANITE;
//	model->iMatOrder[1] = BASALT;
	// It might be better so save M in model and use only mass fractions in fM
	model->fM[0] = 0.3*M;
	model->fM[1] = 0.7*M;
//#endif

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
	fprintf(stderr,"Initializing model:\n");
	fprintf(stderr,"Mtot=%g ucore=%g\n",M,ucore);
	fprintf(stderr,"iMatOrder[%i, %i]\n",model->iMatOrder[0],model->iMatOrder[1]);
	fprintf(stderr,"fM[%g, %g]\n",model->fM[0],model->fM[1]);

	for (i=0; i<model->nMat; i++)
	{
		/*
		** Initialize one material.
		** i=0: Granite
		** i=1: Iron
		** i=2: Basalt
		*/
		model->tillMat[i] = tillInitMaterial(i, model->dKpcUnit, model->dMsolUnit, 100, 100, 100.0, 1200.0, 1);

		// Debug information

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

		/* Generate the look up table needed for tillColdULookup(). */

		tillInitLookup(model->tillMat[i]);
	}

    /* model->uFixed = uFixed/model->dErgPerGmUnit; */
    model->uc = ucore;

    model->nTableMax = 10000; 
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
	return(tillPressure(model->tillMat[iMat],rho,u)/(rho*rho));
}

/*
** Calculate dudrho to solve for the equilibrium model.
*/
double drhodr(MODEL *model, int iMat, double r,double rho,double M,double u)
{
    double dPdrho,dPdu;

	dPdrho=tilldPdrho(model->tillMat[iMat], rho, u); // dP/drho at u=const.
	dPdu = tilldPdu(model->tillMat[iMat], rho, u);; // dP/du at rho=const.

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
     * Check, if the given initial copnditions are physical.
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
 * Write the lookup table to ballic.model.
 */ 
void modelWriteToFile(MODEL *model)
{
	FILE *fp;
	int i;

	fp = fopen("ballic.model","w");
	assert(fp != NULL);
	
	fprintf(fp,"#R  rho  M  u  mat tillPressure  tillTempRhoU\n");
	for (i=0; i<model->nTable;i++)
	{
		fprintf(fp,"%g %g %g %g %i %g %g\n",model->r[i],model->rho[i],model->M[i],model->u[i],model->mat[i],
			tillPressure(model->tillMat[model->mat[i]], model->rho[i], model->u[i]),
			tillTempRhoU(model->tillMat[model->mat[i]], model->rho[i], model->u[i]));
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
    double **M_array;
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

	// Do a nRho x nU grid
	nRho = 1000;
	nU = 1000;

    M_array = modelMatrixAlloc(nRho, nU);

    // Initialize model (M and uc are not used!)
	model = modelInit(mTot, 0.0);

    nLayer = 2;

	// Open output file
    fp = fopen("solve_model_array_multi.txt","w");
    assert(fp != NULL);

	M = mTot;
	/*
	 * First estimate the maximum possible radius.
	 */
	R = cbrt(3.0*M/(4.0*M_PI*model->tillMat[nLayer-1]->rho0));
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

#ifdef _OPENMP
    fprintf(stderr,"Code was compiled with OpenMP.\n");
//    omp_set_num_threads(16);
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
#endif
#pragma omp for private(i, j, rho, u, M, R, us, nThreads)
    for (i=0; i<nRho; i++)
	{
		rho = rhomin + (rhomax-rhomin)/nRho*i;

		for (j=0; j<nU; j++)
		{
			u = umin + (umax-umin)/nU*j;
			// Solve the model for rho, M and u
            M = modelSolveTwoComponent(model, 0, rho, u, dr, &R, &us);
//            printf("i= %i j= %i M= %g\n", i, j, M);
            M_array[i][j] = M;
        }
    }
}
    fprintf(stderr, "All models solved.\n");

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
			if (M_array[i][j] > 0.0 && fabs(M_array[i][j]-mTot) < 0.01*mTot && fabs(us-u_desired) < 5e-0)
			{
				fprintf(fp,"%3i",1);
			} else if (M_array[i][j] < 0.0) {
				fprintf(fp,"%3i",-1);
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
}
