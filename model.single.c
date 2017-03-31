/*
 * This file provides all functions needed to generate a single component model
 * of a planet.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <assert.h>
#include "tipsy.h"
#include "tillotson/tillotson.h"
#include "model.single.h"

MODEL *modelInit(double ucore, int iMat) {
    /* Initialize the model */
	MODEL *model;
    
    model = malloc(sizeof(MODEL));
    assert(model != NULL);

    model->dKpcUnit = 2.06701e-13;
    model->dMsolUnit = 4.80438e-08;
	model->tillMat = malloc(sizeof(TILLMATERIAL));

	/*
	** Initialize one material.
	** i=0: Granite
	** i=1: Iron
	** i=2: Basalt
	** i=3: Ice
	*/
	model->tillMat = tillInitMaterial(iMat, model->dKpcUnit, model->dMsolUnit, 100, 100, 50.0, 50.0, 1);

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

    /* model->uFixed = uFixed/model->dErgPerGmUnit; */
    model->uc = ucore;

    model->nTableMax = 10000; 
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
** dudrho depends on the internal energy profile that we choose!
*/
double dudrho(MODEL *model,double rho,double u) {
#ifdef BALLIC_U_POLYTROPIC
	/*
	** We assume a polytropic internal energy profile!
	*/
	// Not implemented yet
	assert(0);
	return(Pressure(model,rho,u)/(rho*rho));
#else
#ifdef BALLIC_U_ISENTROPIC
	/*
	** We assume an isentropic internal energy profile!
	*/
	return(tillPressure(model->tillMat,rho,u)/(rho*rho));
#else
  	fprintf(stderr,"No thermal profile defined when compiled!\n");
	assert(0);
#endif
#endif
}

/*
** Calculate dudrho to solve for the equilibrium model.
*/
double drhodr(MODEL *model,double r,double rho,double M,double u) {
    double dPdrho,dPdu;

	dPdrho=tilldPdrho(model->tillMat, rho, u); // dP/drho at u=const.
	dPdu = tilldPdu(model->tillMat, rho, u);; // dP/du at rho=const.

	/*
	** drho/dr = -G*M*rho/(dPdrho+dPdu*dudrho)
	*/
	assert(r >= 0.0);
	if (r > 0.0) {
	// We assume G=1
		return(-M*rho/(r*r*(dPdrho + dPdu*dudrho(model,rho,u))));
	}
	else {
		return(0.0);
	}
}

double dudr(MODEL *model,double r,double rho,double M,double u) {
	return(dudrho(model,rho,u)*drhodr(model,r,rho,M,u));
}

/*
** This derivative is independent of the model and only involves geometry.
*/
double dMdr(double r,double rho) {
	assert(r >= 0.0);
	return(4.0*M_PI*r*r*rho);
}

const double fact = 1.0;

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
** This function solves the model as an initial value problem with rho_initial = rho and 
** M_initial = 0 at r = 0. This function returns the mass when rho == model->tillMat[i]->rho0.
*/
double midPtRK(MODEL *model,int bSetModel,double rho,double h,double *pR) {
    FILE *fp;
    double M = 0.0;
    double r = 0.0;
	double u = model->uc;
    double k1rho,k1M,k1u,k2rho,k2M,k2u,x;
    int i;

    if (bSetModel) {
		i = 0;
		model->rho[i] = rho;
		model->M[i] = M;
		model->u[i] = u;
		model->r[i] = r;
		fp = fopen("ballic.model","w");
		assert(fp != NULL);
		fprintf(fp,"%g %g %g %g %g\n",r,rho,M,u,CalcGrav(r,M));
		++i;
	}

    while (rho > fact*model->tillMat->rho0) {
		/*
		** Midpoint Runga-Kutta (2nd order).
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

		if (bSetModel) {
			model->rho[i] = rho;
			model->M[i] = M;
			model->u[i] = u;
			model->r[i] = r;
			fprintf(fp,"%g %g %g %g %g\n",r,rho,M,u,CalcGrav(r,M));
			++i;
		}
	}
    /*
    ** Now do a linear interpolation to rho == fact*rho0.
    */
    x = (fact*model->tillMat->rho0 - rho)/k2rho;
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

		fprintf(fp,"%g %g %g %g %g\n",r,rho,M,u,CalcGrav(r,M));
		fclose(fp);
		++i;
		model->nTable = i;
		model->dr = h;
	}

    *pR = r;
    return(M);
    }


double modelSolve(MODEL *model,double M) {
    const int nStepsMax = 10000;
    int bSetModel;
    double rmax;
    double dr,R;
    double a,Ma,b,Mb,c,Mc;

    /*
    ** First estimate the maximum possible radius.
    */
    R = cbrt(3.0*M/(4.0*M_PI*model->tillMat->rho0));
    dr = R/nStepsMax;
    a = 1.01*model->tillMat->rho0; /* starts with 1% larger central density */
    Ma = midPtRK(model,bSetModel=0,a,dr,&R);
    fprintf(stderr,"first Ma:%g R:%g\n",Ma,R);
    b = a;
    Mb = 0.5*M;
    while (Ma > M) {
		b = a;
		Mb = Ma;
		a = 0.5*(model->tillMat->rho0 + a);
		Ma = midPtRK(model,bSetModel=0,a,dr,&R);
	}
    while (Mb < M) {
		b = 2.0*b;
	   	Mb = midPtRK(model,bSetModel=0,b,dr,&R);	
		fprintf(stderr,"first Mb:%g R:%g\n",Mb,R);
	}

	// (CR) Debug
	fprintf(stderr,"Root bracketed.\n");

    /*
    ** Root bracketed by (a,b).
    */
    while (Mb-Ma > 1e-10*Mc) {
	c = 0.5*(a + b);
        Mc = midPtRK(model,bSetModel=0,c,dr,&R);	
	if (Mc < M) {
	    a = c;
	    Ma = Mc;
	    }
	else {
	    b = c;
	    Mb = Mc;
	    }
//	fprintf(stderr,"c:%.10g Mc:%.10g R:%.10g\n",c/model->tillMat[0]->rho0,Mc,R);
	}
    /*
    ** Solve it once more setting up the lookup table.
    */
    fprintf(stderr,"rho_core: %g cv: %g uc: %g (in system units)\n",c,model->tillMat->cv,model->uc);
    Mc = midPtRK(model,bSetModel=1,c,dr,&R);
    model->R = R;
    return c;
    }
