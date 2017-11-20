#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <assert.h>
#include "tipsy.h"
#include "tillotson/tillotson.h"
//#include "ballic.h"

#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) > (B) ? (B) : (A))

// Isentropic thermal profile
#define BALLIC_U_ISENTROPIC

#define MATERIAL 0

typedef struct icosa_struct {
    float R[180];
    float v[36];
    } ICOSA;

ICOSA *icosaInit(void) {
    ICOSA *ctx;
    ctx = malloc(sizeof(ICOSA));
    assert(ctx != NULL);
    compute_matrices_(&ctx->R);
    compute_corners_(&ctx->v);
    return(ctx);
    }

void icosaPix2Vec(ICOSA *ctx,int i,int resolution,double *vec) {
    float v[3];
    pixel2vector_(&i,&resolution,&ctx->R,&ctx->v,v);
    vec[0] = v[0];
    vec[1] = v[1];
    vec[2] = v[2];
    }


/* -----------------------------------------------------------------------------
 *
 *  Copyright (C) 1997-2005 Krzysztof M. Gorski, Eric Hivon, 
 *                          Benjamin D. Wandelt, Anthony J. Banday, 
 *                          Matthias Bartelmann, 
 *                          Reza Ansari & Kenneth M. Ganga 
 *
 *
 *  This file is part of HEALPix.
 *
 *  HEALPix is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  HEALPix is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with HEALPix; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *  For more information about HEALPix see http://healpix.jpl.nasa.gov
 *
 *----------------------------------------------------------------------------- */
void pix2vec_ring( long nside, long ipix, double *vec) {
  /*
    c=======================================================================
    c     gives theta and phi corresponding to pixel ipix (RING) 
    c     for a parameter nside
    c=======================================================================
  */
  
  int nl2, nl4, npix, ncap, iring, iphi, ip, ipix1;
  double  fact1, fact2, fodd, hip, fihip;
  double PI=M_PI;
  double z, sz, phi;
  //      PARAMETER (pi     = 3.1415926535897932384626434d0)
  //      parameter (ns_max = 8192) ! 2^13 : largest nside available
  
  int ns_max=8192;
  
  if( nside<1 || nside>ns_max ) {
    fprintf(stderr, "%s (%d): nside out of range: %ld\n", __FILE__, __LINE__, nside);
    exit(0);
  }
  npix = 12*nside*nside;      // ! total number of points
  if( ipix<0 || ipix>npix-1 ) {
    fprintf(stderr, "%s (%d): ipix out of range: %ld\n", __FILE__, __LINE__, ipix);
    exit(0);
  }
  
  ipix1 = ipix + 1; // in {1, npix}
  nl2 = 2*nside;
  nl4 = 4*nside;
  ncap = 2*nside*(nside-1);// ! points in each polar cap, =0 for nside =1
  fact1 = 1.5*nside;
  fact2 = 3.0*nside*nside;
  
  if( ipix1 <= ncap ) {  //! North Polar cap -------------
    
    hip   = ipix1/2.;
    fihip = floor(hip);
    iring = (int)floor( sqrt( hip - sqrt(fihip) ) ) + 1;// ! counted from North pole
    iphi  = ipix1 - 2*iring*(iring - 1);
    
    z = 1. - iring*iring / fact2 ;
    phi   = (1.*iphi - 0.5) * PI/(2.*iring);
  }
  else if( ipix1 <= nl2*(5*nside+1) ) {//then ! Equatorial region ------
    
    ip    = ipix1 - ncap - 1;
    iring = (int)floor( ip / nl4 ) + nside;// ! counted from North pole
    iphi  = (int)fmod(ip,nl4) + 1;
    
    fodd  = 0.5 * (1 + fmod((double)(iring+nside),2));//  ! 1 if iring+nside is odd, 1/2 otherwise
    z = (nl2 - iring) / fact1;
    phi   = (1.*iphi - fodd) * PI /(2.*nside);
  }
  else {//! South Polar cap -----------------------------------
    
    ip    = npix - ipix1 + 1;
    hip   = ip/2.;
/* bug corrige floor instead of 1.* */
    fihip = floor(hip);
    iring = (int)floor( sqrt( hip - sqrt(fihip) ) ) + 1;//     ! counted from South pole
    iphi  = (int)(4.*iring + 1 - (ip - 2.*iring*(iring-1)));
    
    z = -1. + iring*iring / fact2 ;
    phi   = (1.*iphi - 0.5) * PI/(2.*iring);
  }

  sz = sqrt( 1.0 - z*z );
  vec[0] = sz * cos(phi);
  vec[1] = sz * sin(phi);
  vec[2] = z;
  
}

double Packed49[49][3] = {
    {1.0820379E-008,-3.2300459E-006,    0.7790824},
    { 0.0851419,      -0.0604194,       0.3497294},
    {-0.0914288,       0.4137149,       0.6537967},
    {-0.4233704,      -0.0167043,       0.6537949},
    { 0.4161716,       0.0795128,       0.6537953},
    { 0.2606476,      -0.3340458,       0.6537932},
    {-0.1783143,      -0.3843534,       0.6537934},
    { 0.3214161,       0.4881783,       0.5151146},
    {-0.4918671,       0.3793286,       0.4702615},
    {-0.0929470,       0.3254455,       0.2208710},
    {-0.5385916,      -0.3977282,       0.3983726},
    { 0.6544342,      -0.1767253,       0.3839966},
    {-0.0464433,      -0.6884808,       0.3616719},
    { 0.6543453 ,      0.2637903 ,      0.3304789},
    { 0.3839234,      -0.5971428,       0.3209247},
    {-0.7108776,       0.0012718,       0.3187802},
    {-0.1876019,      -0.3352532,       0.1368439},
    { 0.0704146,       0.7351211,       0.2481284},
    {-0.3638534,       0.6695231,       0.1622308},
    { 0.3033485,       0.2022102,       0.0692749},
    { 0.4742969,       0.6067250,       0.1178835},
    {-0.6795831,       0.3744220,       0.0703152},
    {-0.3794767,       0.0482692,       0.0303653},
    { 0.6538247,      -0.4232979,       0.0173632},
    { 0.2332925,      -0.2885923,       0.0015460},
    { 0.7790813, -5.7994478E-013,      -0.0012951},
    {-0.5429419,      -0.5585797,      -0.0131201},
    {-0.1452212,      -0.7628716,      -0.0625065},
    {-0.7541588,      -0.1768841,      -0.0832219},
    { 0.2928920,      -0.7156702,      -0.0948672},
    {-0.0815266,       0.7567586,      -0.1662504},
    { 0.6534047,       0.3539813,      -0.2339385},
    {-0.4662671,       0.5642231,      -0.2668646},
    { 0.3250845,       0.6429856,      -0.2964101},
    {-0.6822678,       0.1837015,      -0.3282282},
    { 0.4930282,      -0.4578927,      -0.3927173},
    {-0.3428409,      -0.5690840,      -0.4069064},
    { 0.6530941,      -0.0471244,      -0.4221572},
    {-0.1036092,       0.2867325,      -0.2191358},
    { 0.0858176,      -0.5795982,      -0.5134887},
    {-0.5513357,      -0.1947856,      -0.5148367},
    { 0.0032634,       0.5253046,      -0.5753380},
    {-0.1660916,      -0.1851457,      -0.2781839},
    { 0.3945018,       0.3205518,      -0.5904102},
    {-0.3693146,       0.2921726,      -0.6206539},
    { 0.2268184,       0.0121013,      -0.3221586},
    { 0.2750635,      -0.2203113,      -0.6948182},
    {-0.1633231,      -0.2729136,      -0.7112054},
    { 0.0133638,       0.1279994,      -0.7683794}
    };


typedef struct model_ctx {
	/* Material coefficients from the Tillotson EOS. */
	TILLMATERIAL **tillMat;
	int nMat;
	int *iMatOrder;	// Here we store in which order the materials are
	double *fM;		// Masses of each layer

	/*
	** Some unit conversion factors.
	*/
	double dKpcUnit;
	double dMsolUnit;

	/*
	** The lookup table for the equilibrium model.
	*/
	int nTableMax;
	int nTable;
	double uc; /* u at r = 0 */
	double *M;
	double *rho;
	double *u;
	double *r;
	int *mat;

	double dr;
	double R;
	} MODEL;


MODEL *modelInit(double M,double ucore) {
	int i;
    /* Initialize the model */
	MODEL *model;
    
    model = malloc(sizeof(MODEL));
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
    model->M = malloc(model->nTableMax*sizeof(double));
    assert(model->M != NULL);
    model->rho = malloc(model->nTableMax*sizeof(double));
    assert(model->rho != NULL);
    model->u = malloc(model->nTableMax*sizeof(double));
    assert(model->u != NULL);
    model->r = malloc(model->nTableMax*sizeof(double));
    assert(model->r != NULL);
    model->mat = malloc(model->nTableMax*sizeof(int));
    assert(model->mat != NULL);
    model->dr =  0.0;
    model->nTable = 0;
    
	return(model);
    }

double drhodr(MODEL *model,int iMat,double r,double rho,double M,double u);
double dudr(MODEL *model,int iMat,double r,double rho,double M,double u);

/*
 * dudrho depends on the internal energy profile that we choose!
 */
double dudrho(MODEL *model,int iMat,double rho,double u)
{
	return(eosPressure(model->tillMat[iMat],rho,u)/(rho*rho));
}

/*
 * Calculate dudrho to solve for the equilibrium model.
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
 * This derivative is independent of the model and only involves geometry.
 */
double dMdr(double r,double rho) {
	assert(r >= 0.0);
	return(4.0*M_PI*r*r*rho);
}

/*
** Solve for rho2 and u2 using the b.c. P1=P2 and T1=T2.
*/
void modelSolveBC(MODEL *model, double *prho, double *pu, int iMat1, int iMat2) {
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
** This function integrates the ODEs for a two component model with b.c.
** rho_initial=rho, u_initial=u until rho(r=R)=rho0. It returns the total
** mass of the model. The parameter h sets the stepsize for the RK2 algorithm
** and for bSetModel=1 the results are saved.
*/
double modelSolveTwoComponent(MODEL *model,int bSetModel,double rho,double u,double h,double *pR)
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

	if (bSetModel) {
		model->rho[i] = rho;
		model->M[i] = M;
		model->u[i] = u;
		model->r[i] = r;
		model->mat[i] = iMat;
		++i;
	}

	/*
	** First integrate the core until M == Mcore.
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

    printf("\n");
    printf("rho= %g u= %g dConstGamma= %g (dConstGamma-1.0)= %g P= %g T= %g\n", rho_atm, u_atm, model->tillMat[0]->dConstGamma, model->tillMat[0]->dConstGamma-1.0, (model->tillMat[0]->dConstGamma-1.0)*rho_atm*u_atm, u_atm/model->tillMat[0]->cv);
    printf("\n");

    printf("mantle:   rho= %g u= %g P= %g T= %g\n", rho1, u1, Ps, Ts);
    printf("envelope: rho= %g u= %g P= %g T= %g\n", rho_atm, u_atm, (model->tillMat[0]->dConstGamma-1.0)*rho_atm*u_atm, u_atm/model->tillMat[0]->cv);
//    exit(1);
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

    modelSolveBCAtmosphere(model, model->tillMat[iMat], rho, u, &rho, &u);

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
        printf("%i\n", i);
		fprintf(fp,"%g %g %g %g %i %g %g\n",model->r[i],model->rho[i],model->M[i],model->u[i],model->mat[i],
			eosPressure(model->tillMat[model->mat[i]], model->rho[i], model->u[i]),
			eosTempRhoU(model->tillMat[model->mat[i]], model->rho[i], model->u[i]));
//		fprintf(fp,"%15.7E %15.7E %15.7E %15.7E %i\n",model->r[i],model->rho[i],model->M[i],model->u[i],model->mat[i]);
	}
	fclose(fp);
}

/*
** This function calls modelSolveComponent for every material for a given
** density rhoc and internal energy uc in the code and returns the
** total mass of the resulting planet. If bSetModel=1 the results are saved in
** the lookup table and both nTable and dr are set.5
*/
double modelSolveAll(MODEL *model,int bSetModel,double rhoc,double uc,double h,double *pR)
{
    FILE *fp;
	// Set inital values rhoc, uc and M
	double rho=rhoc;
	double u = uc;
	double M = 0.0;
	// Radius of the model
	double R = 0.0;
	int Index = 0;

	int i;

	fprintf(stderr,"***********************************************************\n");
	fprintf(stderr,"modelSolveAll:\n");
	fprintf(stderr,"modelSolveAll: bSetModel=%i rhoc=%g uc=%g h=%g\n",bSetModel,rhoc,uc,h);
	fprintf(stderr,"***********************************************************\n");
	fprintf(stderr,"\n");

	M = modelSolveTwoComponent(model,bSetModel,rho,u,h,&R);
/*
	for (i=0; i<model->nMat; i++)
	{
		// Solve for component i
		if (i == model->nMat-1)
		{
			fprintf(stderr,"modelSolveAll: Component:%i Index:%i h:%g rho:%g u:%g M1:%g M2:%g R:%g\n",model->iMatOrder[i], Index, h, rho, u, M,M+model->fM[i], R);
			// Enforce rho(r=R)=rho0 for the last layer
			modelSolveComponent(model,model->iMatOrder[i], bSetModel, 1, &Index, h, &rho, &u, &M,M+model->fM[i], &R);
		} else {
			fprintf(stderr,"modelSolveAll: Component:%i Index:%i h:%g rho:%g u:%g M1:%g M2:%g R:%g\n",model->iMatOrder[i], Index, h, rho, u, M,M+model->fM[i], R);
			modelSolveComponent(model,model->iMatOrder[i], bSetModel, 0, &Index, h, &rho, &u, &M,M+model->fM[i], &R);
			// Determine rho2 and u2 for the next material
			modelSolveBC(model, &rho, &u, model->iMatOrder[i], model->iMatOrder[i+1]);
			// tillSolveBC(model->tillMat[model->iMatOrder[i]],model->tillMat[model->iMatOrder[i]],rho,u,&rho,&u);

		}
	}
*/
	fprintf(stderr,"***********************************************************\n");
	fprintf(stderr,"modelSolveAll: Done. Index:%i h:%g rho:%g u:%g M:%g R:%g\n", Index, h, rho, u, M, R);
	fprintf(stderr,"***********************************************************\n");

	if (bSetModel) {
		// Set nTable and dr
		model->nTable = Index;
		model->dr = h;
		// Write the lookup table to the file ballic.model
		modelWriteToFile(model);
	}
	
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
    R = cbrt(3.0*M/(4.0*M_PI*model->tillMat[model->iMatOrder[1]]->rho0)); // Use lower density material for max radius
    dr = R/nStepsMax;
    a = 2.0*model->tillMat[model->iMatOrder[0]]->rho0; /* starts with 100% larger central density */

	// (CR) Debug
	fprintf(stderr,"R: %g a: %g\n",R,a);

	Ma = modelSolveTwoComponent(model,bSetModel=0,a,model->uc,dr,&R);
//    Ma = modelSolveAll(model,bSetModel=0,a,model->uc,dr,&R);
    fprintf(stderr,"first Ma:%g R:%g\n",Ma,R);
    b = a;
    Mb = 0.5*M;
    while (Ma > M) {
		b = a;
		Mb = Ma;
		a = 0.5*(model->tillMat[model->iMatOrder[0]]->rho0 + a);
		Ma = modelSolveTwoComponent(model,bSetModel=0,a,model->uc,dr,&R);
//		Ma = modelSolveAll(model,bSetModel=0,a,model->uc,dr,&R);
	}
    while (Mb < M) {
		b = 2.0*b;
	   	Mb = modelSolveTwoComponent(model,bSetModel=0,b,model->uc,dr,&R);	
//	   	Mb = modelSolveAll(model,bSetModel=0,b,model->uc,dr,&R);	
		fprintf(stderr,"first Mb:%g R:%g\n",Mb,R);
	}

	// (CR) Debug
	fprintf(stderr,"Root bracketed.\n");

    /*
    ** Root bracketed by (a,b).
    */
    while (Mb-Ma > 1e-10*Mc) {
		c = 0.5*(a + b);
        Mc = modelSolveTwoComponent(model,bSetModel=0,c,model->uc,dr,&R);	
//		Mc = modelSolveAll(model,bSetModel=0,c,model->uc,dr,&R);	
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
    fprintf(stderr,"rho_core: %g cv: %g uc: %g (in system units)\n",c,model->tillMat[MATERIAL]->cv,model->uc);
    Mc = modelSolveTwoComponent(model,bSetModel=1,c,model->uc,dr,&R);
	// This is only needed for modelSolveTwoComponent
	modelWriteToFile(model);
//    Mc = modelSolveAll(model,bSetModel=1,c,model->uc,dr,&R);
    model->R = R;
    return c;
    }

/*
** This code depends on the arrays being ordered in r but does not change
** if M, rho or u are not monotonic. So no changes needed for multi
** component models.
*/
double MLookup(MODEL *model,double r) {
    double x,xi,dr;
    int i;

    i = model->nTable-1;
    if (r >= model->r[i]) return(model->M[i]*(1.0 + log(r-model->r[i]+1)));
    x = r/model->dr;
    xi = floor(x);
	//(CR) Debugging
//	fprintf(stderr,"r=%g dr=%g x=%g xi=%g\n",r,model->dr,x,xi);
    assert(xi >= 0.0);
    x -= xi;
    i = (int)xi;
    if (i < 0) {
		fprintf(stderr,"ERROR r:%.14g x:%.14g xi:%.14g i:%d\n",r,x,xi,i);
	}
    assert(i >= 0);
    if (i < model->nTable-2) {
		return(model->M[i]*(1.0-x) + model->M[i+1]*x);
	}
    if (i == model->nTable-2) {
		dr = model->r[i+1] - model->r[i];
		x = r/dr;
		xi = floor(x);
		x -= xi;
		return(model->M[i]*(1.0-x) + model->M[i+1]*x);
	}
    else {
		i = model->nTable - 1;
		return(model->M[i]*(1.0 + log(r-model->r[i]+1)));
	}
    }

double rhoLookup(MODEL *model,double r) {
    double x,xi,dr;
    int i;

    i = model->nTable-1;
    if (r >= model->r[i]) return(model->rho[i]*exp(-(r-model->r[i])));
    x = r/model->dr;
    xi = floor(x);
    assert(xi >= 0.0);
    x -= xi;
    i = (int)xi;
    if (i < 0) {
		fprintf(stderr,"ERROR r:%.14g x:%.14g xi:%.14g i:%d\n",r,x,xi,i);
	}
    assert(i >= 0);
    if (i < model->nTable-2) {
		return(model->rho[i]*(1.0-x) + model->rho[i+1]*x);
	}
    if (i == model->nTable-2) {
		dr = model->r[i+1] - model->r[i];
		x = r/dr;
		xi = floor(x);
		x -= xi;
		return(model->rho[i]*(1.0-x) + model->rho[i+1]*x);
	}
    else {
		i = model->nTable - 1;
		return(model->rho[i]*exp(-(r-model->r[i])));
	}
    }

double uLookup(MODEL *model,double r) {
    double x,xi,dr;
    int i;

    i = model->nTable-1;
    if (r >= model->r[i]) return(model->u[i]*exp(-(r-model->r[i])));
    x = r/model->dr;
    xi = floor(x);
    assert(xi >= 0.0);
    x -= xi;
    i = (int)xi;
    if (i < 0) {
		fprintf(stderr,"ERROR r:%.14g x:%.14g xi:%.14g i:%d\n",r,x,xi,i);
	}
    assert(i >= 0);
    if (i < model->nTable-2) {
		return(model->u[i]*(1.0-x) + model->u[i+1]*x);
	}
    if (i == model->nTable-2) {
		dr = model->r[i+1] - model->r[i];
		x = r/dr;
		xi = floor(x);
		x -= xi;
		return(model->u[i]*(1.0-x) + model->u[i+1]*x);
	}
    else {
		i = model->nTable - 1;
		return(model->u[i]*exp(-(r-model->r[i])));
	}
    }

double Fzero(MODEL *model,int bIcosa,double r,double ri,double m,int ns) {
	long npix = (bIcosa)?(40*ns*(ns-1)+12):(12*ns*ns);
	return(MLookup(model,r)-MLookup(model,ri)-npix*m);
	}

double rShell(MODEL *model,double m,double ri) {
    double a = ri;
    double b = 1.0;
    double c,Mc,Ma;

    Ma = MLookup(model,a);
    Mc = 0.0;
    while (m > (MLookup(model,b)-Ma)) b *= 2.0;
    while (fabs(m-(Mc-Ma)) > 1e-10*m) {
		c = 0.5*(a + b);
		Mc = MLookup(model,c);
//		fprintf(stderr,"c:%.7g Mc:%.7g\n",c,Mc);
		if (m > (Mc-Ma)) a = c;
		else b = c;
	}
	return c;
	}

double rShell2(MODEL *model,int bIcosa,double m,double ri,int ns) {
    double a = ri;
    double b = 1.0;
    double c;
    double z;

    z = 1.0;
    while (Fzero(model,bIcosa,b,ri,m,ns) < 0) b *= 2.0;
    while ((b-a)/c > 1e-10) {
		c = 0.5*(a + b);
		z = Fzero(model,bIcosa,c,ri,m,ns);
		if (z < 0) a = c;
		else b = c;
//		printf("c:%.14g M(c):%.14g\n",c,MLookup(model,c));
	}
	return c;
	}

/*
** New version of ballic that solves the problem of distributing the particles separately for the core and the mantle.
*/
void main(int argc, char **argv) {
    const int bCentral = 1;
    int bIcosa = 0;
    const int bRandomRotate = 1;
    TCTX out;
    long ns,npix,ipix,na,nb;
    struct gas_particle gp;
    double r[3];
    double rhoCenter;
    double ri,ro,rs,roa,rob,ros,rta,rtb,rts;
    double m,mTot,l1,l2,nsf;
    double x,y,ang1,ang2,ang3;
    int j,iShell,nDesired,nReached,nLast;
    double theta,phi;
    float rr[3];
    ICOSA *ctx;
    int nShell,nMaxShell;
    double *rsShell;
    long *nsShell;
    long *isShell;
    double *xyz;
    MODEL *model;
    /*
     * Parameters to solve the model.
     */
    double rhocore, ucore, R, us;
    int bSetModel;
    int iter;
    int nSmooth;
    double d,rho,u,eta,w0,dPdrho;
    FILE *fpi,*fpo;
    int iRet,i;
	double mCore;
    /*
     * The number of components in the model. Note that this is not the same as
     * model->nMat, which is the total number of materials in the library.
     */
    int nLayer;
	/* A first guess for the number of particles of a material. */
	int nDesiredMat;
	int nReachedLastMat = 0;
	int iShellLastMat = 0;
	int iMat;
	double rMat; /* Final Radius of one component. */

    const int nStepsMax = 10000;

#if 0
    if (argc != 5) {
		fprintf(stderr,"Usage: ballic <nDesired> <TotalMass> <rhocore> <ucore> >myball.std\n");
		exit(1);
	}
	
	/* Get command line parameters. */
	nDesired = atoi(argv[1]);
	mTot = atof(argv[2]);
	rhocore = atof(argv[3]);
    ucore = atof(argv[4]);
#endif

    /*
     * Hard coded a model with three layers: core, mantle and envelope.
     */
    nLayer = 3;

	fprintf(stderr,"Model initialized.\n"); // CR
    /*
     * Solve the model for given M, rhoc and uc.
     */
    rhoCenter = rhocore;

#if 0
    /*
     * Model with M=62.366 and three components:
     * Core:     Iron
     * Mantle:   Granite
     * Envelope: Ideal gas
     *
     * 589   176   4.1780000E+01   8.5200000E+00   6.2366582E+01   6.8733496E-02
     */
    rhocore = 4.1780000E+01;
    ucore = 8.5200000E+00;
    mTot = 62.366;
	nDesired = 100000;
#endif

    /*
     * Uranus: M=13.5 ME (844.393) with three components:
     * Core:     Granite
     * Mantle:   Ice
     * Envelope: Ideal gas
     *
     * 951   641   4.42650000e+01   6.96150000e+01   8.44393180e+02   2.90598500e-02
     */
    rhocore = 4.42650000e+01;
    ucore = 6.96150000e+01;
    mTot = 844.393;
	nDesired = 100000;

    model = modelInit(mTot,ucore);

	R = cbrt(3.0*mTot/(4.0*M_PI*model->tillMat[model->iMatOrder[nLayer-2]]->rho0));
    mTot = modelSolveThreeComponent(model, bSetModel=1, rhocore, ucore, R/nStepsMax, &R, &us);
    modelWriteToFile(model);
	fprintf(stderr,"Model solved.\n");	// CR

#if 0
	/*
	** Desired number of particles in the core (Nc=fc*N).
	*/
	mCore = model->fM[0];
	nDesiredCore = mCore/mTot*nDesired;
#endif
	//m = mTot/nDesired;   /* a first guess at the particle mass */

    /*
    ** Initialize icosahedral parameters.
    */
    ctx = icosaInit();

#if (0)
    /*
    ** Test the icosahedral grid.
    */
    ns = 4;
    npix = 40*ns*(ns-1) + 12;
	for (ipix=0;ipix<nipix;++ipix) {
		icosaPix2Vec(ctx,ipix,ns,r);
		printf("%d %g %g %g\n",i,r[0],r[1],r[2]);
	}
#endif
    /*
    ** Initialize the array of shell radii and resolutions.
    */
    nMaxShell = 1000;
    nShell = 0;
    rsShell = malloc(nMaxShell*sizeof(double));
    assert(rsShell != NULL);
    nsShell = malloc(nMaxShell*sizeof(long));
    assert(nsShell != NULL);
    isShell = malloc(nMaxShell*sizeof(long));
    assert(nsShell != NULL);
	
	TipsyInitialize(&out,0,NULL);

	fprintf(stderr,"\n");
	fprintf(stderr,"*******************************************************\n");
	fprintf(stderr,"Distributing shells.\n");
	fprintf(stderr,"nMat=%i nDesired=%i mTot=%g\n",model->nMat,nDesired,mTot);
	fprintf(stderr,"*******************************************************\n");
	fprintf(stderr,"\n");
	for (iMat=0;iMat<model->nMat;iMat++)
	{
		fprintf(stderr,"Material %i:\n",iMat);
		fprintf(stderr,"tillMat=%i fM=%g nDesiredMat=%g m=%g (estimate)\n",
		model->iMatOrder[iMat],
		model->fM[iMat],
		model->fM[iMat]/mTot*nDesired,
		model->fM[iMat]/(model->fM[iMat]/mTot*nDesired));
	}

	/*
	** Distribute the particles in shells for each material separately.
	*/
//	for (iMat=0;iMat<model->nMat;iMat++)

	for (iMat=0;iMat<nLayer;iMat++)
	{
		/*
         * Careful, iMat is the number in which the materials are ordered in iMatOrder.
         */
		nDesiredMat = model->fM[iMat]/mTot*nDesired+nReachedLastMat;
		m = model->fM[iMat]/nDesiredMat;   /* a first guess at the particle mass */
		
		fprintf(stderr,"\n");
		fprintf(stderr,"Material %i: nDesiredMat=%i M=%g m=%g\n",iMat, nDesiredMat, model->fM[iMat],m);

		/*
		** Phase 1: We first determine the number of particles per shell such that
		** the ratio of radial to tangential lengths of their volumes is as close
		** to 1 as possible. This fixes the total number of particles in the model.
		*/
		fprintf(stderr,"PHASE 1: ODE approach (iMat=%i)\n",iMat);
		for (iter=0;iter<2;++iter) {
			if (iMat == 0)
			{
				/* Initialize nReached, ro and iShell for material 0. */
				nReached = 0;
				/* Set starting radius */
				if (bCentral) {
					ro = rShell(model,m,0.0);
				}
				else {
					ro = 0.0;
				}
				iShell = 0;
				rMat = 0.0;
			} else {
				/* For the other components we start from rMat of the previous material. */
				ro = rMat;
				nReached = nReachedLastMat;
				iShell = iShellLastMat;
			}

//			(CR) Debug
//			assert(iter == 0);

			for (;;) {
				ri = ro;

				na = 1;
				roa = rShell2(model,bIcosa,m,ri,na);
				npix = (bIcosa)?(40*na*(na-1)+12):(12*na*na);
			    l1 = roa-ri;
			    l2 = sqrt(M_PI/npix)*(roa+ri);
			    rta = l1/l2;
			    nb = 16;

			    do {
					nb *= 2;
					rob = rShell2(model,bIcosa,m,ri,nb);
					npix = (bIcosa)?(40*nb*(nb-1)+12):(12*nb*nb);
					l1 = rob-ri;
					l2 = sqrt(M_PI/npix)*(rob+ri);
					rtb = l1/l2;
				} while (rtb < 1.0);

				while (nb - na > 1) {
					ns = (na+nb)/2;
					ros = rShell2(model,bIcosa,m,ri,ns);
					npix = (bIcosa)?(40*ns*(ns-1)+12):(12*ns*ns);

//					(CR) Debugging
//					fprintf(stderr,"npix=%i ros=%g ri=%g ns=%i\n",npix,ros,ri,ns);
					l1 = ros-ri;
					l2 = sqrt(M_PI/npix)*(ros+ri);
					rts = l1/l2;
/*					fprintf(stderr,"ns:%d rts:%g\n",ns,rts); */
					if (rts < 1.0) {
						na = ns;
						roa = ros;
						rta = rts;
					}
					else {
						nb = ns;
						rob = ros;
						rtb = rts;
					}
				}
/*
				if (iShell >= 5) {
					fprintf(stderr,"na:%d rta:%g\n",na,rta);
					fprintf(stderr,"nb:%d rtb:%g\n",nb,rtb);
				}
*/
				/*
				** if the two possible ratios differ by less that 1% then we favour
				** the higher resolution spherical grid (nb).
				*/
				if (1/rta+0.01 < rtb) {
					ro = roa;
					ns = na;
					rts = rta;
				}
				else {
					ro = rob;
					ns = nb;
					rts = rtb;
				}

				npix = (bIcosa)?(40*ns*(ns-1)+12):(12*ns*ns);
				if (iShell == nMaxShell) {
					nMaxShell *= 2;
					rsShell = realloc(rsShell,nMaxShell*sizeof(double));
					assert(rsShell != NULL);
					nsShell = realloc(nsShell,nMaxShell*sizeof(long));
					assert(nsShell != NULL);
				}
				nsShell[iShell] = ns;
/*				fprintf(stderr,"nReached:%d npix:%d\n",nReached,npix);*/
				if ((nReached + npix) < nDesiredMat) {
					nReached += npix;
					fprintf(stderr,"iShell:%d ns:%d radial/tangential:%g\n",iShell,ns,rts);
					++iShell;
				}
				else {
					nShell = iShell;
					break;
				}
			}  /* end of iShell loop */
			ns = nsShell[iShell-1];
			npix = (bIcosa)?(40*ns*(ns-1)+12):(12*ns*ns);	
			if (nDesiredMat - nReached > npix/2) {
				nReached += npix;
				nsShell[nShell] = ns;
				fprintf(stderr,"iShell:%d ns:%d radial/tangential:??? (added)\n",nShell,ns);
				nShell++;
			}

			fprintf(stderr,"nReached:%d old mass:%.7g new mass:%.7g\n",
				nReached,m,model->fM[iMat]/nReached);
//			fprintf(stderr,"nReached:%d old mass:%.7g new mass:%.7g\n",
//				nReached,m,mTot/nReached);
			/* Update particles mass for the reached number of particles. */
			m = model->fM[iMat]/(nReached-nReachedLastMat);
			fprintf(stderr,"PHASE 1 update m: nReached=%i nReachedLastMat=%i dReached=%i iShell=%i iShellLastMat=%i nDesired=%i M=%g m=%g\n",nReached,nReachedLastMat,nReached-nReachedLastMat,iShell,iShellLastMat,nDesiredMat,model->fM[iMat],m);
//			m = mTot/nReached;
			nDesiredMat = nReached+1;
		}
//assert(0);
		fprintf(stderr,"PHASE 1 done: nReached=%i nReachedLastMat=%i iShell=%i iShellLastMat=%i nDesired=%i M=%g m=%g\n",nReached,nReachedLastMat,iShell,iShellLastMat,nDesiredMat,model->fM[iMat],m);
		/* Save how many particles we distributed for the last material. */
//		nReachedLastMat = nReached;
//		iShellLastMat = iShell;

		/*
		** Phase 2: With the numbers of particles in each of the shells now fixed
		** we continue by recalculating the radii of the shells based on the updated
		** particle mass.
		*/
		fprintf(stderr,"PHASE 2 (iMat=%i)\n",iMat);
		// (CR) For iMat > 0 we have to set ro properly...
		if (iMat == 0)
		{
			if (bCentral) {
				ro = rShell(model,m,0.0);
			}
			else {
				ro = 0.0;
			}
		} else {
			/* For the other components we start from rMat of the previous material. */
			ro = rMat;
		}
		
		for (iShell=iShellLastMat;iShell<nShell;++iShell) {
			ri = ro;
			ns = nsShell[iShell];
			ro = rShell2(model,bIcosa,m,ri,ns);
			npix = (bIcosa)?(40*ns*(ns-1)+12):(12*ns*ns);
			l1 = ro-ri;
			l2 = sqrt(M_PI/npix)*(ro+ri);
			rts = l1/l2;

			if (bIcosa) {
				d = 1.0;
				if (iShell == 0) d = 2.5;
				if (iShell == 1) d = 1.0;
				if (iShell == 2) d = 3.0;
			}
			else {
				d = 1.9;
				if (iShell == 0) d = 2.0;
				if (iShell == 1) d = 1.5;
			}

			rs = rsShell[iShell] = pow(0.5*(pow(ri,d) + pow(ro,d)),1/d);

			rho = rhoLookup(model,rs);

			u = uLookup(model,rs); /* We also have to look up u from a table */

			// Set rMat to ro
			rMat = ro;
			/* This is old code used to debug the grid. */
#if 0
			eta = rho/model->tillMat[MATERIAL]->rho0;
			/* This was the old code using a constant internal energy uFixed.
			w0 = model->uFixed/(model->par.u0*eta*eta) + 1.0;
			dPdrho = (model->par.a + (model->par.b/w0)*(3 - 2/w0))*model->uFixed + 
					(model->par.A + 2*model->par.B*(eta - 1))/model->par.rho0;

			fprintf(stderr,"iShell:%d r:%g M:%g rho:%g ns:%d radial/tangential:%g dr:%g <? Jeans:%g Gamma:%g\n",iShell,rs,MLookup(model,rs),rho,ns,rts,ro-ri,sqrt(dPdrho/rho),Gamma(model,rho,model->uFixed));
        	*/	
			w0 = u/(model->tillMat[MATERIAL]->u0*eta*eta) + 1.0;
			dPdrho = (model->tillMat[MATERIAL]->a + (model->tillMat[MATERIAL]->b/w0)*(3 - 2/w0))*u + 
					(model->tillMat[MATERIAL]->A + 2*model->tillMat[MATERIAL]->B*(eta - 1))/model->tillMat[MATERIAL]->rho0;

//  		fprintf(stderr,"iShell:%d r:%g M:%g rho:%g u:%g ns:%d radial/tangential:%g dr:%g <? Jeans:%g Gamma:%g\n",iShell,rs,MLookup(model,rs),rho,u,ns,rts,ro-ri,sqrt(dPdrho/rho),Gamma(model,rho,u));
#endif
		}
		fprintf(stderr,"PHASE 2 done: nReached=%i nReachedLastMat=%i iShell=%i iShellLastMat=%i nDesired=%i M=%g m=%g rMat=%g ro=%g\n",nReached,nReachedLastMat,iShell,iShellLastMat,nDesiredMat,model->fM[iMat],m,rMat,ro);
		/*
		** Now generate the coordinates of all the particles in each shell as they are on the unit
		** sphere. This simplifies the later adjusting of the radii of the shells (unit sphere 
		** coordinates are independent of this.
		*/

		/*
		** Phase 3: With the masses and numbers of the particles fixed and the 
		** as well as their angular positions, we now attempt to use an SPH density
		** estimate to calculate the mean radial pressure force on the shell 
		** and converge on the radius of the shell such that it is balanced by
		** gravity on the shell (analytic).
		*/
		if (bCentral) {
			/*
			** First adjust the density of the central particle.
			*/
		}

		/*
		** Now output the particles.
		*/
		for (j=0;j<3;++j) gp.vel[j] = 0.0;
		gp.phi = 0.0;
#if 0
		/*
		** If we have a central particle, the particle mass has to be adjusted.
		*/
		if (bCentral) {
			fprintf(stderr,"\n");
			fprintf(stderr,"mTot= %g m= %g",mTot,m);
			m = mTot/(nReached+1);
			fprintf(stderr,"  m= %g nReached= %i nReachedLastMat= %i\n",m,nReached, nReachedLastMat);
			fprintf(stderr,"Adjusted particle mass accounting for the central particle.\n");
			fprintf(stderr,"\n");
		}
#endif
		gp.mass = m;

		/* Save the material */
		gp.metals = model->iMatOrder[iMat];
		nLast = nReached;

		/* Start writing from nReachedLastMat otherwise we write the inner components more than once. */
		nReached = nReachedLastMat;
//		nReached = 0;

//*****************************************************************************************************************************************
//		(CR) 19.1.16: Writing a central particle that was not originally included in the calculation could cause the mass to deviate from M (this was fixed on 17.4.2017 in the sense that the total mass of the SPH model is consistent now).
//*****************************************************************************************************************************************
		if (bCentral && iMat==0) {
			/* Write the central particle (only for iMat=0). */
			for (j=0;j<3;++j) gp.pos[j] = 0.0;
			gp.temp = uLookup(model, 0);
			// Dont forget to set the material for the central particle
			gp.metals = model->iMatOrder[iMat];
			TipsyAddGas(out,&gp);
		}

		/* Start writing from iShellLastMat otherwise we write the inner components more than once. */
//		for (iShell=0;iShell<nShell;++iShell) {
		for (iShell=iShellLastMat;iShell<nShell;++iShell) {
			rs = rsShell[iShell];
			ns = nsShell[iShell];
			npix = (bIcosa)?(40*ns*(ns-1)+12):(12*ns*ns);
			nReached += npix;
			ang1 = 2.0*M_PI*rand()/(RAND_MAX+1.0);
			ang2 = 2.0*M_PI*rand()/(RAND_MAX+1.0);
			ang3 = 2.0*M_PI*rand()/(RAND_MAX+1.0);
			for (ipix = 0;ipix < npix;++ipix) {
				if (bIcosa) icosaPix2Vec(ctx,ipix,ns,r);
				else pix2vec_ring(ns,ipix,r);
				if (bRandomRotate) {
					y = r[1]*cos(ang1) - r[2]*sin(ang1);
					r[2] = r[1]*sin(ang1) + r[2]*cos(ang1);
			
					x = r[0]*cos(ang2) - r[2]*sin(ang2);
					r[2] = r[0]*sin(ang2) + r[2]*cos(ang2);
					r[0] = x;

					r[1] = y*cos(ang3) - r[2]*sin(ang3);
					r[2] = y*sin(ang3) + r[2]*cos(ang3);
				}

				for (j=0;j<3;++j) gp.pos[j] = rs*r[j];
	    
//				rho = rhoLookup(model,rs);
				gp.temp = uLookup(model,rs);
				TipsyAddGas(out,&gp);
			}
		}

		/* Save how many particles we distributed for the last material. */
		nReachedLastMat = nReached;
		iShellLastMat = iShell;
	} /* iMat */
 	fprintf(stderr,"\n");
	fprintf(stderr,"Writing %d particles. Model R:%g Last Shell r:%g\n",nReached,model->R,rsShell[nShell-1]);

	/* Write all particles to ballic.std */
	TipsyWriteAll(out,0.0,"ballic.std");
	TipsyFinish(out);
	}
