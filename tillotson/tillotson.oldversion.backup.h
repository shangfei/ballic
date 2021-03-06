/*
 ** The header file for the Tillotson EOS library.
 */
#ifndef TILLOTSON_HINCLUDED
#define TILLOTSON_HINCLUDED

#include "interpol/coeff.h"
#include "interpol/interpol.h"
#include "nr/nrcubicspline.h"

#define GRANITE 1
#define IRON 2

#define TILL_N_MATERIAL_MAX 2
/* Degree of the spline function we use for interpolation. */
#define TILL_SPLINE_DEGREE 3

/* We write the look up table as a 1D array where a(i,j)=a(i*N+j) */
#define INDEX(i, j) (((i)*material->nTableV) + (j))

typedef struct till_lookup_entry
{
	double u;
	double rho;
    double u1;		// du/drho
    /*
    ** The following 2 variables are the second derivatives of the above
    ** variables (u and u1) with respect to v (which is the value of the
    ** constant entropy curve (adiabat) at rho_0). Both of these are obtained
    ** by fitting splines to u and u1 runs in the v axis.
    */
    double udv2; 	// d2u/dv2
    double u1dv2;	// d2/dv2(du/drho)
	// Only for debugging
	double udrho2;
} TILL_LOOKUP_ENTRY;

typedef TILL_LOOKUP_ENTRY* TILL_LOOKUP;

typedef struct tillMaterial
{
	int iMaterial;	/* What material is it? */
//	int nTableMax;	/* Max. number of entries in the look up table */
	int nTableRho;	/* Number of entries in the look up table in rho */
	int nTableV;	/* and in v */
	int n;			/* number of steps from rho to zero */
	double rhomax;	/* Max value for the lookup table */
	double vmax;	/* Max value for the lookup table */
	double iExpV;
	/* Unit convertion factors */
	double dKpcUnit;
	double dMsolUnit;
	
	double dGasConst;
	double dErgPerGmUnit;
	double dGmPerCcUnit;
	double dSecUnit;

	/* A material is defined by 10 Tillotson parameters */
	double a;
	double b;
	double A;
	double B;
	double rho0;
	double u0;
	double us;
	double us2;
	double alpha;
	double beta;

	 /* the specific heat capacity (assumed constant) */
	double cv;

	/* The cold curve */	
	TILL_LOOKUP_ENTRY *cold;
//	double delta;
	double drho;
	double dv;

	/* A look up table for u(rho) along an isentrope */
//	TILL_LOOKUP *Lookup;	// this is an array of pointers
	TILL_LOOKUP_ENTRY *Lookup;	// this is an array of pointers
} TILLMATERIAL;

//TILLMATERIAL *tillInitMaterial(int iMaterial, double dKpcUnit, double dMsolUnit, int nTableMax, double rhomax, double vmax);
TILLMATERIAL *tillInitMaterial(int iMaterial, double dKpcUnit, double dMsolUnit, int nTableRho, int nTableV, double rhomax, double vmax, int iExpV);
void tillFinalizeMaterial(TILLMATERIAL *material);
double tilldPdrho(TILLMATERIAL *material, double rho, double u);
double tillSoundSpeed2old(TILLMATERIAL *material, double rho, double u);
double tillPressureSoundold(TILLMATERIAL *material, double rho, double u, double *pcSound);
double tillPressureSound(TILLMATERIAL *material, double rho, double u, double *pcSound);
double tillPressure(TILLMATERIAL *material, double rho, double u);
double tilldPdrho(TILLMATERIAL *material, double rho, double u);
double tilldPdu(TILLMATERIAL *material, double rho, double u);
double tilldTdrho(TILLMATERIAL *material, double rho, double u);
double tilldTdu(TILLMATERIAL *material, double rho, double u);
double tillTempRhoU(TILLMATERIAL *material, double rho, double u);
double tillSoundSpeed(TILLMATERIAL *material, double rho, double u);
double tilldudrho(TILLMATERIAL *material, double rho, double u);
void tillInitColdCurve(TILLMATERIAL *material);
void tillInitLookup(TILLMATERIAL *material);
TILL_LOOKUP_ENTRY *tillSolveIsentrope(TILLMATERIAL *material, double v);
/* Use bsstep.c from the Numerical Recipes */
TILL_LOOKUP_ENTRY *tillSolveIsentropeBS(TILLMATERIAL *material, double v);
/* Stuff for the cubic spline interpolator */
void tillInitSplines(TILLMATERIAL *material);
void tillInitSplineRho(TILLMATERIAL *material);
void tillInitSplinev(TILLMATERIAL *material);
void tillInitSplineU(TILLMATERIAL *material);
void tillInitSplineU1(TILLMATERIAL *material);
double tillCubicInt(TILLMATERIAL *material, double rho, double v);
void cubicint(double u[2],double dudrho[2], double dudv[2], double dudvdrho[2], double rho[2], double rhoint, double *intvalues);
// Just for debugging
double tillSplineIntrho(TILLMATERIAL *material, double rho, int iv);
double tillSplineIntv(TILLMATERIAL *material, double v, int irho);
double tillSplineIntU(TILLMATERIAL *material, double v, int irho);
double tillSplineIntU1(TILLMATERIAL *material, double v, int irho);

float tillFindUonIsentrope(TILLMATERIAL *material,float v,float rho);
/* Used for the root finder */
float denergy(TILLMATERIAL *material,float v,float rho,float u);
float tillFindEntropyCurve(TILLMATERIAL *material,float rho,float u,int iOrder);
double tillLookupU(TILLMATERIAL *material,double rho1,double u1,double rho2,int iOrder);
double tillColdULookup(TILLMATERIAL *material,double rho);
double tillCalcU(TILLMATERIAL *material,double rho1,double u1,double rho2);

/* Defines for the Numerical Recipes routines */

/*
** Bulirsch-Stoer method to integrate ODEs.
**
** y[]:			dependent variable
** dydx[]:		its first derivative at the starting value x
** nv:			number of variables y1,...,yn
** xx:			
** htry:		step size (the algorithm can use a smaller value if needed)
** eps:			required accuracy
** yscal[]:		vector to scale the error
** hdid:		return actual step size
** hnext:		return estimated next step size
** derivs:		function to calculate the right hand side derivatives
*/
void bsstep(float y[], float dydx[], int nv, float *xx, float htry, float eps,
	float yscal[], float *hdid, float *hnext,
	void (*derivs)(float, float [], float []));



#endif

