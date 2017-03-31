/*
 * The header file for model.single.c.
 */
#ifndef MODELSINGLE_HINCLUDED
#define MODELSINGLE_HINCLUDED
#include "tipsy.h"
#include "tillotson/tillotson.h"
#include "grid.h"

#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) > (B) ? (B) : (A))

// Isentropic thermal profile
#define BALLIC_U_ISENTROPIC

typedef struct model_ctx {
	/* Material coefficients from the Tillotson EOS. */
	TILLMATERIAL *tillMat;
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
	double dr;
	double R;
	} MODEL;

/* Initialize the model. */
MODEL *modelInit(double ucore, int iMat);

double dudrho(MODEL *model,double rho,double u);
double drhodr(MODEL *model,double r,double rho,double M,double u);
double dudr(MODEL *model,double r,double rho,double M,double u);
double dMdr(double r,double rho);
double CalcGrav(double r, double M);
double midPtRK(MODEL *model,int bSetModel,double rho,double h,double *pR);
double modelSolve(MODEL *model,double M);
#endif // MODELSINGLE_HINCLUDED
