/*
 ** The header file for modelsolve.
 */
#ifndef MODELSOLVE_HINCLUDED
#define MODELSOLVE_HINCLUDED

#include "tillotson/tillotson.h"

typedef struct model_ctx
{
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
	double *M;
	double *rho;
	double *u;
	double *r;
	double dr;
    double Mtot;
    double us;
	double R;
} MODEL;

/*
 * Functions to initialize and finalize a model.
 */
MODEL *modelInit(int iMat);
void modelFinalize(MODEL *model);

/*
 * Functions needed to calculate the equilibrium model.
 */
double dudrho(MODEL *model, double rho, double u);
double dudm(MODEL *model, double M, double r, double rho, double u);
double drhodm(MODEL *model, double M, double r, double rho, double u);
double drdm(double r, double rho);
int midPtRKIn(MODEL *model, double Mtot, double R, double rho0, double us, double h, double M_mid, double *pR, double *pRho, double *pU);
int midPtRKOut(MODEL *model, int bSetModel, double rhoc, double uc, double h, double M_mid, double *pR, double *pRho, double *pU);

# if 0
double dudr(MODEL *model, double r, double rho, double M, double u);
double drhodr(MODEL *model, double r, double rho, double M, double u);
double dMdr(double r, double rho);
double midPtRKIn(MODEL *model, double rho, double M, double u, double h, double *pRho, double *pU, double *pR);
double midPtRKOut(MODEL *model, int bSetModel, double rho, double u, double h, double *pR);
//double midPtRK(MODEL *model,int bSetModel,double rho,double u,double h,double *pR);
double modelSolve(MODEL *model, double M, double us);
#endif
#endif // MODELSOLVE_HINCLUDED
