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
	double *z;
	double dz;
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
double drhodz(MODEL *model, double r, double rho, double M, double u);
double dudz(MODEL *model, double r, double rho, double M, double u);
double drhodz(MODEL *model, double r, double rho, double M, double u);
double dMdz(MODEL *model, double r, double rho);
double midPtRKIn(MODEL *model, int bSetModel, double rho, double u, double h, double *pR);
double midPtRKOut(MODEL *model, int bSetModel, double rho, double u, double h, double *pR);
//double midPtRK(MODEL *model,int bSetModel,double rho,double u,double h,double *pR);
double modelSolve(MODEL *model, double M, double us);

#endif // MODELSOLVE_HINCLUDED
