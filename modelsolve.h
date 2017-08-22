/*
 ** The header file for modelsolve.
 */
#ifndef MODELSOLVE_HINCLUDED
#define MODELSOLVE_HINCLUDED

#include "tillotson/tillotson.h"
#include "nr/nrutil.h"

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
void modelGetNewInitialGuess(double dR, double drhoc, double duc, double ddR, double ddrho, double ddu, double *pDelta);
void modelSolveToMatchingPoint(MODEL *model, double rhoc, double uc, double M, double R, double rhos, double us, double M_mid, double dm, double *dR, double *drho, double *du);
void modelSolveSingle(MODEL *model, double Mtot, double rhos, double us);
# if 0
double dudr(MODEL *model, double r, double rho, double M, double u);
double drhodr(MODEL *model, double r, double rho, double M, double u);
double dMdr(double r, double rho);
double midPtRKIn(MODEL *model, double rho, double M, double u, double h, double *pRho, double *pU, double *pR);
double midPtRKOut(MODEL *model, int bSetModel, double rho, double u, double h, double *pR);
//double midPtRK(MODEL *model,int bSetModel,double rho,double u,double h,double *pR);
double modelSolve(MODEL *model, double M, double us);
#endif

/*
 * Functions needed to solve a linear sytem of equations using the Gauss algorithm from Numerical Recipes.
 */
/*
 * Gauss-Jordan elimination to solve systems of linear equations.
 *
 * a[1..n][1..n]:   pointer to the matrix containing the coefficients on the 
 *                  left hand side of the equation.
 * n:               number of columns / rows of the matrix a (its a square matrix)
 * b[1..n][1..m]:   pointer to the matrix containing the coefficients on the
 *                  right hand side of the equation.
 * m:               number of right hand side vectors b (usually m=1).
 *
 * Note: The original version uses float which was changed to double.
 */
void gaussj(double **a, int n, double **b, int m);
//void gaussj(float **a, int n, float **b, int m);

#endif // MODELSOLVE_HINCLUDED
