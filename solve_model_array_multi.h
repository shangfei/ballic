/*
 ** The header file for modelsolve.
 */
#ifndef SOLVE_MODEL_ARRAY_MULTI_HINCLUDED
#define SOLVE_MODEL_ARRAY_MULTI_HINCLUDED

#include "tillotson/tillotson.h"
#include "nr/nrutil.h"

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

/*
 * Functions to initialize and finalize a model.
 */
MODEL *modelInit(double M,double ucore);
void modelFinalize(MODEL *model);

/*
 * Functions needed to calculate the equilibrium model.
 */
double drhodr(MODEL *model,int iMat,double r,double rho,double M,double u);
double dudr(MODEL *model,int iMat,double r,double rho,double M,double u);
double drhodr(MODEL *model, int iMat, double r,double rho,double M,double u);
double dudr(MODEL *model,int iMat,double r,double rho,double M,double u);
double dMdr(double r,double rho);

void modelSolveBC(MODEL *model, double *prho, double *pu, int iMat1, int iMat2);
//void modelSolveComponent(MODEL *model,int iMat,int bSetModel,int bLastLayer,int *pIndex,double h,double *prho1,double *pu1,double *pM1,double M2,double *pR);
double modelSolveTwoComponent(MODEL *model,int bSetModel,double rho,double u,double h,double *pR,double *us);
void modelWriteToFile(MODEL *model);

#endif // SOLVE_MODEL_ARRAY_MULTI_HINCLUDED
