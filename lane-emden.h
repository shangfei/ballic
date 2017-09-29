/*
 * The header file for lane-emden.c.
 */
#ifndef LANEEMDEN_HINCLUDED
#define LANEEMDEN_HINCLUDED

#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) > (B) ? (B) : (A))

typedef struct model_ctx {
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
	double uc;      /* u at r = 0 */
    double rhoc;    /* rho at r = 0 */
    double gamma;   /* polytropic index */
    double n;
    double K;
    double alpha;
    /*
     * Dimensionless variables.
     */
    double *w1;
	double *z;
    /*
     * The physical quantities of the model.
     */
	double *M;
	double *rho;
	double *u;
	double *r;
	double dr;
	double R;
	} MODEL;

/* Initialize the model. */
MODEL *modelInit(double rhoc, double ucore, double gamma);

double Pressure(MODEL *model, double rho, double u);
double dw1dz(MODEL *model, double z, double w1, double w2);
double dw2dz(MODEL *model, double z, double w1, double w2);
double dmudz(MODEL *model, double z, double w1, double w2);
double CalcGrav(double r, double M);
double RhoTheta(MODEL *model, double w1);
double UTheta(MODEL *model, double w1);
double RZ(MODEL *model, double z);
double MMu(MODEL *model, double mu);
double LESolver(MODEL *model, int bSetModel, double h, double *pR);
#endif // LANEEMDEN_HINCLUDED
