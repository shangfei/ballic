/*
 ** Copyright (c) 2014-2016 Joachim Stadel and Christian Reinhardt.
 **
 ** ballic provides a low noise particle representation of equilibrium
 ** models.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <assert.h>
#include "tipsy.h"
#include "tillotson/tillotson.h"
#include "grid.h"
#include "model.single.h"
#include "ballic.single.h"

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


void main(int argc, char **argv) {
	/*
	** Do we want a central particle: yes (1) or no (0).
	*/
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
	// Model
    MODEL *model;
    double ucore;
	int iMat;
    int iter;
    int nSmooth;
    double d,rho,u,eta,w0,dPdrho;
    FILE *fpi,*fpo;
    int iRet,i;
	/*
	** These variables are used to find the optimal softening.
	*/
	double l1max = 0.0;
	double l2max = 0.0;

    if (argc != 5) {
	fprintf(stderr,"Usage: ballic <nDesired> <TotalMass> <ucore> <iMat> >myball.std\n");
	exit(1);
	}
    nDesired = atoi(argv[1]);
    mTot = atof(argv[2]);
    ucore = atof(argv[3]);
    iMat = atoi(argv[4]);

	/* Assure that the parameters make sense. */

    assert(nDesired > 0);
    assert(mTot > 0.0);
    assert(ucore > 0.0);
    assert(iMat >= 0);

    model = modelInit(ucore, iMat);
    rhoCenter = modelSolve(model,mTot);

    m = mTot/nDesired;   /* a first guess at the particle mass */
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
    /*
    ** Phase 1: We first determine the number of particles per shell such that
    ** the ratio of radial to tangential lengths of their volumes is as close
    ** to 1 as possible. This fixes the total number of particles in the model.
    */
    fprintf(stderr,"PHASE 1: ODE approach\n");
    for (iter=0;iter<2;++iter) {
	nReached = 0;
 	if (bCentral) {
	    ro = rShell(model,m,0.0);
	    }
	else {
	    ro = 0.0;
	    }
	iShell = 0;
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
		l1 = ros-ri;
		l2 = sqrt(M_PI/npix)*(ros+ri);
		rts = l1/l2;
/*		fprintf(stderr,"ns:%d rts:%g\n",ns,rts); */
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
/*	fprintf(stderr,"nReached:%d npix:%d\n",nReached,npix);*/
	    if ((nReached + npix) < nDesired) {
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
	if (nDesired - nReached > npix/2) {
	    nReached += npix;
	    nsShell[nShell] = ns;
	    fprintf(stderr,"iShell:%d ns:%d radial/tangential:??? (added)\n",nShell,ns);
	    nShell++;
	    }
		fprintf(stderr,"nReached:%d old mass:%.7g new mass:%.7g\n",
	    nReached,m,mTot/nReached);
		m = mTot/nReached;
		nDesired = nReached+1;
	}
    /*
    ** Phase 2: With the numbers of particles in each of the shells now fixed
    ** we continue by recalculating the radii of the shells based on the updated
    ** particle mass.
    */
    fprintf(stderr,"PHASE 2\n");
    if (bCentral) {
	ro = rShell(model,m,0.0);
	}
    else {
	ro = 0.0;
	}
    for (iShell=0;iShell<nShell;++iShell) {

	ri = ro;
	ns = nsShell[iShell];
	ro = rShell2(model,bIcosa,m,ri,ns);
	npix = (bIcosa)?(40*ns*(ns-1)+12):(12*ns*ns);
	l1 = ro-ri;
	l2 = sqrt(M_PI/npix)*(ro+ri);
	rts = l1/l2;
	/*
	** Here we calculate the optimal softening for gravity using the
	** largest bin size. l1: radial size, l2: tangental size
	*/
	if (l1 > l1max) {
			l1max = l1;
	}

	if (l2 > l2max) {
			l2max = l2;
	}

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

	eta = rho/model->tillMat->rho0;
	    /* This was the old code using a constant internal energy uFixed.
	    w0 = model->uFixed/(model->par.u0*eta*eta) + 1.0;
		dPdrho = (model->par.a + (model->par.b/w0)*(3 - 2/w0))*model->uFixed + 
		(model->par.A + 2*model->par.B*(eta - 1))/model->par.rho0;

		fprintf(stderr,"iShell:%d r:%g M:%g rho:%g ns:%d radial/tangential:%g dr:%g <? Jeans:%g Gamma:%g\n",iShell,rs,MLookup(model,rs),rho,ns,rts,ro-ri,sqrt(dPdrho/rho),Gamma(model,rho,model->uFixed));
        */	
    	w0 = u/(model->tillMat->u0*eta*eta) + 1.0;
        dPdrho = (model->tillMat->a + (model->tillMat->b/w0)*(3 - 2/w0))*u + 
        (model->tillMat->A + 2*model->tillMat->B*(eta - 1))/model->tillMat->rho0;

//        fprintf(stderr,"iShell:%d r:%g M:%g rho:%g u:%g ns:%d radial/tangential:%g dr:%g <? Jeans:%g Gamma:%g\n",iShell,rs,MLookup(model,rs),rho,u,ns,rts,ro-ri,sqrt(dPdrho/rho),Gamma(model,rho,u));
        }
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
    TipsyInitialize(&out,0,NULL);

    for (j=0;j<3;++j) gp.vel[j] = 0.0;
    gp.phi = 0.0;
	/*
	** As a first guess we use max(l1max,l2max) for the softening.
	*/
	fprintf(stderr,"l1max: %g l2max: %g epsilon: %g\n",l1max,l2max,MAX(l1max,l2max));
    gp.hsmooth = MAX(l1max,l2max);  /* is actually eps for gasoline */

	/*
	** If we have a central particle, the particle mass has to be adjusted.
	*/
	if (bCentral) {
		m = mTot/(nReached+1);
	}
    gp.mass = m;

	fprintf(stderr,"hsmooth=%g\n",gp.hsmooth);

	gp.hsmooth=0.001;
    
	//gp.temp = model->uFixed;   /* Christian's version of gasoline uses thermal energy instead of temperature as input! */
    nLast = nReached;
    nReached = 0;
    if (bCentral) {
		for (j=0;j<3;++j) gp.pos[j] = 0.0;
		gp.temp = uLookup(model, 0);
		// Dont forget to set the material for the central particle
		gp.metals = iMat;
		TipsyAddGas(out,&gp);
	}
    for (iShell=0;iShell<nShell;++iShell) {
	rs = rsShell[iShell];
	ns = nsShell[iShell];
	npix = (bIcosa)?(40*ns*(ns-1)+12):(12*ns*ns);
	nReached += npix;
	ang1 = 2.0*M_PI*rand()/(RAND_MAX+1.0);
	ang2 = 2.0*M_PI*rand()/(RAND_MAX+1.0);
	ang3 = 2.0*M_PI*rand()/(RAND_MAX+1.0);

	/* Experiment with grav. softening. */
	if (iShell < nShell-1)
	{
		// Inner shell
		gp.hsmooth = 0.001;
	} else {
		// Most outer shell
		gp.hsmooth = 0.01;
	}

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
	    
//	    rho = rhoLookup(model,rs);
	    gp.temp = uLookup(model,rs);
		// Save Material
		gp.metals = iMat;
	    TipsyAddGas(out,&gp);		
	    }
	}
    fprintf(stderr,"Writing %d particles. Model R:%g Last Shell r:%g\n",nReached,model->R,rsShell[nShell-1]);
    TipsyWriteAll(out,0.0,"ballic.std");
    TipsyFinish(out);    
#if 0
    /* Smooth with gather doesn work with my version! */
    system("sleep 1; smooth -s 128 density <ballic.std; sleep 1");

    fpi = fopen("smooth.den","r");
    assert(fpi != NULL);
    fpo = fopen("ballic.den","w");
    assert(fpo != NULL);
    iRet = fscanf(fpi,"%d",&nReached);
    assert(iRet == 1);
    if (bCentral) {
	iRet = fscanf(fpi,"%lf",&rho);
	assert(iRet == 1);
	fprintf(fpo,"0.0 %g\n",rho);
	nReached -= 1;
	}
    for (iShell=0;iShell<nShell;++iShell) {
	rs = rsShell[iShell];
	ns = nsShell[iShell];
	npix = (bIcosa)?(40*ns*(ns-1)+12):(12*ns*ns);
	nReached -= npix;
	for (ipix = 0;ipix < npix;++ipix) {
	    iRet = fscanf(fpi,"%lf",&rho);
	    assert(iRet == 1);
	    fprintf(fpo,"%g %g\n",rs,rho);
	    }
	}
    fclose(fpi);
    fclose(fpo);
    assert(nReached == 0);
#endif
    }
