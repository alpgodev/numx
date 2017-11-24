/*****************************************************************************
 *
 *  LARGE SCALE Monte-Carlo Value-at-Risk
 *
 ****************************************************************************/
#include <math.h>
#include <float.h>
#include <stdlib.h> 
#include <stdio.h>
#include <string.h>
#include "var_module.h"
/*=======================================================================
 *
 *  monteCarloFactor.c									      
 *
 *  This function computes the k-factors
 *
 *	- matrix calibration to ensure SDP matrix
 *    Use sdlscor function
 * 
 *	- spectral decomposition and computation of the k-factors
 *	  Use DSYEVR function (Lapack)
 *
 ------------------------------------------------------------------------*/
int monteCarloFactor(
				/* inputs */
				int *p, int *k, double *epsSDLS, double *corr,
				/* outputs */
				double *F, double *D, double *U, int *info)

/*------------------------------------------------------------------------
 * Inputs
 *	  p       :	integer, number of assets  (p >= k)
 *    k       : integer, number of factors (k >= 1)
 *    epsSDLS : double, filtering level (>0)
 *	  corr    :	double p*p-array, correlation matrix
 *
 * Outputs
 *    F       : double p*k-array, factor matrix
 *    D       : double p-array, eigenvalues
 *    U       : double p*p-array, eigenvectors
 *	  info    : integer, diagnostic argument
 *
 *	Call
 *    SDLSCOR, DSYEVR, YM, PMX, PVX
 --------------------------------------------------------------------------*/
{
	int i=0, j=0;
	int m=0, il=0, iu=0, liwork=0, ldwork=0; 
	int *iwork, *isuppz;
	
	double epsBFGS = 1.E-12; /* BFGS precision */
	double a=0.0,sum=0.0;
	double vl=0.0, vu=0.0, abstol=0.0;
	double *corr1, *corr2, *dwork, *vp;

	/* test parameters */

	/* copy cov matrix */
	corr1 = (double *) malloc((*p)*(*p)*sizeof(double));
	// ym_( p, p, corr, corr1 );

	/* matrix calibration (SDLS)
       semi-definite least square optimization 

	   min || X - Corr ||
		s.t.
			X(i,i) = 1.0   (equality constraint)
			X >= alpha*Id  ( 0.0 < alpha < 1.0 )     
	*/
	dwork  = (double *) malloc(((*p)*(7*(*p)+40))*sizeof(double)); 
	iwork  = (int *) malloc((14*(*p)+2)*sizeof(int));
	sdlscor_( p, corr, &epsBFGS, epsSDLS, iwork, dwork, corr1, info );
	free(dwork);free(iwork);
	if (*info < 0){free(corr1); return 5;}

	/* spectral decomposition */
    /* DSYEVR (cf. LAPACK) 
	   DSYEVR computes selected eigenvalues and, optionally, eigenvectors
       of a real symmetric matrix.  Eigenvalues and eigenvectors can be
       selected by specifying either a range of values or a range of
       indices for the desired eigenvalues */
	vl = 0.0;
	vu = 0.0;
	il = 1;				/* first eigenvalue */
	iu = (*p);			/* last eigenvalue */
	abstol = 0.0;		/* the absolute error tolerance for the eigenvalues */
    liwork = 10*(*p);	/* iwork size */
    ldwork = 26*(*p);   /* dwork size */
	isuppz = (int *) malloc(2*(*p)*sizeof(double));
	dwork  = (double *) malloc(ldwork*sizeof(double)); 
	iwork  = (int *) malloc(liwork*sizeof(int));
	vp = (double *) malloc((*p)*sizeof(double));

    dsyevr_( "V", "A", "U", p, corr1, p,
             &vl, &vu, &il, &iu, &abstol, &m, vp, U, p,
			 isuppz, dwork, &ldwork, iwork, &liwork, info );
	
	/* sort eigenvalues vector in decrease order */
    avod_(p, vp, D);
	
	/* free memory */
	free(isuppz);free(iwork);free(dwork);
	free(corr1);free(vp);

	/* computation of k-factors
     * copy U -> F
     * ensures the uniqueness of the eigenvectors
     */
	/* E(:,j) = U(:,j)*sqrt(D(j)) for j=1,...,k */
	for (j=0; j<(*k); j++) {
		a = pow(D[j],0.5);
		sum = 0.0;
		for (i=0; i<(*p); i++){
			F[i+j*(*p)] = a*U[i - (*p) + (*p)*(*p) - j*(*p)];
			sum = sum + F[i+j*(*p)];
		}
		if (sum < 0) {
    	    for (i=0; i<(*p); i++){
    	        F[i+j*(*p)] = -F[i+j*(*p)];
            }
        }
	} 
	return 0;
}
