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
 *  mcVaR.c									      
 *
 *	Monte Carlo VaR Simulation for Energy Contracts by "Direct Approach"
 *
 *------------------------------------------------------------------------
 *	
 *	 Copyright (c) 2012 NumX
 * 	 515, route de Savoie, 
 *	 38114 Allemont, France
 *	 All rights reserved.
 *
 *	 This software is the confidential and proprietary information
 *	 of NumX. You shall not disclose such Confidential
 *	 Information and shall use it only in accordance with the terms
 *	 of the licence agreement you entered into with NumX.
 *
 *   Author: Vernaz Yann
 *
 ------------------------------------------------------------------------*/
int monteCarloVaR(
				/* inputs */
				int *p, int *k, int *N, double *alpha, double *epsSDLS,
				double *w, double *corr, double *vol, double *rho,
				/* outputs */
				double *VaR, double *X, double *Y, int *info)

/*------------------------------------------------------------------------
 * Inputs
 *	  p       :	integer, number of assets  (p >= k) pmax=5000
 *    k       : integer, number of factors (k >= 1) kmax=1000
 *    N       : integer, number of Monte-Carlo simulation (scenarios) Nmax=1000000
 *    alpha   : double, probability level (0 < alpha < 1) 
 *    epsSDLS : double, filtering level (>0)
 *    w       : double p-array, portfolio weights
 *	  corr    :	double p*p-array, correlation matrix
 *    vol     : double p-array, volatilities
 *	  rho     :	double p-array, expected returns
 *    model   : simulation model, =1 Gaussian N(0,1) 
 *                                =2 Gaussian N(0,e) 
 *                                =3 non-Gaussian 
 *
 * Outputs
 *	  VaR	  : double, Gaussian Value-at-Risk at probability level alpha
 *    X       : double N-array, simulated portfolio returns
 *    Y       : double N*p-array, simulated assets returns
 *	  info	  : integer, diagnostic argument
 *
 *	Call
 *    monteCarloFactor, monteCarloPortfolioVaR, 
 *    YM, IMXRN, EXPNVA, PMX, PVX
 --------------------------------------------------------------------------*/
{
	int i=0,j=0;
	double a=0.0,sum=0.0;
	double epsBFGS=1.E-12;		/* BFGS precision		 */
	double mean=0.0,std=1.0;	/* N(0,1) Gaussian model */
	double *D, *U, *F, *E, *Rnd;

	/* test input parameters */
	int pmax = 5000;
	int kmax = 1000;
	int Nmax = 1000000;
	
	if (*p > pmax) {*info = 0; return 31;}
	if (*k > kmax) {*info = 0; return 32;}
	if (*N > Nmax) {*info = 0; return 33;}
	if ((*alpha < 0.0)||(*alpha>1)) {*info = 0; return 34;} /* test if 0 < alpha < 1 */
	if (*epsSDLS < 0.0) {*info = 0; return 35;} 
	if (*p < *k) {*info = 0; return 32;} /* test if nb assets < nb factors */

		
	/* step 1 - compute k-factors */
	/* F : factor matrix (p*k-array)        */
    /* D : eigenvalues vector (p-array)     */
    /* U : eigenvectors matrix (p*p-array   */
	/* corr = U'*D*U                        */
	F = (double *) malloc((*p)*(*k)*sizeof(double));
	D = (double *) malloc((*p)*sizeof(double));
	U = (double *) malloc((*p)*(*p)*sizeof(double));
	monteCarloFactor(p, k, epsSDLS, corr, F, D, U, info);
	free(D);free(U);
	if (*info < 0){
		free(F);*VaR=a;
		return 5;
	}
	
	/* step 2 - computation of standard shocks */
	E = (double *) malloc((*p)*(*k)*sizeof(double));
	for (j=0; j<(*k); j++){
		for (i=0; i<(*p); i++){
			a = vol[i];
			E[i+j*(*p)] = a*F[i+j*(*p)];
		}
	}
	free(F);

	/* step 3 - randomize, Y(N,p) = Rnd(N,k)*E'(p,k) */
	/* Gaussian model */
	Rnd = (double *) malloc((*N)*(*k)*sizeof(double));
	imxrn_(N, k, &mean, &std, Rnd );

	/* Y(N,p) = Rnd(N,k)*E'(p,k) */
	pmmt_( N, k, p, Rnd, E, Y );
	free(E);free(Rnd);

	/* Y(N,p) = Y(N,p) + rho(p) */	
	for (j=0; j<(*p); j++){
		a = rho[j];
		for (i=0; i<(*N); i++){
			Y[i+j*(*N)] = Y[i+j*(*N)] + a;
		}
	}

	/* step 4 - compute portfolio ex-post Value-at-Risk */
	monteCarloPortfolioVaR(p, N, alpha, w, Y, &a, X, info);
	*VaR=a;
	if (*info < 0){ 
		return 5;
	}
	return 0;
}
