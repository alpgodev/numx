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
 *  mcPortfolioVaR.c									      
 *
 *	Step 1 : construct portfolio
 * 
 *	Step 2 : computation of the Gaussian VaR
 *           Use EXPNVA (Gaussian VaR)
 *
 ------------------------------------------------------------------------*/
int monteCarloPortfolioVaR(
				/* inputs */
				int *p, int *N, double *alpha,
				double *w, double *Y,
				/* outputs */
				double *VaR, double *X, int *info)

/*------------------------------------------------------------------------
 * Inputs
 *	  p       :	integer, number of assets  (p >= k)
 *    N       : integer, number of Monte-Carlo simulation (scenarios) 
 *    alpha   : double, probability level (0 < alpha < 1) 
 *	  w       :	double p-array, portfolio weights
 *    Y	      : double N*p-array, simlated assets returns
 *
 * Outputs
 *	  VaR	  : double, Gaussian Value-at-Risk
 *    X       : double N-array, simulates portfolio returns
 *	  info	  : integer, diagnostic argument
 *
 *	Call
 *    pmv, expnva
 --------------------------------------------------------------------------*/
{
	double a=0.0;

	/* test parameters */

	/* step 1 - weighted portfolio */
	pmv_( N, p, Y, w, X ); /* X(i,j) = X(i,j)*w(j) */

	/* step 2 - compute ex-post Value-at-Risk */
	expnva_(N, X, alpha, &a, info);		/* Ex-post normal (Gaussian) Value-at-Risk */
	//expcvar_(N, X, alpha, &a, info);	/* Ex-post normal (i.e. Gaussian) conditional Value-at-Risk */
	//expmva_(N, X, alpha, &a, info );	/* Ex-post Cornish-Fisher (or modified) Value-at-Risk */
	
	*VaR=a;
	if (*info < 0){ return 5;}
	return 0;
}
