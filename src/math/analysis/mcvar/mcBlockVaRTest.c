#include <math.h>
#include <float.h>
#include <stdlib.h> 
#include <stdio.h>
#include <string.h>
#include "var_module.h"
/*=======================================================================
 *
 *  monteCarloBlockVaRTest.c									      
 *
 *	Monte Carlo VaR Simulation for Energy Contracts by "Submarket Approach"
 *
 *  Step 1 : Divide the market in m submarket that is (R=[R,R2,…,Rm]).
 *
 *  Step 2 : For each submarket
 *
 *  Step 3 : Generation of Shocks on the factors
 *			– Compute the correlation matrix between the local market factors 
 *            and computation of its Choleski decomposition
 *			– Randomize the factors N(0,Gf) to generate N scenarios
 *			– Generate the shocks for the returns
 * 
 *  Step 4 : Computation of the Gaussian VaR(alpha%)
 *
 ------------------------------------------------------------------------*/
int monteCarloBlockVaRTest(
				/* inputs */
				int *p, int *k, int *nblock, int *block, 
				int *N, double *alpha, double *epsSDLS,
				double *w, double *corr, double *vol, double *rho,
				/* outputs */
				double *VaR, double *X, double *Y, double *F, double *G, int *info)

/*------------------------------------------------------------------------
 * Inputs
 *	  p       :	integer, number of assets  (p >= k)
 *    k       : integer, number of factors (k >= 1)
 *    nblock  : integer, number of submarket (m >= 1)
 *    block   : integer nblock-array, submarket indices [100, ...,400]
 *              sum(block[i], i=1,..,nblock) = p 
 *    N       : integer, number of Monte-Carlo simulation (scenarios) (N >= 1) 
 *    alpha   : double, probability level (0 < alpha < 1) 
 *    epsSDLS : double, filtering level (>0)
 *    w       : double p-array, portfolio weights
 *	  corr    :	double p*p-array, correlation matrix
 *    vol     : double p-array, volatility shocks
 *	  rho     :	double p-array, expected returns
 *
 * Outputs
 *	  VaR     : double, Gaussian Value-at-Risk
 *    X       : double N-array, simulated portfolio returns
 *    Y       : double N*p-array, simulated assets returns
 *    F       : double p*k*nblock-array, matrix of factors 
 *    G       : double m*m-array, cross-correlation matrix
 *	  info    : integer, diagnostic argument
 *
 *	Call
 *    monteCarloVaR, monteCarloFactor, monteCarloPortfolioVaR, 
 *    YM, IMXRN, EXPNVA, PMX, PVX, CHOL, PM, PMMT
 --------------------------------------------------------------------------*/
{
	int i=0,j=0,l=0;
	double a=0.0, sum=0.0;
	double eps=0.0;
	double epsBFGS=1.E-12;		/* BFGS precision		 */
	double mean=0, std=1.0;		/* N(0,1) Gaussian model */ 
	double *FF, *H, *cholG, *E, *Rnd; 
	double *Fi, *Fj, *Di, *D, *Ui, *Zij;
	double *invDi, *invDj;
	double *corrk, *corrij;
	
	/* workspaces */
	int *iwork;
	double *dwork;

	int nbk=0,m=0;
	int nbki=0,nbkj=0;
	int rin=1,cin=1,rout=1,cout=1;
	int nbl=0,nbc=0,nbli=0,nbci=0,nblj=0,nbcj=0;

	/* test parameters */
	/* if number of submarket = 1 -> "Direct Method" */
	if ((*nblock) == 1) {
		monteCarloVaR(p, k, N, alpha, epsSDLS, w, corr, vol, rho, VaR, X, Y, info);
		return 0;
	}

	/* test sum(block[i]) = p */
	for (i=0; i<(*nblock); i++){
		j += block[i];
	}
	if ((j < (*p)) || (j > (*p))) {
		*info = -9000;
		return 5;
	}

	/* nblock*k */
	m = (*k)*(*nblock);

	/* factors eigenvalues D(m) */
	D  = (double *) malloc(m*sizeof(double));

	/* compute factor matrix F(p, nblock*k) 
	
	        | F1  0 ... 0 |                                      
	        | 0  F2 ... 0 |
	    F = | .  .  ... . |  Fk, k-ith block factor
	        | .  .  ... . |
	        | 0  0 ... Fk |
    */
	FF = (double *) malloc((*p)*m*sizeof(double));
	a = 0.0;
	imx_( p, &m, F, &a );	/* zero-matrix initialization */
	imx_( p, &m, FF, &a );	/* zero-matrix initialization */

	/* compute cross-correlation G(m,m), m = nblock*k 
	
	        |   Id       Z(1,2) ... Z(1,nblock) |                                      
	        | Z(2,1)      Id    ... Z(2,nblock) |
	    G = |   .         .     ...      .      |  
	        |   .         .     ...      .      |
	        | Z(nblock,1) .     ...      Id     |

		Z(i,j) the cross correlation between block i,j

		Z(i,j) = Fi'*Cov(i,j)*Fj
    */
	a = 1.0;imdx_( &m, G, &a ); /* matrix identity initialization */
	
	nbk = 1;
	/* construct factor matrix F */
	for (i=0; i<(*nblock); i++){

		/* dim of block i */
		nbki = block[i];

		/* extract correlation block j */
		corrk = (double *) malloc(nbki*nbki*sizeof(double));
		ymp_(p, p, corr, &nbki, &nbki, &nbk, &nbk, corrk, info);
		if (*info < 0){
			return 5;
		}

		/* computation of factor i */
		Fi = (double *) malloc(nbki*(*k)*sizeof(double));
		Fj = (double *) malloc(nbki*(*k)*sizeof(double));
		Di = (double *) malloc(nbki*sizeof(double));
		Ui = (double *) malloc(nbki*nbki*sizeof(double));
		monteCarloFactorBlock(&nbki, k, epsSDLS, corrk, Fi, Di, Ui, info);

		/* scale eigenvectors to have Euclidian norm egal to root of eigenvalues */
		/* copy Di -> D (eigenvalues): inv(sqrt(lambda[i])) */
		for (j=0; j<(*k); j++) {
			a = pow(Di[j],0.5);
			D[j+cout-1] = 1.0/a;
			for (l=0; l<nbki; l++){
				Fj[l+j*nbki] = a*Fi[l+j*nbki];
			}
		} 
		free(corrk);free(Di);free(Ui);
		if (*info < 0){
			return 5;
		}
		
		/* copy Fj -> F */
		ympmp_(&nbki, k, &nbki, k, Fj, &rin, &cin, p, &m, F, &rout, &cout, info);
		free(Fj);

		/* copy Fi -> FF */
		ympmp_(&nbki, k, &nbki, k, Fi, &rin, &cin, p, &m, FF, &rout, &cout, info);
		free(Fi);
		if (*info < 0){
			return 5;
		}

		/* iterator */
		nbk += nbki;
		rout = nbk;
		cout += (*k);
	}
		
	/* construct cross-correlation matrix */
	Zij	  = (double *) malloc((*k)*(*k)*sizeof(double));
	invDi = (double *) malloc((*k)*(*k)*sizeof(double));
	invDj = (double *) malloc((*k)*(*k)*sizeof(double));
	
	nbli=1;nbci=1;
	nblj=1;nbcj=1; 
	nbl=1;nbc=1;
	nbk=1;
	for (i=0; i<(*nblock)-1; i++){

		/* dim of block i */
		nbki = block[i];

		nbcj = i*(*k)+1;

		/* extract factor i */
		Fi = (double *) malloc(nbki*(*k)*sizeof(double));
		ymp_(p, &m, FF, &nbki, k, &nbli, &nbci, Fi, info);
		
		nbli+=nbki; /* i-th row of F */
		nblj =nbli;
		nbc  =nbli; 

		for (j=i+1; j<(*nblock); j++){
		
			/* dim of block j */
			nbkj = block[j];
		
			/* index of the j-th column */
			nbcj += (*k);

			/* extract factor j */
			Fj = (double *) malloc(nbkj*(*k)*sizeof(double));
			ymp_(p, &m, FF, &nbkj, k, &nblj, &nbcj, Fj, info);
			
			/* extract correlation block i,j */
			corrij = (double *) malloc(nbki*nbkj*sizeof(double));
			ymp_(p, p, corr, &nbki, &nbkj, &nbl, &nbc, corrij, info);
			if (*info < 0){
				return 5;
			}

			/* dwork(block[i],k) = Cov(block[i],block[j])*Fj(block[j],k) */
			dwork = (double *) malloc(nbki*(*k)*sizeof(double));
			pm_( &nbki, &nbkj, k, corrij, Fj, dwork);
			free(corrij);

			/* Zij(k,k) = Fi'(block[i],k)*dwork(block[i],k) */
			pmtm_( &nbki, k, k, Fi, dwork, Zij );
			free(dwork);

			/* cross-correlation between block i and block j	*/
			/* Zij(k,k) = inv(sqrt(Di))*Zij*inv(sqrt(Dj))		*/
			a = 0.0;
			imx_( k, k, invDi, &a ); /* zero-matrix initialization */
			imx_( k, k, invDj, &a ); /* zero-matrix initialization */
			for (l=0; l<(*k); l++){
				invDi[l+l*(*k)] = D[l+nbci-1];
				invDj[l+l*(*k)] = D[l+nbcj-1];
			}
			dwork = (double *) malloc((*k)*(*k)*sizeof(double));
			pm_( k, k, k, invDi, Zij, dwork);	/* dwork = inv(sqrt(Di))* ZZij  */
			pm_( k, k, k, dwork, invDj, Zij);   /* Zij   = dwork*inv(sqrt(Dj))	*/
			free(dwork);

			/* construct block j of G(m,m) = [0 Zij; Zij' 0] */
			ympmp_(k, k, k, k, Zij, &rin, &cin, &m, &m, G, &nbci, &nbcj, info);
			dwork = (double *) malloc((*k)*(*k)*sizeof(double));
			rm_(k, k, Zij, dwork);
			ympmp_(k, k, k, k, dwork, &rin, &cin, &m, &m, G, &nbcj, &nbci, info);
			free(dwork);
			if (*info < 0){
				return 5;
			}
			free(Fj);
			nbc +=nbkj; 
			nblj+=nbkj;
		}
		nbci+=(*k);
		nbl += nbki;
		free(Fi);
	}
	free(FF);free(D);free(Zij);free(invDi);free(invDj);

	/* computation of standard shocks */
	/* E(p,m) = F(p,m).*vol(p) */
	E = (double *) malloc((*p)*m*sizeof(double));
	for (j=0; j<m; j++) {
		for (i=0; i<(*p); i++){
			a = vol[i];
			E[i+j*(*p)] = a*F[i+j*(*p)];
		}
	}

	/* cross-correlation G(m,m) calibraion (SDLS) */
	H      = (double *) malloc(m*m*sizeof(double));
	dwork  = (double *) malloc((m*(7*m+40))*sizeof(double)); 
	iwork  = (int *) malloc((14*m+2)*sizeof(int));
	ym_( &m, &m, G, H );
	sdlscor_( &m, H, &epsBFGS, epsSDLS, iwork, dwork, G, info );
	free(H);free(dwork);free(iwork);
	if (*info < 0){
		return 5;
	}

	/* compute choleski of cross-correlation matrix: chol(G) */
	cholG = (double *) malloc(m*m*sizeof(double));
	dwork = (double *) malloc(m*m*sizeof(double));
	chol_(&m, G, dwork, cholG, info);
	free(dwork);
	if (*info < 0){
		return 5;
	}
	H = (double *) malloc(m*m*sizeof(double));
	rm_( &m, &m, cholG, H);
	free(cholG);

	/* Gaussian model - Rnd(N,m) */
	Rnd = (double *) malloc((*N)*m*sizeof(double));
	imxrn_(N, &m, &mean, &std, Rnd);

	/* Q(N,m) = Rnd(N,m)*chol[G(m,m)] */
	dwork = (double *) malloc((*N)*m*sizeof(double));
	pm_( N, &m, &m, Rnd, H, dwork);
	free(Rnd);free(H);

	/* Y(N,p) = Q(N,m)*E'(p,m) */
	pmmt_( N, &m, p, dwork, E, Y );
	free(E);free(dwork);

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
