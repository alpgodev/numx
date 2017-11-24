#include <math.h>
#include <float.h>
#include <stdlib.h> 
#include <stdio.h>
#include <string.h>
#include "var_module.h"
/*=======================================================================
 *
 *  monteCarloBlockXVaR.c									      
 *
 *	Monte Carlo VaR Simulation for Energy Contracts by "Submarket Approach"
 *  with a specified number of factors (different) for each sub-market block
 *
 ------------------------------------------------------------------------*/
int monteCarloBlockXVaR(
				/* inputs */
				int *p, int *nblock, int *block, int *k,
				int *N, double *alpha, double *epsSDLS,
				double *w, double *corr, double *vol, double *rho,
				/* outputs */
				double *VaR, double *X, double *Y, int *info)

/*------------------------------------------------------------------------
 * Inputs
 *    p       :	integer, number of assets  (p >= k)
 *    nblock  : integer, number of submarket (m >= 1)
 *    block   : integer nblock-array, submarket indices [100, ...,400]
 *              sum(block[i], i=1,..,nblock) = p 
 *    k       : integer nblock-array, number of factors (k >= 1) for each block
 *    N       : integer, number of Monte-Carlo simulation (scenarios) (N >= 1) 
 *    alpha   : double, probability level (0 < alpha < 1) 
 *    epsSDLS : double, filtering level (>0)
 *    w       : double p-array, portfolio weights
 *    corr    :	double p*p-array, correlation matrix
 *    vol     : double p-array, volatility shocks
 *    rho     :	double p-array, expected returns
 *
 * Outputs
 *    VaR     : double, Gaussian Value-at-Risk
 *    X       : double N-array, simulated portfolio returns
 *    Y       : double N*p-array, simulated assets returns
 *    info    : integer, diagnostic argument
 *
 * Call
 *    monteCarloVaR, monteCarloFactor, monteCarloPortfolioVaR, 
 *    YM, IMXRN, EXPNVA, PMX, PVX, CHOL, PM, PMMT
 --------------------------------------------------------------------------*/
{
	int i=0,j=0,l=0;
	double a=0.0, sum=0.0;
	double epsBFGS=1.E-12;		/* BFGS precision		 */
	double mean=0, std=1.0;		/* N(0,1) Gaussian model */ 
	double *F, *FF, *G, *H, *cholG, *E, *Rnd; 
	double *Fi, *Fj, *Di, *D, *Ui, *Zij;
	double *invDi, *invDj;
	double *corrk, *corrij;
	
	/* workspaces */
	int *iwork;
	double *dwork;

	int nbFactors=0;
	int sumblock=0,sumki=0,ki=0,kj=0;
	int nbk=0,m=0;
	int nbki=0,nbkj=0;
	int rin=1,cin=1,rout=1,cout=1;
	int nbl=0,nbc=0,nbli=0,nbci=0,nblj=0,nbcj=0;

	//FILE *fpt;
	//fpt = fopen("output.txt","w");	

	/* test parameters */
	/* if number of submarket = 1 -> "Direct Method" */
	if ((*nblock) == 1) {
		nbFactors = k[0];
		monteCarloVaR(p, &nbFactors, N, alpha, epsSDLS, w, corr, vol, rho, VaR, X, Y, info);
		return 0;
	}

	/* test sum(block[i]) = p */
	for (i=0; i<(*nblock); i++){
		sumblock += block[i];
	}
	if ((sumblock < (*p)) || (sumblock > (*p))) {
		*info = -9000;
		return 5;
	}

	/* test if k[i] > 0 */
	for (i=0; i<(*nblock); i++){
		if (k[i] < 1) {
			*info = -9011;
			return 5;
		}
	}

	/* test if size block[i] < nb factors[i] */
	for (i=0; i<(*nblock); i++){
		if (block[i] < k[i]) {
			*info = -9010;
			return 5;
		}
	}

	/* sum(k) */
	m = 0;
	for (i=0; i<(*nblock); i++){
		m += k[i];
	}

	/* factors eigenvalues D(m) */
	D  = (double *) malloc(m*sizeof(double));

	/* factor matrix F(p, sum(k)) 
	
	        | F1  0 ... 0 |                                      
	        | 0  F2 ... 0 |
	    F = | .  .  ... . |  Fk, k-ith block factor
	        | .  .  ... . |
	        | 0  0 ... Fk |
    */
	F  = (double *) malloc((*p)*m*sizeof(double));
	FF = (double *) malloc((*p)*m*sizeof(double));
	a = 0.0;
	imx_( p, &m, F, &a );  /* zero-matrix initialization */
	imx_( p, &m, FF, &a ); /* zero-matrix initialization */

	/* cross-correlation G(m,m), m = sum(k) 
	
	        |   Id       Z(1,2) ... Z(1,nblock) |                                      
	        | Z(2,1)      Id    ... Z(2,nblock) |
	    G = |   .         .     ...      .      |  
	        |   .         .     ...      .      |
	        | Z(nblock,1) .     ...      Id     |

		Z(i,j) the cross correlation between block i,j

		Z(i,j) = inv(sqrt(Di))*Fi'*Corr(i,j)*Fj*inv(sqrt(Dj))
    */
	G = (double *) malloc(m*m*sizeof(double));
	a = 1.0;imdx_( &m, G, &a ); /* matrix identity initialization */
	
	nbk = 1;
	/* construct factor matrix F */
	for (i=0; i<(*nblock); i++){

		/* dim of block i */
		nbki = block[i];

		/* nb factors of block i */
		ki = k[i];

		/* extract correlation block j */
		corrk = (double *) malloc(nbki*nbki*sizeof(double));
		ymp_(p, p, corr, &nbki, &nbki, &nbk, &nbk, corrk, info);
		if (*info < 0){return 5;}

		/* computation of factor i */
		Fi = (double *) malloc(nbki*ki*sizeof(double));
		Di = (double *) malloc(nbki*sizeof(double));
		Ui = (double *) malloc(nbki*nbki*sizeof(double));
		monteCarloFactorBlock(&nbki, &ki, epsSDLS, corrk, Fi, Di, Ui, info);
		if (*info < 0){return 5;}

		/* scale eigenvectors to have Euclidian norm egal to root of eigenvalues */
		/* copy Di -> D (eigenvalues): inv(sqrt(lambda[i])) */
		Fj = (double *) malloc(nbki*ki*sizeof(double));
		for (j=0; j<ki; j++) {
			a = pow(Di[j],0.5);
			D[j+cout-1] = 1.0/a;
			for (l=0; l<nbki; l++){
				Fj[l+j*nbki] = a*Fi[l+j*nbki];
			}
		}
		free(corrk);free(Di);free(Ui);
		if (*info < 0){return 5;}
		
		/* copy Fj -> F */
		ympmp_(&nbki, &ki, &nbki, &ki, Fj, &rin, &cin, p, &m, F, &rout, &cout, info);
		free(Fj);
		if (*info < 0){return 5;}

		/* copy Fi -> FF */
		ympmp_(&nbki, &ki, &nbki, &ki, Fi, &rin, &cin, p, &m, FF, &rout, &cout, info);
		free(Fi);
		if (*info < 0){return 5;}

		/* iterator */
		nbk += nbki;
		rout = nbk;
		cout += ki;
	}
		
	/* construct cross-correlation matrix */
	nbli=1;nbci=1;
	nblj=1;nbcj=1; 
	nbl=1;nbc=1;
	nbk=1;
	sumki = 1;
	for (i=0; i<(*nblock)-1; i++){

		nbki  = block[i]; /* dim of block i           */
		ki    = k[i];     /* nb factors of block i    */
		nbcj  = sumki;    /* index of the j-th column */

		/* extract factor i */
		Fi = (double *) malloc(nbki*ki*sizeof(double));
		ymp_(p, &m, FF, &nbki, &ki, &nbli, &nbci, Fi, info);
		//fprintf(fpt,"\n i=%i,nbki=%i,ki=%i,nbli=%i,nbci=%i",i,nbki,ki,nbli,nbci);
		//fprintf(fpt,"\n");
		if (*info < 0){
		//	fprintf(fpt,"\nError in YMP (0).\n");
		//	fclose(fpt);
			return 5;}
		
		nbli+=nbki; /* i-th row of F */
		nblj =nbli; /* j-th row of F */
		nbc  =nbli; 

		for (j=i+1; j<(*nblock); j++){
		
			nbkj = block[j]; /* dim of block j           */
			kj   = k[j];     /* nb factors of block j    */
			nbcj += k[j-1];      /* index of the j-th column */

			/* extract factor j */
			Fj = (double *) malloc(nbkj*kj*sizeof(double));
			ymp_(p, &m, FF, &nbkj, &kj, &nblj, &nbcj, Fj, info);
//			fprintf(fpt,"\n i=%i,j=%i,nbkj=%i,kj=%i,nblj=%i,nbcj=%i",i,j,nbkj,kj,nblj,nbcj);
//			fprintf(fpt,"\n");
			if (*info < 0){
//				fprintf(fpt,"\nError in YMP (1).\n");
//				fprintf(fpt,"\n i=%i,j=%i,nbkj=%i,kj=%i,nblj=%i,nbcj=%i",i,j,nbkj,kj,nblj,nbcj);
//				fclose(fpt);
				return 5;}
			
			/* extract correlation block[i] vs. block[j] */
			corrij = (double *) malloc(nbki*nbkj*sizeof(double));
			ymp_(p, p, corr, &nbki, &nbkj, &nbl, &nbc, corrij, info);
			if (*info < 0){
//				fprintf(fpt,"\nError in YMP (2).\n");
//				fclose(fpt);
				return 5;}

			/* dwork(block[i],kj) = Corr(block[i],block[j])*Fj(block[j],kj) */
			dwork = (double *) malloc(nbki*kj*sizeof(double));
			pm_( &nbki, &nbkj, &kj, corrij, Fj, dwork);
			free(corrij);

			/* Zij(ki,kj) = Fi'(block[i],ki)*dwork(block[i],kj) */
			Zij = (double *) malloc(ki*kj*sizeof(double));
			pmtm_( &nbki, &ki, &kj, Fi, dwork, Zij );
			free(dwork);

			/* cross-correlation between block i and block j        	*/
			/* Zij(ki,kj) = inv(sqrt(Di))*Zij(ki,kj)*inv(sqrt(Dj))		*/
			invDi = (double *) malloc(ki*ki*sizeof(double));
			invDj = (double *) malloc(kj*kj*sizeof(double));
			a = 0.0;
			imx_( &ki, &ki, invDi, &a ); /* zero-matrix initialization */
			imx_( &kj, &kj, invDj, &a ); /* zero-matrix initialization */
			for (l=0; l<ki; l++){
				invDi[l+l*ki] = D[l+nbci-1];
			}
			for (l=0; l<kj; l++){
				invDj[l+l*kj] = D[l+nbcj-1];
			}
			dwork = (double *) malloc(ki*kj*sizeof(double));
			pm_( &ki, &ki, &kj, invDi, Zij, dwork);	/* dwork(ki,kj) = inv(sqrt(Di))* Zij(ki,kj)  */
			pm_( &ki, &kj, &kj, dwork, invDj, Zij);	/* Zij(ki,kj)   = dwork(ki,kj)*inv(sqrt(Dj)) */
			free(dwork);free(invDi);free(invDj);

			/* construct block i,j of G(m,m) = [0 Zij; Zij' 0] */
			ympmp_(&ki, &kj, &ki, &kj, Zij, &rin, &cin, &m, &m, G, &nbci, &nbcj, info);
			if (*info < 0){
//				fprintf(fpt,"\nError in YMPYMP (1).\n");
//				fclose(fpt);
				return 5;}
			
			dwork = (double *) malloc(kj*ki*sizeof(double));
			rm_(&ki, &kj, Zij, dwork); /* transpose Zij -> Zji */
			free(Zij);
			ympmp_(&kj, &ki, &kj, &ki, dwork, &rin, &cin, &m, &m, G, &nbcj, &nbci, info);
			free(dwork);
			if (*info < 0){
//				fprintf(fpt,"\nError in YMPYMP (2).\n");
//				fclose(fpt);
				return 5;}

			free(Fj);
			nbc +=nbkj; 
			nblj+=nbkj;
		}
		sumki += ki;   /* sum of ki */ 
		nbci += ki;
		nbl  += nbki;
		free(Fi);
	}
	free(FF);free(D);

	/* computation of standard shocks */
	/* E(p,m) = F(p,m).*vol(p) */
	E = (double *) malloc((*p)*m*sizeof(double));
	for (j=0; j<m; j++) {
		for (i=0; i<(*p); i++){
			a = vol[i];
			E[i+j*(*p)] = a*F[i+j*(*p)];
		}
	}
	free(F);

	/* cross-correlation G(m,m) calibraion (SDLS) */
	H      = (double *) malloc(m*m*sizeof(double));
	dwork  = (double *) malloc((m*(7*m+40))*sizeof(double)); 
	iwork  = (int *) malloc((14*m+2)*sizeof(int));
	ym_( &m, &m, G, H ); /* copy a vectorized matrix in a vectorized matrix */
	sdlscor_( &m, H, &epsBFGS, epsSDLS, iwork, dwork, G, info );
	free(H);free(dwork);free(iwork);
	if (*info < 0){return 5;}

	/* compute choleski of cross-correlation matrix: chol(G) */
	cholG = (double *) malloc(m*m*sizeof(double));
	dwork = (double *) malloc(m*m*sizeof(double));
	chol_(&m, G, dwork, cholG, info);
	free(G);free(dwork);
	if (*info < 0){return 5;}
	H = (double *) malloc(m*m*sizeof(double));
	rm_( &m, &m, cholG, H); /* transpose chol(G) */
	free(cholG);

	/* Gaussian model - Rnd(N,m) */
	Rnd = (double *) malloc((*N)*m*sizeof(double));
	imxrn_(N, &m, &mean, &std, Rnd); /* gererate matrix X(i,j)=Normal(mean,std) */

	/* Q(N,m) = Rnd(N,m)*chol[G(m,m)] */
	dwork = (double *) malloc((*N)*m*sizeof(double));
	pm_( N, &m, &m, Rnd, H, dwork); /* Rnd(N,m)*cholG'(m,m) */
	free(Rnd);free(H);

	/* Y(N,p) = Q(N,m)*E'(p,m) */
	pmmt_( N, &m, p, dwork, E, Y ); /* Y(N,p)=Rnd(N,m)*cholG'(m,m)*E'(p,m) */
	free(E);free(dwork);

	/* add the mean: Y(N,p) = Y(N,p) + rho(p) */
	for (j=0; j<(*p); j++){
		a = rho[j];
		for (i=0; i<(*N); i++){
			Y[i+j*(*N)] = Y[i+j*(*N)] + a;
		}
	}

	/* step 4 - compute portfolio ex-post Value-at-Risk */
	monteCarloPortfolioVaR(p, N, alpha, w, Y, &a, X, info);
	*VaR=a;
	if (*info < 0){ return 5;}
//	fclose(fpt);
	return 0;
}
