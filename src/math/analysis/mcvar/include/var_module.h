/*****************************************************************************
 *
 *  LARGE SCALE MONTE-CARLO Value-at-Risk
 *
 *	var_module.h
 *
 ****************************************************************************/
/* C header functions are declared here */   
#ifdef WINDOWS 
#include <math.h>
#endif /*WINDOWS*/

extern void transposeDoubleMatrix(double *, int, int);
extern void transposeIntMatrix(int *, int, int);

extern void transposeDoubleMatrix(double *, int, int);
extern void transposeIntMatrixOld(int *, int, int);
extern void transposeDoubleMatrixOld(double *, int, int);
extern int  lowestindex(int, int, int);
extern int  next(int, int, int);
extern int checkSysInfo();
extern int getReturnCode(int);

/* basic macros */
#define SQR(x)    ((x)*(x))
#define MAX(x, y)  ((x) > (y) ? (x) : (y))
#define MIN(x, y)  ((x) < (y) ? (x) : (y))


/* constants for algorithm */
#define MAX_ITER_PLANE 40
#define MAX_LAMBDA2 1e-8  /* plane search stopping crit. */
#define DIV_ALPHA 5
#define MIN_ALPHA 1e-8    /* max. of about 20 line search iterations */
#define MAX_ITER_LINEAR 1000
#define MAX_ITER_OUT 1000

#ifdef nounderscores
#define ddot_ ddot
#define dcopy_ dcopy
#define daxpy_ daxpy
#define dscal_ dscal
#define dgemv_ dgemv
#define dsyevr_ dsyevr
#define dsyr_ dsyr
#define dsyrk_ dsyrk
#define dposvx_ dposvx;
#define dgelss_ dgelss;
#define gwjms_ gwjms;
#define jms_ jms;
#define omcsqj_ omcsqj;
#define omcdmct_ omcdmct;
#define gwomcsqj_  gwomcsqj;
#define pmv_ pmv;
#define ivx_ ivx;
#define yv_ yv;
#define ym_ ym;
#define sev_ sev;
#define svvx_ svvx;
#define sv_ sv;
#define yvp_ yvp;
#define yvpir_ yvpir;
#define yvip_ yvip;
#define yvi_ yvi;
#define pmv_ pmv;
#define xv_ xv;
#define nv_ nv;
#define pvx_ pvx;
#define pmtv_ pmtv;
#define sevi_ sevi;
#define dgesvd_ dgesvd;
#define jmc_ jmc;
#define pmmt_ pmmt;
#define rm_ rm;
#define rn_ rn
#define dsyev_ dsyev;
#define pmx_ pmx;
#define pmx2_ pmx2;
#define pvx2_ pvx2;
#define mcm_ mcm;
#define katavep_ katavep;
#define dm_  dm;
#define ivi_ ivi;
#define ovtmcv_ ovtmcv;
#define pmrm_ pmrm;
#define pmv2_ pmv2;
#define yvlm_ yvlm;
#define nvinf_ nvinf;
#define nm1_ nm1;
#define nminf_ nminf;
#define ndvl2_ ndvl2;
#define prmm_ prmm;
#define zmtm_ zmtm;
#define pvxinv_ pvxinv;
#define pms2_ pms2;
#define slls_ slls;
#define svd_ svd;
#define ymp_ ymp;
#define yvpir2_ yvpir2;
#define pmvx_ pmvx;
#define concatmath_ concatmath ;
#define concatmatv_ concatmatv;
#define addmcz_ addmcz;
#define addmlz_ addmlz;
#define imdx_ imdx;
#define ympmp_ ympmp;
#define evmin_ evmin;
#define evmax_ evmax;
#define expnva_ expnva;
#define eigenvd_ eigenvd;
#define imxrn_ imxrn;
#define sdlscor_ sdlscor;
#define chol_ chol;
#define imx_ imx;
#define sdls_ sdls;
#define snorm_ snorm;

#endif 

void yvpir2_();
void slls_();
void svd_();

/* utmat */
void ymp_();
void pvxinv_();
void gwjms_();
void gwomcsqj_();
void jms_();
void omcsqj_();
void omcdmct_();
void pmv_();
void ivx_();
void yv_();
void svvx_();
void sv_();
void xv_();
void yvp_();
void yvpir_();
void yvip_();
void yvi_();
void ym_();
void pm_();
void pmv_();
void pmx_();
void pmx2_();
void pmtv_ ();
void pvx_();
void pvx2_();
void nv_();
void sevi_();
void sev_();
void jmc_();
void pmmt_();
void pmrm_();
void rm_();
void mcm_();
void katavep_();
void pmv2_();
void yvlm_();
void dm_();
void ivi_();
void ovtmcv_();
void nvinf_();
void nm1_();
void nminf_();
void ndvl2_();
void prmm_();
void zmtm_();
double dinvnr_();
void nv2_();
void imx_();
void pmvv_();
void svvx2_();
void pms2_();
void pmvx_();
void concatmath_();
void concatmatv_();
void addmcz_();
void addmlz_();
void imdx_();
void ympmp_();
void imprv_();
void imprm_();
void evmin_();
void evmax_();
void testsdp_();
void pmtm_();
void eigenvd_();
void imxrn_();
void chol_();
double snorm_();

/* statistics */
void expnva_();

/* random */
double rn_();

/* optim */
void sdlscor_();
void sdls_();

/* BLAS 1 */
double ddot_();
void dcopy_();
void daxpy_();
void dscal_();

/* BLAS 2 */
void dgemv_();
void dsyr_();

/* BLAS 3 */
void dsyrk_();

/* LAPACK */
void dposvx_();
void dgelss_();
void dgesvd_();
void dsyev_();
void dsyevr_();

double dsumdiv(
 int n, 
 double *x, 
 double *y
);

void dzero(
 int n, 
 double *p
);

void dupge(
 int n, 
 double *A
);
void dgapdev(
 int m, 
 int L, 
 int *N, 
 double *u, 
 double *z, 
 double *pgap, 
 double *pdev);

int printMatrix(char * text, int nl, int nc, double * matrix);

int printVector(char * text, int n, double * vector);

/* Monte Carlo VaR Simulation for Energy Contracts by "Direct Approach" */
int monteCarloVaR(
				/* inputs */
				int *p, int *k, int *N, double *alpha, double *epsSDLS,
				double *w, double *corr, double *vol, double *rho,
				/* outputs */
				double *VaR, double *X, double *Y, int *info);

/* Monte Carlo VaR Simulation for Energy Contracts by "Submarket Approach" */
int monteCarloBlockVaR(
				/* inputs */
				int *p, int *k, int *nblock, int *block, 
				int *N, double *alpha, double *epsSDLS,
				double *w, double *corr, double *vol, double *rho,
				/* outputs */
				double *VaR, double *X, double *Y, int *info);

int monteCarloBlockXVaR(
				/* inputs */
				int *p, int *nblock, int *block, int *k,
				int *N, double *alpha, double *epsSDLS,
				double *w, double *corr, double *vol, double *rho,
				/* outputs */
				double *VaR, double *X, double *Y, int *info);

int monteCarloBlockSimul(
				/* inputs */
				int *p, int *k, int *nblock, int *block, 
				int *N, double *corr, double *vol, double *rho,
				/* outputs */
				double *explainVariance, double *Y, int *info);

int monteCarloBlockVaRTest(
				/* inputs */
				int *p, int *k, int *nblock, int *block, 
				int *N, double *alpha, double *epsSDLS,
				double *w, double *corr, double *vol, double *rho,
				/* outputs */
				double *VaR, double *X, double *Y, double *F, double *G, int *info);

/* k-factor decomposition */
int monteCarloFactor(
				/* inputs */
				int *p, int *k, double *epsSDLS, double *corr,
				/* outputs */
				double *F, double *D, double *U, int *info);

int monteCarloFactorBlock(
				/* inputs */
				int *p, int *k, double *epsSDLS, double *corr,
				/* outputs */
				double *F, double *D, double *U, int *info);

/* portfolio returns and Gaussian Value-at-Risk */
int monteCarloPortfolioVaR(
				/* inputs */
				int *p, int *N, double *alpha,
				double *w, double *Y,
				/* outputs */
				double *VaR, double *X, int *info);
