/*=======================================================================

	Java-C-Fortran stubs

=======================================================================*/

/* f77 compilers tend to append the underscore. If the f77 name already
   has one underscore, g77 will add *two* underscores at the end. Define
   F77_ADD_UNDERSCORE in the Makefile. */

#ifndef F77_ADD_UNDERSCORE
# define F77_ADD_UNDERSCORE 1
#endif

#if F77_ADD_UNDERSCORE
# define F77_FUNCTION(f) f##_
# define F77_FUNCTION2(f) f##__
#else
# define F77_FUNCTION(f) f
# define F77_FUNCTION2(f) f
#endif

/* #define int jint  #define double jdouble */

/* Fortran header functions are declared here */

extern int F77_FUNCTION(func1) ();
extern int F77_FUNCTION(func2) (const char *, int *, float *, unsigned long);
extern int F77_FUNCTION2(f_routine) (int *, const char *, unsigned long);

/* Index tracking allocation */
void allocit_();
void allocitrfr_();

/* Maximum draw-down allocation */
void allocmdd_();

/* Risk-return allocation */
void allocmv_();
void allocmvrfr_();

/* Risk-return allocation with transaction costs */
void allocmvtc_();

/* Risk Budgetting allocation */
void allocrb_();
void allocrbrfr_();

/* Risk-return with Skewness constraint allocation */
void allocskew_();

/* Maximum Sharpe Ratio allocation */
void allocsr_();
void allocsrrfr_();

/* Asset Allocation subject to TE constraint */
void allocte_();

/* Asset Allocation - test feasibility */
void alloctest_();
void alloctestrfr_();

/* Maximum Value at Risk allocation */
void allocvar_();

/* Risk budgeting with volatility constraint */
void allocvol_();

/* APT model with constraints */
void aptcst_();

void avoc_();

void avodi_();

/* Bayes-Stein shrinkage estimator (expected returns) */
void bayesstein_();

/* Non-Linear Programming (BFGS) */
void bfgsbox_();
void bfgsboxs_();

void bmn_();

void buildir_();
void buildfeas_();
void buildpar_();
void buildrbit_();
void buildthgen_();

/* SDLS epsilon calibration (mean, median) */
void calepsmean_();
void calepsmedian_();

/* epsilon (eigenvalue) calibration */
void calepsvar_();

/* Calibration of APT model */
void calapt_();

/* cluster analysis */
void cluster_();

void checkfeasit_();
void checklinbox_();
void checkrbcstr_();
void checksdls_();
void checksdlseq_();
void checksdlsineq_();

/* co-kurtosis */
void cokurt_();

/* co-skewness */
void coskew_();

/* exponentially weighted correlation matrix */
void corexp_();

/* */
void corfact_();

/* correlation matrix */
void corm_();

/* exponentially weighted covariance matrix */
void covexp_();

/* covariance matrix with different values lengths */
void covl_();

/* cov. and corr. matrix with missing data */
void covlack_();
void corlack_();

/* covariance matrix */
void covm_();

/* cov. correction with filtering level */
void covfiltering_();
void covfiltering2_();
void covfiltering3_();

/* Constraint modelling */
void cstsec_();

/* eigenvalues */
void dgeev_();

/* diversification risk eliminated */
void diversification_();

/* Down side risk */
void downside_();

void epsinit_();
void epstrack_();

/* eigenvalues */
void eigenvd_();

/* exante concentration */
void exacc_();

/* exante Normal conditional Value at Risk */
void exacvar_();

/* exante intra-portfolio correlation */
void exaipc_();

/* exante modified Value at Risk */
void examva_();

/* exante modified Shrape ratio */
void examsh_();

/* exante normal Shortfall */
void exanshortfall_();

/* exante Normal Value at Risk */
void exanva_();

/* exante omega ratio */
void exaomega_();

/* exante return */
void exaret_();

/* exante information ratio */
void exarir_();

/* exante kurtosis */
void exarku_();

/* exante skewness */
void exarsk_();

/* exante variance */
void exarva_();

/* exante volatility */
void exarvo_();

/* exante Sharpe ratio */
void exasra_();

/* exante STARR ratio */
void exastarr_();

/* exante tracking error */
void exater_();

/* exante value at risk */
void exavri_();

/* expost correlation coeficient */
void expcoefcorr_();

/* expost Normal conditional VaR */
void expcvar_();

/* expost down number */
void expdnb_();

/* expost down capture */
void expdownc_();

/* expost down percent */
void expdpercent_();

/* expost average gain */
void expgmean_();

/* expost gain percent */
void expgpercent_();

/* expost information ratio */
void expira_();

/* expost kurtosis */
void expkur_();

/* expost average loss */
void explmean_();

/* expost maximum loss */
void expmaxloss_();

/* expost modified Sharpe ratio */
void expmsh_();

/* expost modified VaR */
void expmva_();

/* expost non parametric VaR */
void expnpvar_();

/* expost non parametric CVaR */
void expnpcvar_();

/* expost normal shortfall */
void expnshortfall_();

/* expost Normal VaR */
void expnva_();

/* expost omega ratio */
void expomega_();

/* expost return */
void expret_();

/* expost skewness */
void expske_();

/* expost sortino ratio */
void expsor_();

/* expost sharpe ratio */
void expsra_();

/* expost STARR ratio */
void expstarr_();

/* expost tracking error */
void expter_();

/* expost up capture */
void expupc_();

/* expost up number */
void expupnb_();

/* expost up percent */
void expuppercent_();

/* expost variance */
void expvar_();

/* expost volatility */
void expvol_();

/* expost exponentially weighted volatility */
void expvolexp_();

/* expost value at risk */
void expvri_();

/*  maximal element of a vector      */
void evmax_();

/*Fills holes in a matrix of doubles
	(replacing the unknown value by linear interpolation)
	and identifies unborn and dead assets*/
void fml_();

/* future values */
void futval_();

/* random generators */
double genbet_();
double genchi_();
double genexp_();
double genf_();
double gengam_();
double gennor_();
double genunf_();

/* density function of a multivariate normal random variable */
void getmnormalpdf_();

/* density function of a standard normal random variable */
double getnormalpdf_();

void greps_();

/* workspace (getwork) */
void gwallocit_();
void gwallocitrfr_();
void gwallocmdd_();
void gwallocmv_();
void gwallocmvrfr_();
void gwallocmvtc_();
void gwallocrb_();
void gwallocrbrfr_();
void gwallocskew_();
void gwallocsr_();
void gwallocsrrfr_();
void gwallocte_();
void gwalloctest_();
void gwalloctestrfr_();
void gwallocvar_();
void gwallocvol_();

void gwaptcst_();

void gwbayesstein_();

void gwbfgsbox_();
void gwbfgsboxs_();

void gwcalapt_();

void gwcalepsmedian_();

void gwcalepsvar_();

void gwchttrapca_();

void gwcluster_();

void gwcokurt_();

void gwcoskew_();

void gwcorfact_();

void gwcorm_();

void gwcorlack_();

void gwcovexp_();

void gwcovfiltering_();
void gwcovfiltering2_();
void gwcovfiltering3_();

void gwcovl_();

void gwcovlack_();

void gwcovm_();

void gweigenvd_();

void gwepsinit_();

void gwexamsh_();

void gwexamva_();

void gwexaomega_();

void gwexarir_();

void gwexarku_();

void gwexarsk_();

void gwexater_();

void gwexpcorrcoef_();

void gwexpdownside_();

void gwexpira_();

void gwexpnpvar_();

void gwexpnpcvar_();

void gwexpter_();

void gwfactevol_();

void gwfutval_();

void gwimplret_();

void gwjamesstein_();

void gwkb_();

void gwkatave_();

void gwkatoblock_();

void gwkatoblockx_();

void gwmcbase_();

void gwmcbcor_();

void gwmcbhis_();

void gwmcfact_();

void gwmcluster_();

void gwmultivol_();
void gwmultivolrfr_();

void gwmultireg_();

void gwmve_();

void gwmvgarch_();

void gwnavem_();

void gwnls_();

void gwnlsb_();

void gwnlsx_();

void gwnpreg_();

void gwpastval_();

void gwpca_();

void gwpegbmx_();

void gwqp_();

void gwrcho_();

void gwriskbudgetit_();

void gwrnls_();

void gwrnlsb_();

void gwrnlsx_();

void gwrollingpies_();

void gwschurpca_();

void gwschuro_();

void gwsdls_();

void gwsdlsce_();

void gwsdlsci_();

void gwsdlscor_();

void gwsdlsg_();

void gwsdlsgen_();

void gwsdlstrace_();

void gwshrinkagecov_();

void gwsmbox_();

void gwtfact_();

void gwtrapca_();

/* logarithmic returns */
void hlretm_();

/* */
void initfeas_();

/* */
void interpzcp_();

/* implied returns */
void implret_();

/* random generator */
void init_seed_linux_();

/* */
void ivx_();

/* James-Stein shrinkage estimator (expected returns) */
void jamesstein_();

/* semi-Kato decomposition */
void katave_();

/* Kato matrix correction (by block) */
void katoblock_();
void katoblockx_();

/* Gaussian kernel */
void kerneld_();
void kernela_();

/* Log-returns with lack of data */
void logrlack_();

/* Median */
void median_();

/* Monte Carlo */
void mcbase_();

/* Monte Carlo */
void mcbcor_();

/* Monte Carlo */
void mcbhis_();

/* Monte Carlo */
void mcfact_();

/* Arithmetic returns */
void mhoret_();

/* k-th moment (k=1,2,3,4) */
void moment_();

/* Multivariate regression */
void multireg_();

/* risk budgeting asset allocation */
void multivol_();
void multivolrfr_();

/* Efficient frontier */
void mve_();

/* GARCH(1,1) cov. matrix */
void mvgarch_();

/* EM algorithm for missing data */
void navem_();

void ndm_();

/* NLS solvers */
void nls_();
void nlsb_();
void nlsx_();

/* vector norms */
void nm_();
void nminf_();
void nm1_();

void npreg_();

/* past values */
void pastval_();

/* Simulation of a Brownian process (Wiener process) */
void pebm_();

/* Simulation of GARCH(p,q) process */
void pegarch_();

/* Simulation of a Geometric Brownian Motion */
void pegbm_();

/* Simulation of a multivariate Brownian process */
void pegbmx_();

void pemrbm_();

/* Simulation of a Standard Brownian Motion */
void pesbm_();

/* Principal Component Analysys (PCA) */
void pca_();

void pm_();

void pmc_();

void pmmt_();

void pmtm_();

void pmtv_();

void pmv_();

void pmx_();

/* Option vanille pricing */
void prbinomial_();
void prvanilla_();

void pvx_();

/* Linear Quadratic Programming */
void qp_();

/* Performance contribution */
void retcontrib_();

/* arithmetic returns with lack data */
void retlack_();

/* Robust Choleski */
void rcho_();

/* Risk budgeting with tracking error constraint */
void riskbudgetit_();

void rm_();

double rn_();

/* Robust NLS solvers */
void rnls_();
void rnlsb_();
void rnlsx_();

/* Rolling pies */
void rollingpies_();

/* historical weights */
void rolpies_();

/* Schur decomposition */
void schurpca_();
void schuro_();

/* SDLS solvers */
void sdls_();
void sdlsce_();
void sdlsci_();
void sdlscor_();
void sdlsg_();
void sdlsgen_();
void sdlstrace_();

/* expost down side (semi-volatility) */
void semvol_();

/* Ledoit-Wolf shrinkage cov. matrix */
void shrinkagecov_();

void simext_();

void sm_();

void smoothbox_();

/* Tracking error contribution */
void tecontrib_();

void tfact_();

void tm_();

void trapca_();

/* skewness allocation utilities */
void utskew1_();
void utskew2_();
void utskew3_();

/* Value-at-Risk contribution */
void varcontrib_();

/* Volatility contribution */
void volcontrib_();

void xm_();

void xv_();

void wveq_();

/* ex-post analysis with width */
void wexpret_();
void wdownside_();
void wexpcorrcoef_();
void wexpcvar_();
void wexpdownc_();
void wexpira_();
void wexpkur_();
void wexpmaxloss_();
void wexpmsh_();
void wexpmva_();
void wexpnpvar_();
void wexpnshortfall_();
void wexpnva_();
void wexpske_();
void wexpsor_();
void wexpsra_();
void wexpter_();
void wexpvar_();
void wexpvol_();
void wexpvolexp_();
void wexpvri_();
void wkurtos_();
void wsemvar_();
void wskewn_();
