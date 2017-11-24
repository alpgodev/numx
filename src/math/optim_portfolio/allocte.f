c=======================================================================
c
c     subroutine ALLOCTE                                    
c
c     Computing optimal Mean-Variance with Tracking Error constraint
c
c     Max[ rho*w - rhob]          (maximize relative target return)
c      s.t. 
c     Tracking Error(w) <= TEmax  (tracking error constraint)
c     C*w <= b                    (linear constraints) 
c     Cinf <= w <= Csup           (lower/upper bounds)
c
c        w     : portfolio weights 
c        Q     : covariance matrix 
c        rho   : assets performance 
c        rhob  : benchmark performance
c        TEmax : maximum tracking error
c
c-----------------------------------------------------------------------
      SUBROUTINE allocte ( n, cov, rho, covb, varb, rhob,
     &                     neq, nin, ccst, bcst, cinf, csup,
     &                     temax, maxdic, epsdic,
     &                     iwork, dwork,
     &                     wopt, delopt, optte, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c        n      : number of asset(s) (> 1)                       integer
c        cov    : n*n-array, covariance matrix (n*n)              double  
c        rho    : n-array, expected return(s) (n)                 double
c        covb   : n-array, covariance asset vs.index (n)          double  
c        varb   : variance of index/benchmark                     double
c        rhob   : expected return of index/benchmark              double
c        neq    : number equality constraints                    integer
c        nin    : number inequality constraints                  integer
c        ccst   : matrix of linear constraints (n*(neq+nin))      double
c        bcst   : vector of linear constraints (neq+nin)          double
c        cinf   : lower bounds (n)                                double
c        csup   : upper bounds (n)                                double
c        temax  : tracking error constraint (>0)                  double
c        maxdic : maximum iteration (dichotomy)                  integer
c        epsdic : precision (dichotomy)                           double
c
c     WORKSPACE 
c        iwork  : 16*n + 2*nin + neq + 1                         integer 
c        dwork  : n*(7*n + neq + nin + 42 )+ 2*neq + 3*nin + 4    double
c
c     OUTPUT 
c        wopt   : optimal portfolio vector (nasset)              double
c        info   : = 0 successful exit                           integer
c                 = -107 : tracking error constraint too small
c                 = 106  : max. iteration reached - solution non-optimal
c
c     CALL   
c        OPTE    : computing optimization (cf. UTALLOC.F)
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, neq, nin, maxdic, info
      DOUBLE PRECISION rhob, varb, temax, epsdic, epskat, delopt, optte
      DOUBLE PRECISION cov(*), rho(*), cinf(*), csup(*), ccst(*), 
     &                 bcst(*), wopt(*), covb(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER piwo, pdwo
c
c     external subroutines
      EXTERNAL opte
c
c-----------------------------------------------------------------------
c
c     initializations 
      info = 0
c
c     integer workspace : iwork
c     -------------------------
      piwo = 1
c     piwo  : pointer for internal workspaces of OPTE
c             needs ( 16*n + neq + 2*nin + 1 )
c
c     Total size of iwork array = ( 16*n + 2*nin + neq + 1 )
c
c     double work space : dwork
c     -------------------------
      pdwo = 1 
c     pdwo   : pointer for internal workspace of OPTE
c              needs  n*(7*n + neq + nin + 42 )  + 2*neq + 3*nin + 4 
c     
c     Total size   = n*(7*n + nin + neq + 42) + 2*neq + 3*nin + 4
c
c---------------------------------------------------------------------------
c
c     optimization (see UTALLOC.F)
      CALL opte ( n, cov, rho, covb, varb, rhob,
     &             neq, nin, ccst, bcst, cinf, csup,
     &             temax, maxdic, epsdic,
     &             iwork(piwo), dwork(pdwo),
     &             wopt, delopt, optte, info )
      RETURN
      END
