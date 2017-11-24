c=======================================================================
c
c     subroutine SKGENOPT
c
c     Computing the function and gradient of 
c     Theta(dual) = 0.5*[wopt'*Q*wopt] + wopt'*plin + cste
c
c-----------------------------------------------------------------------
      subroutine skgenopt ( n, qlin, rho, omega0, Q, ndual, dual,
     &                      neq, nin, ccst, bcst, cinf, csup,
     &                      iwork, dwork,
     &                      wopt, info)
c-----------------------------------------------------------------------
c
c     INPUT 
c       n       : number of values                               integer
c       qlin    : vector(n) linear part of the objective funtcion double
c       rho     : mean returns (n)                                double 
c       omega0  : vector(n) second order portfolio                double
c       Q     : semi definite matrix (n*n)                        double
c       ndual   : number of variables                            integer
c       dual   : dual solution lambda (ndual)                     double
c       neq     : number equality constraints                    integer
c       nin     : number inequality constraints                  integer
c       ccst    : matrix of constraints (n*(neq+nin))             double
c       bcst    : vector initial of constraints (neq+nin)         double
c       cinf    : lower bound (n)                                 double
c       csup    : upper bound (n)                                 double     
c
c     WORKSPACE 
c       iwork   : 3*n + 2*nin + neq + 1                           integer
c       dwork   : 2*n*n + 8*n + neq + 3*nin -array                 double
c
c
c     OUTPUT 
c       wopt    : optimal portfolio (n)                            double
c       info    : diagnostic argument                             integer
c
c     CALL   
c
c       buildthgen, qp
c
c-----------------------------------------------------------------------
c
      implicit none

c     Input arguments
      integer n, info, ndual, neq, nin
      double precision qlin(*), rho(*), Q(*), ccst(*), bcst(*)
      double precision dual(*), omega0(*), cinf(*), csup(*)
c     Output arguments      
      double precision wopt(*)
c     Workspace
      double precision dwork(*) 
      integer pdw, pdcov, pdlin, pdlagr
      
      integer iwork(*)
      integer piw
      
c     Local variables
      double precision var , opone, scal, lin, perf    
c     Integer workspace 
      piw = 1
c     needs 3*n + 2*nin + neq + 1
c       so 3*n + 2*nin + neq + 1
c     Total size is 3*n + 2*nin + neq + 1
            
c     Double workspace 
      pdcov = 1
c     needs n*n 
c        so n*n
      pdlin = pdcov + n*n
c     needs n
c        so n*n + n
      pdlagr = pdlin + n
c     needs n + neq + nin
c        so n*n + 2*n +neq +nin   
      pdw = pdlagr + n + neq + nin 
c     needs 3*n for buildthgen
c           n*n + 6*n + 2*nin for quapro
c      
c     Total size is 2*n*n + 8*n + neq + 3*nin      
      

       call buildthgen ( n, qlin, rho, omega0, Q, ndual, dual,
     &                   dwork(pdw), dwork(pdcov), dwork(pdlin))
     
       call qp ( n, dwork(pdcov), dwork(pdlin),
     &           neq, nin, ccst, bcst,
     &           cinf, csup, iwork(piw), dwork(pdw), dwork(pdlagr),
     &           wopt, info)
c
      return
      end
