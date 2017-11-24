c=======================================================================
c
c     subroutine SKTHGENOPT
c
c     Computing the function and gradient of 
c
c     Theta(dual) = 0.5*[wopt'*Q*wopt] + wopt'*plin + cste
c
c-----------------------------------------------------------------------
      SUBROUTINE skthgenopt ( n, Q, kapplin, rho, omega, omega0,
     &                        ndual, ydual, delta, mu, dist,
     &                        dwork, funct, grad )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n       : number of values                               integer
c       Q       : semi definite matrix (n*n)                      double
c       kapplin : vector(n) linear part of the objective funtcion double
c       rho     : mean returns (n)                                double 
c       omega   : vector(n) wights of the portfolio               double
c       omega0  : vector(n) second order portfolio                double
c       ndual   : number of variables                            integer
c       ydual   : dual solution lambda (ndual)                    double
c       delta   : maximal distance (with respect to the L2 norm) 
c                 between omega0 and the optimal portfolio        double
c       mu      : value of the target return                      double
c       dist    : distance (with respect to the L2 norm)
c                 between omega0 and omega                        double      
c
c     WORKSPACE 
c       dwork   : n                                               double
c
c     OUTPUT 
c       funct   : dual function value                             double
c       grad    : (2)-array gradient value                        double
c
c     CALL   
c       OVTMCV, XV, SVVX
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info, ndual
      DOUBLE PRECISION funct, cste, delta, mu, dist, Q(*), kapplin(*), 
     &                 omega(*), rho(*), ydual(*), omega0(*), grad(*)
c
c     workspaces
      DOUBLE PRECISION dwork(*)     
c      
c     local variables
      INTEGER pdw
      DOUBLE PRECISION var, scal, lin, perf, NEGONE   
      PARAMETER (NEGONE = -1.E0) 
c
c     external functions
      EXTERNAL OVTMCV, XV, SVVX      
c
c-----------------------------------------------------------------------     
c
c     pointers for double precision work space : dwork
c     ------------------------------------------------      
c     Double workspace 
      pdw = 1
c     
c     Total size of iwork array = n   
c
c----------------------------------------------------------------
c          
      CALL OVTMCV(n, Q, omega, var)  ! w'*Q*w 
      CALL XV(n, omega, kapplin, lin)
      funct = 0.5*var + lin + ydual(1)*mu + ydual(2)*dist 
      funct = -funct
      CALL XV(n, rho, omega, perf)
      grad(1) = -(mu - perf)
      CALL SVVX (n, omega, omega0, NEGONE, dwork(pdw))
      CALL XV(n, dwork(pdw), dwork(pdw), scal)
      grad(2)= -(scal - delta*delta)
      RETURN
      END
