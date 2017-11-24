c=======================================================================
c     subroutine EXPCOEFCORR                                      
c
c     Correlation coefficient
c     (ex-post correlation coef.)
c
c                        Cov(X,Y)
c     corr. coef. = -------------------
c                   sqrt(Var(X)*Var(Y))
c 
c----------------------------------------------------------------------
      SUBROUTINE EXPCOEFCORR ( n, x, y, dwork, coef, info)
c----------------------------------------------------------------------
c
c     INPUT 
c        n     : size of x (n > 1)                              integer
c        x     : x's value(s) (n)                                double
c        y     : y's value(s) (n)                                double
c
c     WORKSPACE 
c        dwork : 4*n + 6                                         double
c
c     OUTPUT 
c        coef  : ex-post correlation coefficient                 double
c        info  : diagnostic argument                            integer
c
c     CALL   
c        COVMVM  : covariance
c        VOLAT : volatility (standard deviation)
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info
      DOUBLE PRECISION coef, x(*), y(*)
c
c     workspaces
      DOUBLE PRECISION dwork(*)     
c
c     local variables
      INTEGER p, pdz, pdw, pdrho, pdcov
      DOUBLE PRECISION volx, voly, cov, EPS
      PARAMETER ( EPS = 1.E-30 )
c
c     external subroutines
      EXTERNAL COVMVM, VOLAT
c
c-----------------------------------------------------------------------
c
c     initialization
      coef = 0.0
      info = 0
      p    = 2
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdz = 1
c     pdz : pointer for Z = [X,Y], ( 2*n )
      pdw = pdz + ( 2*n )
c     pdw : pointer for COVMVM dwork, ( 2*n )
      pdrho = pdw + ( 2*n )
c     pdrho : pointer for RHO vector, ( 2 )      
      pdcov = pdrho + ( 2 )
c     pdcov : pointer for COV matrix, ( 2*2 )      
c
c     Total size of dwork array = 4*n + 6 
c
c----------------------------------------------------------------------
c
c     construct Z = [X, Y]
      CALL YV ( n, x, dwork(pdz) )
      CALL YV ( n, y, dwork(pdz+n) )
c
c     covariance
      CALL COVMVM ( n, p, dwork(pdz), dwork(pdw), dwork(pdrho), 
     &            dwork(pdcov), info)
      cov = dwork(pdcov + 2)
c
c     volatility
      CALL VOLAT ( n, x, volx, info )
      CALL VOLAT ( n, y, voly, info )
c
c     correlation coefficent 
      IF ((volx .LT. EPS) .OR. (voly .LT. EPS)) THEN
        coef = 0.0
      ELSE
        coef = cov/(volx*voly)
      ENDIF  
c
      RETURN
      END

