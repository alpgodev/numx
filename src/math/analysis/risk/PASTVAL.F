c=======================================================================
c
c     subroutine PASTVAL                                  
c
c     Past values of a portfolio
c
c-----------------------------------------------------------------------
      SUBROUTINE PASTVAL ( ndate, nasset, hdata, w, dwork, pvalue, info)
c-----------------------------------------------------------------------
c
c     INPUT 
c       ndate  : number of dates/values                          integer
c       nasset : number of assets                                integer
c       hdata  : asset(s) values (ndate*nasset)                   double
c       w      : weights (nasset)                                 double
c
c     WORKSPACE 
c       dwork  : nasset                                           double
c
c     OUTPUT 
c       pvalue : portfolio past value (ndate)                     double
c       info   : diagnostic argument                             integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER ndate, nasset, info
      DOUBLE PRECISION hdata(ndate,*), w(*), pvalue(*)
c
c     workspaces
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER i, j
      DOUBLE PRECISION sum, lastValue, EPS
      PARAMETER ( EPS = 1.E-15 )
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0

c     computing number of assets at the date
      DO i = 1,nasset
         lastValue = hdata(ndate,i)
         IF (lastValue .LT. EPS ) THEN
            info = -1
            RETURN
         ELSE
            dwork(i) = w(i) / lastValue
         ENDIF
      ENDDO
c
c     computing history
      DO j = 1,ndate
         sum = 0.0
         DO i = 1,nasset
            sum = sum + dwork(i)*hdata(j,i)
         ENDDO
         pvalue(j) = sum
      ENDDO
c
      RETURN
      END
