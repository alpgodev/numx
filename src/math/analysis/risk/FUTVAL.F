c=======================================================================
c
c     subroutine FUTVAL                                     
c
c     Future values of a portfolio
c
c-----------------------------------------------------------------------
      SUBROUTINE FUTVAL (ndate, nasset, hdata, w, dwork, fvalue, info)
c-----------------------------------------------------------------------
c
c     INPUT 
c            ndate  : number of value(s)                         integer
c            nasset : portfolio size                             integer
c            hdata  : assets values (ndate*nasset)                double
c            w      : weights (nasset)                            double
c
c     WORKSPACE 
c            dwork  : vector (nasset)                             double
c
c     OUTPUT 
c            fvalue : future value of potfolio  (ndate)           double
c            info   : = 0 successful exit                        integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER ndate, nasset, info
      DOUBLE PRECISION hdata(ndate,*), w(*), dwork(*), fvalue(*)
c
c     local variables
      INTEGER i, j
      DOUBLE PRECISION sum, v, EPS
      PARAMETER ( EPS = 1.E-15 )
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
c      
      DO i= 1,nasset
         v = hdata(1,i)
         IF ( v .LT. EPS ) THEN
            info = -1
            RETURN
         ELSE
            dwork(i) = w(i) / v
         ENDIF
      ENDDO
c
c     computing history
      DO j = 1,ndate
         sum = 0.
         DO i = 1,nasset
            sum = sum + dwork(i) * hdata(j,i)
         ENDDO
         fvalue(j) = sum
      ENDDO
c
      RETURN
      END
