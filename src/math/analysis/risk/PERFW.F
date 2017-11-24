c=======================================================================
c
c     subroutine PERFW                                    
c
c     This function computes the based-100 performances values 
c     of a portfolio.
c
c-----------------------------------------------------------------------
      SUBROUTINE PERFW (nasset, ndate, hdata, w, dwork, fvalue, info)
c-----------------------------------------------------------------------
c
c     INPUT 
c            nasset : portfolio size                             integer
c            ndate  : number of value(s)                         integer
c            hdata  : assets values (ndate*nasset)                double
c            w      : weights (ndate*nasset)                      double
c
c     WORKSPACE 
c            dwork  : vector (nasset)                             double
c
c     OUTPUT 
c            pvalue : past value of potfolio  (ndate)             double
c            info   : = 0 successful exit                        integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER ndate, nasset, info
      DOUBLE PRECISION hdata(ndate,*), w(ndate,*), dwork(*), fvalue(*)
c
c     local variables
      INTEGER i, j
      DOUBLE PRECISION sum, v, eps
      PARAMETER ( eps = 1.e-50 )
c
c     external subroutines
      EXTERNAL PASTVAL
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
c
      
      CALL PASTVAL ( nasset, ndate, hdata, w, dwork, pvalue, info)
c      
      DO i= 1,nasset
         v = hdata(1,i)
         IF ( v .LE. eps ) THEN
            info = -1
            RETURN
         ELSE
            dwork(i) = w(i) / v
         ENDIF
      ENDDO
c
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
c
c=======================================================================
