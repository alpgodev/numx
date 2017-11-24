c=======================================================================
c
c     subroutine HLRETM
c
c     This function computes the log h-returns from matrix of values
c
c           r(t+h) = log[price(t+h)/price(t)] for t = 1,...,n-h
c
c-----------------------------------------------------------------------
      SUBROUTINE hlretm ( n, p, X, h, R, info )
c-----------------------------------------------------------------------
c
c     INPUT :
c            n      : number of prices (n > h)                   integer
c            p      : number of assets (p >= 1)                  integer
c            X      : assets price (n*p)                          double
c            h      : horizon (h > 0)                            integer
c
c     OUTPUT :
c            R      : log-returns ( n-h * p)                      double
c            info   : diagnostic argument                        integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, p, h, info
      DOUBLE PRECISION X(n,*), R(*)
c
c     local variables
      INTEGER ipri, jpri, k, nl
      DOUBLE PRECISION price, myzero
      PARAMETER ( myzero = 1.E-15 )
c     
c     intrinsic functions
      INTRINSIC log         
c
c-----------------------------------------------------------------------
c
c     initializations
      info = 0
      nl   = n - h
      k    = 1
c
c     test the number of points
      IF (n .le. h) THEN
        info = -2
        RETURN
      ENDIF
c
c     log return(s)      
      DO jpri = 1,p
         DO ipri = 1,nl
            price = X(ipri,jpri)
            IF (price .lt. myzero) THEN
               info = -103
               RETURN
            ENDIF
            R(k) = log( X(ipri + h,jpri) / price )
            k = k + 1
         ENDDO
      ENDDO
      RETURN
      END

