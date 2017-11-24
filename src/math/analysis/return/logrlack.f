c=======================================================================
c
c     subroutine LOGRLACK
c
c     Log h-returns for prices matrix
c
c           r(t+h) = log[price(t+h)/price(t)] for t = 1,...,n-h
c
c     Take into account of missing data (lack of data), 
c     with historical data structure
c
c     X(n,p) = | x(1,1), ..., x(1,p) |  
c              | x(2,1), ..., -1000  |
c              | x(3,1), ..., x(3,p) |
c              |   .   , ...,   .    |
c              |   .   , ...,   .    |
c              | x(n,1), ..., x(n,p) |
c
c     where x(2,p) = -1000 denote a missing return
c
c-----------------------------------------------------------------------
      SUBROUTINE logrlack ( n, p, X, h, R, info )
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
      DOUBLE PRECISION missing, price, eps
      PARAMETER ( eps = 1.E-30 )
c     
c     intrinsic functions
      INTRINSIC log         
c
c-----------------------------------------------------------------------
c
c     initializations
      info = 0
      missing = -1000
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
            IF ((price .LT. eps) .AND. (price .GT. missing)) THEN
               info = -103
               RETURN
            ENDIF
            IF ((X(ipri +h,jpri) .GT. missing).AND.(price .GT. missing))
     &      THEN
                R(k) = log( X(ipri + h,jpri) / price )
            ELSE
                R(k) = -1000
            ENDIF    
            k = k + 1
         ENDDO
      ENDDO
c
      RETURN
      END
