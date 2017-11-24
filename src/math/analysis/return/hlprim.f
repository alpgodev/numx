c=======================================================================
c
c     subroutine HLPRIM
c
c     This function computes the matrix of values with log-returns
c
c         price(t+1) = price(t)*exp[x(t)] fot t=1,...,n-1
c
c-----------------------------------------------------------------------
      SUBROUTINE hlprim ( n, p, logret, inipri, m, mprice, info )
c-----------------------------------------------------------------------
c
c     INPUT :
c            n      : number of log-return(s) (n >= 1)           integer
c            p      : number of asset(s) (p >= 1)                integer
c            logret : log-return (n*p)                            double
c            inipri : initial price(s) (p)                        double
c
c     OUTPUT :
c            m      : number of price(s)                         integer
c            mprice : price(s) (n*p)                              double
c            info   : diagnostic argument                        integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o 
      INTEGER n, p, m, info
      DOUBLE PRECISION mprice(n + 1,*), logret(n,*), inipri(*)
c
c     local variables
      INTEGER i, j
c     
c     intrinsic functions
      INTRINSIC exp       
c
c-----------------------------------------------------------------------
c
c     initialization 
      info = 0
      m    = n + 1
c
      DO j = 1,p
         mprice(1,j) = inipri(j)
         DO i = 1,n
            mprice(i+1,j) = mprice(i,j) * exp (logret(i,j))
         ENDDO
      ENDDO
c
      RETURN
      END

