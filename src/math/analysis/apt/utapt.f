c=======================================================================
c
c     APT Utilities
c
c        RAPT    : Y = alpha + beta*X
c        RAPTV   : computes assets values as a function of values
c                  of factors with the APT regression and volatility
c        RHAPT   : computes historical assets values as a function
c                  of values of factors with the APT regression
c        RHAPTV  : computes historical assets values as a function
c                  of values of factors with the APT regression
c                  and volatility
c
c-----------------------------------------------------------------------
c
c     subroutine RAPT
c
c     computes assets values as a function of values of factors
c     with the APT regression
c
c-----------------------------------------------------------------------
      SUBROUTINE rapt ( p, q, X, alpha, beta, Y )
c-----------------------------------------------------------------------
c 
c     INPUT 
c            p     : number of assets                           integer
c            q     : number of factors                          integer
c            X     : factor(s) values (q)                        double
c            alpha : alpha coefficient(s) (p)                    double
c            beta  : beta  coefficient(s) (q*p)                  double
c            
c     OUTPUT 
c            Y     : asset(s) values (nasset)                    double
c                        
c     CALL   
c            PMTV  : M'*V = vector ( M matrix(n*m), V vector(n), M'*V vector(m) )
c            SV    : sum of 2 vectors
c            
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER p, q 
      DOUBLE PRECISION alpha(*), beta(*), X(*), Y(*)
c
c     external subroutines
      EXTERNAL PMTV, SV      
c
c------------------------------------------------------------------------
c
c     Y = beta*X
      CALL PMTV ( q, p, beta, X, Y )
c
c     Y = alpha + beta*X
      CALL SV ( p, Y, alpha, Y )
c
      RETURN
      END
c
c=======================================================================
c
c     subroutine RAPTV
c
c     computes assets values as a function of values of factors
c     with the APT regression and volatility
c
c-----------------------------------------------------------------------
      SUBROUTINE raptv ( p, q, X, alpha, beta, vol, Y )
c-----------------------------------------------------------------------
c 
c     INPUT 
c            p     : number of assets                           integer
c            q     : number of factors                          integer
c            X     : factor(s) values (q)                        double
c            alpha : alpha coefficient(s) (p)                    double
c            beta  : beta  coefficient(s) (q*p)                  double
c            vol   : regression error volatility  (p)            double
c            
c     OUTPUT 
c            Y     : asset(s) values (p)                         double
c                        
c     CALL    
c            PMTV  : M'*V = vector ( M matrix(n*m), V vector(n), M'*V vector(m) )
c            rn    : generates standard Normal N(0,1)
c            
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER p, q
      DOUBLE PRECISION alpha(*), beta(*), vol(*), X(*), Y(*)
c
c     local variables
      INTEGER i
c
c     external functions            
      DOUBLE PRECISION rn
      EXTERNAL rn
c
c     external subroutines
      EXTERNAL PMTV      
c
c-----------------------------------------------------------------------
c
c     Y = beta*X
      CALL PMTV ( q, p, beta, X, Y )
c
c     Y = alpha + beta*X + alpha + vol*N(0,1)
      DO i = 1,p
         Y(i) = Y(i) + alpha(i) + rn()*vol(i)
      ENDDO
c
      RETURN
      END
c
c=======================================================================
c
c     subroutine RHAPT
c
c     computes historical assets values as a function of values
c     of factors with the APT regression
c
c-----------------------------------------------------------------------
      SUBROUTINE rhapt ( n, p, q, X, alpha, beta, Y )
c-----------------------------------------------------------------------
c 
c     INPUT 
c            n      : number of value(s) (n >= 1)                integer
c            p      : number of asset(s)                         integer
c            q      : number of factors                          integer
c            X      : factor(s) values (n*q)                      double
c            alpha  : alpha coefficient(s) (p)                    double
c            beta   : beta coefficent(s) (q*p)                    double
c            
c     OUTPUT 
c            Y      : asset(s) values (n*p)                       double
c                        
c     CALL    
c            PM      : product of 2 vectorized full matrices
c                     (M matrix n*m, N matrix m*c -> M*N matrix n*c)
c            SVX     : sum of a vector with a scalar
c            
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, p, q 
      DOUBLE PRECISION alpha(*), beta(*), X(*), Y(n,*)
c
c     local variables
      INTEGER i
c
c     external subroutines
      EXTERNAL PM, SVX            
c
c-----------------------------------------------------------------------
c
c     Y = beta*X
      CALL PM ( n, q, p, X, beta, Y )
c
c     Y = alpha + beta*X
      DO i = 1,p
         CALL SVX ( n, Y(1,i), alpha(i), Y(1,i) )
      ENDDO
c
      RETURN
      END
c
c=======================================================================
c
c     subroutine RHAPTV                                      
c
c     computes historical assets values as a function of values
c     of factors with the APT regression and volatility
c
c-----------------------------------------------------------------------
      SUBROUTINE rhaptv ( n, p, q, X, alpha, beta, vol, Y )
c-----------------------------------------------------------------------
c 
c     INPUT 
c            n     : number of value(s) (n>=1)                  integer 
c            p     : number of assets                           integer
c            q     : number of factors                          integer
c            X     : factor(s) values (n*q)                      double
c            alpha : alpha coefficient(s) (p)                    double
c            beta  : beta  coefficient(s) (q*p)                  double
c            vol   : regression error volatility  (p)            double
c            
c     OUTPUT 
c            Y     : asset(s) values (n*p)                       double
c                        
c     CALL   
c            PM    : product of 2 vectorized full matrices
c                    (M matrix n*m, N matrix m*c -> M*N matrix n*c )
c            SVX   : sum of a vector with a scalar
c            BCMN  : generates a Normal random on the columns of a matrix
c            
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, p, q 
      DOUBLE PRECISION alpha(*), beta(*), vol(*), X(*), Y(n,*)
c
c     local variables
      INTEGER i
c
c     external subroutines
      EXTERNAL PM, SVX, BCMN           
c
c-----------------------------------------------------------------------
c
c     Y = beta*X
      CALL PM ( n, q, p, X, beta, Y )
c
c     Y = alpha + beta*X
      DO i = 1,p
         CALL SVX ( n, Y(1,i), alpha(i), Y(1,i) )
      ENDDO
c
c     Y = alpha + beta*X + vol*N(0,1)
      CALL BCMN ( n, p, Y, vol, Y )
c
      RETURN
      END


