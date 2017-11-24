c======================================================================
c                                                                      
c     subroutine EXARKU                                       
c                                                                      
c     Ex-ante (robust) Kurtosis (excess)                               
c                                                                                                                       
c-----------------------------------------------------------------------
      SUBROUTINE EXARKU ( n, p, w, x, cov, exakur, info)
c-----------------------------------------------------------------------
c
c     INPUT 
c       n      : number of value(s)                              integer
c       p      : number of asset(s)                              integer
c       w      : weights (p)                                      double
c       x      : values (n*p)                                     double
c       cov    : covariance matrix (p*p)                          double
c
c     OUTPUT 
c       exakur : kurtosis                                         double
c       info   : diagnostic argument                             integer
c
c     CALL 
c       EXARVO : ex-ante volatility
c       EXARET : ex-ante return
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     argument i/o
      INTEGER n, p, info
      DOUBLE PRECISION exakur
      DOUBLE PRECISION w(*), x(*), cov(*)
c
c     local variables
      INTEGER k
      DOUBLE PRECISION robvar, m1, m2, m3, m4
c
c     external subroutines
      EXTERNAL EXARVA, MOMENT
c
c-----------------------------------------------------------------------
c
c     initialization
      info   = 0
      exakur = 0.0    
c
c-----------------------------------------------------------------------
c
c     test number of values >= 4
      IF (n .LT. 4) THEN
        info = -2
        RETURN
      ENDIF  
c
c     ex-ante variance
      CALL EXARVA (p, cov, w, robvar, info)     
      IF (info .LT. 0) RETURN
c
c     1,2,3-th moment
      k = 1
      CALL MOMENT( n, p, k, w, x, m1, info) 
      IF (info .lt. 0) RETURN  
      k = 2
      CALL MOMENT( n, p, k, w, x, m2, info)  
      IF (info .lt. 0) RETURN  
      k = 3
      CALL MOMENT( n, p, k, w, x, m3, info) 
      IF (info .lt. 0) RETURN  
      k = 4
      CALL MOMENT( n, p, k, w, x, m4, info)   
      IF (info .lt. 0) RETURN  
c
c     ex-ante kurtosis  
      exakur = m4 - 3*(m1**4) + 6*(m1**2)*m2 - 4*m1*m3
      exakur = ((n*(n*n-2*n+3.))/((n-1.)*(n-2.)*(n-3.)))
     &         *(exakur /(robvar**2))
     &         -3.0*(n-1.)*(2*n-3.)/(n*(n-2)*(n-3.)) - 3.0
c
      RETURN
      END
