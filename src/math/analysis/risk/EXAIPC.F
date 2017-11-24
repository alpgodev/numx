c=======================================================================
c     subroutine EXAIPC                                      
c
c     Intra portfolio correlation coefficient 
c     (ex-ante intra-portfolio correlation coef.)
c     
c                                     n
c                                    ---                              
c                                    \   w(i)*w(j)*corr(i,j) 
c                                    / 
c                                    ---
c                                   i,j=1
c                                   i<>j
c     intra-portfolio corr. coef. = ------------------------
c                                        n
c                                       ---                              
c                                       \   w(i)*w(j) 
c                                       / 
c                                       ---
c                                      i,j=1
c                                      i<>j 
c
c     Where w(i) is the fraction invested in asset i, w(j) is the fraction invested in asset j, 
c     corr(i,j) is the correlation between assets i and j, and n is the number of different assets.
c
c     A measure of -1 means that the assets within the portfolio perform perfectly oppositely: 
c     whenever one asset goes up, the other goes down. A measure of 0 means that the assets 
c     fluctuate independently, i.e. that the performance of one asset cannot be used to predict 
c     the performance of the others. A measure of 1, on the other hand, means that whenever 
c     one asset goes up, so do the others in the portfolio. To eliminate diversifiable risk completely, 
c     one needs an intra-portfolio correlation of -1.
c
c----------------------------------------------------------------------
      SUBROUTINE EXAIPC ( n, corr, w, coef, info)
c----------------------------------------------------------------------
c
c     INPUT 
c       n     : portfolio size (n>=1)                          integer
c       corr  : correlation matrix (n*n)                        double
c       w     : weights (n)                                     double
c
c     OUTPUT 
c       coef  : intra-portfolio correlation coefficient         double
c       info  : diagnostic argument                            integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info
      DOUBLE PRECISION coef, w(*), corr(n,*)
c
c     local variables
      INTEGER i, j
      DOUBLE PRECISION x, y, EPS
      PARAMETER ( EPS = 1.E-30 )
c
c-----------------------------------------------------------------------
c
c     initialization
      coef = 0.0
      info = 0
c
c----------------------------------------------------------------------
c
      x = 0.0
      y = 0.0
      DO i = 1,n
        DO j = 1,n
            IF (i .NE. j) THEN
                x = x + w(i)*w(j)*corr(i,j)
                y = y + w(i)*w(j)
            ENDIF    
        ENDDO
      ENDDO
      IF (y .GT. EPS) THEN
        coef = x/y
      ELSE
        coef = 1.0
      ENDIF      
c
      RETURN
      END

