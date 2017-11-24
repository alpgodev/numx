c=======================================================================
c
c     subroutine NORMINV                                     
c
c     Inverse distribution of N(m,s).       
c
c----------------------------------------------------------------------
      SUBROUTINE NORMINV(p, mean, std, zc, info)
c----------------------------------------------------------------------
c
c     INPUT 
c            p    : probability                                 double
c            mean : mean                                        double
c            std  : std                                         double    
c    
c     OUTPUT :
c            zc   : p-quantile                                  double
c            info : diagnostic argument                         integer 
c  
c    CALL
c          dinvnr : N(0,1) distribution inverse
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o      
      INTEGER info
      DOUBLE PRECISION p, mean, std, zc
c
c     local variables      
      DOUBLE PRECISION q
c
c     external functions
      DOUBLE PRECISION dinvnr
      EXTERNAL dinvnr
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
      zc   = 0.0
c
      IF ((p .LT. 0.).or.(p .GT. 1.)) THEN
         info = -102
         RETURN
      ENDIF
c
c     probability 1-p
      q = 1-p
c
c     p-quantile of N(0,1)
      zc = dinvnr( p, q )
c
c     p-quantile of N(m,s)      
      zc = std*zc + mean
c
      RETURN 
      END
