c=======================================================================
c
c     subroutine EXARIR                                      
c
c     Ex-ante information ratio
c
c-----------------------------------------------------------------------
      SUBROUTINE EXARIR (ndate, nasset, covemp,
     &                    rmean, covben, retben, varben, weight, dwork,
     &                    infrat, info)
c-----------------------------------------------------------------------
c
c     INPUT 
c            nasset : number of asset(s)                         integer
c            ndate  : estimation period                          integer
c            weight : weight(s) (nasset)                          double
c            covemp : covariance matrix (nasset*nasset)           double
c            covben : covariance of benchmark (nasset)            double
c            varben : variance of benchmark                       double
c            rmean  : asset(s) mean return(s) (nasset)            double
c            retben : benchmark returns (ndate)                   double
c
c     WORKSPACE 
c            dwork  : nasset*(2*nasset + 3) + 2                   double
c
c     OUTPUT 
c            infrat : information ratio                           double
c            info   : diagnostic argument                        integer
c
c     CALL   
c           EXATER  : computes ex-ante robust tracking error
c           EXARET  : computes ex-ante return
c           MV      : computing the mean of a vector
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER nasset, ndate,  info
      DOUBLE PRECISION varben, infrat
      DOUBLE PRECISION covemp(*), weight(*), covben(*), rmean(*), 
     &                 retben(*)
c
c     workspaces
      DOUBLE PRECISION dwork(*)
c
c     local variables
      DOUBLE PRECISION tracke, exantr, mrben, EPS
      PARAMETER ( EPS = 1.E-30 )
c
c     external subroutines
      EXTERNAL EXARET, EXATER, MV
c
c-----------------------------------------------------------------------
c
c     initialisation
      infrat = 0
c 
c     computing ex-ante tracking error
      CALL EXATER ( nasset, covemp, covben, varben, weight, dwork,
     &              tracke, info)   
      IF (info .LT. 0) RETURN
c
c     test if Tracking Error < EPS
      IF (tracke .LT. EPS) THEN
         info = -105
         RETURN
      ENDIF
c
c     ex-ante return
      CALL EXARET ( nasset, rmean, weight, exantr, info )
      IF (info .LT. 0) RETURN
c
c     benchmark mean
      CALL MV ( ndate, retben, mrben )
c
c     computing information ratio
      infrat = ( exantr - mrben ) / tracke
c
      RETURN
      END
