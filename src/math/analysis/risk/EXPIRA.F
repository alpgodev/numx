c=======================================================================
c
c     subroutine EXPIRA                                      
c
c     Ex-post Information ratio
c
c-----------------------------------------------------------------------
      SUBROUTINE EXPIRA ( ndate, retpor, retben, dwork, infrat, info)
c-----------------------------------------------------------------------
c
c     INPUT :
c            ndate  : estimation period                          integer
c            retpor : returns portfolio during the
c                     estimation period vector (ndate)            double
c            retben : returns benchmark during the
c                     estimation period vector (ndate)            double
c     WORKSPACE :
c            dwork  : vector ( ndate )                            double
c
c     OUTPUT :
c            infrat : information ratio                           double
c            info   : diagnostic argument                        integer
c
c     CALL   :
c            EXPTER : ex-post tracking error
c            EXPRET : ex-post return
c            MV     : mean of a vector
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER ndate, info
      DOUBLE PRECISION infrat
      DOUBLE PRECISION retpor(*), retben(*), dwork(*)
c
c     local variables
      DOUBLE PRECISION myzero, tracke, expr, mrben
c
c     external functions
      DOUBLE PRECISION spmpar
      EXTERNAL spmpar
c
c     external subroutines
      EXTERNAL EXPRET, MV
c
c-----------------------------------------------------------------------
c
c     machine precision
      myzero = spmpar(2)*(1.E+4)
c
c     ex-post tracking error
      CALL EXPTER ( ndate, retpor, retben, dwork, tracke, info)
      IF (info .lt. 0) RETURN
c  
c     test if Tracking Error is too small
      IF (tracke.le.myzero) THEN
         info = -1
         RETURN
      ENDIF
c
c     ex-post mean return
      CALL EXPRET ( ndate, retpor, expr, info )
      IF (info .lt. 0) RETURN
c
c     benchmark ex-post mean return
      CALL MV ( ndate, retben, mrben )
c
c     information ratio
      infrat = ( expr - mrben ) / tracke
c
      RETURN
      END

