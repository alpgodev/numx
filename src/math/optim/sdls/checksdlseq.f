c=======================================================================
c
c     subroutine  CHECKSDLS                                   
c
c     tests SDLS equalities constraints
c
c-----------------------------------------------------------------------
      SUBROUTINE checksdlseq ( n, cmat, mcte, acte, bcte,
     &                         dwork, info)
c-----------------------------------------------------------------------
c
c     INPUT 
c       n     : matrix size                                      integer
c       cmat  : optimized matrix (n*n)                            double
c       mcte  : number of equal constraints                      integer
c       acte  : mcte Ai symmetric matrices(n*n) of
c               equal constraints (mcte*n*n)                      double
c       bcte  : vector value of equal constraints (mcte)          double
c
c     WORKSPACE
c       dwork : n*n                                               double
c
c     OUTPUT 
c       info  : diagnostic argument                              integer
c
c     CALL   
c       PMC :
c       TM  :
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o arguments
      INTEGER n, mcte, info
      DOUBLE PRECISION cmat(n,*), acte(*), bcte(*)
c
c     workspaces
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER i, sm, ksae, pdmac
      DOUBLE PRECISION trace, diff, EPS        
      PARAMETER (EPS =1.E-5 )
c
c-----------------------------------------------------------------------
c     
      pdmac = 1 
c     Total size of dwork array = n*n 
c      
      sm = n*n
      ksae = 1
      DO i = 1,mcte
        CALL PMC (n, acte(ksae), cmat, dwork(pdmac))
        CALL TM(n, dwork(pdmac), trace)
        diff = trace - bcte(i)
        diff = ABS(diff)
        IF (diff .GT. EPS) THEN
            info = -1109
            RETURN
        ENDIF
        ksae = ksae + sm
      ENDDO
      RETURN
      END
