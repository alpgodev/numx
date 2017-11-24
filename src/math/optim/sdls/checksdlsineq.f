c======================================================================
c
c     subroutine  CHECKSDLSINEQ
c
c     Test SDLS inequalities constraints
c
c----------------------------------------------------------------------
      SUBROUTINE checksdlsineq ( nmat, cmat, mcti, acti, lcti, ucti,
     &                           dwork, info)
c----------------------------------------------------------------------
c
c     INPUT 
c            nmat   : dimension of the matrix                    integer
c            cmat   : symmetric matrix(nmat*nmat) to optimize     double
c            mcti   : number of inequal constraints              integer
c            acti   : mcti Bi symmetric matrices(nmat*nmat) of
c                     inequal constraints  vector(mcti*nmat*nmat) double
c            lcti   : lower bounds of inequal constraints
c                     vector(mcti)                                double
c            ucti   : upper bounds of inequal constraints
c                     vector(mcti)                                double
c
c     WORKSPACE:
c            dwork  : double workspace (nmat*nmat)                double
c
c     OUTPUT 
c            info   : diagnostic argument                         integer
c
c     CALL   
c        PMC : squared matrix product
c        TM  : matrix trace
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o arguments
      INTEGER nmat, mcti, info
      DOUBLE PRECISION cmat(nmat,*), acti(*), lcti(*), ucti(*)
c
c     local variables
      INTEGER i, sm, ksai
      DOUBLE PRECISION EPS, trace
c
c     Workspace
      DOUBLE PRECISION dwork(*)
      INTEGER pdmac
c
      PARAMETER(EPS = 1.E-5 )
c
c-----------------------------------------------------------------------
c     
      pdmac = 1
c     pointer on a matrix needs nmat * nmat 
c     double workspace nmat * nmat
c      
      sm = nmat*nmat
      ksai = 1
      DO i=1,mcti
        CALL PMC (nmat, acti(ksai), cmat, dwork(pdmac))
        CALL TM (nmat, dwork(pdmac), trace)
        IF ((trace .lt. lcti(i)-EPS) .or. 
     &    (trace .gt. ucti(i)+EPS)) THEN
            info = -1110
            RETURN
        ENDIF
        ksai = ksai + sm
      ENDDO
c
      RETURN
      END
