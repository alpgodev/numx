c=======================================================================
c
c     subroutine  CHECKSDLS
c
c     Computing the constraints matrices for SDLSCI
c
c-----------------------------------------------------------------------
      SUBROUTINE checksdls ( nmat, cmat, mcte, acte, bcte, mcti, acti,
     &                       lcti, ucti, dwork, info)
c-----------------------------------------------------------------------
c
c     INPUT 
c            nmat   : dimension of the matrix                    integer
c            cmat   : symmetric matrix(nmat*nmat) to optimize     double
c            mcte   : number of equal constraints                integer
c            acte   : mcte Ai symmetric matrices(nmat*nmat) of
c                     equal constraints, vector(mcte*nmat*nmat)   double
c            bcte   : value of equal constraints, vector(mcte)    double
c            mcti   : number of inequal constraints              integer
c            acti   : mcti Bi symmetric matrices(nmat*nmat) of
c                     inequal constraints  vector(mcti*nmat*nmat) double
c            lcti   : lower bounds of inequal constraints
c                     vector(mcti)                                double
c            ucti   : upper bounds of inequal constraints
c                     vector(mcti)                                double
c
c     WORKSPACE:
c            dwork  : double workspace needs nmat*nmat            double
c
c     OUTPUT 
c            info   : diagnostic argument                         integer
c
c     CALL   
c       checksdlseq
c       checksdlsineq
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o arguments
      INTEGER nmat, mcte, mcti, info
      DOUBLE PRECISION cmat(nmat,*), acte(*), bcte(*), acti(*), lcti(*)
      DOUBLE PRECISION ucti(*)
c
c     Workspace
      double precision dwork(*)
c
c-----------------------------------------------------------------------
c
      CALL checksdlseq ( nmat, cmat, mcte, acte, bcte,
     &                   dwork, info )
      IF (info .LT. 0) RETURN
c
      CALL checksdlsineq ( nmat, cmat, mcti, acti, lcti, ucti,
     &                     dwork, info)
c
      RETURN
      END
