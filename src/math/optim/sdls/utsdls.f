c=======================================================================
c
c     SDLS Utility Functions
c
c-----------------------------------------------------------------------
c
c        ADJOINT : computing the symmetric matrix sum of Ai*y(i)
c
c        ADJOINTG: computing the symmetric matrix
c                  sum(Ai*lambda(i)) + sum(Bj*mu(j)) - sum(Bj*gamma(j))
c
c        OPERATA : computing the operator_A : vector(m) result of
c                  scalar product of m Ai constraints matrices with
c                  matrix X
c
c        CSDLSCE : computing the equal constraints matrices for SDLSCE
c
c        CSDLSCI : computing the constraints matrices for SDLSCI
c
c=======================================================================
c
c     subroutine ADJOINT
c
c     Computes symmetric matrix sum of Ai*y(i)
c
c-----------------------------------------------------------------------
      SUBROUTINE adjoint ( mct, nmat, vecty, matai, matadj )
c-----------------------------------------------------------------------
c
c     INPUT 
c            mct    : number of constraints                      integer
c            nmat   : dimension of the matrices X and Ai         integer
c            vecty  : y vector(mct)                               double
c            matai  : mct vectorized symmetric matrices Ai
c                     vector( mct * (nmat*(nmat+1)/2) )           double
c
c     OUTPUT 
c            matadj : sum(Ai.y(i)) vectorized symmetric matrix
c                     vector(nmat*(nmat+1)/2)                     double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o arguments
      INTEGER mct, nmat
      DOUBLE PRECISION vecty(*), matai(*), matadj(*)
c
c     local variables
      INTEGER im, in, pai, ssm
      DOUBLE PRECISION sum
c
c     initialization
      pai  = 0
      ssm  = (nmat*(nmat+1)/2)
c
      DO in=1,ssm
         pai = 0
         sum = 0.
         DO im=1,mct
            sum = sum + matai(pai+in)*vecty(im)
            pai = pai + ssm
         ENDDO
         matadj(in) = sum
      ENDDO
c
      RETURN
      END
c
c=======================================================================
c
c     subroutine ADJOINTG
c
c     Computing symmetric matrix
c     sum(Ai*lambda(i)) + sum(Bj*mu(j)) - sum(Bj*gamma(j))
c
c-----------------------------------------------------------------------
      SUBROUTINE adjointg ( mcte, mcti, nmat, ydual, acte, acti,
     &                      matadj )
c-----------------------------------------------------------------------
c
c     INPUT 
c            mcte   : number of constraints                      integer
c            mcti   : number of inequal constraints              integer
c            nmat   : dimension of the matrices X and Ai         integer
c            ydual  : dual solution        vector(mcte + 2*mcti)  double
c            acte   : mcte Ai symmetric matrices(nmat*nmat)
c                     of constraints
c                     vectorized  vector( mct*(nmat*(nmat+1)/2) ) double
c            acti   : mcti Bj symmetric matrices(nmat*nmat) of
c                     inequal constraints  vector(mcti*nmat*nmat)
c                     vector( mcti*(nmat*(nmat+1)/2) )            double
c
c     OUTPUT 
c            matadj : sum(Ai.y(i)) vectorized symmetric matrix
c                     vector(nmat*(nmat+1)/2)                     double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o arguments
      INTEGER mcte, mcti, nmat
      DOUBLE PRECISION ydual(*), acte(*), acti(*), matadj(*)
c
c     local variables
      INTEGER im, in, pact, ssm, plam, pgam, pmu
      DOUBLE PRECISION sum
c
c-----------------------------------------------------------------------
c
c     pointers for dual vector ydual
      plam = 1
c     plam   : pointer for lambda Lagrange coefficients, so mcte more
      pmu = plam + mcte
c     pmu    : pointer for mu Lagrange coefficients, so mcti more
      pgam = pmu + mcti
c     pgam   : pointer for gam Lagrange coefficients, so mcti more
c
      pact = 0
      ssm  = (nmat*(nmat+1)/2)
c
      DO in = 1,ssm
         sum  = 0.
         pact = 0
         DO im = 1,mcte
            sum = sum + acte(pact+in)*ydual(plam+im-1)
            pact = pact + ssm
         ENDDO
         pact = 0
         DO im = 1,mcti
            sum = sum + acti(pact+in)*(ydual(pmu+im-1)-ydual(pgam+im-1))
            pact = pact + ssm
         ENDDO
         matadj(in) = sum
      ENDDO
c
      RETURN
      END
c
c=======================================================================
c
c     subroutine OPERATA
c
c     Computes operator_A : vector(m) result of scalar product of
c     m Ai constraints matrices with matrix X
c
c-----------------------------------------------------------------------
      SUBROUTINE operata ( mct, nmat, matx, matai, vect )
c-----------------------------------------------------------------------
c
c     INPUT 
c            mct    : number of constraints                      integer
c            nmat   : dimension of the matrices X and Ai         integer
c            matx   : vectorized symmetric matrix X
c                     vector(nmat*(nmat+1)/2)                     double
c            matai  : mct vectorized symmetric matrices Ai
c                     vector( mct * (nmat*(nmat+1)/2) )           double
c
c     OUTPUT 
c            vect   : Ai.X vector(mct)                            double
c
c     CALL   
c            XMS    : computing scalar product of two vectorized
c                     symmetric matrices M1.M2
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o arguments
      INTEGER mct, nmat
      DOUBLE PRECISION matx(*), matai(*), vect(*)
c
c     local variables
      INTEGER im, pai, ssm
c
c     initialization
      pai = 1
      ssm = (nmat*(nmat+1)/2)
      DO im = 1,mct
c
c        computing scalar product between matrix Ai and matrix X
         CALL XMS ( nmat, matx, matai(pai), vect(im) )
         pai = pai + ssm
      ENDDO
c
      RETURN
      END
c
c=======================================================================
c
c     subroutine  CSDLSCE                                    
c
c     Equal constraints matrices for SDLSCE
c
c-----------------------------------------------------------------------
      SUBROUTINE csdlsce ( nmat, cmat, imatce, matce, mcte, acte, bcte )
c-----------------------------------------------------------------------
c
c     INPUT 
c            nmat   : dimension of the matrix                    integer
c            cmat   : symmetric matrix(nmat*nmat) to optimize     double
c            imatce : symmetric matrix(nmat*nmat)                integer
c                     equal constraints indicators
c                     if imatce(i,j) = 1 : xmat(i,j) = cmat(i,j)
c                     if imatce(i,j) = 2 : xmat(i,j) = matce(i,j)
c                     only the lower part is used
c            matce  : symmetric matrix(nmat*nmat)                 double
c                     equal constraint value if imatce(i,j) = 2
c                     only the lower part is used
c     OUTPUT 
c            mcte   : number of equal constraints                integer
c            acte   : mcte Ai symmetric matrices(nmat*nmat) of
c                     equal constraints, vector(mcte*nmat*nmat)   double
c            bcte   : value of equal constraints, vector(mcte)    double
c
c     CALL   
c        IMX     : initialization at a scalar of a matrix
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o arguments
      INTEGER nmat, mcte
      INTEGER imatce(nmat,*)
      DOUBLE PRECISION cmat(nmat,*), matce(nmat,*), acte(*), bcte(*)
c
c     local variables
      INTEGER i, j, ksa, sm, kb, ice
      DOUBLE PRECISION dzero
      PARAMETER ( dzero=0.E0 )
c
c-----------------------------------------------------------------------
c
      sm  = nmat*nmat
      ksa = 1
      kb  = 0
      DO i = 1,nmat
         DO j = 1,i
            ice = imatce(i,j)
            IF ( ice .NE. 0 ) THEN
               kb = kb + 1
               IF (ice .EQ. 1) bcte(kb) = cmat(i,j)
               IF (ice .EQ. 2) bcte(kb) = matce(i,j)
               CALL IMX ( nmat, nmat, acte(ksa), dzero )
               acte(ksa+i+(j-1)*nmat-1) = 1
               ksa = ksa + sm
            ENDIF
         ENDDO
      ENDDO
c
      mcte = kb
c
      RETURN
      END
c
c=======================================================================
c
c     subroutine  CSDLSCI
c
c     Computing the constraints matrices for SDLSCI
c
c-----------------------------------------------------------------------
      SUBROUTINE csdlsci ( nmat, cmat, imatc, matc, matlci,
     &                     mcte, acte, bcte, mcti, acti, lcti, ucti )
c-----------------------------------------------------------------------
c
c     INPUT 
c            nmat   : dimension of the matrix                    integer
c            cmat   : symmetric matrix(nmat*nmat) to optimize     double
c            imatc  : symmetric matrix(nmat*nmat)                integer
c                     constraints indicators (see method)
c                     only the lower part is used
c            matc   : symmetric matrix(nmat*nmat)                 double
c                     equal constraint value if imatc(i,j) = 2
c                     inequal constraint tolerance if imatc(i,j) = 3
c                     inequal constraint up value if imatc(i,j) = 4
c                     only the lower part is used
c            matlci : symmetric matrix(nmat*nmat)                 double
c                     inequal constraint low value if imatc(i,j) = 4
c                     only the lower part is used
c                     only the lower part is used
c     OUTPUT 
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
c     CALL   
c        IMX     : initialization at a scalar of a matrix
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o arguments
      INTEGER nmat, mcte, mcti
      INTEGER imatc(nmat,*)
      DOUBLE PRECISION cmat(nmat,*), matc(nmat,*), matlci(nmat,*), 
     &                 acte(*), bcte(*), acti(*), lcti(*), ucti(*)
c
c     local variables
      INTEGER i, j, sm, ik, icte, icti, ksae, ksai, icc
      DOUBLE PRECISION dzero
      PARAMETER ( dzero = 0.E0 )
c
c-----------------------------------------------------------------------
c
      sm = nmat*nmat
      icte = 0
      icti = 0
c
      ksae = 1
      ksai = 1
      DO i = 1,nmat
         DO j = 1,i
            ik  = i+(j-1)*nmat
            icc = imatc(i,j)
c
            IF ( icc .EQ. 1 ) THEN
               icte = icte + 1
               bcte(icte) = cmat(i,j)
               CALL IMX ( nmat, nmat, acte(ksae), dzero )
               acte(ksae+ik-1) = 1
               ksae = ksae + sm
            ENDIF
c
            IF ( icc .EQ. 2 ) THEN
               icte = icte + 1
               bcte(icte) = matc(i,j)
               CALL IMX ( nmat, nmat, acte(ksae), dzero )
               acte(ksae+ik-1) = 1
               ksae = ksae + sm
            ENDIF
c
            IF ( icc .EQ. 3 ) THEN
               icti = icti + 1
               ucti(icti) = cmat(i,j) + matc(i,j)
               lcti(icti) = cmat(i,j) - matc(i,j)
               CALL IMX ( nmat, nmat, acti(ksai), dzero )
               acti(ksai+ik-1) = 1
               ksai = ksai + sm
            ENDIF
c
            IF ( icc .EQ. 4 ) THEN
               icti = icti + 1
               ucti(icti) = matc(i,j)
               lcti(icti) = matlci(i,j)
               CALL IMX ( nmat, nmat, acti(ksai), dzero )
               acti(ksai+ik-1) = 1
               ksai = ksai + sm
            ENDIF
c
         ENDDO
      ENDDO
c
      mcte = icte
      mcti = icti
c
      RETURN
      END
