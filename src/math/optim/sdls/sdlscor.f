c=======================================================================
c
c     subroutine  SDLSCOR
c
c     Semi-Definite Least Square optimization
c     for a correlation matrix with M >= alpha     ( 0.0 =<alpha < 1.0 )  
c
c-----------------------------------------------------------------------
      SUBROUTINE sdlscor ( nmat, cmat, epsbfg, alpha, iwork, dwork,
     &                     xmat, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c            nmat   : dimension of the matrix                    integer
c            cmat   : symmetric matrix(nmat*nmat) to optimize     double
c            epsbfg : BFGS precision stop test                    double
c            alpha  : alpha  ( 0.0 =<alpha < 1.0 )                double
c
c     WORKSPACE 
c            iwork  : ( 14*nmat + 2 )                            integer
c            dwork  : ( nmat*(7*nmat + 40) )                      double       
c
c     OUTPUT 
c            xmat   : matrix solution (nmat*nmat)                 double
c            info   : diagnostic argument                        integer
c
c     CALL   
c        SEMDX   : computing the sum of the diagonal elements of a
c                  matrix with a scalar : B = A + x*I
c                  ( A square matrix(n*n), x scalar,
c                    gives B square matrix(n*n) )
c        PMX     : computing M*X = matrix
c                  ( M matrix(n*m), X scalar, gives M*X matrix(n*m) )
c        YM      : copy a vectorized matrix in a vectorized matrix
c        CMCMS    : converting a vector of full square matrix
c                   in vector of symmetric matrix
c        BFGSSDLS : quasi-Newton optimizer with BFGS method
c                   ( short call version for SDLS )
c        SIMSCOR  : subroutine used by BFGS to compute the value of the
c                   dual function and its gradient
c                   specific to obtain a correlation matrix
c        SOLSCOR  : computing the solution of SDLS optimization from the
c                   dual solution   x = pK( C + Ai.y(i) ) 
c                   specific to obtain a correlation matrix
c        CMSMC    : converting a vector of symmetric matrix
c                   in vector of full square matrix
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
      EXTERNAL simscor, gesterr
c
c     i/o arguments
      INTEGER nmat, info
      DOUBLE PRECISION epsbfg, alpha, cmat(*), xmat(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER i, ssm, mct, totitr, totsim, piw, pinfo, pdw, pdcv, pdxs,
     &        pdyd, pdbinf, pdbsup, pdgrad, pdcs
      DOUBLE PRECISION eps, infini, opalph, invopa, funct
      PARAMETER ( eps = 1.e-15, infini = 1.d20 )
c
c-----------------------------------------------------------------------
c
c     initializations
      info = 0
c
c     test if alpha in [0,1]      
      IF ( (alpha .LT. 0.).or.(alpha .GT. (1. - eps)) ) THEN
         info = -2001
         RETURN
      ENDIF
c
      ssm    = nmat*(nmat+1)/2
      mct    = nmat
      opalph = (1.-alpha)
      invopa = 1./opalph
c
c     pointers for integer work space  : iwork
c     ----------------------------------------
      piw   = 1
c     piw   : pointer for local workspaces,
c             BFGSSDLS needs (2*mct+1)
c                      + communication with SIMSCOR (12*nmat+1),
c                      so ((2*mct+1)+(12*nmat+1))
c             SOLSCOR needs (12*nmat)
c             so, union of spaces : ((2*mct+1)+(12*nmat+1))
c             for correlation mct = nmat,
c             so : (14*nmat+2)
c
c     Total size of iwork array = (14*nmat+2)
c
      pinfo  = piw + (2*mct+1)
c     pinfo  : communication pointer with SIMCOR for info
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdcv   = 1
c     pdcv   : pointer for changed variable matrix C(nmat*nmat)
      pdxs   = pdcv + (nmat*nmat)
c     pdxs   : pointer for the vectorized symmetric matrix X,
c              so ssm = (nmat*(nmat+1)/2) more
      pdyd   = pdxs + ssm
c     pdyd   : pointer for initial value of ydual, so mct more
      pdbinf = pdyd + mct
c     pdbinf : pointer for inferior bounds binf, so mct more
      pdbsup = pdbinf + mct
c     pdbsup : pointer for superior bounds bsup, so mct more
      pdgrad = pdbsup + mct
c     pdgrad : pointer for gradient at solution, so mct more
      pdw    = pdgrad + mct
c     pdw    : pointer for local workspaces,
c             BFGSSDLS needs ( mct*(mct+11)/2 )
c                      + communication with SIMSCOR ( mct*(5*mct+30) ),
c                      for correlation mct = nmat
c                      so ( nmat*(nmat+11)/2 )+ ( nmat*(5*nmat+30) )
c             SOLSCOR needs (nmat*(nmat+1)/2 + nmat*(3*nmat+27))
c             so, union of spaces: (nmat*(nmat+11)/2)+(nmat*(5*nmat+30 )
c
c     Total size of dwork array = (nmat*nmat) + ssm + 4*mct
c                           + ( mct*(mct+11)/2 ) + ( mct*(5*mct+30) )
c     for correlation mct = nmat
c        = (nmat*nmat) + (nmat*(nmat+1)/2) + 4*nmat
c          + ( nmat*(mnat+11)/2 ) + ( nmat*(5*nmat+30) )
c        = nmat*( nmat + ((nmat+1)/2)+((nmat+11)/2)+(5*nmat+30)+4 )
c        = nmat * ( nmat + (2*nmat+12)/2) + (5*nmat+30) + 4 )
c        = nmat * ( nmat + (nmat+6) + (5*nmat+30) + 4 )
c
c     = nmat * (7*nmat+40)
c
      pdcs  = pdw +  mct*(mct+11)/2
c     pdcs  :  communication pointer with SIMCOR
c              for the vectorized symmetric matrix C,
c
c-----------------------------------------------------------------------
c
c     chgt. variable CA = (C - alpha*I) / (1-alpha)), if alpha > epsilon
      IF (alpha .GT. eps) THEN
         CALL SEMDX ( nmat, cmat, -alpha, dwork(pdcv) )
         CALL PMX ( nmat, nmat, dwork(pdcv), invopa, dwork(pdcv) )
      ELSE
         CALL YM ( nmat, nmat, cmat, dwork(pdcv) )
      ENDIF
c
c     initialization of ydual
      DO i = 1,mct
         dwork(pdyd+i-1) = 0.
      ENDDO
c
c     initializations of bounds
      DO i = 1,mct
         dwork(pdbinf+i-1) = -infini
         dwork(pdbsup+i-1) = infini
      ENDDO
c
c     construction of iwork vector (communication with simscor)
      iwork(pinfo) = 0
c
c     conversion of the full matrix C in vectorized symmetric matrix
c     in dwork vector (communication with simscor)
      CALL CMCMS ( nmat, dwork(pdcv), dwork(pdcs) )
c
c     non-linear optimization (BFGS)
c
c     (the second parameter is not used here, we can replace GESTERR
c      by an other name of subroutine)
      CALL bfgssdls ( simscor, gesterr, mct, dwork(pdyd), epsbfg,
     &                dwork(pdbinf), dwork(pdbsup),
     &                iwork(piw), dwork(pdw),
     &                funct, dwork(pdgrad), totitr, totsim,
     &                info)
c
c     gestion of SIMSCOR errors
      info = iwork(pinfo)
c
c     construction of the matrix solution from dual solution
      CALL solscor ( nmat, dwork(pdcs), dwork(pdyd),
     &               iwork(piw), dwork(pdw),
     &               dwork(pdxs), info )
c
c     conversion of the vectorized symmetric matrix X in full matrix
      CALL CMSMC ( nmat, dwork(pdxs), dwork(pdcv) )
c
c     changing variable X = Z + alpha*I if alpha > epsilon
      IF (alpha .GT. eps) THEN
         CALL PMX ( nmat, nmat, dwork(pdcv), opalph, dwork(pdcv) )
         CALL SEMDX ( nmat, dwork(pdcv), alpha, xmat )
      ELSE
         CALL YM ( nmat, nmat, dwork(pdcv), xmat )
      ENDIF
c
c     restore corr(i,i) = 1.0
      DO i = 1,nmat
        xmat(i + (i-1)*nmat) = 1.0
      ENDDO
c
      RETURN
      END
