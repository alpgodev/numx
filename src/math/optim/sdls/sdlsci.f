c=======================================================================
c
c     subroutine  SDLSCI                                     
c
c     Semi-Definite Least Square optimization
c     with equal and inequal unitary constraints
c
c-----------------------------------------------------------------------
      SUBROUTINE sdlsci ( nmat, cmat, imatc, matc, matlci,
     &                    epsbfg, alpha, iwork, dwork,
     &                    xmat, info )
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
c            epsbfg : BFGS precision stop test                    double
c            alpha  : alpha  ( 0.0 =<alpha )                      double
c
c     WORKSPACE 
c            iwork  : vector(2*mcte+4*mcti+12*nmat+5)            integer
c            dwork  : vector                                      double
c                     (  (mcti+mcte)*(nmat*nmat)
c                         + (mcti) + (mcti) + mcti*(nmat*nmat)
c                         + (mcte+2*mcti)*((mcte+2*mcti+27)/2)
c                         + (2*mcte+2*mcti+5)*(nmat*(nmat+1)/2)
c                         + (nmat*(4*nmat+27))  )
c
c     OUTPUT 
c            xmat   : symmetric matrix solution (nmat*nmat)       double
c            info   : = 0 successful exit                        integer
c
c     CALL   
c        CSDLSCI : computing the constraints matrices for SDLSCI
c        SDLSG   : semi-Definite Least Square optimization
c                  general version with equal and inequal constraints
c                  with M >= alpha >= 0.0
c                  ( the constraints matrices Ai
c                    and Bi may be non-symmetric )
c
c     METHOD 
c            if imatc(i,j) = 1, xmat(i,j) = cmat(i,j)
c            if imatc(i,j) = 2, xmat(i,j) = matc(i,j)
c            if imatc(i,j) = 3,
c               cmat(i,j)-matc(i,j) < xmat(i,j) < cmat(i,j)+matc(i,j)
c            if imatc(i,j) = 4, matlci(i,j) < xmat(i,j) < matc(i,j)
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o arguments
      INTEGER nmat, info, imatc(nmat,*)
      DOUBLE PRECISION epsbfg, alpha, cmat(*), matc(*),matlci(*),xmat(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER mcte, mcti, i, j, pdae, pdbe, pdai, pdli, pdui, piw, pdw,
     &        n, mcteq
c
c-----------------------------------------------------------------------
c
c     initializations
      info = 0
c
c     computing the number of equal and inequal constraints
      mcteq = 0
      mcte = 0
      mcti = 0
      DO i = 1,nmat
         DO j = 1,i
            IF ( imatc(i,j) .EQ. 1 ) THEN 
                mcte = mcte + 1
                mcteq = mcteq + 1
            ENDIF    
            IF ( imatc(i,j) .EQ. 2 ) mcte = mcte + 1
            IF ( imatc(i,j) .EQ. 3 ) mcti = mcti + 1
            IF ( imatc(i,j) .EQ. 4 ) mcti = mcti + 1
         ENDDO
      ENDDO
c
c     pointers for integer work space : iwork
c     ---------------------------------------
      piw = 1
c     piw   : pointer for local workspaces,
c             SDLSG needs (2*mcte+4*mcti+12*nmat+5)
c
c     Total size of iwork array : (2*mcte+4*mcti+12*nmat+5)
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdbe  = 1
c     pdbe  : pointer for vector b(mcte)
      pdae  = pdbe + (mcte)
c     pdae  : pointer for the mcte Ai matrix(nmat*nmat),
c             so mcte*(nmat*nmat) more
      pdai  = pdae + mcte*(nmat*nmat)
c     pdai  : pointer for the mcti Bi matrix(nmat*nmat),
c             so mcti*(nmat*nmat) more
      pdli  = pdai + mcti*(nmat*nmat)
c     pdli  : pointer for the low vector(mcti)
      pdui  = pdli + (mcti)
c     pdui  : pointer for the up vector(mcti)
      pdw   = pdui + (mcti)
c     pdw   : pointer for local workspaces,
c             SDLSG needs (  (mcte+2*mcti)*((mcte+2*mcti+25)/2)
c                          + (2*mcte+2*mcti+5)*(nmat*(nmat+1)/2)
c                          + (nmat*(4*nmat+27))  )
c
c     Total size of dwork array : (mcte)  + mcte*(nmat*nmat)
c                               + (mcti) + (mcti) + mcti*(nmat*nmat)
c                               + (mcte+2*mcti)*((mcte+2*mcti+25)/2)
c                               + (2*mcte+2*mcti+5)*(nmat*(nmat+1)/2)
c                               + (nmat*(4*nmat+27))  )
c
c                         = (mcte+2*mcti) + (mcti+mcte)*(nmat*nmat)
c                         + (mcti) + (mcti) + mcti*(nmat*nmat)
c                         + (mcte+2*mcti)*((mcte+2*mcti+25)/2)
c                         + (2*mcte+2*mcti+5)*(nmat*(nmat+1)/2)
c                         + (nmat*(4*nmat+27))  )
c
c                         = (mcti+mcte)*(nmat*nmat)
c                         + (mcti) + (mcti) + mcti*(nmat*nmat)
c                         + (mcte+2*mcti)*((mcte+2*mcti+27)/2)
c                         + (2*mcte+2*mcti+5)*(nmat*(nmat+1)/2)
c                         + (nmat*(4*nmat+27))  )
c
c-----------------------------------------------------------------------
c
c     case with eq. constraints mcteq = n*(n+1)/2
      n = nmat*(nmat + 1)/2
      IF (mcteq .EQ. n) THEN
        CALL  YM ( nmat, nmat, cmat, xmat )
        RETURN
      ENDIF
c
c     Ai and Bi constraints matrices (cf. utsdls.f)
      CALL csdlsci ( nmat, cmat, imatc, matc, matlci,
     &               mcte,  dwork(pdae), dwork(pdbe),
     &               mcti, dwork(pdai), dwork(pdli), dwork(pdui) )
c
c     SDLS optimization
      CALL sdlsg ( nmat, cmat, mcte,  dwork(pdae), dwork(pdbe),
     &             mcti, dwork(pdai), dwork(pdli), dwork(pdui),
     &             epsbfg, alpha, iwork(piw), dwork(pdw),
     &             xmat, info )
c
      RETURN
      END
