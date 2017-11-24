c=======================================================================
c
c     subroutine ALLOCMDD                                    
c
c     Asset allocation strategy s.t. Maximum Drawdown constraint
c
c-----------------------------------------------------------------------
      SUBROUTINE allocmdd ( ndate, nasset, radata, cov, rmean, neq, nin,
     &                      ccst, bcst, cinf, csup,
     &                      delta, lambda,
     &                      iwork, dwork, wopt, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       ndate  : number of dates                            integer
c       nasset : number of risky assets                     integer
c       radata : returns assets  matrix(ndate*nasset)        double
c       cov    : covariance  matrix(nasset*nasset)           double
c       rmean  : mean returns  vector(nasset)                double
c       neq    : number of equality constraints             integer
c       nin    : number of inequality constraints           integer
c       ccst   : matrix of constraints (nasset*(neq+nin))    double
c       bcst   : vector of constraints  (neq+nin)            double
c       cinf   : lower bound (nasset)                        double
c       csup   : upper bound (nasset)                        double
c       delta  : maximum drawdown                            double
c       lambda : risk aversion                               double
c
c     WORKSPACE 
c       iwork  : 13*nasset+2*ndate+2*nin+neq+1    integer
c       dwork  :                                       double
c                    (   nasset*(6*nasset+ndate+nin+neq+29)
c                        + 4*(ndate+nin) + 2*neq             )
c
c     OUTPUT 
c        wopt   : optimal portfolio (nasset)                 double
c        info   : diagnostic argument                        integer
c
c     CALL   
c        YM      : copy a vectorized matrix in a vectorized matrix
c        PMX     : computing M*X = matrix
c                  ( M matrix(n*m), X scalar, gives M*X matrix(n*m) )
c        YV      : copy a vector in a vector
c        IVX     : initialization at a scalar of a vector
c        QP      : quadratic solver (vectorized version)
c        TESTSDP : test if covariance matrix is SDP
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER ndate, nasset, neq, nin, info
      DOUBLE PRECISION delta, lambda
      DOUBLE PRECISION radata(ndate,*), cov(*), rmean(*), 
     &                 ccst(nasset,*), bcst(*), cinf(*), csup(*), 
     &                 wopt(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER ncst, neqtot, nintot, ncstot, i, j, k, icc, 
     &        piwork, pdwork, piwk, pgk, psp, pss, psa, pdwk, piwq, 
     &        plin, pdwq, plagr, pbcs, pccs, pro, pcov, pisdp, pdsdp
      DOUBLE PRECISION sum
c
c     external subroutines
      EXTERNAL qp, YM, PMX, YV, IVX, testsdp
c
c-----------------------------------------------------------------------
c
c     initialisations
      info = 0 
      ncst = neq + nin
      neqtot = neq
      nintot = nin + ndate
      ncstot = neqtot + nintot
c
c     pointers for integer work space : iwork
c     ---------------------------------------
      piwork = 1
      piwk = piwork
c     piwk  : pointer for KATAVEP who needs (12*nasset)
      pgk = piwk + (12*nasset)
c     pgk   : pointer Kato groups, so (nasset) more
c     the call of KATAVEP needs (13*nasset)
      pisdp  = piwk + (12 * nasset)
c     pisdp  : pointer for TESTSDP (nasset)  
c
c     pointers for QUAPRO who uses a part of the same space
      piwq = piwork
c     piwq  : pointer for QUAPRO who needs
c             (  3*nasset + 2*nintot + neqtot + 1 ),
c             so (3*nasset + 2*(nin+ndate) + neq + 1 ) more
c     pinext = piwq + ( 3*nasset + 2*ndate + 2*nin + neq + 1 )
c     pinext : pointer for the next iwork array
c
c     together KATAVEP/TESTSDP and QUAPRO need (union of space)
c     ( 14*nasset + 2*nin + neq + 7 )
c
c     Total size of iwork array = ( 14*nasset + 2*ndate + 2*nin + neq + 1 )
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pro = 1
c     pro  : pointer for robust cov. matrix, so (nasset*nasset) more
      pcov = pro + (nasset*nasset)
c     pcov : pointer for covariance matrix, so (nasset*nasset) more
      pdwork = pcov + (nasset*nasset)
c
      psp = pdwork
c     psp  : pointer for vector spectr of KATAVEP, so (nasset) more
      pss = psp + (nasset)
c     pss  : pointer for vector eigsor of KATAVEP, so (nasset) more
      psa = pss + (nasset)
c     psa  : pointer for vector spectr of KATAVEP, so (nasset) more
      pdwk = psa + (nasset)
c     pdwk : pointer for KATAVEP workspace, needs nasset*(4*nasset+26)
c     so KATAVEP needs nasset*(4*nasset+29)
      pdsdp = pdwk + nasset*(4*nasset + 26)
c     pdsdp : pointer for TESTSDP (nasset*(2*nasset + 7))  
c     
c     pointers for QUAPRO who uses a part of the same space
c     than KATAVEP
      plin = pdwork
c     plin  : pointer for linear part vector, so ( nasset ) more
      pdwq = plin + ( nasset )
c     pdwq  : pointer for QUAPRO workspace who needs
c             (nasset*nasset + 6*nasset + 2*nintot),
c             so ( nasset*(nasset+6) + 2*nintot ) more
      plagr = pdwq + ( nasset*(nasset+6) + 2*nintot )
c     plagr  : pointer for Lagrange multipliers vector
c              ( nasset + nintot + neqtot ),
c              so ( nasset + nintot + neqtot ) more
      pbcs = plagr +  ( nasset + nintot + neqtot )
c     pbcs : pointer for the constraints vector ( nintot + neqtot ),
c              so ( nintot + neqtot ) more
      pccs = pbcs +  ( nintot + neqtot )
c     pccs : pointer for the constraints matrix
c              nasset*( nintot + neqtot ),
c              so nasset*( nintot + neqtot ) more
c     so QUAPRO needs :  nasset
c                      + nasset*(nasset+6) + 2*(nintot)
c                      + nasset + nintot + neqtot
c                      + nintot + neqtot
c                      + nasset*( nintot + neqtot ) 
c
c                      = nasset*(1+(nasset+6)+1 + nintot+neqtot)
c                      + 2*nintot + nintot + neqtot + nintot + neqtot
c
c                      = nasset*( 2*nasset + nintot + neqtot + 8 )
c                      + 4*nintot + 2*neqtot
c
c                      = nasset*( 2*nasset + ndate + nin + neq + 8 )
c                      + 4*(ndate+nin) + 2*neq
c
c     together KATAVEP/TESTSDP and QUAPRO need (union of space)
c                        nasset*( 6*nasset + ndate + nin + neq + 36 )
c                      + 4*(ndate+nin) + 2*neq
c
c     Total size of dwork array
c                      = nasset*nasset
c                      + nasset*nasset
c                      + nasset*( 6*nasset + ndate + nin + neq + 36 )
c                      + 4*(ndate+nin) + 2*neq
c
c                      = nasset*( 8*nasset + ndate + nin + neq + 36 )
c                      + 4*(ndate+nin) + 2*neq
c
c-----------------------------------------------------------------------
c
      CALL YM ( nasset, nasset, cov, dwork(pro) )
c
c     covariance matrix SDP ?
      CALL testsdp (nasset,dwork(pro),iwork(pisdp),dwork(pdsdp),info)
      IF (info .NE. 0) THEN
         info = -108
         RETURN
      ENDIF        
c
c     covariance matrix * lambda
      CALL PMX ( nasset, nasset, dwork(pro), lambda, dwork(pcov) )
c
c     linear part
      DO i=1,nasset
         dwork(plin+i-1) = - rmean(i)
      ENDDO
c
c     copy vector of constraints and computes new constraints
      CALL YV ( ncst, bcst, dwork(pbcs) )
      CALL IVX ( ndate, dwork(pbcs+ncst), delta )
c
c     copy matrix of constraints 
      CALL YV ( nasset*ncst, ccst, dwork(pccs) )
c     
c     drawdown constraints
      icc = pccs + nasset*ncst
      DO i=1,nasset
         DO j=1,ndate
            sum = 0.0
            DO k=1,ndate
               sum = sum + ( radata(j,i) - radata(k,i) )
            ENDDO
            dwork(icc) = sum / ndate
            icc = icc + 1
         ENDDO
      ENDDO        
c
c     quadratic solver
      CALL qp ( nasset, dwork(pcov), dwork(plin), neqtot, nintot,
     &          dwork(pccs), dwork(pbcs), cinf, csup,
     &          iwork(piwq), dwork(pdwq), dwork(plagr), wopt,
     &          info )
c
      RETURN
      END
