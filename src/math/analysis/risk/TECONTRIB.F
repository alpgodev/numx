c=======================================================================
c
c     subroutine TECONTRIB                                  
c
c     Tracking Error Contribution (marginal risk contribution)
c
c     TE(i) = w(i)*[w'*Cov] / (Total Tracking Error)
c
c-----------------------------------------------------------------------
      SUBROUTINE TECONTRIB( n, cov, covb, varb, w, iwork, dwork,
     &                      mte, mpte, info)
c-----------------------------------------------------------------------
c
c     INPUT 
c       n      : number of asset(s)                              integer
c       w      : weight(s) (n)                                    double 
c       cov    : covariance matrix (n*n)                          double
c       covb   : covariance of benchmark (n)                      double
c       varb   : variance of benchmark                            double
c
c     WORKSPACE 
c       iwork  : 0                                              integer
c       dwork  : n*(6*n + 32) + 3                                double
c
c     OUTPUT 
c       mte   : marginal Tracking Error contribution (n)         double
c       mpte  : marginal Tracking Error contribution in % (n)    double
c       info  : diagnostic argument                             integer
c
c     CALL   
c       EXATER : ex-ante Tracking Error
c       PMV    : Matrix * Vector = Vector
c       PVX    : Vector * Scalar = Vector
c
c-----------------------------------------------------------------------
c      
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info
      DOUBLE PRECISION varb
      DOUBLE PRECISION w(*), cov(*), covb(*), mte(*), mpte(*)
c
c     workspaces      
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER i, piwork, pdwork, piwt, piwk, pigk, pdwt, pdwk,
     &        pdro, pdsp, pdss, pdsa  
      DOUBLE PRECISION  EPS, te, invte, tmp, sum
      PARAMETER (EPS = 1.E-30)
c
c     external subroutines
      EXTERNAL EXATER, IVX, PMV, PVX
c
c-----------------------------------------------------------------------
c
c     initialisation
      info = 0
      CALL IVX ( n, mte, 0.0 )
      CALL IVX ( n, mpte, 0.0 )
c      
c     pointers for integer work space : iwork
c     ---------------------------------------
      piwork = 1
      piwt   = piwork 
c     piwt  : pointer for EXATER workspace ( 13*n )      
      piwk   = piwork
c     piwk  : pointer for KATAVEP workspace ( 12*n )
      pigk   = piwk + ( 12*n )
c     pigk  : pointer Kato groups, so (p) more
c
c     Total size of iwork array = ( 13*n )
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdwork = 1  
      pdwt   = pdwork 
c     pdwt  : pointer for EXATER workspace ( n*(6*n + 32) + 3 ) 
      pdwk   = pdwork
c     pdwk  : pointer for KATAVEP workspace ( n*(4*n + 26) )       
      pdsp   = pdwk + ( n*(4*n + 26) )
c     pdsp   : pointer for vector spectr of KATAVEP, so ( n ) more
      pdss   = pdsp + ( n )
c     pdss   : pointer for vector eigsor of KATAVEP, so ( n ) more
      pdsa   = pdss + ( n )
c     pdsa   : pointer for vector spectr of KATAVEP, so ( n ) more
      pdro   = pdwk + ( n )
c     pdro   : pointer for the robust quadratic part, so ( n*n ) more  
c
c     Total size of dwork array = n*(6*n + 32) + 3
c
c--------------------------------------------------------------------        
c
c     absolute portfolio Tracking Error EXATER (n, cov, covben, varben, w, dwork, tracke, info)
      CALL EXATER ( n, cov, covb, varb, w, dwork(pdwt), te, info)
      IF (info .LT. 0) RETURN
c
c     if TE is too smal then contribution = w(i)
      IF (ABS(te) .LT. EPS) THEN
        DO i = 1,n
            mpte(i) = w(i)
        ENDDO
        info = -1
        RETURN
      ENDIF  
c
c     compute w'*cov      
      CALL PMV ( n, n, cov, w, mte )
c
c     compute (w'*cov)(i) - covb(i)
      DO i = 1,n
        mte(i) = mte(i) - covb(i)
      ENDDO
c
c     compute ((w'cov)(i) - covb(i))/te
      invte = 1./te
      CALL PVX ( n, mte, invte, mte )
c
c     compute (varb - w'*covb)/te
      tmp = 0.0
      DO i = 1,n
        tmp = tmp - w(i)*covb(i)
      ENDDO  
      tmp = invte*(tmp + varb)
c
c     marginal contribution w(i)*[((w'cov)(i) - covb(i))/te + (varb - w'*covb)/te]
      DO i = 1,n
        mte(i) = w(i)*(mte(i) + tmp)
      ENDDO  
      sum = 0.
      DO i = 1, n 
        sum = sum + mte(i)  
      ENDDO
      sum = 1./sum      
c
c     marginal contribution in %, 100*TE*MarginalTE(i)   
      CALL PVX ( n, mte, invte, mpte )
c
c     cas TE = 0%
      
c
c      open(unit=1,file='TE.txt',status='unknown')
c      write(1,*) te
c      write(1,*) '-----------'
c      write(1,*) tmp
c      write(1,*) '-----------'
c      do i = 1,n
c		 write(1,*)  mte(i)
c	end do   
c	close(unit=1)    
c
      RETURN
      END
