c=======================================================================
c
c     subroutine KATOBLOCK
c
c     This function implements matrix correction by matrix sub-block.
c     Input - Covariance Matrix and block description.  
c
c------------------------------------------------------------------------
      subroutine katoblock ( n, cov, nclass, class, eps,
     &                       iwork, dwork, covcal, info )
c------------------------------------------------------------------------
c
c     INPUT 
c            n      : matrix size                                integer
c            cov    : covariance matrix matrix (n by n)           double
c            nclass : number of block                            integer
c            class  : class definition (n)                       integer
c            eps    : sensibility parameter                       double
c
c     WORKSPACE 
c            iwork  : ( 14*n + 3 )                               integer 
c            dwork  : ( 5*(n*(n + 1)/2) + n*(8*n + 30))           double
c     OUTPUT 
c            covcal : calibrated matrix (n*n)                     double
c            info   : = 0 successful exit                        integer
c
c     CALL   
c            YM      : copy a vectorized matrix in a vectorized matrix
c            YMCPI   : copy a part of a square matrix
c                      in a smaller square matrix with index
c            TM      : computing the trace of a vectorized full matrix
c            KATAVEP : Calibrates a symmetric matrix
c                      with positive eigenvalues ( > eps kato )
c                      and average Kato eigenvalues
c            YMCPIR  : copy a square matrix
c                      in a part of a bigger square matrix with index
c            EVMIN   : computing the minimum of a vector
c            EMDMIN  : computing the minimum of the diagonal of a matrix
c            SDLS    : Semi-Definite Least Square optimization
c                      with M >= alpha >= 0.0
c                    ( the constraints matrices Ai may be non-symmetric )
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, nclass, info
      INTEGER class(*)
      DOUBLE PRECISION eps
      DOUBLE PRECISION cov(n,*), covcal(n,*)
c 
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER i, j, npk, nbkgrp, piw, pig, pind, pdw, pdv1, pdv2, pdv3,
     &        pdm1, pdm2, pdcov, mct
      DOUBLE PRECISION epskato, myzero, infini, som, trace, alpha, 
     &                 mineig, epsbfg
      DOUBLE PRECISION aict(1), bct(1)
      PARAMETER ( myzero = 1.0d-30, infini = 1.0d50 )
c
c     external subroutines
      EXTERNAL YM, YMCPI, TM, katavep, YMCPIR, EVMIN, EMDMIN, sdls
c     
c     intrinsic functions
      INTRINSIC MIN      
c
c-----------------------------------------------------------------------
c
c     pointers for integer work space : iwork
c     ---------------------------------------
      pig = 1
c     pig   : pointer for vector ( n )
      pind = pig + n
c     pind  : pointer for vector ( n )
      piw = pind + n
c     piw   : pointer for work space used by KATAVEP and SDLS
c        KATAVEP  needs : 12*n 
c        SDLS needs     : 12*n + 3 
c        so common workspace is ( 12*n + 3 )
c
c     Total size of iwork array = ( 14*n + 3 ) 
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdv1 = 1
c     pdv1  : pointer for vector( n )
      pdv2 = pdv1 + n
c     pdv2  : pointer for vector( n )
      pdv3 = pdv2 + n
c     pdv3  : pointer for vector( n )
      pdm1 = pdv3 + n
c     pdm1  : pointer for matrix( n*n )
      pdm2 = pdm1 + ( n*n )
c     pdm2  : pointer for matrix( n*n )
      pdcov = pdm2 + ( n*n )
c     pdmcov: pointer for covariance matrix( n*n )
      pdw = pdcov + ( n*n )
c     pdw   : pointer for work space used by KATAVEP and SDLS
c      KATAVEP needs : n*(4*n + 26)
c      SDLS needs    : 5*(n*(n + 1)/2) + (n*(5*n + 27))
c      so common workspace is ( 5*(n*(n + 1)/2) + (n*(5*n + 27)))
c
c     Total size of dwork array
c             = ( 5*(n*(n + 1)/2) + n*(8*n + 30)
c
c------------------------------------------------------------------------
c
c     initialization
      info = 0
      alpha = infini
      som   = 0
c
c     copy: cov -> dwork(pdcov)
      CALL YM ( n, n, cov, dwork(pdcov) )
c
c     if nb class is 1 exit (no computation)
      IF (nclass .EQ. 1) RETURN
c
c     calibrate the covariance matrix by block
      DO i = 1,nclass
c
c        dimension of block i 
         npk = 0
         DO j = 1,n
            IF (class(j).eq.i) THEN
               npk = npk + 1
               iwork(pind + npk - 1) = j
            ENDIF
         ENDDO
         som = som + npk
         IF (npk.gt.1) THEN
c         
c           extract block i
            CALL YMCPI ( n, cov, npk, iwork(pind),
     &                   dwork(pdm1), info )
c           YMCPI : info = 0 allways by construction
c
c           Trace(sub-block i)
            CALL TM ( npk, dwork(pdm1), trace )
c
c           Kato sensibility: epskato = Trace(block)/(npk*eps) 
            IF ( eps .GT. myzero ) THEN 
                epskato = trace/( npk*eps )
            ELSE
                epskato = trace/( npk*myzero )   
            ENDIF
c
c           Kato correction of the sub-block
            CALL  katavep ( npk, dwork(pdm1), epskato, iwork(piw),
     &                      dwork(pdw),
     &                      dwork(pdm2), dwork(pdv1),
     &                      dwork(pdv2), dwork(pdv3),
     &                      nbkgrp, iwork(pig), info)
            IF (info .LT. 0) RETURN
c
c           SDLS parameter: alpha = min(block eigenvalues)
            CALL EVMIN ( npk, dwork(pdv3), mineig )
            alpha = MIN( alpha, mineig )
c
c           sub-block -> cov. matrix
            CALL YMCPIR (npk, dwork(pdm2), n, iwork(pind),
     &                   dwork(pdcov), info)
c           YMCPIR : info = 0 allways by construction
         ENDIF
      ENDDO
c
      IF (som .NE. n) THEN
         info = -2007
         RETURN
      ENDIF
c
c     precision of SDLS constraints: epsbfg
      CALL EMDMIN ( n, dwork(pdcov), epsbfg )
      epsbfg = epsbfg / n
c
c     SDLS solver 
      mct = 0
      CALL sdls( n, dwork(pdcov), mct, aict, bct, epsbfg, alpha,
     &            iwork(piw), dwork(pdw), covcal, info )
c      
      RETURN
      END
