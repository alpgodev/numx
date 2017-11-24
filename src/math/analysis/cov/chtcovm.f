c=======================================================================
c
c     subroutine CHTCOVM
c
c     covariance matrix with poor historic
c    (computes the covariances independantly with the available
c     historic and correct the matrix with SDLS)
c
c-----------------------------------------------------------------------
      SUBROUTINE chtcovm ( ndate, n, nintro, dasset, age, unborn,
     &                     hole, iwork, dwork, mean, cov, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c            ndate  : lagest estimation period                   integer
c            n      : nb of assets                               integer
c            nintro : nb of days after the birth of an asset 
c                     before introducing it in the covariance 
c                     matrix                                     integer
c            dasset : values (returns, spreads...) of the assets
c                     during the estimation period
c                     (with value = hole if unknown)
c                     vectorized matrix (ndate*nasset)            double
c            age    : age of the assets (number of dates since
c                     which their quotations are known in dasset)
c                     0<=age(i)<=ndate
c                     vector(nasset)                             integer
c            unborn : stand for unborn quotations                 double
c            hole   : stands for unknown value                    double
c
c     WORKSPACE 
c            iwork  : 14*n + 3                                   integer
c            dwork  : n*(8*n + 40) + 1                            double
c
c     OUTPUT 
c            mean   : mean of the assets (n)                      double
c            cov    : covariance matrix (n*n)                     double
c            info   : diagnostic argument                        integer
c
c     CALL   
c        CHTMCM  : mean of assets with a poor historic (with unknown values)
c        SDLS    : Semi-Definite Least-Square optimization
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o arguments
      INTEGER ndate, n, nintro, info, age(*)
      DOUBLE PRECISION dasset(ndate,*), cov(n,*), mean(*), 
     &                 unborn, hole
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER i, j, k, nestim, ind1st, pisdls, pdsdls, pdempcov, pdaict,
     &        pdbct
      DOUBLE PRECISION mcte, alpha, epsbfg, som
      PARAMETER ( epsbfg = 1.E-5 )
c
c-----------------------------------------------------------------------
c
c     initializations
      info = 0
c      
c     pointers for integer work space : iwork
c     ---------------------------------------
c     pisdls : pointer for SDLS iwork ( 14*n + 3 ) 
      pisdls = 1      
c
c     Total size of iwork: 14*n + 3      
c
c     pointers for double precision work space : dwork
c     ------------------------------------------------         
      pdsdls = 1
c     pdsdls : pointer for SDLS dwork n*(6*n + 40)
c
      pdempcov = pdsdls + ( n*(6*n + 40) )
c     pdempcov : pointer for covariance matrix ( n*n )     
c
      pdaict = pdempcov + ( n*n )
c     pdaict : pointer for mct Ai constraint matrices (n*n)
c
      pdbct = pdaict + ( n*n ) 
c     pdbct : pointer for constraints ( 1 )  
c
c     Total size of dwork = ( n*(6*n + 40))    
c                         + n*n
c                         + n*n
c                         + 1
c                         = ( n*(8*n + 40)) + 1
c
c-----------------------------------------------------------------------
c
c     computing means of rents
c     by construction, unknown values of this historic are
c     unborn (historic was 
c        1) interpolated to fill holes
c        2) built with BLDHIST.F (to omit dead values)
c        before calling CHTCOVM)
      CALL chtmcm ( ndate, n, dasset, unborn, age, mean, info )
      IF (info .LT. 0) RETURN
c
c     computing empiric covariance matrix
      do i = 1,n
         do j=1,i
c
c           for each couple of assets (i,j) s.t. i<=j
c           compute the estimation period = min(age(i),age(j))
            if (age(i) .le. age(j)) then
               nestim = age(i)
            else
               nestim = age(j)
            end if
c           if asset i or j has a too short historic (age<nintro)
c           then empcov(i,j) = unknown = hole
            if (nestim .lt. nintro) then
               dwork(pdempcov+(j-1)*n + i - 1) = hole
               info = 1
            else
c
c              if nestim too big, nestim:=ndate
               if (nestim .gt. ndate) nestim = ndate
c 
c              compute the 1st index of the estimation period
               ind1st = ndate - nestim + 1
c              compute the empirical covariance bewteen those 2 assets
c              with nestim historical values
               som = 0
               do k = ind1st,ndate
                  som = som + (dasset(k,i)-mean(i))
     &                       *(dasset(k,j)-mean(j))
               end do
                dwork(pdempcov+(j-1)*n + i - 1) = som/(nestim-1)
            end if
c
            if (i .ne. j) then
c              empcov is symmetric : no use to make this computation 
c              twice
               dwork(pdempcov+(i-1)*n + j - 1) 
     &              = dwork(pdempcov+(j-1)*n + i - 1)
            end if
c
         end do
      end do
c
c     covariance matrix correction by SDLS without constraints
c     min. eigenvalue = 1.e-10
      mcte  = 0.
      alpha = 1.e-10
      CALL sdls ( n, dwork(pdempcov), mcte, dwork(pdaict), dwork(pdbct),
     &            epsbfg, alpha,
     &            iwork(pisdls),dwork(pdsdls),
     &            cov, info)
c
      RETURN
      END

