c=======================================================================
c     subroutine EXPCOEFCORRS                                      
c
c     This function computes the Spearman rank correlation coefficient 
c     (ex-post correlation coef.)
c                        
c                                 n
c                         6      ---
c     corr. coef = 1 - --------  \   D(i)^2
c                      n(n*n-1)  / 
c                                ---
c                                i=1
c
c     where D(i) = X(i) - Y(i)
c
c----------------------------------------------------------------------
      SUBROUTINE EXPCOEFCORRS ( n, x, y, iwork, dwork, coef, info)
c----------------------------------------------------------------------
c
c     INPUT 
c        n     : size of x (n > 1)                              integer
c        x     : x's value(s) (n)                                double
c        y     : y's value(s) (n)                                double
c
c     WORKSPACE 
c        iwork : n                                              integer
c        dwork : 4*n + 6                                         double
c
c     OUTPUT 
c        coef  : Spearman rank correlation coefficient           double
c        info  : diagnostic argument                            integer
c
c     CALL   
c        AVOCI : Sort a vector in increase order with index
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info
      DOUBLE PRECISION coef, x(*), y(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)     
c
c     local variables
      INTEGER i, ind
      INTEGER pindex, pdx, pdxrank, pdyrank
      DOUBLE PRECISION sum, tmp
c
c     external subroutines
      EXTERNAL AVOCI
c
c-----------------------------------------------------------------------
c
c     initialization
      coef = 0.0
      info = 0
c
c     pointers for integer work space : iwork
c     ---------------------------------------
      pindex = 1
c     pindex : pointer for index vector (n)
c
c     Total size of iwork array = n 
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdx     = 1
c     pdx     : pointer for local vector (n)
      pdxrank = pdx + ( n )
c     pdxrank : pointer for X rank (n)
      pdyrank = pdxrank + ( n )
c     pdyrank : pointer for Y rank (n)      
c
c     Total size of dwork array = 3*n 
c
c----------------------------------------------------------------------
c
c     rank of vectors x and y
      CALL AVOCI (n, x, dwork(pdx), iwork(pindex))
      DO i = 1,n
        ind = dwork(pindex + i - 1) 
        dwork(pdxrank + ind - 1) = i
      ENDDO
      CALL AVOCI (n, y, dwork(pdx), iwork(pindex))
      DO i = 1,n
        ind = iwork(pindex + i - 1) 
        dwork(pdyrank + ind - 1) = i
      ENDDO
c
c     Spearman rank correlation coefficient
      sum = 0.0
      DO i = 1,n
        tmp = dwork(pdxrank + i -1) - dwork(pdyrank + i - 1)
        sum = sum + tmp*tmp
      ENDDO
      tmp = sum /(n*(n*n-1))
      coef = 1.0 - 6.0*tmp
c
      RETURN
      END

