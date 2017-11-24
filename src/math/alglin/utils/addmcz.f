c=======================================================================
c
c     subroutine addmcz                                 
c
c     This function adds columns of zeros to a matrix
c
c-----------------------------------------------------------------------
        SUBROUTINE addmcz(nl, nc, mat, nczeros, matout, info)
c-----------------------------------------------------------------------
c
c     INPUT 
c       nl      : number of lines in the matrix mat                              integer
c       nc      : number of column in the matrix mat                              integer
c       mat     : matrix constituing the upper block of matout         (nl, nc)-double 
c       nczeros : number of zero lines to be addded                            integer
c
c     OUTPUT
c       matout  : vertical concatenation of matup and matdown (nl, nc + nczeros)-double 
c                (ncdown should be equal to ncup)
c       info    : diagnostic argument                             integer   
c 
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c      
c     i/o arguments          
      INTEGER nl, nc, nczeros,info
      DOUBLE PRECISION mat(nl, * ), matout(nl,*)

c     local variables            
      INTEGER ncnew, size, idxone, i, j
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0) 
c       
      ncnew  = nc + nczeros
      size   = nl * ncnew
      idxone = 1
      DO j = 1,nl
        DO i = 1,nczeros
            matout(j,nc+i) = ZERO
        ENDDO
      ENDDO
c
c     copy a part of a matrix in a part of a matrix      
      CALL YMPMP ( nl, nc, nl, nc, mat, idxone, idxone, nl, ncnew,
     &             matout, idxone, idxone, info )
                     
      RETURN
      END
