c=======================================================================
c
c     subroutine addmlz                                  
c
c     This function adds lines of zeros to a matrix.
c
c-----------------------------------------------------------------------
        SUBROUTINE addmlz(nl, nc, mat, nlzeros, matout, info)
c-----------------------------------------------------------------------
c
c     INPUT 
c       nl      : number of lines in the matrix mat              integer
c       nc      : number of column in the matrix mat             integer
c       mat     : matrix (nl, nc)                                 double 
c       nlzeros : number of zero lines to be added               integer
c
c     OUTPUT
c       matout  : matrix (nl + nlzeros, nc)                       double                              
c       info    : diagnostic argument                            integer
c                 -1 ncup != ncdown   
c 
c
c------------------------------------------------------------------------
c
      IMPLICIT NONE
c      
c     i/o arguments           
      INTEGER nl, nc, nlzeros,info
      DOUBLE PRECISION mat(nl, * ), matout((nlzeros+nl),*)
c
c     local variables
      INTEGER nlnew, size, idxone, i, j 
      DOUBLE PRECISION ZERO 
      PARAMETER (ZERO=0.0)
c       
      nlnew  = nl + nlzeros
      size   = nlnew * nc
      idxone = 1
      DO j = 1,nlzeros
        DO i = 1,nc
            matout(nl+j,i) = ZERO
        ENDDO
      ENDDO
c
c     copy a part of a matrix in a part of a matrix            
      CALL YMPMP ( nl, nc, nl, nc, mat, idxone, idxone, nlnew, nc,
     &             matout, idxone, idxone, info )          
      RETURN
      END
