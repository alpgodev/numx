c=======================================================================
c
c     subroutine CALAPT                                      
c
c     APT (arbitrage Pricing Theory) model -> alpha and beta coefficients
c     for a given set of factor(s) X
c
c            Y = alpha + X*beta'
c 
c     where Y is the matrix (n*p) of asset(s) value(s), X the matrix (n*q) 
c     of factor(s) value(s) and beta=[beta(1), ...,beta(q)]
c
c-----------------------------------------------------------------------
      SUBROUTINE calapt ( n, p, q, Y, X, iwork, dwork, alpha, beta,info)
c-----------------------------------------------------------------------
c
c     INPUT 
c            n      : number of value(s) (n > 0)                 integer
c            p      : number of asset(s) (p > 0)                 integer
c            q      : number of factors  (q > 0)                 integer
c            Y      : asset(s) values (n*p)                       double
c            X      : factor(s) values (n*q)                      double
c
c     WORKSPACE 
c            iwork  : q + 1                                      integer
c            dwork  : (q + 1)*(n + 2*p + q + 2)                   double
c
c     OUTPUT 
c            alpha  : alpha coefficient(s) (p)                    double
c            beta   : beta coefficient(s) (q*p)                   double
c            info   : diagnostic argument                        integer
c
c     CALL    
c            IVX    : initialization at a scalar of a vector
c            YV     : copy a vector in a vector
c            PRMM   : computing M'*M 
c            PMTM   : product of 2 vectorized full matrices
c                    (M matrix n*m, N matrix n*c, gives M'*N matrix m*c)
c            DGETRF : LU factorization (cf. LAPACK)
c            DGETRI : inverse matrix (cf. LAPACK)
c            PM     : product of 2 vectorized full matrices
c                     (M matrix n*m, N matrix m*c, gives M*N matrix n*c)
c            YLMV   : copy a row of a vectorized matrix in a vector
c            YMP    : copy a part of a vectorized matrix in a vectorized matrix
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, p, q, info
      DOUBLE PRECISION Y(*), X(*), alpha(*), beta(q,*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER nfp1, piw, pf, pdw, pftf, pftr, px
      DOUBLE PRECISION DUN
      PARAMETER ( DUN = 1.D0 )
c
c-----------------------------------------------------------------------
c
c     initializations
      info = 0 
      nfp1 = q + 1
c
c     pointers for integer work space : iwork
c     ---------------------------------------
      piw = 1
c     piw   : pointer for DGETRF, DGETRI  who need  nfp1
c
c     Total size of iwork array = nfp1 = q + 1
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pf = 1
c     pf   : pointer for matrix F, so ( n*nfp1 ) more
      pftf = pf + ( n*nfp1 )
c     pftf : pointer for matrix F'*F, so ( nfp1*nfp1 ) more
      pftr = pftf + ( nfp1*nfp1 )
c     pftr : pointer for matrix F'*R, so ( nfp1*p ) more
      px = pftr + ( nfp1*p )
c     px   : pointer for matrix X, so ( nfp1*p ) more
      pdw = px + ( nfp1*p )
c     pdw  : DGETRI, ( nfp1 )
c
c     Total size of dwork array = (n*nfp1) + (nfp1*nfp1) + (nfp1*p)
c                               + (nfp1*p) + nfp1
c                               = nfp1*( n + 2*p + nfp1 + 1 )
c                                 with nfp1 = (q + 1)
c
c     Total size of dwork array = (q + 1)*(n + 2*p + q + 2)
c
c-----------------------------------------------------------------------
c
c     construction of matrix F (1. in first column, rfdata in the rest)
      CALL IVX ( n, dwork(pf), DUN )
      CALL YV ( q*n, X, dwork(pf + n) )
c
c     F'*F
      CALL PRMM ( n, nfp1, dwork(pf), dwork(pftf) )
c
c     F'*R
      CALL PMTM ( n, nfp1, p, dwork(pf), Y, dwork(pftr) )
c     
c     DGETRF (cf. LAPACK) computes an LU factorization of a general 
c     M-by-N matrix A using partial pivoting with row interchanges.
c
c     The factorization has the form
c         A = P * L * U
c     where P is a permutation matrix, L is lower triangular with unit
c     diagonal elements (lower trapezoidal if m > n), and U is upper
c     triangular (upper trapezoidal if m < n).
      CALL DGETRF ( nfp1, nfp1, dwork(pftf), nfp1, iwork(piw), info )
      IF (info .LT. 0) RETURN
c
c     DGETRI (cf. LAPACK) computes the inverse of a matrix using 
c     the LU factorization computed by DGETRF.
c
c     This method inverts U and then computes inv(A) by solving the system
c     inv(A)*L = inv(U) for inv(A).
      CALL DGETRI ( nfp1, dwork(pftf), nfp1, iwork(piw), dwork(pdw),
     &              nfp1, info )
      IF (info .LT. 0) RETURN
c
c     [alpha, beta] := inv(F'*F) * (F'*R)
      CALL PM ( nfp1, nfp1, p, dwork(pftf), dwork(pftr), dwork(px) )
c
c     copy [alpha, beta] from X
      CALL YLMV ( nfp1, p, dwork(px), 1, alpha, info )
      CALL YMP( nfp1, p, dwork(px), q, p, 2, 1, beta, info )
c
      RETURN
      END

