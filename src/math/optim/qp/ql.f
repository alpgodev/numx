c=======================================================================
c
c     subroutine QL                                      
c
c     Quadratic solver (vectorized version)
c
c     MINIMIZE        .5*X'*Q*X + P'*X
c     SUBJECT TO      A(I)*X   =  B(I)    ,  I=1,...,nceg
c                     A(J)*X  <=  B(J)    ,  J=1,...,ncineg
c                     Cinf  <=  X  <=  Csup
c
c-----------------------------------------------------------------------	
      SUBROUTINE ql ( N, Q, P, nceg, ncineg, A, B, Cinf, Csup,
     &                iwork, dwork,
     &                lagr, wopt, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       N      : number of variables                             integer
c       Q      : quadratic matrix (N*N)                           double
c       P      : linear part vector (N)                           double
c       A      : matrix of constraints (N*(nceg+ncineg))          double
c       B      : vector of constraints (nceg+ncineg)              double
c       Cinf   : lower bounds vector (N)                          double
c       Csup   : upper bounds vector (N)                          double
c       nceg   : number of equality constraints                  integer
c       ncineg : number of inequality constraints                integer
c
c     WORKSPACE 
c       iwork  : N                                               integer
c       dwork  : (3*N*N)/2 + 10*N + 2*(nceg + ncineg + 1)
c                 + 2*N*(nceg + ncineg)                           double 
c
c     OUTPUT 
c       lagr   : Lagrange multipliers (2*N+ncineg+nceg)           double
c       wopt   : optimal solution vector (N)                      double
c       info   : diagnostic argument                             integer
c
c-------------------------------------------------------------------------------
c     USED FUNCTION
c
c     From rne/optim/ directory: 
c     ql0001.f
c 
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER N, nceg, ncineg, info
      DOUBLE PRECISION Q(N,*), wopt(*), P(*), lagr(*),
     &                 Cinf(*), Csup(*), A(N,*), B(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER M, ME, MMAX,NMAX, MNN, IOUT, IFAIL, IPRINT, LWAR, LIWAR,
     &        pd1, pd2, pdwork
      DOUBLE PRECISION EPS, ONE
      PARAMETER (ONE = -1.E0)
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
      M    = nceg + ncineg    ! total number of constraints
      ME   = nceg             ! number of equality constraints
      MMAX = M                ! row dimension of A (constraint matrix)
      NMAX = N                ! row dimension of Q (quadratic part)
      MNN  = M + N + N        ! dimension of Lagrange multipliers
      IOUT   = 1              ! desired output
      IPRINT = 1              ! output control (=0 no output)
      EPS =  2.2204e-16       ! termination accuracy
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pd1 = 1
c     pd1 : temporary constraint matrix (N*M)
      pd2  = pd1 + ( N*M )
c     pd2 :  temporary constraint matrix (N*M)
      pdwork = pd2 + ( N*M )
c
c     Total size of dwork array = (3/2)*N*N + 10*N 
c                               + 2*(nceg + ncineg + 1)
c                               + 2*N*(nceg + ncineg)
c
      LWAR   = (3*NMAX*NMAX)/2 + 10*NMAX + 2*(MMAX + 1)
      LIWAR  = N
c----------------------------------------------------------------------------
c
c     Cholesky decomposition internally computed
      iwork(1)=1
c     
c     transform matrix constraint A -> -A'
      CALL RM ( N, M, A, dwork(pd1) )
      CALL PMX ( M, N, dwork(pd1), ONE, dwork(pd2) ) 
c
c     QL Solver
c
c     MINIMIZE        .5*X'*C*X + D'*X
c     SUBJECT TO      A(J)*X  +  B(J)   =  0 ,  J=1,...,ME
c                     A(J)*X  +  B(J)  >=  0 ,  J=ME+1,...,M
c                     XL  <=  X  <=  XU
c
      CALL ql0001(M, ME, MMAX, N, NMAX, MNN,
     &            Q, P, dwork(pd2), B, Cinf, Csup,
     &            wopt, lagr, IOUT, IFAIL, IPRINT,
     &            dwork(pdwork), LWAR, iwork, LIWAR, EPS)
c
c     error management
      IF (IFAIL .EQ. 0)  info = 0     ! optimality conditions satisfied
      IF (IFAIL .EQ. 1)  info =  1007 ! max. nb. iterations exceeded
      IF (IFAIL .EQ. 2)  info = -1008 ! accuracy insufficient
      IF (IFAIL .EQ. 3)  info = -1    ! division by zero
      IF (IFAIL .EQ. 5)  info = -1009 ! lenght of work. array too short
      IF (IFAIL .GT. 10) info = -1010 ! inconsistent constraints
c
      RETURN
      END
