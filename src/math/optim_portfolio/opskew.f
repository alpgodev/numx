c=======================================================================
c
c     Utilities for ALLOCSKEW
c
c-----------------------------------------------------------------------
      SUBROUTINE utskew1 ( n, w, iwork, dwork, Q, info )
c-----------------------------------------------------------------------
c     SDLS input matrix
c     Q = | W    w |
c         | w'   1 |
c      
c-----------------------------------------------------------------------
c
c     INPUT 
c       n       : portfolio size                                 integer
c       w       : benchmark portfolio weights (n)                integer
c
c     WORKSPACE 
c       iwork   : n                                              integer
c       dwork   : n*n                                             double
c
c     OUTPUT 
c       Q       : SDLS input matrix (n+1 by n+1)                  double
c       info    : diagnostic argument                            integer
c
c     CALL   
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info
      DOUBLE PRECISION w(*), Q(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER pindex, pdW
      INTEGER i, p
c
c     external subroutines
      EXTERNAL PV1, YMCPIR
c
c-----------------------------------------------------------------------
c
c     initializations 
      info = 0
c
c     pointers for integer work space  : iwork
c     -----------------------------------------
      pindex = 1
c     pindex : pointer for index vector (n)
c
c     Total size of iwork array = n
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdW = 1
c     pdW : pointer for matrix W (n*n)                
c
c     Total size of dwork array = n*n
c
c-----------------------------------------------------------------------
c                        
c     computes matrix W(p,p) = w(p)*w(p)'
      CALL PV1 ( n, w, dwork(pdW) )
c      
c     Q = | W   w | 
c         | w'  1 |
      p = n + 1
      DO i = 1,n
        iwork(pindex + i - 1) = i
      ENDDO
      CALL YMCPIR ( n, dwork(pdW), p, iwork(pindex), Q, info )      
      DO i = 1,n
        Q(n*i + i)     = w(i)
        Q((n+1)*n + i) = w(i)
      ENDDO
      Q(p*p) = 1.0
c
      RETURN
      END
c
c-----------------------------------------------------------------------
      SUBROUTINE utskew2( n, neq, ccst, bcst, nbEqCst, Ceq, bEq, info )
c-----------------------------------------------------------------------
c     SDLS equalities constraints
c-----------------------------------------------------------------------
c
c     INPUT
c       n       : portfolio size (number of assets)              integer
c       neq     : number equality constraints                    integer
c       ccst    : matrix of constraints (n*(neq+nin))             double
c       bcst    : vector initial of constraints (neq+nin)         double
c
c     WORKSPACE 
c       dwork   :                                                 double
c
c     OUTPUT 
c       nbEqCst : number of equalities constraints               integer          
c       Ceq     : symmetric matrices of equal constraints 
c                 (nbEqCst*(n+1)*(n+1))                           double
c       bEq     : vector of equal constraints (nbEqCst)           double
c       info    : diagnostic argument                            integer
c
c     CALL   
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, neq, nbEqCst, info
      DOUBLE PRECISION ccst(*), bcst(*), Ceq(*), bEq(*)
c
c     local variables
      INTEGER i, j, p
      DOUBLE PRECISION ZERO
      PARAMETER ( ZERO = 0.E+0 )
c
c     external subroutines
      EXTERNAL IVX, IMX
c
c-----------------------------------------------------------------------
c
c     initializations 
      info = 0      ! diagnostic argument
      p    = n + 1  ! problem size (number of assets + 1)
c
c     number of eq. constraints
      nbEqCst = neq + 1
c
c-----------------------------------------------------------------------
c
c     bEq, Ceq initialization
      CALL IVX ( nbEqCst, bEq, ZERO )
      CALL IVX ( nbEqCst*p*p, Ceq, ZERO )
c
c     SDLS equalities constraints (u=1)
      bEq(1)   = 1.0
      Ceq(p*p) = 1.0
c
c     general equalities constraints
c     case neq > 0      
      IF (neq .GT. 0) THEN
        DO i = 1,neq
            bEq(1+i) = bcst(i)
            DO j = 1,n
                Ceq(p*p*i + n*j + j) = 0.5*ccst(neq*(i-1) + j)
                Ceq(p*p*i + p*n + j) = 0.5*ccst(neq*(i-1) + j)
            ENDDO
        ENDDO        
      ENDIF       
c
c
      RETURN
      END
c      
c-----------------------------------------------------------------------
      SUBROUTINE utskew3 ( n, p, x, w, rho, cov, neq, nin, ccst, bcst,
     &                     cinf, csup, mu, sigma, skew,
     &                     iwork, dwork,
     &                     nbIneqCst, Cineq, bLowerIneq, bUpperIneq, info)
c-----------------------------------------------------------------------
c     SDLS inequalities constraints
c-----------------------------------------------------------------------
c
c     INPUT
c       n          : number of values (>3 )                      integer
c       p          : portfolio size (number of assets)           integer
c       x          : returns values (n*p)                         double
c       w          : benchmark portfolio weights (p)              double   
c       rho        : expected returns (p)                         double
c       cov        : covariance matrix (n*n)                      double
c       neq        : number equality constraints                  integer
c       nin        : number inequality constraints               integer
c       ccst       : matrix of constraints (p*(neq+nin))          double
c       bcst       : vector initial of constraints (neq+nin)      double
c       cinf       : lower bound (p)                              double
c       csup       : upper bound (p)                              double
c       mu         : performance target                           double
c       sigma      : volatility target                            double
c       skew       : skewness constraints                         double
c
c     WORKSPACE 
c       iwork      : p                                           integer
c       dwork      : p*(16*p + n + 32)                            double
c
c     OUTPUT 
c       nbIneqCst  : number of equalities constraints            integer          
c       Cineq      : symmetric matrices of inequal constraints 
c                   (nbIneqCst*(p+1)*(p+1))                       double
c       bLowerIneq : lower bounds of inequal cst. (nbIneqCst)     double
c       bUpperIneq : upper bounds of inequal cst. (nbIneqCst)     double
c       info       : diagnostic argument                         integer
c
c     CALL   
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, p, neq, nin, nbIneqCst, info
      DOUBLE PRECISION mu, sigma, skew, x(*), w(*), rho(*), cov(*), 
     &                 ccst(*), bcst(*), cinf(*), csup(*),
     &                 Cineq(*), bLowerIneq(*), bUpperIneq(*)
c
c     workspaces
      DOUBLE PRECISION dwork(*)
      INTEGER iwork(*) 
c
c     local variables
      INTEGER pindex, pdm3, pdtansork, pdW, pdH, pdHk, pdWH, pdtmpH,
     &        pdexarsk
      INTEGER i, j, k, m
      DOUBLE PRECISION mtarget, m3, vol, vol3, skewb, tmp, ZERO, sigma3
      PARAMETER ( ZERO = 0.E+0 )
c
c     external subroutines
      EXTERNAL IMX, PMX, SM, YM, PMC, TM, TANSORK, MOMENT3, EXARVO, 
     &         EXARSK, YMCPIR
c
c-----------------------------------------------------------------------
c
c     initializations 
      info = 0      ! diagnostic argument
      m    = p + 1  ! problem size (number of assets + 1)     
c
c     pointers for integer work space  : iwork
c     -----------------------------------------
      pindex = 1
c     pindex : pointer for index vector (p)
c
c     Total size of iwork array = p  
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdexarsk = 1
c     pdexarsk : pointer for EXARSK (p*(5*p + 1))      
      pdW = pdexarsk + ( p*(5*p + 1 ))
c     pdW : pointer for matrix W (p*p)
      pdtansork = pdW + ( p*p )
c     pdtansor : pointer for TANSORK (p)
      pdHk = pdtansork + ( p )
c     pdHk : pointer for Hk matrix (p*p)
      pdtmpH = pdHk + ( p*p )
c     pdtmpH : pointer for tamporary matrix (p*p)
      pdH = pdtmpH + ( p*p )
c     pdH : pointer for matrix H (p*p)
      pdWH = pdH + ( p*p )
c     pdWH : pointer for matrix WH (p*p)                                              
      pdm3 = pdWH + ( p*p )
c     pdm3 : pointer for MOMENT3 (p*(5*p + 1))
c
c     Total size of dwork array = p*(15*p + 3)
c
c-----------------------------------------------------------------------
c
c     number of inequalities:
c              cst (nin)  
c            + risk budget (1) + target return (1) 
c            + 3th moment constraint (1) + lower/upper bounds (p)
      nbIneqCst = 1 + 1 + 1 + p + nin
c      
c     computes matrix W(p,p) = w(p)*w(p)'
      CALL PV1 ( p, w, dwork(pdW) )
c
c     bLowerIneq, BUpperIneq and Cineq initialization
      CALL IVX ( nbIneqCst, bLowerIneq, ZERO )
      CALL IVX ( nbIneqCst, bUpperIneq, ZERO )
      CALL IVX ( m*m*nbIneqCst, Cineq, ZERO )
c
c     risk budget constraint <W, Cov> <= sigma2
      bLowerIneq(1) = 0.0
      bUpperIneq(1) = sigma*sigma
      DO i = 1,p
        iwork(pindex + i - 1) = i
      ENDDO
      CALL YMCPIR ( p, cov, m, iwork(pindex), Cineq, info )
c
c     SDLS target return constraint (mu)
      bLowerIneq(2) = mu
      bUpperIneq(2) = 1.E+5
      DO i = 1,p
        Cineq(m*m + m*p + i) = rho(i)
      ENDDO
c
c     SDLS 3th moment constraint (skew)
      CALL MOMENT3 ( n, p, w, x, dwork(pdm3), m3, info)
c
c     ex-ante volatility (benchmark)      
      CALL EXARVO (p, cov, w, vol, info)
      IF (info .LT. 0) RETURN
c
c     ex-ante skewness (benchmark)
      CALL EXARSK ( n, p, w, x, cov, dwork(pdexarsk), skewb, info)
      IF (info .LT. 0) RETURN
c
c      momentTarget = (skewTarget + skewb)*(vol**3)
      vol3   = vol*vol*vol      
      sigma3 = sigma*sigma*sigma
c      mtarget = skew*sigma3
      mtarget = skew*sigma3 + skewb*vol3
c          
      bLowerIneq(3) = mtarget + m3
      bUpperIneq(3) = 1.E+8
c 
      CALL IVX ( p*p, dwork(pdtmpH), ZERO ) ! initialization
      DO k = 1,p
c
c       computes H[k]
        CALL TENSORK ( n, p, k, x, dwork(pdtansork), dwork(pdHk), info)
        IF (info .LT. 0) RETURN
c
c       computes matrix Wref*H[k]
        CALL PMC ( p, dwork(pdW), dwork(pdHk), dwork(pdWH) )  
c
c       computes <Wref, H[h]> := Trace(Wref*H[k]) 
        CALL TM ( p, dwork(pdWH), tmp )
c
c       constraint
        Cineq(m*m + m*m + m*k)     = 0.5*tmp
        Cineq(m*m + m*m + m*p + k) = 0.5*tmp
c
c       computes w(k)*H(k)
        CALL PMX ( p, p, dwork(pdHk), w(k), dwork(pdWH) )
c
c       computes H = sum(w(k)*H(k))
        CALL SM ( p, p, dwork(pdWH), dwork(pdtmpH), dwork(pdH) )
c
c       copy H -> tmpH
        CALL YM ( p, p, dwork(pdH), dwork(pdtmpH) )
      ENDDO
c      
      DO i = 1,p
        DO j = 1,p
            Cineq(m*m + m*m + m*(i-1) + j) = 
     &            dwork(pdH + p*(i-1) + j)
        ENDDO
      ENDDO
c
c     SDLS lower/upper constraints (Cinf/Csup)
      DO i = 1,p
        bLowerIneq(3 + i) = cinf(i)
        bUpperIneq(3 + i) = csup(i)
        Cineq(m*m + m*m + m*m + m*m*(i-1) + m*p + i) = 1.0
      ENDDO
c
c     general inequalities constraints (ccst)
c     case nin > 0      
      IF (nin .GT. 0) THEN
        DO i = 1,nin
            bLowerIneq(p + 1 + 1 + 1 + i) = bcst(neq + i) 
            bUpperIneq(p + 1 + 1 + 1 + i) = 1.0E+15
            DO j = 1,p
                Cineq(m*m + m*m + m*m + m*m*p + m*p + j) 
     &        = ccst(neq*p + nin*(i-1) + j)
            ENDDO
        ENDDO        
      ENDIF  
c
      RETURN
      END

