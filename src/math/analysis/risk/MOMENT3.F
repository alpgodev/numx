c=======================================================================
c
c     subroutine MOMENT3.F - Estimation of 3-th moment
c
c              n
c             ---         
c     h(w) =  \  w(i)*w(j)*w(k)*Hijk  
c             / 
c             ---
c             t=1
c
c                     n
c                    ---         
c          = < w*w', \  w(k)*H[k] >  
c                    / 
c                    ---
c                    k=1
c
c-----------------------------------------------------------------------
      SUBROUTINE MOMENT3 ( n, p, w, x, dwork, m3, info)
c-----------------------------------------------------------------------
c
c     INPUT 
c       n      : number of values                                integer
c       p      : number of asset(s)                              integer
c       w      : weights (p)                                      double
c       x      : values (n*p)                                     double
c
c     WORKSPACES
c       dwork  : p*(5*p + 1)                                      double
c
c     OUTPUT 
c       m3     : 3-th moment                                      double
c       info   : diagnostic argument                             integer
c
c     CALL 
c       TANSORK
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, p, info
      DOUBLE PRECISION m3, w(*), x(*)
c
c     workspaces
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER k, pdW, pdtansor, pdHk, pdtmpH, pdH, pdWH
      DOUBLE PRECISION ZERO
      PARAMETER ( ZERO = 0.0)
c
c-----------------------------------------------------------------------
c
c     initializations 
      info = 0
      m3   = 0.0
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdW = 1
c     pdW : pointer for matrix W (p*p)
      pdtansor = pdW + ( p*p )
c     pdtansor : pointer for TANSORK (p)
      pdHk = pdtansor + ( p )
c     pdHk : pointer for Hk matrix (p*p)
      pdtmpH = pdHk + ( p*p )
c     pdtmpH : pointer for tamporary matrix (p*p)
      pdH = pdtmpH + ( p*p )
c     pdH : pointer for matrix H (p*p)
      pdWH = pdH + ( p*p )
c     pdWH : pointer for matrix WH (p*p)                                              
c
c     Total size of dwork array = p*(5*p + 1) 
c
c-----------------------------------------------------------------------
c                           
c     computes matrix W(p,p) = w(p)*w(p)'
      CALL PV1 ( p, w, dwork(pdW) )
c
c     computes sum of tansor
c          p
c         ---         
c         \   w(k)*H[k]  
c         / 
c         ---
c         k=1
c
      CALL IMX ( p, p, dwork(pdtmpH), ZERO ) ! initialization
      DO k = 1,p
        CALL TENSORK ( n, p, k, x, dwork(pdtansor), dwork(pdHk), info)
        IF (info .LT. 0) RETURN
c
c       computes w(k)*H(k)
        CALL PMX ( p, p, dwork(pdHk), w(k), dwork(pdHk) )     
c
c       computes H = sum(w(k)*H(k))
        CALL SM ( p, p, dwork(pdHk), dwork(pdtmpH), dwork(pdH) )  
c
c       copy H -> tmpH
        CALL YM ( p, p, dwork(pdH), dwork(pdtmpH) )   
      ENDDO
c
c     computes matrix W*H
      CALL PMC ( p, dwork(pdW), dwork(pdH), dwork(pdWH) )
c
c     3-th moment is egal to Trace(W*H) 
      CALL TM ( p, dwork(pdWH), m3 )
c 
      RETURN
      END
