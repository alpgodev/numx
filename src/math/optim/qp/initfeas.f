c=======================================================================
c
c     subroutine initfeas                                   
c
c     test QP problem feasibility
c
c-----------------------------------------------------------------------
      SUBROUTINE initfeas(n, cinf, csup, neq, nin, ccst, bcst,
     &                    iwork, dwork, wopt, info)
c-----------------------------------------------------------------------
c
c     INPUT 
c       n      : problem size                                    integer
c       cinf   : lower bounds (n)                                 double
c       csup   : upper bounds (n)                                 double
c       neq    : number of equality constraints                  integer
c       nin    : number of inequality constraints                integer
c       ccst   : constraints matrix (n*(neq + nin + 1))           double
c       bcst   : constraints vector (neq + nin + 1)               double
c
c     WORKSPACE 
c       iwork  : 3*n + 2*nin + neq + 1                           integer
c       dwork  : n*(neq+nin+2*n+15)+2*neq+4*nin+13                double 
c
c     OUTPUT 
c       wopt   : optimal portfolio (n)                            double
c       info   : diagnostic argument                             integer
c
c     CALL
c       buildfeas, qp, YV
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, neq, nin, info
      DOUBLE PRECISION cinf(*), csup(*), ccst(*), bcst(*), wopt(*)
c
c     workspaces 
      DOUBLE PRECISION dwork(*)
      INTEGER iwork(*)
c
c     local variables 
      INTEGER i, n1, infot, piw, pdwopt, pdcov, pdlin, pdccst,
     &        pdcinf, pdcsup, pdlagr, pdw
      DOUBLE PRECISION dmax, temp, EPS, DZERO
      PARAMETER (EPS=1.E-8, DZERO=0.0)
c
c-----------------------------------------------------------------------
c
c     pointers for integer work space  : iwork
c     ---------------------------------------- 
      piw  = 1
c     piw  : pointer for QP (3*n + neq + 2*nin + 1)
c             
c     Total so size of iwork array = 3*n + neq + 2*nin + 1
c     
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      n1 = n + 1 
      pdwopt = 1 
c     needs n1
c        so n1
      pdcov = pdwopt + n1
c     needs n1*n1 
c       so n1*n1+n1
      pdlin = pdcov + n1*n1
c     needs n1 
c       so n1*n1+2*n1
      pdccst = pdlin + n1
c     needs n1*(neq+nin)
c       so n1*(neq+nin+n1+2)
      pdcinf = pdccst + n1*(neq+nin)
c     needs n1
c      so n1*(neq+nin+n1+3)
      pdcsup = pdcinf + n1
c     needs n1
c      so n1*(neq+nin+n1+4)
      pdlagr = pdcsup + n1
c     needs n1+neq+nin   
c       so n1*(neq+nin+n1+5)+neq+nin  
      pdw = pdlagr + n1+neq+nin
c     needs n1*(n1 + 6) + 2*nin
c       so n1*(neq+nin+2*n1+11)+neq+3*nin  = n*(neq+nin+2*(n+1)+11)+neq+3*nin+neq+nin+2*(n+1)+11
c           = n*(neq+nin+2*n+15)+2*neq+4*nin+13
c
c------------------------------------------------------------------------
c
      CALL buildfeas (n, cinf, csup, neq, nin, ccst, n1, dwork(pdcinf),
     &                dwork(pdcsup), dwork(pdccst), dwork(pdcov), 
     &                dwork(pdlin), info)
      IF (info .LT. 0) THEN
        info = -8888
        RETURN
      ENDIF
c
c     quadratic solver
      CALL qp(n1, dwork(pdcov), dwork(pdlin), neq, nin,
     &        dwork(pdccst), bcst, dwork(pdcinf), dwork(pdcsup),
     &        iwork(piw), dwork(pdw), dwork(pdlagr),
     &        dwork(pdwopt), info)
      CALL YV(n, dwork(pdwopt), wopt)
      IF (info .NE. 0) RETURN
      CALL PMTV(n, (neq+nin), ccst, wopt, dwork(pdw))
      dmax = 0.0
      DO i = 1,(neq+nin)
        temp = dwork(pdw+i-1)-bcst(i)
        IF (temp .GT. dmax) dmax = temp
      ENDDO
      IF (dmax .GT. (DZERO+EPS)) THEN
            info = 1111
            RETURN
      ENDIF
      CALL checklinbox (n, wopt,
     &                  neq, nin, ccst, bcst, cinf, csup,
     &                  dwork(pdw), infot)
      IF (infot .LT. 0) THEN
        info = infot
        RETURN
      ENDIF
      RETURN
      END
