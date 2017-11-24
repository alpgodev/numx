c=======================================================================
c
c     subroutine buildfeas                                   
c
c     test the problem feasibility
c
c-----------------------------------------------------------------------
      SUBROUTINE buildfeas (n, cinf, csup, neq, nin, ccst, n1, cinf1,
     &                      csup1, ccst1, cov, lin, info)
c-----------------------------------------------------------------------
c
c     INPUT 
c       n      : problem size                                    integer
c       cinf   : lower bounds (n)                                 double
c       csup   : upper bounds (n)                                 double
c       neq    : number of equality constraints                  integer
c       nin    : number of inequality constraints                integer
c       ccst   : constraints matrix (n*(neq + nin))               double
c
c     OUTPUT 
c       n1     : problem size + one variable                     integer
c       cinf1  : lower bounds (n1)                                double
c       csup1  : upper bounds (n1)                                double
c       ccst1  : constraints matrix (n1*(neq + nin))              double
c
c     CALL
c       YV, IMX, IVX, YMPMP
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, neq, nin, n1, info
      DOUBLE PRECISION cinf(*), csup(*), ccst(*), cinf1(*), csup1(*), 
     &                 ccst1(*), cov(*), lin(*)

c     local variables      
      INTEGER ncst, IONE
      DOUBLE PRECISION DZERO, INFINI
      PARAMETER (IONE=1, DZERO=0.0, INFINI=1.E+15)
c      
      CALL YV(n, cinf, cinf1) ! copy cinf1 <- cinf
      CALL YV(n, csup, csup1) ! copy csup1 <- csup
      cinf1(n1) = -1.0
      csup1(n1) = INFINI
      ncst = neq + nin
      CALL IMX(n1, n1, cov, DZERO)
      CALL IVX(n, lin, DZERO)
      lin(n1)= 1.0
      CALL IMX(n1, ncst, ccst1, DZERO)
      CALL YMPMP(n, ncst, n, ncst, ccst, IONE, IONE,
     &           n1, ncst, ccst1, IONE, IONE, info)
      RETURN
      END
