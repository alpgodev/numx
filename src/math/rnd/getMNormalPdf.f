c=======================================================================
c
c     function getMNormalPdf              
c
c     Returns the density function of a multivariate normal random variable
c
c-----------------------------------------------------------------------
c
c     Copyright (c) 2014 NumX
c     All rights reserved.
c 
c     This software is the confidential and proprietary information
c     of NumX. You shall not disclose such Confidential
c     Information and shall use it only in accordance with the terms
c     of the licence agreement you entered into with NumX.
c
c     Author : Yann Vernaz
c
c-----------------------------------------------------------------------
      SUBROUTINE getMNormalPdf(d, x, mu, cov, iwork, dwork, z, info)
c-----------------------------------------------------------------------
c
c     INPUT 
c       d        : dimension (d>0)                               integer
c       x        : data (d)                                       double
c       mu       : mean (d)                                       double
c       cov      : cov. matrix (d*d)                              double
c
c     WORKSPACE
c       iwork    : d                                             integer
c       dwork    : d*(2*d + 7)                                    double
c
c     OUTPUT
c       z        : multivariate normal Pdf                        double
c       info     : diagnostic argument                           integer
c
c----------------------------------------------------------------------- 
c
c     arguments i/o
      INTEGER d, info
      DOUBLE PRECISION z, x(*), mu(*), cov(d,*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables      
      LOGICAL select
      INTEGER i, sdim, lwork, pdcov, pdeig, piw, pip, pdv, pdw
      DOUBLE PRECISION PI, invPI, detCOV, s, EPS
      PARAMETER (PI=3.1415926535897932384626433832795028841971693993751)
      PARAMETER (EPS = 1.E-30) 
c
c     intrinsic functions
      INTRINSIC MIN, MAX, SQRT, EXP
c
c-----------------------------------------------------------------------
c
c     initialisations 
      info = 0
c
c     pointers for integer work space  : iwork
c     ----------------------------------------
      piw = 1
c     piw : pointer for DGEES, so ( d ) more
c
c     Total size of dwork array = ( d ) 
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdcov = 1
c     pdcov : pointer for cov. matrix, so ( d*d ) more
      pdeig = pdcov + ( d*d )
c     pdeig : pointer for eigenvalues, so ( d ) more
      pip   = pdeig + ( d )
c     pip   : pointer for DGEES, so ( d ) more
      pdv   = pip + ( d )
c     pdv   : pointer for DGEES, so ( d*d ) more
      pdw   = pdv + ( d*d )
c     pdw   : pointer for DGEES, so ( 5*d ) more
c
c     Total size of dwork array = d*(2*d + 7) 
c
c----------------------------------------------------------------------
c
c     invPI = (2*PI)^(-d/2)
      invPI = (2.0*PI)**(-d/2)   
c
c     case d=1
      IF (d .EQ. 1) THEN
        variance = MAX(cov(1,1),EPS)
        detCOV   = 1.0/SQRT(variance)
        s        = ((x(1)-mu(1))**2)/variance
        z        = invPI*detCOV*EXP(-0.5*s) 
        z        = MIN(z, 1.0)
        z        = MAX(z, 0.0)
        RETURN
      ENDIF
c
c     |cov|
      CALL YM ( d, d, cov, dwork(pdcov) )
      sdim = 0
      lwork = 5*d
      CALL DGEES( 'V', 'N', select, d, dwork(pdcov), d,
     &            sdim, dwork(pdeig), dwork(pip), dwork(pdv),
     &            d, dwork(pdw), lwork, iwork(piw), info )
      detCOV = 1.0
      DO i = 1,d
        detCOV = detCOV*dwork(pdeig + i - 1)  
      ENDDO
c
c     |cov|^(-1/2)  
      detCOV = 1./SQRT(detCOV)
c
c     cov^(-1)     
      CALL JMS ( d, cov, dwork(pdv), dwork(pdcov), info )
c
c     y =  x - mu
      CALL DV ( d, x, mu, dwork(pdeig) )
c
c     (x-mu)'*invCOV*(x-mu)
      CALL OVTMCV ( d, dwork(pdcov), dwork(pdeig), s )
c
c     EXP(-1/2*[(x - mu)'*(cov^(-1))*(x-mu)])
      s = EXP(-0.5*s)
c
c     multivariate Normal Pdf
      z = invPI*detCOV*s
      z = MIN(z, 1.0)
      z = MAX(z, 0.0)
c
      RETURN
      END
