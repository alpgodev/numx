c=======================================================================
c
c     MVE
c
c     Efficient Frontier (Lagrange Method)				
c
c-----------------------------------------------------------------------
      SUBROUTINE mve( nbValues, retA, retB, retC, volA, volB, volC,
     &                minExpectedReturn, maxExpectedReturn,
     &                iwork, dwork,
     &                mveRet, mveVol, info)
c-----------------------------------------------------------------------
c
c	Inputs
c		nbValues: integer, number of point(s) (>1)
c		retA: double, mean return of portfolio A
c		retB: double, mean return of portfolio B
c		retC: double, mean return of portfolio C
c		volA: double, volatility of portfolio A (>=0)
c		volB: double, volatility of portfolio B (>=0)
c		volC: double, volatility of portfolio C (>=0)
c		minExpectedReturn: double, minimum return of universe
c		maxExpectedReturn: double, maximum return of universe
c
c	Outputs 
c		mveReturns: double array of dimension nbValues, vector of mean returns
c		mveVolatilities: double array of dimension nbValues, vector of volatilities
c		info: integer array of dimension 1, diagnostic argument
c
c	Call
c		nls 
c
c======================================================================= 
c			
      IMPLICIT NONE
c
c     i/o arguments
      INTEGER nbValues, info
      DOUBLE PRECISION retA, retB, retC, volA, volB, volC
      DOUBLE PRECISION minExpectedReturn, maxExpectedReturn
      DOUBLE PRECISION mveRet(*), mveVol(*)
      EXTERNAL simext
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER i, pdxinit, pdgrid, pdudat, pdxopt, pdgrad, 
     &        pdobs, pdnls, piudat, pinls
      
      INTEGER mxopt, nfct, spadim, liudat, ldudat
      DOUBLE PRECISION epsbfg, eps, funct, step, tmp 
      PARAMETER ( eps = 1.E-15 )    
c
c-----------------------------------------------------------------------
c
c     initializations
      info   = 0
      mxopt  = 3
      nfct   = 3
      spadim = 1
      epsbfg = 1.E-8
      liudat = 0
      ldudat = 0
      funct = 0.
c
c     Double workspaces 
c	
      pdxinit = 1
c     needs mxopt 	
      pdobs = pdxinit + mxopt
c     needs nfct so mxopt + nfct 
      pdgrid = pdobs + nfct
c     needs nfct so mxopt + nfct + nfct*spadim
      pdudat = pdgrid + nfct*spadim
c     needs 1 so mxopt + nfct + nfct*spadim + 1
      pdxopt = pdudat +1
c     needs mxopt so 2*mxopt + nfct + nfct*spadim + 1
      pdgrad = pdxopt + mxopt
c     needs mxopt so 3*mxopt + nfct + nfct*spadim + 1
      pdnls = pdgrad + mxopt 
c     needs ( mxopt*(mxopt+17)/2 + nfct*(mxopt+3+spadim) + ldudat ) so 
c     so mxopt*(mxopt+23)/2 + nfct*(mxopt+4+2*spadim) + 1
c  
c     Integer workspaces 
      piudat = 1
c     needs 1 so 1 
      pinls = piudat + 1
c     needs 2*mxopt + 6 + liudat  
c     so 2*mxopt + 6 + 1 = 2* mxopt + 7 
c
      dwork(pdxinit)     = 0.1
      dwork(pdxinit + 1) = 0.1
      dwork(pdxinit + 2) = 0.1
c
c     volatilities    
      dwork(pdgrid)   = retA
      dwork(pdgrid+1) = retB
      dwork(pdgrid+2) = retC
c
c     performances (returns)
      dwork(pdobs)     = volA
      dwork(pdobs + 1) = volB
      dwork(pdobs + 2) = volC
c
c     non-linear least-square (NLS)
      CALL nls ( simext, mxopt, nfct, dwork(pdxinit), dwork(pdobs),
     &           spadim, dwork(pdgrid), epsbfg, liudat, iwork(piudat),
     &           ldudat, dwork(pdudat), iwork(pinls), dwork(pdnls),
     &           dwork(pdxopt), funct, dwork(pdgrad), info)
      IF (info .LT. 0) RETURN
c
      step = (maxExpectedReturn - minExpectedReturn) / (nbValues - 1)
c
      DO i = 1,nbValues
        mveRet(i)= minExpectedReturn + i * step
      ENDDO
      DO i = 1,nbValues
        tmp = mveRet(i)
        mveVol(i) = dwork(pdxopt) + dwork(pdxopt+1)*tmp +
     &  dwork(pdxopt+2)*tmp*tmp
      ENDDO
      DO i = 1,nbValues
        IF (mveVol(i) .LT. 0.) mveVol(i) = 0.
      ENDDO
c
      RETURN
      END
