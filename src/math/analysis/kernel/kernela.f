c=======================================================================
c
c     subroutine KERNELA                                      
c
c     Density Estimation by Non Parametric Gaussian Kernel method 
c     with iterative plug-in bandwidth selection. This version only uses 
c     the gauss kernel and estimates only the density itself and not its 
c     derivatives.
c
c-----------------------------------------------------------------------
      SUBROUTINE kernela ( N, M, X, Z, F, H, info)
c-----------------------------------------------------------------------
c
c     INPUT 
c            N      : number of histogram points (n>1)           integer
c            M      : number of gird points      (m>0)           integer
c            X(N)   : sorted returns/histogram (n)                double
c            Z(M)   : sorted output gird (m)                      double
c
c     OUTPUT 
c            F(M)   : estimated density (m)                       double
c            H      : estimated iterative plugin bandwitdth       double
c            info   : diagnostic argument                        integer                                                                      
c                                                                       
c-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER info
      DIMENSION X(N),Z(M),F(M)
c
c
      info = 0      
      XN=DBLE(N)
      XIQR=X(NINT(.75*XN))-X(NINT(.25*XN)+1)
      PI=3.141592645D0
      RT2PI = DSQRT(2.D0*PI)
      ITER=5
C-------  Estimate inflation constant C
          H2=(.920*XIQR)/(XN**(1.D0/7.D0))
          H3=(.912*XIQR)/(XN**(1.D0/9.D0))
          S2 = 0.D0
          S3 = 0.D0
          DO 20 I = 1, N-1
           DO 30 J = I+1, N
                     D2 = ((X(I) - X(J))/H2)**2
                     D3 = ((X(I) - X(J))/H3)**2
                     IF(D2.GT.50.AND.D3.GT.60) GOTO 20
                     E2 = DEXP(-0.5D0*D2)
                     E3 = DEXP(-0.5D0*D3)
                     S2 = S2 + (D2**2 - 6.D0*D2 + 3.D0)*E2
               S3 = S3 + (D3**3 - 15.D0*D3**2 + 45.D0*D3 - 15.D0)*E3
30          CONTINUE
20         CONTINUE
          RHAT2 = (2.D0*S2)/((XN**2)*(H2**5)*RT2PI)
          RHAT2 = RHAT2 + 3.D0/(RT2PI*XN*(H2**5))
          RHAT3 = (-2.D0*S3)/((XN**2)*(H3**7)*RT2PI)
          RHAT3 = RHAT3 + 15.D0/(RT2PI*XN*(H3**7))
          CO1 = 1.357D0*(RHAT2/RHAT3)**(1.D0/7.D0)
C-
C-------  Compute constant of asymptotic formula
C-
       CONST=1.D0/(2.D0*DSQRT(PI))
       A=1.132795764/RHAT3**(1.D0/7.D0)*XN**(-1.D0/2.D0)
C-
C------  Loop over iterations
C-
       DO 100 IT=1,ITER
C-                                                                     
C-------  Estimate functional
C
        S=0.D0                                                    
        DO 40 I = 1, N-1
           DO 50 J = I+1, N
                     D2 = ((X(I) - X(J))/A)**2
                     IF(D2.GT.50) GOTO 40
                     E2 = DEXP(-0.5D0*D2)
                     S = S + (D2**2 - 6.D0*D2 + 3.D0)*E2
50          CONTINUE
40       CONTINUE
        R2 = (2.D0*S)/((XN**2)*(A**5)*RT2PI)
        R2 = R2 + 3.D0/(RT2PI*XN*(A**5))
C-                                                                     
C------- Estimate bandwidth by asymptotic formula
C-
         H=(CONST/(R2*XN))**(.2D0)                                    
         A=CO1*H**(5./7.)
100     CONTINUE
C-
C------- Estimate density with plugin bandwidth
C-             
      JBEGIN=1
      JEND=1                                             
      DO 200 I=1,M                                                      
         S=0.D0    
         DO 210 J=JBEGIN,JEND
                T=(Z(I)-X(J))/H
                IF(T.GT.5.0.AND.JBEGIN.LT.N) THEN
                   JBEGIN=JBEGIN+1
                   GOTO 210
                END IF
                S=S+DEXP(-T*T/2.)
210             CONTINUE         
         DO 220 JEND=J,N
                T=(Z(I)-X(JEND))/H
                IF(T.LT.-5.0) GOTO 230
                S=S+DEXP(-T*T/2.)
220             CONTINUE
230         F(I)=S/(XN*H*RT2PI)
            JEND=JEND-1
C-                         
200      CONTINUE                                             
      RETURN                                                            
      END                                                               
