c=======================================================================
c
c     subroutine LOKERN                                      
c
c     GENERAL SUBROUTINE FOR KERNEL SMOOTHING:
c     COMPUTATION OF ITERATIVE PLUG-IN ALGORITHM FOR LOCAL BANDWIDTH
c     SELECTION FOR KERNELS WITH (NUE,KORD) = (0,2),(0,4),(1,3) OR (2,4).
c
c     THE RAW DATA SHOULD BE GIVEN BY THE POINTS (T(1),X(1)),...,(T(N),X(N))
c
c     THE RESULTING ESTIMATOR OF THE NUE-TH DERIVATIVE OF THE
c     REGRESSION CURVE IS GIVEN THROUGH THE POINTS
c     (TT(1),Y(1)),...,(TT(M),Y(M))
c
c     THE LOCAL PLUG-IN BANDWIDTH ARRAY IS GIVEN BY BAN(1),...,BAN(M)
c
c-----------------------------------------------------------------------
       SUBROUTINE LOKERN(T, X, N, TT, M, IHOM, NUE, KORD, IRND,
     .                   ISMO, M1, TL, TU, S, SIG, WN, W1, WM, BAN, Y)
c-----------------------------------------------------------------------
c
c  INPUT    T(N)         INPUT GRID (T(1)<T(2)<...<T(N))
c  INPUT    X(N)         DATA
c  INPUT    N            LENGTH OF X
c
C  INPUT    TT(M)        OUTPUT GRID
C  INPUT    M            LENGTH OF TT
C
C  INPUT    IHOM         HOMOSZEDASTICY OF VARIANCE
C                        0: HOMOSZEDASTIC ERROR VARIABLES
C                        <> 0: IF THE VARIANCE SHOULD ESTIMATED AS
C                              SMOOTH FUNCTION.
C                        ****** DEFAULT VALUE: IHOM=0
C
C  INPUT    NUE          ORDER OF DERIVATIVE (0-4)
C                        ****** DEFAULT VALUE: NUE=0
C
C  INPUT    KORD         ORDER OF KERNEL (<=6), FOR ISMO=0  ONLY
C                        NUE=0, KORD=2 OR KORD=4
C                        OR NUE=1, KORD=3 OR NUE=2, KORD=4 ARE ALLOWED
C                        ****** DEFAULT VALUE: KORD=NUE+2
C
C  INPUT    IRND         0: IF RANDOM GRID POINTS T MAY OCCUR
C                        <>0 ELSE  (ONLY NECESSARY IF S SHOULD BE
C                            COMPUTED)
C                        ****** DEFAULT VALUE IRND=0
C
C  INPUT    ISMO         0:ESTIMATING THE OPTIMAL LOCAL BANDWIDTH
C                        <>0 USING LOCAL INPUT BANDWIDTH-ARRAY IN BAN
C                        ****** DEFAULT VALUE ISMO=0
C
C  INPUT    M1           >=10, LENGTH OF W1, LARGE VALUES WILL INCREASE
C                        THE ACCURACY OF THE INTEGRAL APPROXIMATION
C                        ****** DEFAULT VALUE: M1=400
c
c IN/OUTPUT TL/TU        LOWER/UPPER BOUND FOR INTEGRAL APPROXIMATION
c                        AND VARIANCE ESTIMATION (IF SIG=0 AND IHOM=0),
c                        IF TU<=TL, [TL,TU] ARE COMPUTED AS ABOUT
C                        THE 87% MIDDLE PART OF [T(1),T(N)]
C                        ****** DEFAULT VALUES: TL=1.0, TU=0.0
C
C IN/OUTPUT S(0:N)       IF S(N)<=S(0) THIS ARRAY IS COMPUTED AS
C                        MIDPOINTS OF T, FOR NON-RANDOM DESIGN AND AS
C                        SMOOTHED QUANTILES FOR RANDOM DESIGN
C                        ****** DEFAULT VALUES: S(0)=1.0, S(N)=0.0
C                               AND THE OTHER S(I) CAN BE UNDEFINED
C
C IN/OUTPUT SIG          RESIDUAL VARIANCE, ESTIMATED FOR SIG=0 OR
C                        IHOM<>0, ELSE GIVEN BY INPUT
C                        ****** DEFAULT VALUE: SIG=-1.0
C
C IN/OUTPUT BAN(M)       LOCAL PLUG-IN BANDWIDTH ARRAY
C                        ****** WILL BE SET IN SUBROUTINE FOR ISMO=0
C
C WORK     WN(0:N,5)     WORK ARRAY FOR KERNEL SMOOTHING ROUTINE   
C                        ****** WILL BE SET IN SUBROUTINE
C WORK     W1(M1,3)      WORK ARRAY FOR INTEGRAL APPROXIMATION
C                        ****** WILL BE SET IN SUBROUTINE
C WORK     WM(M)         WORK ARRAY FOR INTEGRAL APPROXIMATION
c                        ****** WILL BE SET IN SUBROUTINE
c
c OUTPUT   Y(M)          KERNEL ESTIMATE WITH BOP(=B0 FOR ISMO<>0)
c                        ****** WILL BE SET IN SUBROUTINE
c-----------------------------------------------------------------------
c  USED SUBROUTINES: COFF, RESEST, KERNEL WITH FURTHER SUBROUTINES
c                     WHICH ARE CONTAINED IN THE FILE subs.f
c-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N),T(N),TT(M),Y(M),WN(0:N,5),S(0:N)
      DIMENSION W1(M1,3),BAN(M),WM(M)
      DIMENSION BIAS(2,0:2),VARK(2,0:2),FAK2(2:4)
C-
C-------- 1. INITIALISATIONS AND SOME ERROR-CHECKS
      DATA BIAS/.2,.04762,.4286,.1515,1.33,.6293/
      DATA VARK/.6,1.250,2.143,11.93,35.0,381.6/
      DATA FAK2/4.,36.,576./
      NYG=0
      INPUTS=0
C-------- IF NO ERRORS SHOULD BE WRITTEN ON STANDARD OUTPUT, SET IPRINT > 0
C-------- IF ERRORS AND VERY DETAILED WARNINGS SHOULD BE WRITTEN ON 
C--------           STANDARD OUTPUT, SET IPRINT < 0
      IPRINT=0
C
      IF(NUE.GT.4.OR.NUE.LT.0) THEN
        IF(IPRINT.EQ.0) 
     .     PRINT *,'lokern: Order of derivative not allowed'
        STOP
      END IF
      IF(NUE.GT.2.AND.ISMO.EQ.0) THEN
        IF(IPRINT.EQ.0) 
     .     PRINT *,'lokern: Order of derivative not allowed'
        STOP
      END IF
      IF(N.LE.2) THEN
        IF(IPRINT.EQ.0) 
     .     PRINT *,'lokern: Number of data too small'
        STOP
      END IF
      IF(M.LT.1) THEN
        IF(IPRINT.EQ.0) 
     .    PRINT *,'lokern: No output points'
        STOP
      END IF
      IF(M1.LT.10) THEN
        IF(IPRINT.EQ.0) 
     .     PRINT *,'lokern: Variable M1 is choosen too small'
        STOP
      END IF
      KK=(KORD-NUE)/2
      IF(2*KK+NUE.NE.KORD) THEN
        IF(IPRINT.EQ.0) 
     .    PRINT *,'lokern: Kernel order not allowed, set to ',NUE+2
        KORD=NUE+2
      END IF
      IF(KORD.GT.4.AND.ISMO.EQ.0) THEN
        IF(IPRINT.EQ.0) 
     .    PRINT *,'lokern: Kernel order not allowed, set to ',NUE+2
        KORD=NUE+2
      END IF
      IF(KORD.GT.6.OR.KORD.LE.NUE) THEN
        IF(IPRINT.EQ.0) 
     .    PRINT *,'lokern: Kernel order not allowed, set to ',NUE+2
        KORD=NUE+2
      END IF
      RVAR=SIG
C-
C-------- 2. COMPUTATION OF S-SEQUENCE 
      S0=1.5*T(1)-0.5*T(2)
      SN=1.5*T(N)-0.5*T(N-1)
      IF(S(N).LE.S(0)) THEN 
        INPUTS=1
        DO 20 I=1,N-1
20        S(I)=.5*(T(I)+T(I+1))
        S(0)=S0
        S(N)=SN
        IF(ISMO.NE.0.AND.IRND.NE.0) GOTO 230
      ELSE
        IF(ISMO.NE.0) GOTO 230
      END IF
C-
C-------- 3. COMPUTATION OF MINIMAL, MAXIMAL ALLOWED GLOBAL BANDWIDTH
      BMAX=(SN-S0)*.5
      BMIN=(SN-S0)/DBLE(N)*DBLE(KORD-1)*.6
C-
C-------- 4. WARNINGS IF TT-GRID LARGER THAN T-GRID
      IF(TT(1).LT.S0.AND.TT(M).GT.SN.AND.IPRINT.LT.0) PRINT *,
     .  'lokern: Extrapolation at both boundaries not optimized'      
      IF(TT(1).LT.S0.AND.TT(M).LE.SN.AND.IPRINT.LT.0) PRINT *,
     .   'lokern: Extrapolation at left boundary not optimized'      
      IF(TT(1).GE.S0.AND.TT(M).GT.SN.AND.IPRINT.LT.0) PRINT *,
     .   'lokern: Extrapolation at right boundary not optimized'      
C-
C-------- 5. COMPUTE TL,TU AND THEIR T-GRID AS INNER PART FOR
C            INTEGRAL APPROXIMATION IN THE ITERATIONS
      ITT=0
51    IF (TU.LE.TL) THEN
        TL=.933*S0+.067*SN
        TU=.067*S0+.933*SN
      ITT=ITT+1
      END IF
      TL=MAX(S0,TL)
      TU=MIN(SN,TU)
      IL=1
      IU=N
      WN(1,1)=0.0
      WN(N,1)=0.0
      DO 50 I=1,N
        IF(T(I).LE.TL.OR.T(I).GE.TU) WN(I,1)=0.0
        IF(T(I).GT.TL.AND.T(I).LT.TU) WN(I,1)=1.0
        IF(T(I).LT.TL) IL=I+1
50      IF(T(I).LE.TU) IU=I
      NN=IU-IL+1
      IF(NN.EQ.0.AND.ITT.EQ.0) THEN
        TU=TL-1.0
        GOTO 51
      END IF
      IF(NN.EQ.0.AND.ITT.EQ.1) THEN
        TU=SN
        TL=S0
        GOTO 51
      END IF
C-
C-------- 6. COMPUTE T-GRID FOR INTEGRAL APPROXIMATION 
      DO 60 I=1,M1
         W1(I,2)=1.0
60       W1(I,1)=TL+(TU-TL)*DBLE(I-1)/DBLE(M1-1)
C-
C-------- 7. CALCULATION OF WEIGHT FUNCTION
      ALPHA=1.D0/DBLE(13)
      DO 70 I=IL,IU
        XI=(T(I) - TL)/ALPHA/(TU-TL)
        IF(XI.GT.1) GOTO 71
70      WN(I,1)=(10.0-15*XI+6*XI*XI)*XI*XI*XI
71    DO 72 I=IU,IL,-1
        XI=(TU-T(I))/ALPHA/(TU-TL)
        IF(XI.GT.1) GOTO 73
72      WN(I,1)=(10.0-15*XI+6*XI*XI)*XI*XI*XI
73    DO 74 I=1,M1
        XI=(W1(I,1)-TL)/ALPHA/(TU-TL)
        IF(XI.GT.1) GOTO 75
74      W1(I,2)=(10.0-15*XI+6*XI*XI)*XI*XI*XI
75    DO 76 I=M1,1,-1
        XI=(TU-W1(I,1))/ALPHA/(TU-TL)
        IF(XI.GT.1) GOTO 77
76      W1(I,2)=(10.0-15*XI+6*XI*XI)*XI*XI*XI
77    CONTINUE 
C-
C-------- 8. COMPUTE CONSTANTS FOR ITERATION
      EX=1./DBLE(KORD+KORD+1)
      KK2=(KORD-NUE)
      KK=KK2/2
C-
C-------- 9. ESTIMATING VARIANCE AND SMOOTHED PSEUDORESIDUALS
      IF((SIG.LE..0).AND.(IHOM.EQ.0)) 
     .     CALL RESEST(T(IL),X(IL),NN,WN(IL,2),R2,SIG)
      IF(IHOM.NE.0) THEN 
        CALL RESEST(T,X,N,WN(1,2),SNR,SIG)
        BRES=MAX(BMIN,.2*NN**(-.2)*(S(IU)-S(IL-1)))
        DO 91 I=1,N
          WN(I,3)=T(I)
91        WN(I,2)=WN(I,2)*WN(I,2)
        CALL KERNEL(T,WN(1,2),N,BRES,0,KK2,NYG,S,
     .            WN(1,3),N,WN(1,4))
      ELSE
        CALL COFF(WN(1,4),N,SIG)
      END IF
C-
C-------- 10. ESTIMATE/COMPUTE INTEGRAL CONSTANT
100   VI=0.
      DO 101 I=IL,IU
101      VI=VI+WN(I,1)*N*(S(I)-S(I-1))**2*WN(I,4)       
C- 
C-------- 11. REFINEMENT OF S-SEQUENCE FOR RANDOM DESIGN
      IF(INPUTS.EQ.1.AND.IRND.EQ.0) THEN
        DO 110 I=0,N
          WN(I,5)=DBLE(I)/DBLE(N+1)
          WN(I,2)=(DBLE(I)+.5)/DBLE(N+1)
110       WN(I,3)=WN(I,2)
        EXS=-DBLE(3*KORD+1)/DBLE(6*KORD+3)
        EXSVI=DBLE(KORD)/DBLE(6*KORD+3)
        BS=0.1*(VI/(SN-S0)**2)**EXSVI*N**EXS
        CALL KERNEL(WN(1,5),T,N,BS,0,2,NYG,WN(0,3),WN(0,2),N+1,S(0))
111     ISORT=0
        VI=0.0
        DO 112 I=1,N
        VI=VI+WN(I,1)*N*(S(I)-S(I-1))**2*WN(I,4)       
          IF(S(I).LT.S(I-1)) THEN
            SSI=S(I-1)
            S(I-1)=S(I)
            S(I)=SSI
            ISORT=1
          END IF
112     CONTINUE
        IF(ISORT.EQ.1) GOTO 111
        IF(ISMO.NE.0) GOTO 230
      END IF
      B=BMIN*2.
C-
C-------- 12. COMPUTE INFLATION CONSTANT AND EXPONENT AND LOOP OF ITERATIONS
      CONST=DBLE(2*NUE+1)*FAK2(KORD)*VARK(KK,NUE)*VI
     .       /(DBLE(2*KORD-2*NUE)*BIAS(KK,NUE)**2*DBLE(N))
      FAC=1.1*(1.+(NUE/10.)+0.05*(KORD-NUE-2.))
     .       *N**(2./DBLE((2*KORD+1)*(2*KORD+3)))
      ITENDE=1+2*KORD+KORD*(2*KORD+1)

      DO 120 IT=1,ITENDE
C-
C-------- 13. ESTIMATE DERIVATIVE OF ORDER KORD IN ITERATIONS
        B2=B*FAC
        B2=MAX(B2,BMIN/DBLE(KORD-1)*DBLE(KORD+1))
        B2=MIN(B2,BMAX)
        CALL KERNEL(T,X,N,B2,KORD,KORD+2,NYG,S,W1(1,1),M1,W1(1,3))
C-
C-------- 14. ESTIMATE INTEGRALFUNCTIONAL IN ITERATIONS
        XMY2=.75*(W1(1,2)*W1(1,3)*W1(1,3)+W1(M1,2)*W1(M1,3)*W1(M1,3))
        DO 140 I=2,M1-1
140       XMY2=XMY2+W1(I,2)*W1(I,3)*W1(I,3)
        XMY2=XMY2*(TU-TL)/DBLE(M1)
C-
C-------- 15. FINISH OF ITERATIONS
        B=(CONST/XMY2)**EX
        B=MAX(BMIN,B)
        B=MIN(BMAX,B)

120     CONTINUE
C-------- 16  COMPUTE SMOOTHED FUNCTION WITH GLOBAL PLUG-IN BANDWIDTH
160   CALL KERNEL(T,X,N,B,NUE,KORD,NYG,S,TT,M,Y)
C-
C-------- 17. VARIANCE CHECK 
      IF(IHOM.NE.0) SIG=RVAR
      IF(RVAR.EQ.SIG.OR.R2.LT..88.OR.IHOM.NE.0.OR.NUE.GT.0) GOTO 180
      II=0
      IIL=0
      J=2
      TLL=MAX(TL,TT(1))
      TUU=MIN(TU,TT(M))
      DO 170 I=IL,IU
         IF(T(I).LT.TLL.OR.T(I).GT.TUU) GOTO 170
         II=II+1
         IF(IIL.EQ.0) IIL=I
171       IF(TT(J).LT.T(I)) THEN
            J=J+1
            IF(J.LE.M) GOTO 171
            END IF
       WN(II,3)=X(I)-Y(J)+(Y(J)-Y(J-1))*(TT(J)-T(I))/(TT(J)-TT(J-1))
170    CONTINUE
      CALL RESEST(T(IIL),WN(1,3),II,WN(1,4),SNR,RVAR)
      Q=SIG/RVAR
      CALL COFF(WN(1,4),N,SIG) 
      IF(Q.LE.2.) GOTO 180
      IF(Q.GT.5..AND.R2.GT..95) RVAR=RVAR*.5
      SIG=RVAR
      CALL COFF(WN(1,4),N,SIG)
      GOTO 100
C-
C-------- 18. LOCAL INITIALIZATIONS
180   BVAR=B
      NUEV=0
      KORDV=2
      NYL=1
C-
C-------- 19. COMPUTE INNER BANDWIDTHS
        G1=0.86*(1.+DBLE(KORD-NUE-2)*.05)*B
        G1=G1*DBLE(N)**(4./DBLE(2*KORD+1)/(2*KORD+5))
        G1=MAX(G1,BMIN/DBLE(KORD-1)*DBLE(KORD+1))
        G1=MIN(G1,BMAX)

        G2=1.4*(1.+DBLE(KORD-NUE-2)*0.05)*B
        G2=G2*DBLE(N)**(2./DBLE(2*KORD+1)/DBLE(2*KORD+3))
        G2=MAX(G2,BMIN)
        G2=MIN(G2,BMAX)
C-
C-------- 20. ESTIMATE/COMPUTE INTEGRAL CONSTANT VI LOCALLY
        DO 200 I=1,N
200       WN(I,4)=DBLE(N)*WN(I,4)*(S(I)-S(I-1))
        DO 201 J=1,M
          BAN(J)=BVAR
             WM(j)=tt(j)
           IF(TT(J).LT.S(0)+G1) THEN
             DIST=((TT(J)-G1-S(0))/G1)**2
             BAN(J)=BVAR*(1.0+1.0*DIST)
             BAN(J)=MIN(BAN(J),BMAX)
             WM(J)=TT(J)+.5*DIST*G1
           ELSE IF(TT(J).GT.S(N)-G1) THEN
             DIST=((TT(J)-S(N)+G1)/G1)**2
             BAN(J)=BVAR*(1.0+1.0*DIST)
             BAN(J)=MIN(BAN(J),BMAX)
             WM(J)=TT(J)-.5*DIST*G1
           END IF
201     CONTINUE
        CALL KERNEL(T,WN(1,4),N,BVAR,NUEV,KORDV,NYL,S,WM,M,BAN)
C-
C-------- 21. ESTIMATION OF KORDTH DERIVATIVE LOCALLY
        WSTEP=(TT(M)-TT(1))/DBLE(M1-2)
        DO 210 J=2,M1
           W1(J,2)=TT(1)+DBLE(J-2)*WSTEP
210     W1(J,1)=TT(1)+DBLE(J-1.5)*WSTEP
        W1(1,1)=TT(1)+.5*WSTEP

        CALL KERNEL(T,X,N,G1,KORD,KORD+2,NYG,S,W1(2,2),M1-1,W1(2,3))

        DO 211 J=2,M1
211        W1(J,3)=W1(J,3)*W1(J,3)
        DO 212 J=1,M
           Y(J)=G2
           IF(TT(J).LT.S(0)+G1) THEN           
             Y(J)=G2*(1.0+1.0*((TT(J)-G1-S(0))/G1)**2)
             Y(J)=MIN(Y(J),BMAX)
           ELSE IF(TT(J).GT.S(N)-G1) THEN
             Y(J)=G2*(1.0+1.0*((TT(J)-S(N)+G1)/G1)**2)
             Y(J)=MIN(Y(J),BMAX)
           END IF
212      CONTINUE

       CALL KERNP(W1(2,2),W1(2,3),M1-1,G2,NUEV,KORDV,NYL,
     .            W1(1,1),WM,M,Y)
C-
C-------- 22. FINISH
       DO 220 J=1,M
          XH=BMIN**(2*KORD+1)*ABS(Y(J))*VI/CONST
          XXH=CONST*ABS(BAN(J))/VI/BMAX**(2*KORD+1)
          IF(BAN(J).LT.XH) THEN
             BAN(J)=BMIN
          ELSE
             IF(Y(J).LT.XXH) THEN
               BAN(J)=BMAX
             ELSE
                BAN(J)=(CONST*BAN(J)/Y(J)/VI)**EX
             END IF
          END IF
220    CONTINUE
C-
C-------- 23. COMPUTE SMOOTHED FUNCTION WITH LOCAL PLUG-IN BANDWIDTH
230     DO 231 J=1,M
231        Y(J)=BAN(J)
        CALL KERNEL(T,X,N,B,NUE,KORD,NYL,S,TT,M,Y)

        RETURN
        END
