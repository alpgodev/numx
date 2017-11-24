c=======================================================================
c
c     SGAMMA
c 
c     Standard Gamma distribution generator G(a)
c
c-----------------------------------------------------------------------
c
c     Copyright (c) 2014 NumX
c     All rights reserved.
c 
c     This software is the confidential and proprietary information
c     of NumX. You shall not disclose such Confidential
C     Information and shall use it only in accordance with the terms
c     of the licence agreement you entered into with NumX.
c
c     author: Yann Vernaz
c
c----------------------------------------------------------------------- 
      DOUBLE PRECISION FUNCTION sgamma(a)
c-----------------------------------------------------------------------      
c
c     INPUT :                                                          
c      a : Gamma parameter (>= 1.0)                               double
c
c     Method - for details see:                                                                                                                
c                                                                      
c     AHRENS, J.H. AND DIETER, U.                            
c     GENERATING GAMMA VARIATES BY A                         
c     MODIFIED REJECTION TECHNIQUE.                          
c     COMM. ACM, 25,1 (JAN. 1982), 47 - 54.                  
c                                                                      
c     STEP NUMBERS CORRESPOND TO ALGORITHM 'GD' IN THE ABOVE PAPER     
c                                 (STRAIGHTFORWARD IMPLEMENTATION)     
c                                                                      
c     PARAMETER  0.0 < A < 1.0  !                            
c                                                                      
c     AHRENS, J.H. AND DIETER, U.                            
c     COMPUTER METHODS FOR SAMPLING FROM GAMMA,              
c     BETA, POISSON AND BINOMIAL DISTRIBUTIONS.              
c     COMPUTING, 12 (1974), 223 - 246.                       
c                                                                      
c     (ADAPTED IMPLEMENTATION OF ALGORITHM 'GS' IN THE ABOVE PAPER)
c
c     COEFFICIENTS Q(K) - FOR Q0 = SUM(Q(K)*A**(-K))
c     COEFFICIENTS A(K) - FOR Q = Q0+(T*T/2)*SUM(A(K)*V**K)
c     COEFFICIENTS E(K) - FOR EXP(Q)-1 = SUM(E(K)*Q**K)
c
c-----------------------------------------------------------------------
c
c     scalar arguments
      DOUBLE PRECISION a
c     
c     local scalars
      DOUBLE PRECISION a1,a2,a3,a4,a5,a6,a7,aa,aaa,b,b0,c,d,e,e1,e2,e3,
     +     e4,e5,p,q,q0,
     +     q1,q2,q3,q4,q5,q6,q7,r,s,s2,si,sqrt32,t,u,v,w,x
c     
c     external functions
      DOUBLE PRECISION ranf, sexpo, snorm
      EXTERNAL ranf, sexpo, snorm
c     
c     intrinsic functions 
      INTRINSIC abs, log, exp, sign, sqrt
c     
c     save statement 
      SAVE aa,aaa,s2,s,d,q0,b,si,c,q1,q2,q3,q4,q5,q6,q7,a1,a2,a3,a4,a5,
     +     a6,a7,e1,e2,e3,e4,e5,sqrt32
c     
c     data statements
c     PREVIOUS A PRE-SET TO ZERO - AA IS A', AAA IS A"
c     SQRT32 IS THE SQUAREROOT OF 32 = 5.656854249492380
      DATA q1,q2,q3,q4,q5,q6,q7/.04166669,.02083148,.00801191,.00144121,
     +     -.00007388,.00024511,.00024240/
      DATA a1,a2,a3,a4,a5,a6,a7/.3333333,-.2500030,.2000062,-.1662921,
     +     .1423657,-.1367177,.1233795/
      DATA e1,e2,e3,e4,e5/1.,.4999897,.1668290,.0407753,.0102930/
      DATA aa/0.0/,aaa/0.0/,sqrt32/5.656854/
c
c     executable statements
      IF (a.EQ.aa) GO TO 10
      IF (a.LT.1.0) GO TO 130
c
c     STEP  1:  RECALCULATIONS OF S2,S,D IF A HAS CHANGED
      aa = a
      s2 = a - 0.5
      s = sqrt(s2)
      d = sqrt32 - 12.0*s
c
c     STEP  2:  T=STANDARD NORMAL DEVIATE,
c               X=(S,1/2)-NORMAL DEVIATE.
c               IMMEDIATE ACCEPTANCE (I)
   10 t = snorm()
      x = s + 0.5*t
      sgamma = x*x
      IF (t.GE.0.0) RETURN
c
c     STEP  3:  U= 0,1 -UNIFORM SAMPLE. SQUEEZE ACCEPTANCE (S)
      u = ranf()
      IF (d*u.LE.t*t*t) RETURN
c
c     STEP  4:  RECALCULATIONS OF Q0,B,SI,C IF NECESSARY
      IF (a.EQ.aaa) GO TO 40
      aaa = a
      r = 1.0/a
      q0 = ((((((q7*r+q6)*r+q5)*r+q4)*r+q3)*r+q2)*r+q1)*r
c
c               APPROXIMATION DEPENDING ON SIZE OF PARAMETER A
c               THE CONSTANTS IN THE EXPRESSIONS FOR B, SI AND
c               C WERE ESTABLISHED BY NUMERICAL EXPERIMENTS
      IF (a.LE.3.686) GO TO 30
      IF (a.LE.13.022) GO TO 20
c
c               CASE 3:  A .GT. 13.022
      b = 1.77
      si = .75
      c = .1515/s
      GO TO 40
c
c               CASE 2:  3.686 .LT. A .LE. 13.022
   20 b = 1.654 + .0076*s2
      si = 1.68/s + .275
      c = .062/s + .024
      GO TO 40
c
c               CASE 1:  A .LE. 3.686
   30 b = .463 + s + .178*s2
      si = 1.235
      c = .195/s - .079 + .16*s
c
c     STEP  5:  NO QUOTIENT TEST IF X NOT POSITIVE
   40 IF (x.LE.0.0) GO TO 70
c
c     STEP  6:  CALCULATION OF V AND QUOTIENT Q
      v = t/ (s+s)
      IF (abs(v).LE.0.25) GO TO 50
      q = q0 - s*t + 0.25*t*t + (s2+s2)*log(1.0+v)
      GO TO 60
c
   50 q = q0 + 0.5*t*t* ((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v
c
c     STEP  7:  QUOTIENT ACCEPTANCE (Q)
   60 IF (log(1.0-u).LE.q) RETURN
c
c     STEP  8:  E=STANDARD EXPONENTIAL DEVIATE
c               U= 0,1 -UNIFORM DEVIATE
c               T=(B,SI)-DOUBLE EXPONENTIAL (LAPLACE) SAMPLE
   70 e = sexpo()
      u = ranf()
      u = u + u - 1.0
      t = b + sign(si*e,u)
c
c     STEP  9:  REJECTION IF T .LT. TAU(1) = -.71874483771719
   80 IF (t.LT. (-.7187449)) GO TO 70
c
c     STEP 10:  CALCULATION OF V AND QUOTIENT Q
      v = t/ (s+s)
      IF (abs(v).LE.0.25) GO TO 90
      q = q0 - s*t + 0.25*t*t + (s2+s2)*log(1.0+v)
      GO TO 100
c
   90 q = q0 + 0.5*t*t* ((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v
c
c     STEP 11:  HAT ACCEPTANCE (H) (IF Q NOT POSITIVE GO TO STEP 8)
  100 IF (q.LE.0.0) GO TO 70
      IF (q.LE.0.5) GO TO 110
c
c     through line 125 to handle large Q case
      IF (q.LT.15.0) GO TO 105
c
c     Here Q is large enough that Q = log(exp(Q) - 1.0) (for DOUBLE PRECISION Q)
c     reformulate test at 120 in terms of one EXP, if not too big
c     87.49823 is close to the largest DOUBLE PRECISION which can be
c     exponentiated (87.49823 = log(1.0E38))
      IF ((q+e-0.5*t*t).GT.87.49823) GO TO 125
      IF (c*abs(u).GT.exp(q+e-0.5*t*t)) GO TO 70
      GO TO 125
c
 105  w = exp(q) - 1.0
      GO TO 120
c
  110 w = ((((e5*q+e4)*q+e3)*q+e2)*q+e1)*q
c
c               IF T IS REJECTED, SAMPLE AGAIN AT STEP 8
  120 IF (c*abs(u).GT.w*exp(e-0.5*t*t)) GO TO 70
 125  x = s + 0.5*t
      sgamma = x*x
      RETURN
c
c     ALTERNATE METHOD FOR PARAMETERS A BELOW 1  (.3678794=EXP(-1.))
c
c     changed B to B0 (which was added to declarations for this)
c     in 130 to END to fix rare and subtle bug.
c     Line: '130 aa = 0.0' was removed (unnecessary, wasteful).
c     Reasons: the state of AA only serves to tell the A .GE. 1.0
c     case if certain A-dependant constants need to be recalculated.
c     The A .LT. 1.0 case (here) no longer changes any of these, and
c     the recalculation of B (which used to change with an
c     A .LT. 1.0 call) is governed by the state of AAA anyway.
 130  b0 = 1.0 + .3678794*a
  140 p = b0*ranf()
      IF (p.GE.1.0) GO TO 150
      sgamma = exp(log(p)/a)
      IF (sexpo().LT.sgamma) GO TO 140
      RETURN
c
  150 sgamma = -log((b0-p)/a)
      IF (sexpo().LT. (1.0-a)*log(sgamma)) GO TO 140
      RETURN
      END
