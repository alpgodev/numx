c=======================================================================
c
c     subroutine NPREG                                      
c
c     Regression Estimation by Non Parametric Kernel method 
c       
c     THE RAW DATA SHOULD BE GIVEN BY THE POINTS 
c     (T(1),X(1)),...,(T(N),X(N))
c
c     THE RESULTING ESTIMATOR OF THE NUE-TH DERIVATIVE OF THE
c     REGRESSION CURVE IS GIVEN THROUGH THE POINTS
c     (TT(1),Y(1)),...,(TT(M),Y(M))
c
c-----------------------------------------------------------------------
      SUBROUTINE npreg ( N, M, NUE, T, X, TT, dwork, Y, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       N     : number of points (N > 1)                         integer
c       M     : number of (output) gird points (M > 0)           integer
c       NUE   : derivative order of the estimated regression,    integer    
c               only NUE=0,1 or 2 are possible. 
c       T(N)  : input gird (N)                                    double 
c       X(N)  : data (N)                                          double 
c       TT(M) : output gird (M)                                   double           
c
c     WORKSPACE 
c       dwork : 6*(N + 201)                                       double
c
c     OUTPUT 
c       Y(M)  : estimates (M)                                     double
c       info  : diagnostic argument                              integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o arguments      
      INTEGER NUE, N, M, info
      DOUBLE PRECISION T(*), X(*), TT(*), Y(*)
c
c     workspaces
      DOUBLE PRECISION dwork(*)
c
c     local variables      
      INTEGER pds, pdwn, pdw1
      INTEGER IHOM, IRND, ISMO, M1, KORD
      PARAMETER(IHOM = 0, IRND = 1, ISMO = 0, M1 = 400)
c
c     IHOM, SIG: These variables indicate which variance estimator should be used for IRND=0 
c           or ISMO=0 (default values). They are not used at all for both IRND<>0 and ISMO<>0.
c           For IHOM=0 (default) homoscedastic variance is assumed, i.e. that the variance of 
c           the error variables equals a constant and thus does not depend on the location or 
c           on the regression function. For IHOM=0 and SIG<=0 (both default values), 
c           the variance is estimated nonparametrically and can be found as output of SIG 
c           as long as IRND=0 or ISMO=0 (default values). For IHOM=0 and SIG>0 the algorithm 
c           uses the input value of SIG as variance estimator for IRND=0 or ISMO=0 (default values).
c           For IHOM<>0 a smooth variance function is assumed and estimated by the algorithm 
c           for IRND=0 or ISMO=0 (default values). The variable SIG is not used than.
c     IRND: This variable will specify the weights of the convolution kernel estimator 
c           and the calculation of the S-array. The classical Gasser-Müller form of the weights 
c           is used for IRND<>0. This will be useful for equidistant and nearly equidistant 
c           design points T(1),...,T(N). If the distances T(I)-T(I-1) are very variable, 
c           e.g. if T(1),...,T(N) are ordered variables from a random sample, such a convolution 
c           kernel estimator is inefficient. Than it is useful to perform an easy modification 
c           which is done for IRND=0 (default). For large sample size N this modification can 
c           also handle multiple design points automatically.
c     ISMO  This variable simply indicates if previously specified local or global 
c           bandwidths should be used. Than put ISMO<>0 and specify a global bandwidth 
c           by the variable B in glkern.f and a local bandwidth array BAN of length M in lokern.f.
c           Data-adaptive plugin bandwidths are calculated for ISMO=0 (default).
c     M1    >=10, length of W1, large value will increase the accuracy of the integral 
c           approximation (default M1=400).
c     KORD  This variable specifies the kernel order which is used by the kernel regression 
c           estimator. KORD-NUE has to be an even number, KORD=NUE+2 will be the default choice.
c           Only NUE+2<=KORD<=4 is allowed for automatic bandwidth choice (ISMO=0, default) and 
c           NUE+2<=KORD<=6 is allowed for smoothing with specified bandwidths (ISMO<>0).
c
      DOUBLE PRECISION TL, TU, SIG, B
      PARAMETER (TL = 1.0, TU = 0.0, SIG = -1.0)
c
c     TL, TU These variables specify lower and upper bounds of the region for the global bandwidth
c            estimation step (even in lokern.f). If TL>=TU (default) they are set automatically 
c            and exclude a small part of the boundary region.
c     S(0:N) This array is of methodological interest mainly. It specifies the weights of the 
c            convolution kernel estimator. For S(N)>=S(0) (default) this array is computed 
c            automatically according to the IRND-variable. For IRND=0 (default) the computation 
c            of this array does not depend on the bandwidths, the same array may be used if the 
c            subroutines are called for different bandwidths (ISMO<>0). It is very important to 
c            perform a new initialization for different design points and useful to do so for 
c            different NUE, KORD and variance too.
c            For IRND<>0 the S-array does only depend on T(1),...,T(N).
c     SIG    Residual variance, estimated for SIG=0 or IHOM<>0, else given by input 
c            (default SIG = -1.0).
c     B      Global plug-in bandwith (B can be undefined if ISMO=0).
c
c--------------------------------------------------------------------
c
c     initialization
      info  = 0
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pds = 1
c     pds  : pointer for S array, so ( N + 1 )
      pdwn = pds + ( N + 1 )
c     pdwn : pointer for WN array, so ( 5*(N + 1) )
      pdw1 = pdwn + ( 5*(N + 1) )
c     pdw1 : pointer for W1 array, so (3*M1)
c
c     Total size of dwork array = N + 1 
c                               + 5*(N + 1)
c                               + 3*400
c                               = 6*(N + 201)
c
c--------------------------------------------------------------------
c      
c     kernel order
      KORD = NUE + 2
c
c      DOUBLE PRECISION S(0:N), WN(0:N,5),W1(M1,3)
c
c     S(0)=1.0, S(N)=0.0
      dwork(pds) = 1.0
      dwork(pds + N + 1) = 0.0
c
c     call non parametric regression estimation      
      CALL GLKERN(T, X, N, TT, M, IHOM, NUE, KORD, IRND, ISMO, M1, TL, 
     &            TU, dwork(pds), SIG, dwork(pdwn), dwork(pdw1), B, Y)
c
      RETURN
      END
