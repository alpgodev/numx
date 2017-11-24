c=======================================================================
c
c     subroutine QP
c
c     Quadratic solver
c
c     Minimization of function  f(x) = 1/2*x'*quamat*x + plin'*x
c
c            cinf(i) <= x(i) <= csup(i), i=1,nmat
c            (cmat(1,j),...,cmat(nmat,j))*x  = bvect(j) for j=1,nceg
c            (cmat(1,j),...,cmat(nmat,j))*x <= bvect(j) 
c                                               for j=nceg+1,nceg+ncineg
c
c-----------------------------------------------------------------------
      SUBROUTINE qp ( nmat, quamat, plin, nceg, ncineg, cmat, bvect,
     &                cinf, csup, iwork, dwork,
     &                lagr, wopt, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c            nmat   : size of the matrix                         integer
c            quamat : quadratic matrix (nmat*nmat)                double
c            plin   : linear part vector (nmat)                   double
c            nceg   : number of equality constraints             integer
c            ncineg : number of inequality constraints           integer
c            cmat   : matrix of constraints (nmat*(nceg+ncineg))  double
c            bvect  : vector of constraints (nceg+ncineg)         double
c            cinf   : inferior limit vector (nmat)                double
c            csup   : superior limit vector (nmat)                double
c
c     WORKSPACE 
c            iwork  : (3*nmat + 2*ncineg + nceg + 1)             integer
c            dwork  : (nmat*nmat + 6*nmat + 2*ncineg)             double 
c
c     OUTPUT 
c            lagr   : Lagrange multipliers (nmat+ncineg+nceg)     double
c            wopt   : optimal vector (nmat)                       double
c            info   : diagnostic argument                        integer
c
c-------------------------------------------------------------------------------
c     USED FUNCTION
c
c     From rne/optim/ directory: 
c     anfm01.f anfm03.f anfm05.f anrs01.f auxo01.f dimp03.f dnrm0.f optr03.f 
c     pasr03.f zthz.f anfm02.f anfm04.f anfm06.f anrs02.f desr03.f dipvtf.f 
c     optr01.f opvf03.f plcbas.f 
c
c     From BLAS library: 
c     daxpy.f dcopy.f ddot.f dnrm2.f dscal.f dswap.f idamax.f 
c 
c     From LAPACK library: 
c     dlamch.f   
c
c     From rne/calelm directory: 
c     add.f ddif.f dmmul.f 
c
c     PLCBAS Error Codes
c                   =   1 maximal number of iterations exceeded
c                   =  -1 no inferior bound supplied
c                   =  -2 degenerate point with indefinite cycle
c                   =  -3 may be there is a problem of bound. Too large
c                         gap between 2 successive iterations 
c                   =  -4 incorrect data 
c                   = -11 incompatibility of equality constraints
c                         detected when calling OPTR01.
c                   = -12 no admissible points OPTR01 with the given 
c                         constraints 
c                   = -13 degenerate point with indefinite cycle in
c                         OPTR01.
c                   = -14 maximal number of iterations for OPTR01
c                         reached without finding an admissible point
c                   = -15 cinf constraints no verified
c                   = -16 csup constraints no verified  
c 
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER nmat, nceg, ncineg, info, iter
      DOUBLE PRECISION quamat(nmat,*), wopt(*), plin(*), lagr(*),
     &                 cinf(*), csup(*), cmat(nmat,*), bvect(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER modo, imp, ira, iosort
      DOUBLE PRECISION fopt
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
c
c     modo = 1 : normal call
c     modo = 2 : as modo = 1,may be util when matrix is definite positive
c     modo = 3 : call with initial conditions (compatible with restrictions)
      modo = 1
c
c     parameters for output: no output here
      imp    = 6
      iosort = 6
c
c     ira = 0  : no bounds
c     ira = 1  : only lower bounds
c     ira = 2  : only upper bounds
c     ira = 3  : lower and upper bounds
      ira = 3
c
c     Number maximum of iterations
C     if ITER <= 0, number max=14*(nmat+nceg+ncineg)
      iter = 0
c
      CALL plcbas ( quamat, plin, cmat, bvect, cinf, csup, ira,
     &              nceg, ncineg, wopt, fopt, dwork, iwork, lagr, imp,
     &              iosort, nmat, modo, info, iter )
c------------------------------------------------------------------------
c    INPUTS
c
c     quamat    Matrix of dimension (N,N), the quadratic part. 
c               By default only the triangular inferior part is supplied.
c
c     plin      Vector  N-dimensional of the linear part  
c
c     cmat      Matrix of dimension (N,MI+MD). Contiene los coeficientes
c               de las restricciones  de igualdad y  desigualdad. Los de
c               igualdad se almacenan en las  MI primeras columnas de la
c               matriz  C  y los de desigualdad en las restantes.
c
c     bvect     Vector of dimension  MI+MD. Contiene los coeficientes de
c               los  terminos  independientes  de  las  restricciones de
c               igualdad y desigualdad.
c
c     cinf      Si  IRA= 0  o  2  esta  variable  no sera  utilizada. Si
c               IRA= 1  o  3,   CI  sera  un  vector  de  dimension  N
c               conteniendo las cotas inferiores de  X. Si  X(I) no esta
c               acotado  inferiormente, CI(I)  debera ser  menor  que la
c               raiz cuadrada negativa  de la constante real  mas grande
c               de la maquina.
c
c     csup      Si  IRA= 0  o  1  esta  variable no  sera  utilizada. Si
c               IRA= 2  o  3,  CS   sera  un  vector  de  dimension  N
c               conteniendo las cotas superiores de  X. Si  X(I) no esta
c               acotado  superiormente,  CS(I)  debera  ser  mayor  que
c               la raiz cuadrada de la constante real  mas grande  de la
c               maquina.
c
c     ira       Variable  que  indica  si  existen  restricciones  de
c               acotacion en el problema. Toma los valores:
c                  0  : No existen restricciones de acotacion.
c                  1  : Se  tienen  solo  restricciones  de  acotacion
c                       inferior.
c                  2  : Existen solo restricciones de acotacion superior
c                  3  : Se tienen ambos tipos de restriccion.
c
c     nceg      Numero de restricciones de igualdad del problema.
c
c     ncineg    Numero  de  restricciones  de desigualdad  del  problema
c
c        W      Vector de trabajo de dimension  N*N+6*N+2*MD
c
c        IV     Vector entero  de dimension  3*N+2*MD+MI+1.
c
c     IMP    Indicador del nivel de impresion de salida de resultados
c               Toma los valores:
c                  6  : Sin salida de resultados.
c                  7  : Escribe el motivo de finalizacion del proceso.
c                  8  : Se obtiene informacion referente  a la solucion:
c                       el numero de iteraciones; el  optimo  calculado,
c                       la norma del vector de Kuhn-Tucker, el numero de
c                       restricciones activas y cuales   son    (para
c                       identificar estas restricciones  se  sigue el
c                       siguiente convenio para los valores que se
c                       obtienen:
c
c                       [-N,-1]         : Restriciones de cota inferior
c                       [ 1,N ]         : Restricciones de cota superior
c                       [N+1,N+MI]      : Restricciones de igualdad
c                       [N+MI+1,N+MI+MD]: Restricciones de desigualdad)
c
c                       y los multiplicadores de  Lagrange  asociados  a
c                       ellas  (si IND=0);  el valor  del funcional  ;
c                       la norma del punto calculado (si  IND=-3).
c
c                  9  : Se obtiene  informacion referente  al desarrollo
c                       de las iteraciones:  el numero de  restricciones
c                       activas   y cuales  son  (para identificarlas se
c                       sigue el convenio anterior), la restriccion que
c                       se anade  o  se  elimina  del conjunto activo,
c                       el  tipo  de  direccion   de  descenso calculada
c                       y las iteraciones  con  punto degenerado.
c
c               > 10  : Ademas,  se  obtiene  el punto calculado en cada
c                       iteracion e informacion sobre el proceso llevado
c                       a cabo por OPTR01. (Ver valores posibles de  IMP
c                       en  OPTR01).
c
c        IO     Output chanel number.
c
c        N      Number of problem variables (N > 1).
c
c        MODO   Variable que indica el modo de comienzo del proceso con
c               los valores:
c                  1  : No se facilita un punto inicial,ni restricciones
c                       recomendadas como activas.
c                  2  : Igual que MODO = 1, pero la subrutina trabaja
c                       de forma diferente. Este caso debe utilizarse
c                       cuando hay valores de seguridad para las
c                       variables o las restricciones que son grandes y
c                       que previsiblemente no seran alcanzados.
c                       Tambien puede ser util este MODO cuando la
c                       forma cuadratica es definida positiva.
c                  3  : Se facilita un punto inicial. En este caso
c                       el punto ha de ser admisible para todas las
c                       restricciones.
c
c     OUTPUTS
c
c        X      Vector of dimension N, optimal vector,
c
c        F      Variable  que  contiene  el  valor  del funcional  en el
c               optimo si el proceso  ha finalizado sin problemas.
c
c      LAGR     Vector of double precision of dimension N+MI+MD si IRA es
c               mayor que cero y de dimension MI+MD si IRA=0.En las N
c               primeras componentes contiene los multiplicadores de
c               Lagrange asociados a las restricciones de cota (si IRA
c               es mayor que cero) y en las MI+MD ultimas componentes
c               contiene los multiplicadores asociados a las de igualdad
c               y desigualdad respectivamente.
c
c      INFO     Diagnostic argument
c-------------------------------------------------------------------------     
c     error code
      IF (info .EQ. 1)   info =  1007
      IF (info .EQ. -12) info =  1001
      IF (info .EQ. -1)  info = -1002
      IF (info .EQ. -2)  info = -1003
      IF (info .EQ. -3)  info = -1004
      IF (info .EQ. -4)  info = -1005
      IF (info .EQ. -11) info = -1006
      IF (info .EQ. -13) info = -1003
      IF (info .EQ. -14) info =  1007
c
      RETURN
      END
