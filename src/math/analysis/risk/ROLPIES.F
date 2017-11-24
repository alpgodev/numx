c=======================================================================
c
c     subroutine ROLPIES                                     
c
c     Historical weights (rolling pies)
c
c     w[i] = w[0] * (x[i]/x[0])*(y[0]/y[i])
c
c     w[t] : weight at time t
c     x[t] : asset value at time t
c     y[t] : portfolio value at time t
c
c-----------------------------------------------------------------------
      SUBROUTINE ROLPIES ( n, wini, x, y, mode, w, info)
c-----------------------------------------------------------------------
c
c     INPUT :
c            n     : number of value(s) (n > 1)                 integer
c            wini  : initial weight                               double
c            x     : asset value(s) (n)                           double
c            y     : portfolio value(s) (n)                       double
c            mode  : =1 past values, =2 futur values             integer
c
c     OUTPUT :
c            w     : historical weight(s) (n)                    double
c            info  : = 0 successful exit                        integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, mode, info
      DOUBLE PRECISION wini
      DOUBLE PRECISION x(*), y(*), w(*)
c
c     local variables      
      INTEGER i
      DOUBLE PRECISION invx, myzero 
      PARAMETER ( myzero = 1.E-30 )
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
      w(1) = wini
c
c     test mode value
      IF ((mode .NE. 1) .AND. (mode .NE. 2)) THEN
        info = -106
        RETURN
      ENDIF   
c
c     test initial asset value
      IF (x(1) .lt. myzero ) THEN
        info = -103
        RETURN
      ENDIF
c
c     historical weight(s)
      IF (mode .eq. 1) THEN
        invx = 1./x(n) 
        DO i = 2,n
            IF (y(i) .lt. myzero ) THEN
                info = -103
                RETURN
            ENDIF
            w(i) = w(1) * (invx*x(n - i + 1)) * (y(n)/y(n - i + 1))  
        ENDDO
      ENDIF
c
      IF (mode .eq. 2) THEN
        invx = 1./x(1) 
        DO i = 2,n
            IF (y(i) .lt. myzero ) THEN
                info = -103
                RETURN
            ENDIF    
            w(i) = w(1) * (invx*x(i)) * (y(1)/y(i))  
        ENDDO
      ENDIF
c      
      RETURN
      END
