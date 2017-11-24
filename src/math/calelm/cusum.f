      SUBROUTINE cusum(n,w)
c     Utility fct: cumulated sum
      DOUBLE PRECISION w(*), t
      t = 0.0d0
      DO 10 k=1,n
      w(k) = t+w(k)
      t    = w(k)
   10 CONTINUE
      END
