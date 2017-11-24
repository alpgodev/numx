      SUBROUTINE cupro(n,w)
c     Utility fct: cumulated product
      DOUBLE PRECISION w(*), t
      t = 1.0d0
      DO 10 k=1,n
      w(k) = t*w(k)
      t    = w(k)
   10 CONTINUE
      END
