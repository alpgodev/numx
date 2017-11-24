      subroutine cuproi(n,wr,wi)
c     Utility fct: cumulated product, complex argument
      double precision wr(*),wi(*),tr,ti,wwr
      tr=1.0d0
      ti=0.0d0
      DO 10 k=1,n
c     w(k)=t*w(k)
      wwr=wr(k)
      wr(k)=tr*wwr-ti*wi(k)
      wi(k)=tr*wi(k)+ti*wwr
      tr=wr(k)
      ti=wi(k)
   10 continue
      END
