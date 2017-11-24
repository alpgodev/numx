c=======================================================================
c
c     subroutine CSTSEC                                    
c
c     Sector constraints used by asset allocation methods
c
c-----------------------------------------------------------------------
      SUBROUTINE cstsec ( nasset, nsect, nsubs, asect,
     &                    minsec, equsec, maxsec,
     &                    neq, nin, ccst, bcst )
c-----------------------------------------------------------------------
c
c     INPUT :
c            nasset : size of portfolio                          integer
c            nsect  : number of sectors                          integer
c            nsubs  : number of sub-sectors/sector vector(nsect) integer
c            asect  : assets dispatching by sector
c                     matrix (nasset*nsect)                      integer
c            minsec : sectors min weight matrix (nsect,nssmax)    double
c            equsec : sectors equ weight matrix (nsect,nssmax)    double
c            maxsec : sectors max weight matrix (nsect,nssmax)    double
c                     nssmax : max of nsubs
c
c     OUTPUT :
c            neq    : number of equality constraints             integer
c            nin    : number of inequality constraints           integer
c            ccst   : matrix of constraints (nasset*(nin+neq))    double
c            bcst   : vector of constraints (nin+neq)             double
c                     neq = ( number of equsec(i,j)>=0 )
c                     nin = 2*(somme(nsubs(i)),i=1,nsect)
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER nasset, nsect, neq, nin
      INTEGER nsubs(*), asect(nasset,*)
      DOUBLE PRECISION minsec(nsect,*), equsec(nsect,*), 
     &                 maxsec(nsect,*), bcst(*), ccst(nasset,*)
c
c     local variables
      INTEGER ic, is, jss, ka, nsstot
c
c-----------------------------------------------------------------------
c
c     nb of subsectors
      nsstot = 0
      DO is = 1,nsect
         nsstot = nsstot + nsubs(is)
      ENDDO
c
c     nb of inequality constrainsts
      nin = 2*nsstot
c
c     equality constraints matrix and vector
      neq = 0
      ic = 0
      DO is=1,nsect
         DO jss=1,nsubs(is)
            IF (equsec(is,jss).ge.0.0) THEN
               neq = neq + 1
               ic = ic + 1
               bcst(ic) = equsec(is,jss)
               DO ka=1,nasset
                  IF (asect(ka,is).eq.jss) THEN
                     ccst(ka,ic) = 1.0
                  ELSE
                     ccst(ka,ic) = 0.0
                  ENDIF
               ENDDO
            ENDIF    
         ENDDO
      ENDDO
c
c     inequality constraints matrix and vector
      DO is=1,nsect
         DO jss=1,nsubs(is)
            ic = ic + 1
            bcst(ic) = maxsec(is,jss)
            bcst(ic+nsstot) = -minsec(is,jss)
            DO ka=1,nasset
               IF (asect(ka,is).eq.jss) THEN
                  ccst(ka,ic) = 1.
                  ccst(ka,ic+nsstot) = -1.0
               ELSE
                  ccst(ka,ic) = 0.0
                  ccst(ka,ic+nsstot) = 0.
               ENDIF      
            ENDDO
         ENDDO
      ENDDO
c
      RETURN
      END
