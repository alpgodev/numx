c=======================================================================
c
c     subroutine MCGEN              
c
c     PURPOSE Computes the generic MC
c
c-----------------------------------------------------------------------
c
c     Input:
c     
c     PAYOFF: external routine for the option PAYOFF.
c     MCSIMU:  external routine for the asset diffusion.
c     nopt:  number of options we price with MC,               integer.
c     indpo: nopt-vector, 
c               indpo(i)=0 means no pricing of the i-th option
c               1 means prcing of the i-th option.             integer. 
c     nsimul: number of simulations,                           integer.
c     nasset: number of assets,                                integer.
c     asset: vector of asset quotations at the present time t, double.
c     matur: maturity (in years)vector-nopt,                   double.
c     nstep: number of steps to maturity,                      integer.
c     npari: size of pari vector-7,                            integer.
c     pari: integer user data for mcsim vector-npari,          integer.
c     npard: size of pard vector-(nt+nu+nx+4*nasset+17),       integer.
c     pard: double user data for mcsimu and payoff,            double.
c
c     WORKSPACE:
c           
c     dwork: double workspace size (nstep+1)*nasset+2*nopt.
c
c     Output:
c     
c     vpay:  Payoff simulations matrix(nopt,nsimul),           double.
c     price: Price of the option vector-nopt,                  double. 
c     cinfo: path of error                                     character*80.
c     info: = 0 successful exit                                integer.
c           = -1 nsimul=0
c
c     CALL: 
c
c     GESTERR:  Computes the error paths.
c     MV:       Computes the mean of a vector.
c
c---------------------------------------------------------------------------
c
      subroutine MCGEN(PAYOFF,MCSIMU,nopt,indpo,nsimul,nasset,asset,
     &     matur,nstep,npari,pari,npard,pard,dwork,vpay,price,
     &     cinfo,info)


      implicit none

      external PAYOFF,MCSIMU

      integer nopt,nsimul,nasset,nstep,npari,npard,info
      integer pari(npari,*),indpo(*)
      double precision asset(*),matur(*),pard(npard,*),price(*)
      double precision myzero,vpay(nopt,*),coef
      character*80 cinfo
      
      double precision dwork(*)
      integer pdtr,pdval,pdpri


      integer i,iosort

      parameter (myzero=0.,iosort=6)
#ifdef INOUT     
      integer j
#endif
c
c     seed for random
c
c#ifdef INTELFOR 
c      call TIME(seed)
c      call RANDOM_SEED(seed)
c#else
c      seed = time()
c      call srand(seed)
c#endif  

c
c     double workspace
c
  
      pdtr=1
cpdpri+nopt

c
c     needs (nstep+1)* nasset
c     
      pdval=pdtr+(nstep+1)* nasset

c
c     needs nopt
c
      pdpri=pdval+nopt
c
c     needs nopt
c


c     global size is (nstep+1)*nasset+2*nopt
c
      if (nsimul .eq. 0) then
         info=-1
         cinfo='MCGEN'
         return
      endif


      do i=1,nsimul

c
c     MC simulation 
c

         call MCSIMU(nasset,asset,nstep,matur,npard,pard,npari,pari
     &    ,dwork(pdtr),cinfo,info)
         
         if(info.ne.0) then
            call GESTERR ( cinfo, 'MCSIM/MCGEN/        ' )
            return
         endif

c
c        computing the payoff
c

         call PAYOFF ( nopt,indpo,nasset, nstep+1, dwork(pdtr), matur,
     &        npari,pari, npard, pard,
     &        dwork(pdval), cinfo,info )

         if (info .ne. 0) then
            call GESTERR ( cinfo, 'PAYOFF/MCGEN/        ' )
            return
         endif
c
c Copy the i-th simulation price in vpay i-th column
c
         
 
         call YVCM ( nopt,nsimul,vpay,dwork(pdval), i,info)

         if(info.ne.0) then
            call GESTERR ( cinfo, 'PAYOFF/MCGEN/        ' )
            return
         endif
c
      end do

   

c
c     computing price
c
      call IVX(nopt,price,myzero)
c
c     Sum of the lines of vpay
c
  

      call SLM ( nopt, nsimul, vpay, price)


c
c     Mean of the simulations
c     
     

      coef=1./nsimul
     

      call PVX(nopt,price,coef,dwork(pdpri))
      
      call YV ( nopt, dwork(pdpri), price )

#ifdef INOUT     
      call IMPRV(iosort,' price ',nopt,price)
       write(iosort,*)  'nsimul ',nsimul
       do i=1,nopt
          write(*,*) (vpay(i,j),j=1,nsimul)
       enddo
#endif
       


       return
      end
