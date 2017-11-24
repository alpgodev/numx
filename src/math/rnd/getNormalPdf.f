!---------------------------------------------------------------------------------------------------------------
! getNormalPdf returns the density function of a standard normal random variable
      double precision function getNormalPdf(x)
      double precision :: x
      
      getNormalPdf = 0.398942280401433D0 * exp(-x * x * 0.5D0) 
      !0.398942280401433 = 1/sqrt(2*pi)
      end function
