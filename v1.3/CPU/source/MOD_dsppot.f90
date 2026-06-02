!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!
!     ###############################################################
!     ##                                                           ##
!     ##  module dsppot  --  dispersion interaction scale factors  ##
!     ##                                                           ##
!     ###############################################################
!
!
module dsppot
   implicit none
   real*8 :: dsp2scale !<scale factor for 1-2 dispersion energy interactions
   real*8 :: dsp3scale !<scale factor for 1-3 dispersion energy interactions
   real*8 :: dsp4scale !<scale factor for 1-4 dispersion energy interactions
   real*8 :: dsp5scale !<scale factor for 1-5 dispersion energy interactions
   logical :: use_dcorr !<flag to use long range dispersion correction
   save
end
