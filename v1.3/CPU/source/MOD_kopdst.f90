!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #######################################################################
!     ##                                                                   ##
!     ##  module kopdst  --  forcefield parameters for out-plane distance  ##
!     ##                                                                   ##
!     #######################################################################
!
!
module kopdst
   implicit none
   integer :: maxnopd !<maximum number of out-of-plane distance entries
   parameter (maxnopd=500)
   real*8 :: opds(maxnopd) !<force constant parameters for out-of-plane distance
   character*16 :: kopd(maxnopd) !<string of atom classes for out-of-plane distance
   save
end
