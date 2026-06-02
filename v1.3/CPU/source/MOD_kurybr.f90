!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #######################################################################
!     ##                                                                   ##
!     ##  module kurybr  --  forcefield parameters for Urey-Bradley terms  ##
!     ##                                                                   ##
!     #######################################################################
!
!
!     maxnu   !<maximum number of Urey-Bradley parameter entries
!
!     ucon    !<force constant parameters for Urey-Bradley terms
!     dst13   !<ideal 1-3 distance parameters for Urey-Bradley terms
!     ku      !<string of atom classes for Urey-Bradley terms
!
!
module kurybr
   implicit none
   integer :: maxnu !<maximum number of Urey-Bradley parameter entries
   integer :: maxnups !<maximum number of angle repulsion parameter entries
   integer :: maxnuq !<maximum number of quartic urey parameter entries
   parameter (maxnu=2000)
   parameter (maxnups=500)
   parameter (maxnuq=500)
   real*8 :: ucon(maxnu) !<force constant parameters for Urey-Bradley terms
   real*8 :: dst13(maxnu) !<ideal 1-3 distance parameters for Urey-Bradley terms
   real*8 :: uconps(maxnups) !<force constant parameters for angle repulsion terms
   real*8 :: dst13ps(maxnups) !<ideal 1-3 distance parameters for angle repulsion terms
   real*8 :: uconq(maxnuq) !<force constant parameters for quartic urey terms
   real*8 :: dst13q(maxnuq) !<ideal 1-3 distance parameters for quartic urey terms
   character*12 :: ku(maxnu) !<string of atom classes for Urey-Bradley terms
   character*12 :: kups(maxnups) !<string of atom classes for angle repulsion terms
   character*12 :: kuq(maxnuq) !<string of atom classes for quartic urey terms
   save
end
