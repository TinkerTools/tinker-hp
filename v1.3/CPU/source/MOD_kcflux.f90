!
!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###############################################################
!     ##                                                           ##
!     ##  module kcflux -- charge flux term forcefield parameters  ##
!     ##                                                           ##
!     ###############################################################
!
!
module kcflux
   implicit none
   integer :: maxncfb !<maximum number of bond stretch charge flux entries
   integer :: maxncfa !<maximum number of angle bend charge flux entries
   parameter (maxncfb=2000)
   parameter (maxncfa=2000)
   real*8 :: cflb(maxncfb) !<charge flux over stretching of a bond length
   real*8 :: cfla(2,maxncfa) !<charge flux over bending of a bond angle
   real*8 :: cflab(2,maxncfa) !<charge flux over asymmetric bond within an angle
   character*8 :: kcfb(maxncfb) !<string of atom classes for bond stretch charge flux
   character*12 :: kcfa(maxncfa) !<string of atom classes for angle bend charge flux
   save
end
