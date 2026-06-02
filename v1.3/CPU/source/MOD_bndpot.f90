!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #####################################################################
!     ##                                                                 ##
!     ##  module bndpot  --  specifics of bond stretch functional forms  ##
!     ##                                                                 ##
!     #####################################################################
!
!
module bndpot
   implicit none
   integer :: winbndtyp !<window object for bndtyp
   real*8 :: cbnd !<cubic coefficient in bond stretch potential
   real*8 :: qbnd !<quartic coefficient in bond stretch potential
   real*8 :: bndunit !<convert bond stretch energy to kcal/mole
   character*8 :: bndtyp_default !<default bond stretch potential energy function
   character*8, pointer :: bndtyp(:)  !<type of bond stretch potential energy function
   save
end
