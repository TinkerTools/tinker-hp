!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###################################################################
!     ##                                                               ##
!     ##  module angpot  --  specifics of angle bend functional forms  ##
!     ##                                                               ##
!     ###################################################################
!
!
module angpot
   implicit none
   integer :: winangtyp!<window object corresponding to angtyp
   integer :: idx_ps(245,3)!<indices for the PS water potential
   character*8 :: opbtyp!<type of out-of-plane bend potential energy function
   character*8, pointer :: angtyp(:)!<type of angle bending function for each bond angle
   real*8 :: c5z_ps(245)!<c5z coefficient for the PS water potential
   real*8 :: angunit!<convert angle bending energy to kcal/mole
   real*8 :: stbnunit!<convert stretch-bend energy to kcal/mole
   real*8 :: aaunit!<convert angle-angle energy to kcal/mole
   real*8 :: opbunit!<convert out-of-plane bend energy to kcal/mole
   real*8 :: opdunit!<convert out-of-plane distance energy to kcal/mole
   real*8 :: cang!<cubic coefficient in angle bending potential
   real*8 :: qang!<quartic coefficient in angle bending potential
   real*8 :: pang!<quintic coefficient in angle bending potential
   real*8 :: sang!<sextic coefficient in angle bending potential
   real*8 :: copb!<cubic coefficient in out-of-plane bend potential
   real*8 :: qopb!<quartic coefficient in out-of-plane bend potential
   real*8 :: popb!<quintic coefficient in out-of-plane bend potential
   real*8 :: sopb!<sextic coefficient in out-of-plane bend potential
   real*8 :: copd!<cubic coefficient in out-of-plane distance potential
   real*8 :: qopd!<quartic coefficient in out-of-plane distance potential
   real*8 :: popd!<quintic coefficient in out-of-plane distance potential
   real*8 :: sopd!<sextic coefficient in out-of-plane distance potential
   save
end
