!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ##############################################################
!     ##                                                          ##
!     ##  module angtor  --  angle-torsions in current structure  ##
!     ##                                                          ##
!     ##############################################################
!
!
module angtor
   implicit none
   integer :: nangtor!<total number of angle-torsion interactions
   integer :: nangtorloc!<local number of angle-torsion interactions
   integer, pointer :: iat(:,:)!<torsion and angle numbers used in angle-torsion
   integer, pointer :: nbangtor(:)!<number of angle-torsion interactions before each atom
   real*8, pointer :: kant(:,:)!<1-, 2- and 3-fold angle-torsion force constants
   integer :: winiat!<window object corresponding to iat
   integer :: winnbangtor!<window object corresponding to nbangtor
   integer :: winkant!<window object corresponding to kant
   save
end
