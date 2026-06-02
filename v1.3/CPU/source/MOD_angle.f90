!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ##################################################################
!     ##                                                              ##
!     ##  module angle  --  bond angles within the current structure  ##
!     ##                                                              ##
!     ##################################################################
!
!
module angle
   implicit none
   integer :: nangle!<total number of bond angles in the system
   integer :: nangleloc!<numbers of bond angles in the local domain
   integer :: winiang!<window object corresponding to iang
   integer :: winak!<window object corresponding to ak
   integer :: winanat!<window object corresponding to anat
   integer :: winafld!<window object corresponding to afld
   integer, allocatable :: angleloc(:)!<correspondance between global and local bond angles
   integer, pointer :: iang(:,:)!<numbers of the atoms in each bond angle
   real*8, pointer ::  ak(:)!<harmonic angle force constant (kcal/mole/rad**2)
   real*8, pointer ::  anat(:)!<ideal bond angle or phase shift angle (degrees)
   real*8, pointer ::  afld(:)!<periodicity for Fourier bond angle term
   save
end
