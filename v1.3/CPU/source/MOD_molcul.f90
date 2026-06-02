!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #####################################################################
!     ##                                                                 ##
!     ##  module molcul  --  individual molecules within current system  ##
!     ##                                                                 ##
!     #####################################################################
!
!
module molcul
   implicit none
   integer :: nmol !<total number of separate molecules in the system
   integer :: nmoleloc !<local number of separate molecules in the system
   integer :: winmolcule !<window object corresponding to molcule
   integer :: winkmol !<window object corresponding to kmol
   integer :: winimol !<window object corresponding to imol
   integer :: winmolmass !<window object corresponding to molmass
   integer, pointer :: molcule(:) !<number of the molecule to which each atom belongs
   integer, pointer :: kmol(:) !<contiguous list of the atoms in each molecule
   integer, pointer :: imol(:,:) !<first and last atom of each molecule in the list
   real*8 :: totmass !<total weight of all the molecules in the system
   real*8, pointer :: molmass(:) !<molecular weight for each molecule in the system
   save
end
