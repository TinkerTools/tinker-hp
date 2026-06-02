!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #################################################################
!     ##                                                             ##
!     ##  module group  --  partitioning of system into atom groups  ##
!     ##                                                             ##
!     #################################################################
!
!
module group
   use sizes
   implicit none
   integer :: ngrp !<total number of atom groups in the system
   integer, allocatable :: kgrp(:) !<contiguous list of the atoms in each group
   integer, allocatable :: grplist(:) !<number of the group to which each atom belongs
   integer, allocatable :: igrp(:,:) !<first and last atom of each group in the list
   real*8, allocatable :: grpmass(:) !<total mass of all the atoms in each group
   real*8, allocatable :: wgrp(:,:) !<weight for each set of group-group interactions
   logical :: use_group !<flag to use partitioning of system into groups
   logical :: use_intra !<flag to include only intragroup interactions
   logical :: use_inter !<flag to include only intergroup interactions
   integer :: natgroup !<number of atoms involved in a group
   integer :: nlocatgroup !<number of local atoms involved in a group
   integer :: npolegroup !<number of multipoles involved in a group
   integer :: npolelocgroup !<number of multipoles involved in a group
   integer, allocatable :: globglobgroup(:) !<associated indexes
   integer, allocatable :: loclocgroup(:) !<global-group correspondance
   integer, allocatable :: globgroup(:) !<associated group indexes
   integer, allocatable :: locgroup(:) !<global group - local group correspondance
   integer, allocatable :: poleglobgroup(:) !<local-global correspondance for the multipoles in the group frame
   integer, allocatable :: polelocgroup(:) !<global-local correspondance in the group frame
   integer, allocatable :: domlengroup(:) !<number of atoms in the group per domain
   integer, allocatable :: ipolegroup(:) !<associated multipoles indexes
   integer, allocatable :: domlenpolegroup(:) !<number of multipoles in the group per domain
   integer, allocatable :: bufbegpolegroup(:) !<"bufbegpole" equivalent in the group frame
   integer, allocatable :: bufbeggroup(:) !<"bufbeg" equivalent in the group frame
   integer, allocatable :: pollistgroup(:) !<atoms-multipoles correspondance in the group frame
   real*8, allocatable :: depgroup(:,:) !<polarization energy derivatives associated to the group
   real*8, allocatable :: uindgroup(:,:) !<induced dipoles associated to the group
   real*8, allocatable :: uinpgroup(:,:) !<induced dipoles associated to the group
   real*8 :: vir_group(3,3) !<virial associated to the polarization energy of the group
   real*8 :: epgroup !<polarization energy associated to the group
   save
end
