!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ######################################################################
!     ##                                                                  ##
!     ##  module vdw  --  van der Waals parameters for current structure  ##
!     ##                                                                  ##
!     ######################################################################
!
!
!
module vdw
   implicit none
   integer :: nvdw !<total number van der Waals active sites in the system
   integer :: nvt !<number of distinct vdw types/classes in the system
   integer :: nvdwloc !<local number van der Waals active sites in the system
   integer :: nvdwbloc !<local+neighbors number van der Waals active sites in the system
   integer :: nvdwlocnl!<localnl number van der Waals active sites in the system
   integer :: winjvdw !<window object corresponding to jvdw
   integer :: winivdw !<window object corresponding to ivdw
   integer :: winired !<window object corresponding to ired
   integer :: winivt !<window object corresponding to ivt
   integer :: winjvt !<window object corresponding to jvt
   integer :: winnbvdw !<window object corresponding to nbvdw
   integer, pointer :: jvdw(:) !<type or class index into vdw parameters for each atom
   integer, pointer :: ivdw(:) !<number of the atom for each van der Waals active site
   integer, pointer :: ired(:) !<attached atom from which reduction factor is applied
   integer, pointer :: ivt(:) !<type/class index for each distinct vdw type or class
   integer, pointer :: jvt(:) !<frequency of each vdw type or class in the system
   integer, pointer :: nbvdw(:) !<number of 'vdw' atoms before each atom
   integer :: winradmin !<window object corresponding to radmin
   integer :: winepsilon !<window object corresponding to epsilon
   integer :: winradmin4 !<window object corresponding to radmin4
   integer :: winepsilon4 !<window object corresponding to epsilon4
   integer :: winradhbnd !<window object corresponding to radhbnd
   integer :: winepshbnd !<window object corresponding to epshbnd
   integer :: winkred !<window object corresponding to kred
   integer :: nvdwblocloop
   integer, allocatable :: vdwlocnl(:) !<glob-locnl vdw correspondance
   real*8, pointer :: radmin(:,:) !<minimum energy distance for each atom class pair
   real*8, pointer :: epsilon(:,:) !<well depth parameter for each atom class pair
   real*8, pointer :: radmin4(:,:) !<minimum energy distance for 1-4 interaction pairs
   real*8, pointer :: epsilon4(:,:) !<well depth parameter for 1-4 interaction pairs
   real*8, pointer :: radhbnd(:,:) !<minimum energy distance for hydrogen bonding pairs
   real*8, pointer :: epshbnd(:,:) !<minimum energy distance for hydrogen bonding pairs
   real*8, pointer :: kred(:) !<value of reduction factor parameter for each atom
   save
end
