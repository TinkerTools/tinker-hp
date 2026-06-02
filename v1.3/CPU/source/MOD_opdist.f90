!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ######################################################################
!     ##                                                                  ##
!     ##  module opdist  --  out-of-plane distances in current structure  ##
!     ##                                                                  ##
!     ######################################################################
!
!
module opdist
   implicit none
   integer :: nopdist !<total number of out-of-plane distances in the system
   integer :: nopdistloc !<local number of out-of-plane distances in the system
   integer :: winiopd !<window object corresponding to iopd
   integer :: winnbopdist !<window object corresponding to nbopdist
   integer :: winopdk !<window object corresponding to opdk
   integer, pointer :: iopd(:,:) !<numbers of the atoms in each out-of-plane distance
   integer, pointer :: nbopdist(:) !<number of angle used in out-of-plane distance before each atom
   real*8, pointer ::  opdk(:) !<force constant values for out-of-plane distance
   save
end
