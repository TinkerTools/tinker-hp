!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ######################################################################
!     ##                                                                  ##
!     ##  module opbend  --  out-of-plane bends in the current structure  ##
!     ##                                                                  ##
!     ######################################################################
!
!
module opbend
   implicit none
   integer :: nopbend !<total number of out-of-plane bends in the system
   integer :: nopbendloc !<local number of out-of-plane bends in the system
   integer :: winiopb !<window object corresponding to iopb
   integer :: winnbopbend !<window object corresponding to nbopbend
   integer :: winopbk !<window object corresponding to opbk
   integer, pointer :: iopb(:) !<bond angle numbers used in out-of-plane bending
   integer, pointer :: nbopbend(:) !<number of angle used in out-of-plane bending before each atom
   real*8, pointer ::  opbk(:) !<force constant values for out-of-plane bending
   save
end
