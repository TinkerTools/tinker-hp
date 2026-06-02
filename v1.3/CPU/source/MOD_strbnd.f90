!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #################################################################
!     ##                                                             ##
!     ##  module strbnd  --  stretch-bends in the current structure  ##
!     ##                                                             ##
!     #################################################################
!
!
!
module strbnd
   implicit none
   integer :: nstrbnd !<total number of stretch-bend interactions
   integer :: nstrbndloc !<local number of stretch-bend interactions
   integer, pointer :: isb(:,:) !<angle and bond numbers used in stretch-bend
   integer, pointer :: nbstrbnd(:) !<number of stretch-bend interactions before each atom
   integer :: winisb !<window object corresponding to isb
   integer :: winnbstrbnd !<window object corresponding to nbstrbnd
   integer :: winsbk !<window object corresponding to sbk
   real*8, pointer :: sbk(:,:) !<force constants for stretch-bend terms
   save
end
