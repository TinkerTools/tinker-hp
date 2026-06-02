!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ####################################################################
!     ##                                                                ##
!     ##  module strtor  --  stretch-torsions in the current structure  ##
!     ##                                                                ##
!     ####################################################################
!
!
module strtor
   implicit none
   integer, pointer :: ist(:,:) !<torsion and bond numbers used in stretch-torsion
   integer, pointer :: nbstrtor(:) !<number of stretch-torsion interactions before each atom
   integer :: nstrtor !<total number of stretch-torsion interactions
   integer :: nstrtorloc !<local number of stretch-torsion interactions
   integer :: winist !<window object corresponding to ist
   integer :: winnbstrtor !<window object corresponding to nbstrtor
   integer :: winkst !<window object corresponding to kst
   real*8, pointer :: kst(:,:) !<1-, 2- and 3-fold stretch-torsion force constants
   save
end
