!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #################################################################
!     ##                                                             ##
!     ##  module couple  --  near-neighbor atom connectivity lists   ##
!     ##                                                             ##
!     #################################################################
!
!
!     maxn13   !<maximum number of atoms 1-3 connected to an atom
!     maxn14   !<maximum number of atoms 1-4 connected to an atom
!     maxn15   !<maximum number of atoms 1-5 connected to an atom
!
!     n12      !<number of atoms directly bonded to each atom
!     winn12    !<window object corresponding to n12
!     i12      !<atom numbers of atoms 1-2 connected to each atom
!     wini12    !<window object corresponding to i12
!     n13      !<number of atoms in a 1-3 relation to each atom
!     winn13    !<window object corresponding to n13
!     i13      !<atom numbers of atoms 1-3 connected to each atom
!     wini13    !<window object corresponding to i13
!     n14      !<number of atoms in a 1-4 relation to each atom
!     winn14    !<window object corresponding to n14
!     i14      !<atom numbers of atoms 1-4 connected to each atom
!     wini14    !<window object corresponding to i14
!     n15      !<number of atoms in a 1-5 relation to each atom
!     winn15    !<window object corresponding to n15
!     i15      !<atom numbers of atoms 1-5 connected to each atom
!     wini15    !<window object corresponding to i15
!
!
module couple
   use sizes
   implicit none
   integer :: maxn13 !<maximum number of atoms 1-3 connected to an atom
   integer :: maxn14 !<maximum number of atoms 1-4 connected to an atom
   integer :: maxn15 !<maximum number of atoms 1-5 connected to an atom
   parameter (maxn13=3*maxvalue)
   parameter (maxn14=3*maxvalue)
   parameter (maxn15=3*maxvalue)
   integer, allocatable :: n12(:) !<number of atoms directly bonded to each atom
   integer, allocatable :: i12(:,:) !<atom numbers of atoms 1-2 connected to each atom
   integer, pointer :: n13(:) !<number of atoms in a 1-3 relation to each atom
   integer, pointer :: i13(:,:) !<number of atoms in a 1-3 relation to each atom
   integer, pointer ::  n14(:) !<number of atoms in a 1-4 relation to each atom
   integer, pointer :: i14(:,:) !<atom numbers of atoms 1-4 connected to each atom
   integer, pointer :: n15(:) !<number of atoms in a 1-5 relation to each atom
   integer, pointer :: i15(:,:) !<atom numbers of atoms 1-5 connected to each atom
   integer :: winn13 !<window object corresponding to n13
   integer :: wini13 !<window object corresponding to i13
   integer :: winn14 !<window object corresponding to n14
   integer :: wini14 !<window object corresponding to i14
   integer :: winn15 !<window object corresponding to n15
   integer :: wini15 !<window object corresponding to i15
   save
end
