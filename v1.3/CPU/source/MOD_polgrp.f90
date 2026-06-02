!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #####################################################################
!     ##                                                                 ##
!     ##  module polgrp  --  polarizable site group connectivity lists   ##
!     ##                                                                 ##
!     #####################################################################
!
!
!
module polgrp
   use sizes
   implicit none
   integer :: maxp11 !<maximum number of atoms in a polarization group
   integer :: maxp12 !<maximum number of atoms in groups 1-2 to an atom
   integer :: maxp13 !<maximum number of atoms in groups 1-3 to an atom
   integer :: maxp14 !<maximum number of atoms in groups 1-4 to an atom
   integer :: maxscalp !<maximum number of atoms in scaled polar interaction per atom
   parameter (maxp11=200)
   parameter (maxp12=200)
   parameter (maxp13=200)
   parameter (maxp14=200)
   parameter (maxscalp=480)
   integer :: winnp11 !<window object corresponding to np11
   integer :: winip11 !<window object corresponding to ip11
   integer :: winnp12 !<window object corresponding to np12
   integer :: winip12 !<window object corresponding to ip12
   integer :: winnp13 !<window object corresponding to np13
   integer :: winip13 !<window object corresponding to ip11
   integer :: winnp14 !<window object corresponding to np14
   integer :: winip14 !<window object corresponding to ip14
   integer, pointer :: np11(:) !<number of atoms in polarization group of each atom
   integer, pointer :: ip11(:,:) !<atom numbers of atoms in same group as each atom
   integer, pointer :: np12(:) !<number of atoms in groups 1-2 to each atom
   integer, pointer :: ip12(:,:) !<atom numbers of atoms in groups 1-2 to each atom
   integer, pointer :: np13(:) !<number of atoms in groups 1-3 to each atom
   integer, pointer :: ip13(:,:) !<atom numbers of atoms in groups 1-3 to each atom
   integer, pointer :: np14(:) !<number of atoms in groups 1-4 to each atom
   integer, pointer :: ip14(:,:) !<atom numbers of atoms in groups 1-4 to each atom
   save
end
