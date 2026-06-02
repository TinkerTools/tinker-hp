!
!
!     ###################################################
!     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
!     ##              All Rights Reserved              ##
!     ###################################################
!
!     ######################################################################
!     ##                                                                  ##
!     ##  module divcon  --  specifics of DC-JI/DIIS polarization solver  ##
!     ##                                                                  ##
!     ######################################################################
!
!
!
module divcon
   implicit none
   integer :: maxdim !<max number of atoms per block 
   integer :: km !<number of blocks
   integer :: zdim
   integer :: natprblk !<block size
   integer :: clsttype !<cluster algorithm
   integer :: dcx 
   integer :: dcy
   integer :: dcz
   integer :: precomp !<rely on precomputation
   integer :: nocomdiis !<diis extrapolation without MPI communications
   integer, allocatable :: grplst(:) !<group list 
   integer, allocatable :: atmofst(:) !<atom offset
   integer, allocatable :: kofst(:)  !<where each Zmat begins in the array
   integer, allocatable :: npergrp(:) !<number of atoms per group
   integer, allocatable :: knblist(:) !<neighbor list indexing
   integer, allocatable :: point(:) !<neighbor list indexing
   integer, allocatable :: klst(:,:) !<neighbor list
   real*8, allocatable :: zmat(:) !<Z matrix
   real*8, allocatable :: means(:,:) !<mean values in blocks
   real*8, allocatable :: oldmeans(:,:) !<old mean values
   real*8, allocatable :: ytab(:) !<induced dipole eletric field
   save
end
