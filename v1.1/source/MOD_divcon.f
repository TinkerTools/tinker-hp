c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ######################################################################
c     ##                                                                  ##
c     ##  module divcon  --  specifics of DC-JI/DIIS polarization solver  ##
c     ##                                                                  ##
c     ######################################################################
c
c
c
      module divcon
      implicit none
      integer maxdim,km,zdim,kblks,natprblk 
      integer clsttype
      integer dcx,dcy,dcz,precomp,nocomdiis
      integer, allocatable :: grplst(:)
      integer, allocatable :: atmofst(:)
      integer, allocatable :: kofst(:)
      integer, allocatable :: npergrp(:)
      integer, allocatable :: knblist(:)
      integer, allocatable :: point(:)
      integer, allocatable :: klst(:,:)
      real*8, allocatable :: zmat(:)
      real*8, allocatable :: means(:,:)
      real*8, allocatable :: oldmeans(:,:)
      real*8, allocatable :: ytab(:)
      real*8 kxmin,kymin,kzmin
      real*8 kxmax,kymax,kzmax
      save
      end
