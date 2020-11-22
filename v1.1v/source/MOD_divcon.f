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
      !DIR$ ATTRIBUTES ALIGN:64 :: grplst
      integer, allocatable :: grplst(:)
      !DIR$ ATTRIBUTES ALIGN:64 :: atmofst
      integer, allocatable :: atmofst(:)
      !DIR$ ATTRIBUTES ALIGN:64 :: kofst
      integer, allocatable :: kofst(:)
      !DIR$ ATTRIBUTES ALIGN:64 :: npergrp
      integer, allocatable :: npergrp(:)
      !DIR$ ATTRIBUTES ALIGN:64 :: knblist
      integer, allocatable :: knblist(:)
      !DIR$ ATTRIBUTES ALIGN:64 :: point
      integer, allocatable :: point(:)
      !DIR$ ATTRIBUTES ALIGN:64 :: klst
      integer, allocatable :: klst(:,:)
      !DIR$ ATTRIBUTES ALIGN:64:: zmat
      real*8, allocatable :: zmat(:)
      !DIR$ ATTRIBUTES ALIGN:64:: means
      real*8, allocatable :: means(:,:)
      !DIR$ ATTRIBUTES ALIGN:64:: oldmeans
      real*8, allocatable :: oldmeans(:,:)
      !DIR$ ATTRIBUTES ALIGN:64:: ytab
      real*8, allocatable :: ytab(:)
      real*8 kxmin,kymin,kzmin
      real*8 kxmax,kymax,kzmax
      save
      end
