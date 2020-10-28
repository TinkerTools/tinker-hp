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
#include "tinker_precision.h"
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
      !DIR$ ATTRIBUTES ALIGN:64:: zmat, zmatDn
      real(t_p), allocatable :: zmat(:), zmatDn(:)
      !DIR$ ATTRIBUTES ALIGN:64:: means
      real(t_p), allocatable :: means(:,:)
      !DIR$ ATTRIBUTES ALIGN:64:: oldmeans
      real(t_p), allocatable :: oldmeans(:,:)
      !DIR$ ATTRIBUTES ALIGN:64:: ytab
      real(t_p), allocatable :: ytab(:)
      real(t_p) kxmin,kymin,kzmin
      real(t_p) kxmax,kymax,kzmax
      save
      end
