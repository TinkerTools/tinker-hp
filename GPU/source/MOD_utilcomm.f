c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ######################################################
c     ##                                                  ##
c     ##  utils communication module  -- Tinker-HP        ##
c     ##                                                  ##
c     ######################################################
c
c
c
c     Arrays to be used with MPI to exchange data and tracks requests
c     Some of them serve mainly as a memory pool to minimise
c     reallocation

c
c     do_not_commpole  switch to decide whether or not we should
c                      communicate __poleglob__
c     no_commdir       switch to avoid commfield communications
c     skpPcomm         switch to avoid commdirdir communication routine
c

#include "tinker_macro.h"
      module utilcomm

      logical::  do_not_commpole=.false.
     &          , no_commdir=.false.

      !  Operations  Options
      !      +bytes   +count
      !         765    43210
      !         000    00000
      ! Module parameters
      enum,bind(C)
      enumerator ucBig,ucNeig,ucShort,ucRec,ucDirRec
      enumerator ucNfea
      enumerator ucSend=32,ucRecv=64,ucWait=128
      end enum
      integer    ucComm
      logical    skpPcomm
      character*6 ucOptName(ucNfea)

      ! Requests arrays for async communication
      integer   ,allocatable,target::
     &            reqs_poleglob(:)  , reqr_poleglob(:)
     &          , reqs_dird_a(:)    , reqr_dird_a(:)
     &          , reqs_dirdir(:)    , reqr_dirdir(:)
     &          , reqs_recdir(:)    , reqr_recdir(:)
     &          , reqs_dirrec(:)    , reqr_dirrec(:)
     &          , reqs_recrec(:)    , reqr_recrec(:)
     &          , reqs_recdirsolv(:), reqr_recdirsolv(:)

      ! polarisation MPI Buffer
      real(t_p) ,allocatable,target:: buff_field(:,:,:,:)
      real(t_p) ,allocatable,target:: buffermpi2d(:,:),buffermpi2d1(:,:)
     &          , buffermpimu1(:,:,:),buffermpimu2(:,:,:)
     &          , buffermpi3d (:,:,:),buffermpi3d1(:,:,:)

      ! Control variables
      logical   ,protected ::l_r1,l_r2,l_p1,l_p2,l_p3,l_f0,l_f1

      ! General Pool
      integer   ,allocatable,target:: buffMpi_i1(:),buffMpi_i2(:)
      real(t_p) ,allocatable,target:: buffMpi_r1(:),buffMpi_r2(:)
      real(r_p) ,allocatable,target:: buffMpi_p1(:),buffMpi_p2(:)
     &          , buffMpi_p3(:)
      mdyn_rtyp ,allocatable,target:: buffMpi_f0(:),buffMpi_f1(:)

      type ucPool 
      integer   ,allocatable:: buffMpi_i1(:),buffMpi_i2(:)
      real(t_p) ,allocatable:: buffMpi_r1(:),buffMpi_r2(:)
      real(r_p) ,allocatable:: buffMpi_p1(:),buffMpi_p2(:)
     &          , buffMpi_p3(:)
      mdyn_rtyp ,allocatable:: buffMpi_f0(:),buffMpi_f1(:)
      end type

      type(ucPool),protected,target:: ucP0,ucP1

      parameter (ucComm=ucSend+ucRecv+ucWait
     &          ,skpPcomm=.false.
     &          ,ucOptName=['Big','Neig','Short','Rec','DirRec']
     &          )

      interface commDDd_ext
      module subroutine commDDd_ext_c(extb,dim,opt)
      implicit  none
      real(t_p),intent(inout):: extb(*)
      integer  ,optional     :: dim,opt
      end subroutine
      end interface
      interface commDDd_add
      module subroutine commDDd_add_c(addb,dim,opt)
      implicit  none
      real(t_p),intent(inout):: addb(*)
      integer  ,optional     :: dim,opt
      end subroutine
      end interface
      interface commDDrd
      module subroutine commDDrd_add_c(dirb,recb,dim,opt)
      implicit  none
      real(t_p),intent(inout):: dirb(*)
      real(t_p),intent(inout):: recb(*)
      integer  ,optional     :: dim,opt
      end subroutine
      end interface

      end module
