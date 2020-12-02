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

c     do_not_commpole  switch to decide whether or not we should
c                      communicate __poleglob__

#include "tinker_precision.h"
#include "tinker_types.h"

      module utilcomm

      logical::do_not_commpole=.false.

      ! Requests arrays for async communication
      integer   ,allocatable,target::
     &            reqs_poleglob(:)  , reqr_poleglob(:)
     &          , reqs_dirdir(:)    , reqr_dirdir(:)
     &          , reqs_recdir(:)    , reqr_recdir(:)
     &          , reqs_recdirsolv(:), reqr_recdirsolv(:)

      real(t_p) ,allocatable,target:: buff_field(:,:,:,:)
      ! polarisation MPI Buffer
      real(t_p) ,allocatable,target:: buffermpi2d(:,:),buffermpi2d1(:,:)
     &          , buffermpimu1(:,:,:),buffermpimu2(:,:,:)
     &          , buffermpi3d (:,:,:),buffermpi3d1(:,:,:)

      ! Pool
      real(r_p) ,allocatable,target:: buffMpi_p1(:),buffMpi_p2(:)
      mdyn_rtyp ,allocatable,target:: buffMpi_fr(:)

      end module


