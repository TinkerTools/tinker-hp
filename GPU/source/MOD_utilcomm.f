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
#include "tinker_precision.h"

      module utilcomm

      logical::do_not_commpole=.false.

      integer   ,allocatable,target::reqsendpoleglob(:)
     &          ,reqrecvpoleglob(:)
      real(t_p) ,allocatable,target::buff_field(:,:,:,:)
      real(t_p) ,allocatable,target::buffermpi2d(:,:),buffermpi2d1(:,:)
     &          ,buffermpimu1(:,:,:),buffermpimu2(:,:,:)
     &          ,buffermpi3d(:,:,:),buffermpi3d1(:,:,:)
      real(r_p) ,allocatable,target:: buffMpi_p1(:),buffMpi_p2(:)

      end module


