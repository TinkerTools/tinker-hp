c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##################################################################
c     ##                                                              ##
c     ##  module colvars  --  Colvars/Tinker-HP interface             ##
c     ##                                                              ##
c     ##################################################################
c
#include "tinker_precision.h"
      module colvars
      use iso_c_binding
      implicit none
#ifdef COLVARS
      interface
        subroutine allocate_colvars () bind
     $       ( C, name = "allocate_colvars" )
        end subroutine allocate_colvars
        subroutine compute_colvars_tinker () bind
     $  ( C, name = "compute_colvars" )
        end subroutine compute_colvars_tinker
        subroutine delete_colvars () bind
     $  ( C, name = "delete_colvars" )
        end subroutine delete_colvars
      end interface
#endif
      logical :: use_colvars
      integer :: ncvatoms
      integer  ,allocatable:: cvatoms_ids(:)
      real(r_p),allocatable:: cv_pos(:,:),decv(:,:),decv_tot(:,:)
      real(r_p),target:: dt_sim
      real(r_p),target:: temp_rand

      end module
