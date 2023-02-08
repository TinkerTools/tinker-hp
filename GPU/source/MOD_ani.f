c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##################################################################
c     ##                                                              ##
c     ##  module ani  --                                              ##
c     ##                                                              ##
c     ##################################################################
c
c
c     MLpot   machine learning potential name
c     s_save
c     MLpotcut
c     xyz_c
c
#include "tinker_macro.h"
      module ani
      use, intrinsic :: iso_c_binding

      integer, private:: s_save=0
      integer    naml, namloc
      integer :: ml_embedding_mode=1
      integer, allocatable:: list(:), iaml(:), amloc(:)
     &       , atomic_ani(:)
      logical :: ml_resources_initialized=.FALSE.
      logical, allocatable:: laml(:), grpmllist(:) 
      real(t_p) mlpotcut
      real(t_p) :: mlpotscale=1.0
      real(c_float) cell_a(3,3)
      real(c_float),allocatable,target:: xyz_c(:)
      real(c_float),allocatable:: d2(:), dp(:), aeani(:)
      real(c_float),allocatable:: d_ani(:)
      character*16   MLpot
      character*240  model_file

      interface
        ! Initiate devices and python ressources for ani
        function init_ml_ressources(rank,nn_name,model_file
     &       ,debug)  bind(C)
          import c_int32_t
          integer(c_int32_t) :: init_ml_ressources
          integer(c_int32_t),value :: rank, debug
          character(*) :: nn_name
          character(*) :: model_file
        end function

        function ml_models(m_c, a_en, f_c, ce_a, a_c, 
     &                           l1, l2, dist, dv, istrict,
     &                           nb, nbp, nstrict, dograd
     &                     ) bind (c)
            import c_float,c_double,c_int32_t,c_int64_t
            real(c_float)  m_c(*), ce_a(*)
            real(c_float)  dist(*), dv(*), a_en(*)
            real(c_float) f_c(*)
            integer(c_int32_t) a_c(*)
            integer(c_int32_t) l1(*), l2(*),istrict(*)
            integer(c_int32_t),value :: nb, nbp,nstrict
            integer(c_int32_t),value :: dograd
            integer(c_int32_t) :: ml_models
        end function ml_models
      end interface

      interface
        subroutine set_embedding_weights
        end subroutine
      end interface

      contains

      subroutine realloc_position_buffer(rsize,npairs)
      integer,intent(in):: rsize,npairs

      if (rsize.gt.s_save) then
         if (allocated(xyz_c)) then
!$acc exit data delete(xyz_c,atomic_ani,aeani,d_ani)
            deallocate(xyz_c,atomic_ani,aeani,d_ani)
         end if
         allocate(xyz_c(rsize*3),atomic_ani(rsize),aeani(rsize)
     &           ,d_ani(3*rsize))
!$acc enter data create(xyz_c,atomic_ani,aeani,d_ani)
         s_save = rsize
      end if
      if (npairs.gt.size(d2)) then
         if (allocated(d2)) then
!$acc exit data delete(d2,dp)
            deallocate(d2,dp)
         end if
         allocate(d2(npairs),dp(3*npairs))
!$acc enter data create(d2,dp)
      end if

      end subroutine realloc_position_buffer

      end module
