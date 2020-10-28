c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ####################################################################
c     ##                                                                ##
c     ##  module utilcu  -- CUDA utility functions and params           ##
c     ##                                                                ##
c     ####################################################################
c
#ifdef _CUDA
#define TINKER_CUF
#include "tinker_precision.h"
      module utilcu
        use cudafor
#if defined(MIXED) || defined(SINGLE)
        use libm ,only: f_sign=>copysignf, f_abs=>fabsf
     &           , f_sqrt=>sqrtf, f_floor=>floorf, f_erfc=>erfcf
     &           , f_exp=>expf
        use cudadevice ,only: f_inv=>__frcp_rn
#else
        use libm ,only: f_sign=>copysign, f_abs=>fabs
     &           , f_floor=>floor, f_erfc=> erfc
        use cudadevice ,only: f_inv=>__drcp_rn
#endif
        use sizes ,only: tinkerdebug
        implicit none
        integer  ,parameter  :: all_lanes=Z'ffffffff'
        integer  ,parameter  :: BLOCK_DIM=128
        integer  ,parameter  :: PME_BLOCK_DIM=64
        integer  ,parameter  :: PME_BLOCK_DIM1=32
        integer  ,constant   :: nproc, ndir
        real(t_p),constant :: xcell,ycell,zcell,i_xcell,i_ycell,i_zcell
     &           ,xcell2,ycell2,zcell2,eps_cell,box34
        logical  ,constant   :: octahedron

        contains

        attributes(global) subroutine disp_cell_cu()
        implicit none
        if (threadIdx%x.eq.1) then
           print*,xcell,ycell,zcell
           print*,xcell2,ycell2,zcell2
        end if
        end subroutine

        subroutine copy_cell_cu(xcell_,ycell_,zcell_,
     &             xcell2_,ycell2_,zcell2_,eps_cell_,
     &                            octahedron_,box34_)
        implicit none
        real(t_p),intent(in),value::xcell_,ycell_,zcell_,
     &     xcell2_,ycell2_,zcell2_,eps_cell_,box34_
        logical  ,intent(in),value::octahedron_

        !sets the default stream for all subsequent high-level CUDA Fortran operations
        !ierr = sudaforSetDefaultStream(rec_stream)

        xcell   = xcell_
        ycell   = ycell_
        zcell   = zcell_
        xcell2  = xcell2_
        ycell2  = ycell2_
        zcell2  = zcell2_
        eps_cell= eps_cell_
        i_xcell = 1/xcell_
        i_ycell = 1/ycell_
        i_zcell = 1/zcell_
        ! set octahedron data
        octahedron = octahedron_
        box34   = box34_
        !call disp_cell_cu<<<1,1>>>()
        end subroutine

        subroutine check_launch_kernel(msg)
        implicit none
        character(*),optional,intent(in):: msg
        integer ierrSync

        if (tinkerdebug.gt.0) then
           ierrSync = cudaDeviceSynchronize()
        else
           ierrSync = cudaGetLastError()
        end if

 125    format("Error ",I5," made launching kernel ",/,3x,A)
 126    format("Error ",I5," made launching kernel ",A
     &         ,/,3x,A)

        if (ierrSync.ne.cudaSuccess) then
           if (present(msg)) then
              write(*,126) ierrSync,msg,cudaGetErrorString(ierrSync)
           else
              write(*,125) ierrSync,cudaGetErrorString(ierrSync)
           end if
        end if
        end subroutine

        subroutine copy_data_to_cuda_env(data,flag)
        implicit none
        integer,intent(in),value::data
        integer,intent(in),optional::flag
        nproc = data
        if (present(flag)) then
           if (flag.eq.0) ndir = data
        end if
        end subroutine

      end module
#else
      ! For portability purpose
      subroutine void_utilcu
      end
#endif
