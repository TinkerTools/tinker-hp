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
#include "tinker_types.h"
      module utilcu
#if TINKER_SINGLE_PREC + TINKER_MIXED_PREC
        use libm ,only: f_sign=>copysignf, f_abs=>fabsf
     &           , f_sqrt=>sqrtf, f_floor=>floorf, f_erfc=>erfcf
     &           , f_exp=>expf
        use cudadevice ,only: f_inv=>__frcp_rn
#  ifdef USE_DETERMINISTIC_REDUCTION
     &                 , __tp2ll_rz=>__float2ll_rz
     &                 , __rp2ll_rz=>__double2ll_rz
#  endif
#else
        use libm ,only: f_sign=>copysign, f_abs=>fabs
     &           , f_floor=>floor, f_erfc=> erfc
        use cudadevice ,only: f_inv=>__drcp_rn
#endif
        use cudadevice ,only: __ll_as_tp=>__longlong_as_double
        use cudafor
        use sizes ,only: tinkerdebug
        implicit none
        integer  ,parameter :: all_lanes=Z'ffffffff'
        integer  ,parameter :: BLOCK_DIM=128
        integer  ,parameter :: PME_BLOCK_DIM=64  ! TODO Remove PME_BLOCK_DIM[1]
        integer  ,parameter :: PME_BLOCK_DIM1=32
        integer  ,parameter :: PME_GRID_BDIM=PME_BLOCK_DIM1
        integer  ,parameter :: PME_FPHI_BDIM=PME_BLOCK_DIM
        integer  ,parameter :: VDW_BLOCK_DIM=BLOCK_DIM
        integer  ,parameter :: TRP_BLOCK_DIM=64
        integer  ,constant  :: nproc, ndir
        real(t_p),constant  :: xcell,ycell,zcell
     &           ,i_xcell,i_ycell,i_zcell
     &           ,xcell2,ycell2,zcell2,mincell2,eps_cell,box34
        logical  ,constant  :: octahedron
        logical  ,constant  :: use_virial
        logical  ,constant  :: skipvdw12
        logical  ,constant  :: vcouple
        logical  ,constant  :: balanced_comput
        integer,allocatable,device :: utilcu_buffer(:)
        integer,allocatable,device :: utilcu_buffer_i(:)

        interface dmem_set
          module procedure dmem_set_i4
        end interface

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
        mincell2= min(xcell2_,ycell2_,zcell2_)
        eps_cell= eps_cell_
        i_xcell = 1/xcell_
        i_ycell = 1/ycell_
        i_zcell = 1/zcell_
        ! set octahedron data
        octahedron = octahedron_
        box34   = box34_
        !call disp_cell_cu<<<1,1>>>()
        end subroutine

        subroutine cu_update_vir_switch(use_vir)
        implicit none
        logical,intent(in):: use_vir
        use_virial = use_vir
        end subroutine

        subroutine cu_update_vcouple(vcouple_)
        implicit none
        integer,intent(inout)::vcouple_
        if      (vcouple_.eq.1) then
           vcouple = .true.
        else if (vcouple_.eq.0) then
           vcouple = .false.
        else
           write(*,*) " WARNING !! Unusual value of vcouple",vcouple_
           write(*,*) " Disabling van der Waals coupling method"
           vcouple  = .false.
           vcouple_ = 0
        end if
        end subroutine

        subroutine cu_update_balanced_comput(commdir_)
        implicit none
        logical,intent(in)::commdir_
        balanced_comput=commdir_
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

        subroutine chk_cuAPI_O(i,msg)
        implicit none
        integer,intent(in):: i
        character(*),intent(in):: msg
        integer ierrSync

 126    format(A,' call return error',I5,/,3x,A)

        if (i.ne.cudaSuccess)
     &     write(0,126) msg,i,cudaGetErrorString(ierrSync)
        if (tinkerdebug.gt.0) ierrSync = cudaDeviceSynchronize()
        end subroutine

        subroutine cu_update_skipvdw12(skipvdw12_)
        logical,intent(in):: skipvdw12_
        skipvdw12 = skipvdw12_
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

      subroutine mem_download_i(dst,src,stream)
      integer(4),device:: src
      integer(4),intent(out):: dst
      integer(cuda_stream_kind),intent(in):: stream
      integer ierr

      ierr = cudaMemCpyAsync(dst,src,1,cudaMemcpyDeviceToHost,stream)
      if (ierr.ne.cudasuccess) then
         write(*,*) 'error',ierr,'detected in mem_down_i4'
         write(*,*) cudageterrorstring(ierr)
      end if
      end subroutine

      subroutine prdmem_request(d_array,n,opt,stream)
      implicit none
      integer,allocatable,device:: d_array(:)
      integer n, opt
      integer(cuda_stream_kind) :: stream
      integer ierr

      if (.not.allocated(d_array)) then
         !ierr = cudamalloc(d_array,n)
         allocate( d_array(n),stat=ierr )
         call chk_cuAPI_O(ierr,'cudaMalloc')
      else if (n.gt.size(d_array)) then
         !call chk_cuAPI_O(cudaFree(d_array),'cudaFree')
         !ierr = cudamalloc(d_array,n)
         deallocate(d_array,stat=ierr )
         allocate( d_array(n),stat=ierr )
         call chk_cuAPI_O(ierr,'cudaMalloc')
      end if

      end subroutine

      subroutine dmem_set_i4(dst,val,n,stream,offset)
      integer(4),device    :: dst(*)
      integer(4),intent(in):: val
      integer(cuda_stream_kind),intent(in):: n
      integer(cuda_stream_kind),intent(in):: stream
      integer(cuda_stream_kind),optional:: offset
      integer(cuda_stream_kind) offset_
      integer ierr
 14   format(A,I4,A,/,4X,A)

      offset_ = 1
      if (present(offset)) offset_=offset
      ierr = cudaMemsetAsync(dst(offset_),val,n,stream)
      call chk_cuAPI_O(ierr,'utilcu_cuMemset')
      end subroutine

        attributes(global)
     &      subroutine transpose_z3cl(v_in,v_out,width,sloc,loc)
        implicit none
        integer,parameter:: nd=3,TRP_BLOCK_DIM=64,tdim=3
        integer,value:: width,sloc
        integer,device:: loc(*)
        real(t_p),device:: v_in(width*nd),v_out(nd*width)
        integer xIdx,yIdx,idxi,idxo,offset,i1,xIdx_,yIdx_
        real(t_p),shared:: tile(tdim*TRP_BLOCK_DIM)

        xIdx = (blockIdx%x-1)*TRP_BLOCK_DIM + threadIdx%x
        yIdx = (blockIdx%y-1)*nd            + threadIdx%y

        if( xIdx.le.sloc .and. yIdx.le.nd ) then
           idxi = ( yIdx-1 )*width + xIdx            ! v_in index
           idxo = (threadIdx%x-1)*tdim+threadIdx%y   ! shared memory Index
           tile ( idxo ) = v_in( idxi )
           v_in ( idxi ) = 0
        end if
        call syncthreads

        ! Map transposed index to v_out
        offset = ( blockIdx%x-1)*TRP_BLOCK_DIM*nd
        idxi   = (threadIdx%y-1)*TRP_BLOCK_DIM + threadIdx%x
        yIdx   = (offset+idxi-1)/nd+1
        if (yIdx.le.sloc) then
           yIdx   = loc(yIdx)
           xIdx   = mod(idxi-1,nd)+1
           idxo   = (yIdx-1)*nd + xIdx

           if( xIdx.le.nd .and. yIdx.le.width ) then
              if( nd.ne.tdim ) then
                 ! Map idxi to shared memory index
                 yIdx = (idxi-1)/nd
                 xIdx = mod(idxi-1,nd)+1
                 idxi = yIdx*tdim + xIdx
              end if
              v_out(idxo) = v_out(idxo) + tile(idxi)
           end if
        end if
        end subroutine

        attributes(global)
     &      subroutine transpose_z3ml(v_in,v_out,width,sloc,loc)
        implicit none
        integer,parameter:: nd=3,TRP_BLOCK_DIM=64,tdim=3
        integer,value:: width,sloc
        integer,device:: loc(*)
        real(r_p),device:: v_in(width*nd),v_out(nd*width)
        integer xIdx,yIdx,idxi,idxo,offset,i1,xIdx_,yIdx_
        real(r_p),shared:: tile(tdim*TRP_BLOCK_DIM)

        xIdx = (blockIdx%x-1)*TRP_BLOCK_DIM + threadIdx%x
        yIdx = (blockIdx%y-1)*nd            + threadIdx%y

        if( xIdx.le.sloc .and. yIdx.le.nd ) then
           idxi = ( yIdx-1 )*width + xIdx            ! v_in index
           idxo = (threadIdx%x-1)*tdim+threadIdx%y   ! shared memory Index
           tile ( idxo ) = v_in( idxi )
           !if (blockIdx%x.eq.gridDim%x) then
           !   print*, yIdx,xIdx,idxi,idxo
           !end if
           v_in ( idxi ) = 0
        end if
        call syncthreads

        ! Map transposed index to v_out
        offset = ( blockIdx%x-1)*TRP_BLOCK_DIM*nd
        idxi   = (threadIdx%y-1)*TRP_BLOCK_DIM + threadIdx%x
        !idxo    = offset + idxi
        yIdx   = (offset+idxi-1)/nd+1
        if (yIdx.le.sloc) then
           yIdx   = loc(yIdx)
           xIdx   = mod(idxi-1,nd)+1
           idxo   = (yIdx-1)*nd + xIdx

           if( xIdx.le.nd .and. yIdx.le.width ) then
              if( nd.ne.tdim ) then
                 ! Map idxi to shared memory index
                 yIdx = (idxi-1)/nd
                 xIdx = mod(idxi-1,nd)+1
                 idxi = yIdx*tdim + xIdx
              end if
              !if (blockIdx%x.eq.gridDim%x) then
              !   print*, yIdx,xIdx,idxi,idxo,'in'
              !end if
              v_out(idxo) = v_out(idxo) + tile(idxi)
           end if
        end if
        end subroutine

        attributes(global)
     &      subroutine transpose_z3fl(v_in,v_out,width,sloc,loc)
        implicit none
        integer,parameter:: nd=3,TRP_BLOCK_DIM=64,tdim=3
        integer,value:: width,sloc
        integer,device:: loc(*)
        mdyn_rtyp,device:: v_in(width*nd),v_out(nd*width)
        integer xIdx,yIdx,idxi,idxo,offset,i1,xIdx_,yIdx_
        mdyn_rtyp,shared:: tile(tdim*TRP_BLOCK_DIM)

        xIdx = (blockIdx%x-1)*TRP_BLOCK_DIM + threadIdx%x
        yIdx = (blockIdx%y-1)*nd            + threadIdx%y

        if( xIdx.le.sloc .and. yIdx.le.nd ) then
           idxi = ( yIdx-1 )*width + xIdx            ! v_in index
           idxo = (threadIdx%x-1)*tdim+threadIdx%y   ! shared memory Index
           tile ( idxo ) = v_in( idxi )
           v_in ( idxi ) = 0
        end if
        call syncthreads

        ! Map transposed index to v_out
        offset = ( blockIdx%x-1)*TRP_BLOCK_DIM*nd
        idxi   = (threadIdx%y-1)*TRP_BLOCK_DIM + threadIdx%x
        yIdx   = (offset+idxi-1)/nd+1
        if (yIdx.le.sloc) then
           yIdx   = loc(yIdx)
           xIdx   = mod(idxi-1,nd)+1
           idxo   = (yIdx-1)*nd + xIdx

           if( xIdx.le.nd .and. yIdx.le.width ) then
              if( nd.ne.tdim ) then
                 ! Map idxi to shared memory index
                 yIdx = (idxi-1)/nd
                 xIdx = mod(idxi-1,nd)+1
                 idxi = yIdx*tdim + xIdx
              end if
              v_out(idxo) = v_out(idxo) + tile(idxi)
           end if
        end if
        end subroutine

        attributes(global)
     &      subroutine transpose_z3f(v_in,v_out,width)
        implicit none
        integer,parameter:: nd=3,TRP_BLOCK_DIM=64,tdim=3
        integer,value:: width
        mdyn_rtyp,device:: v_in(width*nd),v_out(nd*width)
        integer xIdx,yIdx,idxi,idxo,offset,i1,xIdx_,yIdx_
        mdyn_rtyp,shared:: tile(tdim*TRP_BLOCK_DIM)

        xIdx = (blockIdx%x-1)*TRP_BLOCK_DIM + threadIdx%x
        yIdx = (blockIdx%y-1)*nd            + threadIdx%y

        if( xIdx.le.width .and. yIdx.le.nd ) then
           idxi = ( yIdx-1 )*width + xIdx            ! v_in index
           idxo = (threadIdx%x-1)*tdim+threadIdx%y   ! shared memory Index
           tile ( idxo ) = v_in( idxi )
           v_in ( idxi ) = 0
        end if
        call syncthreads

        ! Map transposed index to v_out
        offset = ( blockIdx%x-1)*TRP_BLOCK_DIM*nd
        idxi   = (threadIdx%y-1)*TRP_BLOCK_DIM + threadIdx%x
        yIdx   = (offset+idxi-1)/nd+1
        xIdx   = mod(idxi-1,nd)+1
        idxo   = (yIdx-1)*nd + xIdx

        if( xIdx.le.nd .and. yIdx.le.width ) then
           if( nd.ne.tdim ) then
              ! Map idxi to shared memory index
              yIdx = (idxi-1)/nd
              xIdx = mod(idxi-1,nd)+1
              idxi = yIdx*tdim + xIdx
           end if
           v_out(idxo) = v_out(idxo) + tile(idxi)
        end if
        end subroutine

        attributes(global)
     &      subroutine transpose_z3m(v_in,v_out,width)
        implicit none
        integer,parameter:: nd=3,TRP_BLOCK_DIM=64,tdim=3
        integer,value:: width
        real(r_p),device:: v_in(width*nd),v_out(nd*width)
        integer xIdx,yIdx,idxi,idxo,offset,i1,xIdx_,yIdx_
        real(r_p),shared:: tile(tdim*TRP_BLOCK_DIM)

        xIdx = (blockIdx%x-1)*TRP_BLOCK_DIM + threadIdx%x
        yIdx = (blockIdx%y-1)*nd            + threadIdx%y

        if ( xIdx.le.width .and. yIdx.le.nd ) then
           idxi = ( yIdx-1 )*width + xIdx            ! v_in index
           idxo = (threadIdx%x-1)*tdim+threadIdx%y   ! shared memory Index
           tile ( idxo ) = v_in( idxi )
           v_in ( idxi ) = 0
        end if
        call syncthreads

        ! Map transposed index to v_out
        offset = ( blockIdx%x-1)*TRP_BLOCK_DIM*nd
        idxi   = (threadIdx%y-1)*TRP_BLOCK_DIM + threadIdx%x
        yIdx   = (offset+idxi-1)/nd+1
        xIdx   = mod(idxi-1,nd)+1
        idxo   = (yIdx-1)*nd + xIdx

        if ( xIdx.le.nd .and. yIdx.le.width ) then
           if ( nd.ne.tdim ) then
              ! Map idxi to shared memory index
              yIdx = (idxi-1)/nd
              xIdx = mod(idxi-1,nd)+1
              idxi = yIdx*tdim + xIdx
           end if
           v_out(idxo) = v_out(idxo) + tile(idxi)
        end if
        end subroutine

      end module
#else
      ! For portability purpose
      subroutine void_utilcu
      end
#endif
