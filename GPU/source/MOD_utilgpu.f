c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ####################################################################
c     ##                                                                ##
c     ##  module utilgpu  -- GPU utility functions and params           ##
c     ##                                                                ##
c     ####################################################################
c
c     ngpus        numbers of gpus available 
c     devicenum    device attach to MPI hostrank
c     gpu_gangs    attribute #gangs to create inside a kernels
c     gpu_workers  attribute #workers for parallelism mapping
c     gpu_vector   attribute #vector_length in kernels or parallel clauses
c     dir_queue    direct space computation async queue
c     rec_queue    reciprocal space computation async queue
c     def_queue    default async queue
c     Nblock       number of k block kernel comput in direct space
c     maxscaling   maximum number of scaling factor of an atom among interactions 1-2 to 1-5
c     maxscaling1  maximum number of scaling factor of an atom among p_interactions
c     warning      warning message to print in case of trouble
c     ftc          a transformation matrix to convert a multipole object in fraction coordinates to Cartesian
c     ctf          a transformation matrix to convert a multipole object in Cartesian coordinates to fractional
c     qi1 and qi2  are useful to compute ftc and ctf
c     devProp      Properties of the device
c     nSMP         number of streaming Multi-processor on the device
c     nSPcores     number of cores od the device
c
c
#include "tinker_macro.h"
      module utilgpu
      use tinheader ,only: ti_p,re_p
      use boxes     ,only: orthogonal,octahedron,box34
      use cell
      use couple
      use chgpot
      use domdec    ,only: rank,nproc,hostrank
      use mplpot
      use polpot
      use polgrp
      use vdwpot
      use sizes
      use tinMemory
      use tinTypes
      use iso_c_binding,only:c_ptr
#ifdef _OPENACC
      use openacc
      use cudafor
#else
#define cuda_stream_kind 8
#endif
      implicit none
      integer ngpus
      integer :: devicenum=-1
      integer gpu_gangs,ngangs_rec
      integer gpu_workers,gpu_vector
      integer rec_queue,dir_queue,def_queue
      integer(mipk),private:: zero8
      logical mod_gangs
#ifdef _OPENACC
      type(cudaDeviceProp) devProp
      integer(cuda_stream_kind) dir_stream,rec_stream
      integer(cuda_stream_kind) def_stream
      type(cudaEvent) dir_event,rec_event,def_event
      integer,parameter :: zero_flag=0
#else
      integer(mipk) dir_stream,rec_stream,def_stream,dir_event
     &       ,rec_event,def_event
      parameter(dir_stream=0,rec_stream=0,def_stream=0
     &         ,dir_event=0,rec_event=0,def_event=0)
#endif
      integer nSMP,nSPcores,cores_SMP
      integer Ndir_block,Ndir_async_block,count_block
      real(t_p),parameter:: inf  =4*huge(0.0_ti_p)
      real(r_p),parameter:: inf_r=4*huge(0.0_re_p)

      integer,parameter:: maxscaling=64
      integer,parameter:: maxscaling1=256
      integer,parameter:: maxBlock=ishft(1,16)
      integer,parameter:: BLOCK_SIZE=32
      integer,parameter::  WARP_SIZE=32
      integer,parameter::  WARP_SHFT=5
      integer,parameter::RED_BUFF_SIZE=ishft(1,15)

      integer,parameter,dimension(6)::qi1=(/ 1, 2, 3, 1, 1, 2 /)
      integer,parameter,dimension(6)::qi2=(/ 1, 2, 3, 2, 3, 3 /)

      integer   nred_buff(  RED_BUFF_SIZE)
      real(t_p) vred_buff(6*RED_BUFF_SIZE)
      real(t_p) vred_buf1(6*RED_BUFF_SIZE)
      real(r_p) ered_buff(  RED_BUFF_SIZE)
      real(r_p)  lam_buff(  RED_BUFF_SIZE)
      ener_rtyp ered_buf1(  RED_BUFF_SIZE)
      real(t_p),dimension(10,10)::ftc,ctf

      ! WorkSpace buffer to be used across the code
      integer  ,allocatable:: ug_workS_i(:),ug_workS_i1(:)
      real(t_p),allocatable:: ug_workS_r(:),ug_workS_r1(:)
      real(r_p),allocatable:: ug_workS_m(:),ug_workS_m1(:)
      ! -----------------------------------
      ! ug_workS_r   is used within
      !              - tmatxb_pme_core3
      !              - empole1cgpu, epreal1cgpu, emreal3d_cu, epreal3d_cu, emrecip1gpu, eprecip1gpu
      ! -----------------------------------
      ! ug_workS_r1  is used within
      !              - emrecip1gpu, eprecip1gpu
      !

      character*128 warning,sugest_vdw

      parameter( zero8=0 )

#ifdef _OPENACC
      interface
        function acc_malloc_bc(bytes) result(ptr) 
     &                  bind(C,name='acc_malloc')
        import C_PTR
        integer,value,intent(in )::bytes
        type(C_PTR)  ,intent(out)::ptr
        end function
      end interface

      interface
      function acc_get_cuda_stream_c( async ) result(queue)
     &     bind(C,name='acc_get_cuda_stream')
      import c_int
      integer(c_int),value:: async
      integer(c_int),intent(out):: queue
      end function
      end interface

      interface
        subroutine acc_free_bc(ptr) bind(C,name='acc_free')
        import C_PTR
        type(C_PTR),value ::ptr
        end subroutine
      end interface

      interface dev_sort
        module procedure dev_sort_i4
      end interface
#endif

      interface Tinker_shellEnv
        module procedure Tinker_shellEnv_i
        module procedure Tinker_shellEnv_f4
        module procedure Tinker_shellEnv_f8
      end interface

      interface mem_move
         module procedure mem_move_i4
         module procedure mem_move_i8
         module procedure mem_move_r4
         module procedure mem_move_r4_2d
         module procedure mem_move_r8
         module procedure mem_move_r8_2d
      end interface

      interface mem_set
         module procedure mem_set_i4
         module procedure mem_set_i8
         module procedure mem_set_r4
         module procedure mem_set_r4_2d
         module procedure mem_set_r4_3d
         module procedure mem_set_r4_5d
         module procedure mem_set_r8
         module procedure mem_set_r8_2d
         module procedure mem_set_r8_3d
         module procedure mem_set_r8_5d
      end interface

      interface reduce_buffer
        module procedure reduce_buffer_m
      end interface

      interface reduce_energy_virial
        module procedure reduce_energy_virial0
        module procedure reduce_energy_virial1
        module procedure reduce_energy_virial2
      end interface

      interface reduce_energy_action
         module procedure reduce_energy_action0
         module procedure reduce_energy_action1
      end interface

      interface transpose_az3
         module procedure transpose_az_l3
         module procedure transpose_az_f3
         module procedure transpose_az_d3
      end interface
c
!$acc declare create(gpu_gangs,gpu_workers,gpu_vector,devicenum,
!$acc& ngpus,warning,sugest_vdw,ftc,ctf)
c
!$acc declare copyin(qi1,qi2)

      contains

c     ########################################################################
c     ##                                                                    ##
c     ##  subroutine imagevec_acc  --  compute the minimum image distance  ##
c     ##                                                                    ##
c     ########################################################################
c
c     "imagevec3" takes the components of pairwise distances between
c     two points in a periodic box and converts to the components
c     of the minimum image distances. Indice i designs x, y or z
c     direction on a graphics card
c
c     xcell    length of the a-axis of the complete replicated cell
c     ycell    length of the b-axis of the complete replicated cell
c     zcell    length of the c-axis of the complete replicated cell
c     xcell2   half the length of the a-axis of the replicated cell
c     ycell2   half the length of the b-axis of the replicated cell
c     zcell2   half the length of the c-axis of the replicated cell
      pure subroutine imagevec_acc(xpos2,ypos2,zpos2,n,sizepos2) 
!$acc routine vector
      implicit none
      integer,intent(in):: n,sizepos2
      real(t_p), intent(inout):: xpos2(sizepos2)
      real(t_p), intent(inout):: ypos2(sizepos2)
      real(t_p), intent(inout):: zpos2(sizepos2)
      integer j
      real(t_p) cel
      
!$acc loop vector
      do j=1,n
         if (abs(xpos2(j)) .gt. xcell2) then
            cel      = sign(xcell,xpos2(j))
            xpos2(j) = xpos2(j) 
     &               - cel*floor((abs(xpos2(j))+xcell2)/xcell)
         end if
         if (abs(ypos2(j)) .gt. ycell2) then
            cel      = sign(ycell,ypos2(j))
            ypos2(j) = ypos2(j) 
     &               - cel*floor((abs(ypos2(j))+ycell2)/ycell)
         end if
         if (abs(zpos2(j)) .gt. zcell2) then
            cel      = sign(zcell,zpos2(j))
            zpos2(j) = zpos2(j) 
     &               - cel*floor((abs(zpos2(j))+zcell2)/zcell)
         end if
      end do
      end
c
c
c
c     ###################################################################
c     ##                                                               ##
c     ##  subroutine setscale2_acc  --  compute the scaling factors    ##
c     ##                                                               ##
c     ###################################################################
c
c
c     "setscale2" takes the indice of the reference atom, the global indices
c     of the neighbours, the scaletype and returns the scaling factors

      subroutine setscale2_acc(indice,iscal,scal,nscal,
     &                     nn12,nn13,nn14,ntot,
     &                     kvec,nb,scaletype)
!$acc routine vector
      implicit none
      integer, intent(in)   :: indice,nb,nscal
      integer, intent(in)   :: nn12,nn13,nn14,ntot
      integer, intent(in)   :: scaletype
      integer, intent(in)   :: kvec(nscal)
      integer, intent(inout):: iscal(nscal)
      real(t_p) , intent(inout):: scal(nscal)
      integer i,j,k,nnp11
      real(t_p)  scale2,scale3,scale4,scale41,scale5,scale_t
      real(t_p),parameter :: zero=0.0_ti_p, one=1.0_ti_p
c     character(len=1) :: scaletype_temp
c     scaletype_temp = TRANSFER(scaletype,MOLD=scaletype_temp)
cdeb  if(rank.eq.0.and.tinkerdebug) print*, 'setscale2_acc'
      if(nb>nscal-ntot) then
         print*, 'Data collision in scal array\n',
     &        'Fix setscale2_acc by extending scal size'
      endif
c 
c     select scaling coefficients
c
      select case (scaletype)
         case(0)
c           scale2  = c2scale
c           scale3  = c3scale
c           scale4  = c4scale
c           scale5  = c5scale
         case(1)
            scale2  = m2scale
            scale3  = m3scale
            scale4  = m4scale
            scale5  = m5scale
         case(2)
            scale2  = p2scale
            scale3  = p3scale
            scale4  = p4scale
            scale41 = p41scale
            scale5  = p5scale
         case(3)
            scale2  = v2scale
            scale3  = v3scale
            scale4  = v4scale
            scale5  = v5scale
         case default
            scale2  = 1.0_ti_p
            scale3  = 1.0_ti_p
            scale4  = 1.0_ti_p
            scale5  = 1.0_ti_p
      endselect
c
c     get list of atoms indices and there scaling factor
c     give scalinf factors to selected atoms
!$acc loop vector
      do i=1,nn12
         iscal(i) = i12 (i,indice)
          scal(nscal+1-i) = scale2
      end do
!$acc loop vector
      do i=nn12+1,nn13
         iscal(i) = i13 (i-nn12,indice)
          scal(nscal+1-i) = scale3
      end do
!$acc loop vector
      do i=nn13+1,nn14
         iscal(i) = i14 (i-nn13,indice)
          scal(nscal+1-i) = scale4
      end do
!$acc loop vector
      do i=nn14+1,ntot
         iscal(i) = i15 (i-nn14,indice)
          scal(nscal+1-i) = scale5
      end do
!$acc loop vector
      do i=ntot+1,ntot+ntot
         iscal(i) = zero
      end do

!$acc loop vector
      do i=1,nb
         scal(i) = one
      end do
c
c     Find where additional factor for 1-4 intragroup polarization 
c     should be taken into account
      if (scaletype==2) then
         nnp11 = np11(indice)
!$acc loop vector
         do i=1,nnp11
            j = ip11(i,indice)
!$acc loop seq
            do k=nn13+1,nn14
               if (iscal(k).eq.j) then
                  scal(nscal+1-k) = scale4*scale41
               endif
            end do
         end do
      end if
c
c     for all indices in nblist, select the corresponding atom
c     and apply his scaling factor in kvec indexation
c
!$acc loop vector
      do i=1,nb
         j = kvec(i)
!$acc loop seq
         do k=1,ntot
            if(iscal(k).eq.j) then
               iscal(ntot+k) = i
            end if
         end do
      end do
c
c     Replace scaling factor along kvec index
c
!$acc loop vector
      do i=1,ntot
         j=iscal(ntot+i)
         if (j/=0) scal(j)=scal(nscal+1-i)
      end do
      end
c
c
c
      subroutine fillscale_acc(ind,nn12,nn13,nn14,ntot,
     &             scale2,scale3,scale4,scale5,
     &             iscal,fscal,tscal,
     &             scale41)
!$acc routine vector
      implicit none
      integer,value, intent(in) :: ind,nn12,nn13,nn14,ntot
      real(t_p),value,intent(in)::scale2,scale3,scale4,scale5
      real(t_p),value,intent(in),optional:: scale41
      integer, intent(inout) :: iscal(64)
      integer, intent(inout), optional :: tscal(64)
      real(t_p) , intent(inout) :: fscal(64)

      integer i,j,k,nnp11

      if(ntot>64) print*, 'scaling array too short'
c 
c     pack scaling coefficients
c
      if (present(tscal)) then
!$acc loop vector
         do j=1,nn13
            if      (j.le.nn12) then
               iscal(j) = i12 (j, ind)
               fscal(j) = scale2
               tscal(j) = 2
            else
               iscal(j) = i13 (j-nn12,ind)
               fscal(j) = scale3
               tscal(j) = 4
            end if
         end do
!$acc loop vector
         do j=nn13+1,ntot
            if      (j.le.nn14) then
               iscal(j) = i14 (j-nn13,ind)
               fscal(j) = scale4
               tscal(j) = 8
            else
               iscal(j) = i15 (j-nn14,ind)
               fscal(j) = scale5
               tscal(j) = 16
            end if
         end do
      else
!$acc loop vector
         do j=1,nn13
            if      (j.le.nn12) then
               iscal(j) = i12 (j, ind)
               fscal(j) = scale2
            else
               iscal(j) = i13 (j-nn12,ind)
               fscal(j) = scale3
            end if
         end do
!$acc loop vector
         do j=nn13+1,ntot
            if      (j.le.nn14) then
               iscal(j) = i14 (j-nn13,ind)
               fscal(j) = scale4
            else
               iscal(j) = i15 (j-nn14,ind)
               fscal(j) = scale5
            end if
         end do
      end if
c
c     Find where additional factor for 1-4 intragroup polarization 
c     should be taken into account
c
      if (present(scale41)) then
         nnp11 = np11(ind)
!$acc loop vector
         do i=1,nnp11
            j = ip11(i,ind)
!$acc loop seq
            do k=nn13+1,nn14
               if (iscal(k).eq.j) then
                  fscal(k) = scale4*scale41
               endif
            end do
         end do
      end if

      end subroutine

 
c     ####################################################################
c     ##                                                                ##
c     ##  subroutine setscale2p_acc  --  compute the scaling factors    ##
c     ##                                                                ##
c     ####################################################################
c
c
c     "setscale2" takes the indice of the reference atom, the global indices
c     of the neighbours, the scaletype and returns the scaling factors

      subroutine setscale2p_acc(indice,iscal,scal,nscal,
     &                     nnp11,nnp12,nnp13,nnp14,
     &                     kvec,nb,scaletype)
!$acc routine vector
      implicit none
      integer, intent(in)   :: indice,nb,nscal
      integer, intent(in)   :: nnp11,nnp12,nnp13,nnp14
      integer, intent(in)   :: scaletype
      integer, intent(in)   :: kvec(nscal)
      integer, intent(inout):: iscal(nscal)
      real(t_p) , intent(inout):: scal(nscal)
      integer i,j,k
      real(t_p)  scale1,scale2,scale3,scale4
cdeb  if(rank.eq.0.and.tinkerdebug) print*, 'setscale2p_acc'
      if(nb>nscal-nnp14) then
         print*, 'Data collision in scal array\n',
     &        'Fix setscale2p_acc by extending scal size'
      endif
c 
c 
c     select scaling coefficients
c
      select case (scaletype)
         case(0)
            scale1  = d1scale
            scale2  = d2scale
            scale3  = d3scale
            scale4  = d4scale
         case(1)
            scale1  = u1scale
            scale2  = u2scale
            scale3  = u3scale
            scale4  = u4scale
         case default
            scale1  = 1.0_ti_p
            scale2  = 1.0_ti_p
            scale3  = 1.0_ti_p
            scale4  = 1.0_ti_p
      endselect
c
c     get list of atoms indices and there scaling factor
c     give scalinf factors to selected atoms
!$acc loop vector
      do i=1,nnp11
         iscal(i) = ip11 (i,indice)
          scal(nscal+1-i) = scale1
      end do
!$acc loop vector
      do i=nnp11+1,nnp12
         iscal(i) = ip12 (i-nnp11,indice)
          scal(nscal+1-i) = scale2
      end do
!$acc loop vector
      do i=nnp12+1,nnp13
         iscal(i) = ip13 (i-nnp12,indice)
          scal(nscal+1-i) = scale3
      end do
!$acc loop vector
      do i=nnp13+1,nnp14
         iscal(i) = ip14 (i-nnp13,indice)
          scal(nscal+1-i) = scale4
      end do
!$acc loop vector
      do i=nnp14+1,2*nnp14
         iscal(i) = 0
      end do
!$acc loop vector
       do i=1,nb
          scal(i) = 1.0_ti_p
       end do
c
c     for all indices in the lists, select the corresponding atom
c
!$acc loop vector
      do i=1,nb
         j = kvec(i)
!$acc loop seq
         do k=1,nnp14
            if(iscal(k).eq.j) iscal(nnp14+k) = i
         end do
      end do
c
c     Replace scaling factor along kvec index
c
!$acc loop vector
      do i=1,nnp14
         j=iscal(nnp14+i)
         if (j/=0) scal(j)=scal(nscal+1-i)
      end do
      end
c
c
      subroutine fillscalep_acc(ind,nnp11,nnp12,nnp13,nnp14,
     &                          scale1,scale2,scale3,scale4,
     &                          iscal,fscal,tscal)
!$acc routine vector
      implicit none
      integer, intent(in):: ind,nnp11,nnp12,nnp13,nnp14
      real(t_p), intent(in):: scale1,scale2,scale3,scale4
c     character*1, intent(in) :: scaletype
      integer, intent(inout) :: iscal(:)
      integer, intent(inout),optional :: tscal(:)
      real(t_p), intent(inout) :: fscal(:)

      integer i,j,k

      if(nnp14>96) print*, 'scaling array too short'
c 
c     select scaling coefficients
c
c     select case (scaletype)
c        case('d')
c           scale1  = d1scale
c           scale2  = d2scale
c           scale3  = d3scale
c           scale4  = d4scale
c        case('u')
c           scale1  = u1scale
c           scale2  = u2scale
c           scale3  = u3scale
c           scale4  = u4scale
c        case default
c           scale1  = 1.0_ti_p
c           scale2  = 1.0_ti_p
c           scale3  = 1.0_ti_p
c           scale4  = 1.0_ti_p
c     endselect
c
c     get list of atoms indices and there scaling factor
c     give scalinf factors to selected atoms
      if (present(tscal)) then
!$acc loop vector
         do i=1,nnp12
            if (i.le.nnp11) then
               iscal(i) = ip11 (i,ind)
               fscal(i) = scale1
               tscal(i) = 1
            else
               iscal(i) = ip12 (i-nnp11,ind)
               fscal(i) = scale2
               tscal(i) = 2
            end if
         end do
!$acc loop vector
         do i=nnp12+1,nnp14
            if (i.le.nnp13) then
               iscal(i) = ip13 (i-nnp12,ind)
               fscal(i) = scale3
               tscal(i) = 4
            else
               iscal(i) = ip14 (i-nnp13,ind)
               fscal(i) = scale4
               tscal(i) = 8
            end if
         end do
      else
!$acc loop vector
         do i=1,nnp12
            if (i.le.nnp11) then
               iscal(i) = ip11 (i,ind)
               fscal(i) = scale1
            else
               iscal(i) = ip12 (i-nnp11,ind)
               fscal(i) = scale2
            end if
         end do
!$acc loop vector
         do i=nnp12+1,nnp14
            if (i.le.nnp13) then
               iscal(i) = ip13 (i-nnp12,ind)
               fscal(i) = scale3
            else
               iscal(i) = ip14 (i-nnp13,ind)
               fscal(i) = scale4
            end if
         end do
      end if

      end subroutine

      subroutine Tinker_shellEnv_i(envName,envValue,defVal)
      implicit none
      character(*),intent(in):: envName
      integer,intent(in) :: defVal
      integer,intent(out):: envValue
      character*64 value
      integer length,status

      call get_environment_variable('TINKER_'//envName,value,
     &                               length,status)
      envValue = defVal
      if (status.eq.0) read(value,*) envValue
      end subroutine

      subroutine Tinker_shellEnv_f8(envName,envValue,defVal)
      implicit none
      character(*),intent(in):: envName
      real(8),intent(in) :: defVal
      real(8),intent(out):: envValue
      character*64 value
      integer length,status

      call get_environment_variable('TINKER_'//envName,value,
     &                               length,status)
      envValue = defVal
      if (status.eq.0) read(value,*) envValue
      end subroutine

      subroutine Tinker_shellEnv_f4(envName,envValue,defVal)
      implicit none
      character(*),intent(in):: envName
      real(4),intent(in) :: defVal
      real(4),intent(out):: envValue
      character*64 value
      integer length,status

      call get_environment_variable('TINKER_'//envName,value,
     &                               length,status)
      envValue = defVal
      if (status.eq.0) read(value,*) envValue
      end subroutine
c
c     updating the number of gangs to launch on the gpu if authorized
c
      subroutine update_gang(gang_)
      implicit none
      integer gang_
      if (mod_gangs.and.ngpus>0) gpu_gangs = gang_
      end
c
c
      subroutine update_gang_rec(natoms)
      implicit none
      integer natoms
      ngangs_rec = nSMP
      end subroutine
c
c
      subroutine set_warning
      implicit none

      warning ='Warning, system moved too much since last neighbor list'
     &       //'update, try lowering nlupdate '
      sugest_vdw = 'VDW'
!$acc enter data copyin(warning,sugest_vdw)
      end
c
c
c
      subroutine switch_queue
      implicit none

      if (def_queue.eq.rec_queue) then
         def_queue = dir_queue
      else
         def_queue = rec_queue
      end if
      end
c
c
      subroutine reset_def_to_dir
      implicit none

      def_queue = dir_queue
      count_block = 0
      end subroutine
c
c
      subroutine switch_dir_queue
      implicit none

      if (count_block.lt.Ndir_async_block) then
         if (def_queue.eq.dir_queue) then
            def_queue = rec_queue
         else
            def_queue = dir_queue
            count_block = count_block+1
         end if
      else
         def_queue = dir_queue
      end if
      end

      subroutine openacc_abort(message)
      use mpi
      character(*) :: message
      integer :: ierr
#ifdef _OPENACC
      print *, "OpenACC abort : ",message
      call MPI_Abort(MPI_COMM_WORLD, -9, ierr)
#endif
      end subroutine

#ifdef _OPENACC
c
c     Allocate cache memory for Thrust cache library
c
      subroutine thrust_cache_init(n)
      use thrust,only: thrust_alloc_cache_memory
      implicit none
      integer,intent(in)::n
      logical,save:: f_in=.true.
      integer siz

      if (f_in) then
         siz   = max(250000,n*7)
         call thrust_alloc_cache_memory(siz)
         sd_prmem = sd_prmem + siz*sizeof(n)
      end if
      f_in = .false.
      end
c
c     Bind device to MPI Thread
c
      subroutine bind_gpu
      character(len=6) :: local_rank_env
      integer          :: local_rank_env_status, local_rank
      integer device_type, device_start

      !Initialisation OpenAcc
      call get_environment_variable(name="PGI_ACC_DEVICE_NUM",
     &     value=local_rank_env, status=local_rank_env_status)
      if (local_rank_env_status == 0) return

      call get_environment_variable(name="ACC_DEVICE_NUM",
     &     value=local_rank_env, status=local_rank_env_status)
      if (local_rank_env_status == 0) return

      call get_environment_variable(name="ACC_DEVICE_HOST_RANK",
     &     value=local_rank_env, status=local_rank_env_status)
      if (local_rank_env_status == 0) return

      ! Look for CUDA environment
      call get_environment_variable(name="CUDA_VISIBLE_DEVICE",
     &     value=local_rank_env, status=local_rank_env_status)
      if (local_rank_env_status == 0) then
         read(local_rank_env, *) local_rank
         goto 100
      end if

      ! Recuperation du rang local du processus via la variable d'environnement
      ! positionnee par Slurm, l'utilisation de MPI_Comm_rank n'etant pas encore
      ! possible puisque cette routine est utilisee avant l'initialisation de MPI
      call get_environment_variable(name="SLURM_LOCALID",
     &     value=local_rank_env, status=local_rank_env_status)
      if (local_rank_env_status == 0) then
         read(local_rank_env, *) local_rank
         !if (local_rank.eq.0) print*,'SLURM detected'
         goto 100
      end if

      ! OpenMPI
      call get_environment_variable(name="OMPI_COMM_WORLD_LOCAL_RANK",
     &     value=local_rank_env, status=local_rank_env_status)
      if (local_rank_env_status == 0) then
         read(local_rank_env, *) local_rank
         !if (local_rank.eq.0) print*,'OpenMPI detected'
         goto 100
      end if

      ! Intel MPI
      call get_environment_variable(name="MPI_LOCALRANKID",
     &     value=local_rank_env, status=local_rank_env_status)
      if (local_rank_env_status == 0) then
         read(local_rank_env, *) local_rank
         !if (local_rank.eq.0) print*,'Intel MPI detected'
         goto 100
      end if

 34   format( "Warning : Could not bind GPU device before MPI_Init",
     &    ,/, "          Switched to default binding procedure")
      write(0,34)
      return

 100  continue
      ! Look for device ID to bind with using OpenACC

      device_type = acc_get_device_type()
      ngpus = acc_get_num_devices(device_type)

      ! Circle on total number of device present on the host
      ! to prevent errors from `acc_set_device_num`
      if (local_rank.ge.ngpus) local_rank = mod(local_rank,ngpus)

      ! Bind to device ID local_rank
      call acc_set_device_num( local_rank,device_type )
      devicenum=local_rank
      end subroutine bind_gpu
c
c     Initiate device parameters
c
      subroutine selectDevice
      use mpi
      implicit none
      integer ::i,cuda_success=0
      character*64 value,value1,value2
      character*22 name
      integer length,length1,length2,status,status1,status2
      integer device_type
      integer hostnm,getpid
      integer device_start,host_comm,host_rank,ierr

c
c     Query number of device 
c     Link a device to a MPI process
c     If not bind by slurm local ID
c
      host_rank = hostrank
      if (devicenum==-1) then
         device_type = acc_get_device_type()
         device_start= 0
         call get_environment_variable('ACC_DEVICE_NUM',value,
     &                                  length,status)
         call get_environment_variable('PGI_ACC_DEVICE_NUM',value1,
     &                                  length1,status1)
         call get_environment_variable('ACC_DEVICE_HOST_RANK',value2,
     &                                  length2,status2)
         if (status1.eq.0) then
            ngpus = 1
            read (value1,'(i)') devicenum
         else if (status.eq.0) then
            ngpus = 1
            read (value,'(i)') devicenum
         else
            if (status2.eq.0) read(value2,*) device_start
            ngpus     = acc_get_num_devices(device_type)
            call MPI_Comm_split_type(MPI_COMM_WORLD,
     $           MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, host_comm,ierr)
            call MPI_Comm_rank(host_comm,host_rank,ierr)
            call MPI_Comm_free(host_comm,ierr)
            devicenum = mod( device_start+host_rank,ngpus )
         end if

         call acc_set_device_num(devicenum, device_type)
      end if
c
c     Query device(s) propertie(s)
c
      call getDeviceProp
c
c     Display device number and MPI process
c
      status = hostnm(name)
      if (nproc.gt.16.or.nproc.eq.1) return
      write(0,'(a,i2,a,i2,a,i8,a,a)')
     &     "rank ", host_rank,
     &     " attach to device ", devicenum, 
     &     " pid : ",getpid(), 
     &     " on ", name
      end subroutine
c
c     Initiate Device streams
c
      subroutine initDevice_stream
      implicit none
      integer ::i,cuda_success=0
      character*64 value
      integer int_val
      integer length,status,device_type
      integer high_priority,low_priority

      rec_queue  = acc_async_noval
      dir_queue  = acc_async_noval
      def_queue  = acc_async_noval
      rec_stream = acc_get_cuda_stream(acc_async_noval)
      dir_stream = rec_stream
      def_stream = rec_stream

      call Tinker_shellEnv_i("ASYNC_COVER",int_val,0)

      if (int_val.eq.0.or.(nproc.eq.1.and.int_val.eq.1)) return

      if (rec_stream.eq.0) then
 16      format(/,
     &'| WARNING !! cannot recover CUDA stream from OpenACC async',
     &' queue',/,
     &'             Skiping multi queues feature activation')
         if (rank.eq.0) write(*,16)
         return
      end if

      dir_queue = acc_async_noval + 1
 13   format(' ***** Asynchronous Computation Overlapping enable'
     &      ,I5,/)
      if (rank.eq.0) write(*,13) acc_async_noval
c
c     Create stream(s) and attach them to Acc queue
c
      cuda_success = cudaDeviceGetStreamPriorityRange(low_priority,
     &              high_priority)
 
      ! Create main stream for TINKER
      cuda_success = cuda_success +
     &  cudastreamcreate(rec_stream,cudaStreamNonBlocking,high_priority)
      call acc_set_cuda_stream(rec_queue,rec_stream)

      ! Create Secondary stream if necessary
      if (dir_queue.ne.rec_queue) then
         cuda_success = cuda_success + 
     &   cudastreamcreate(dir_stream,cudaStreamNonBlocking,low_priority)
         call acc_set_cuda_stream(dir_queue,dir_stream)
      else
         dir_stream = rec_stream
         def_stream = rec_stream
      end if

      if (cuda_success.ne.0)
     &   print*,'error creating cuda async queue >',cuda_success
c
c     Creating Event for stream syncrhonisation
c
      cuda_success = cuda_success + 
     &       cudaEventCreateWithFlags(dir_event,cudaEventDisableTiming)
      cuda_success = cuda_success + 
     &       cudaEventCreateWithFlags(rec_event,cudaEventDisableTiming)
      cuda_success = cuda_success + 
     &       cudaEventCreateWithFlags(def_event,cudaEventDisableTiming)

      if (cuda_success.ne.0)
     &   print*,'error creating cuda Events  >', cuda_success

      end

      subroutine syncWithStream(stream)
      integer,optional::stream
      integer stream_,istat

      if (present(stream)) then
         stream_= stream
      else
         stream_= rec_stream
      end if
      istat  = cudaStreamSynchronize(stream_)
      if (istat.ne.0) print '(A,I5,A,I14)',' WARNING ! Error',istat,
     &   ' synchronizing with stream ',stream_
      end subroutine

      subroutine comparetonoval(queue)
      implicit none
      integer queue
      integer(cuda_stream_kind)::stream

      stream = acc_get_cuda_stream(queue)
      if (rec_stream.ne.stream) then
  12  format(A,3I14,2I5)
         print 12,'stream are differents'
     &   ,rec_stream,dir_stream,stream,queue,rank
      else
         print*, 'equal stream',rec_stream
      end if
      end subroutine

      subroutine Finalize_async_recover
      implicit none
      integer cuda_success

      if (dir_queue.ne.rec_queue) then
         !cuda_success = cudaStreamDestroy(dir_stream)
c12      format('CUDA Error destroying rec/direct streams',I10)
c        if (cuda_success.ne.0) then
c           print 12, cuda_success
c        end if
 13     format(/,1x,'***** ',
     &        'Asynchronous Computation Overlapping disabled *****')
         if (rank.eq.0) write(*,13)
         dir_stream = rec_stream
         def_stream = rec_stream
         dir_queue  = rec_queue
         def_queue  = rec_queue
      end if

      end subroutine
c
c     A cuda implementation of _pgi_acc_wait_async
c     Synchonisation within CUDA streams doesn't seem to work with OpenACC
c
      subroutine stream_wait_async(wait_s,async_s,event_)
      implicit none
      integer(cuda_stream_kind),intent(in)::wait_s,async_s
      type(cudaEvent),optional::event_
      type(cudaEvent) :: event
      integer cuda_success

      cuda_success = 0
      event = def_event
      if (present(event_)) event=event_
         
      cuda_success = cuda_success +
     &     cudaEventRecord(event,wait_s)
      cuda_success = cuda_success +
     &     cudaStreamWaitEvent(async_s,event,zero_flag)

      if (cuda_success.ne.0)
     &   write(*,*) 'Error in stream_wait_async ! ',cuda_success
      end

      subroutine start_dir_stream_cover
      implicit none
      integer cuda_success
      cuda_success = 0

      cuda_success = cuda_success +
     &     cudaEventRecord(dir_event,rec_stream)
      cuda_success = cuda_success +
     &     cudaStreamWaitEvent(dir_stream,dir_event,zero_flag)

      if (cuda_success.ne.0)
     &   write(*,*) 'Error in start_dir_stream_cover!',cuda_success

      end subroutine

      subroutine end_dir_stream_cover
      implicit none
      integer cuda_success
      cuda_success = 0

      cuda_success = cuda_success +
     &     cudaEventRecord(rec_event,dir_stream)
      cuda_success = cuda_success +
     &     cudaStreamWaitEvent(rec_stream,rec_event,zero_flag)

      if (cuda_success.ne.0)
     &   write(*,*) 'Error in start_dir_stream_cover!',cuda_success

      end subroutine

      subroutine getDeviceProp
      implicit none
      integer:: ierr=0,i

      ierr = CUDAGETDEVICEPROPERTIES(devProp,devicenum)
      if (ierr.ne.0) 
     &   print*, 'error getting device properties ',ierr

      nSMP = devProp%multiProcessorCount
      select case (devProp%major)
         case(2)
            if (devProp%minor.eq.1) then
               cores_SMP = 48
            else
               cores_SMP = 32
            end if
         case(3)
            cores_SMP = 192
         case(5)
            cores_SMP = 128
         case(6)
            if (devProp%minor.eq.1) then
               cores_SMP = 128
            else if (devProp%minor.eq.0) then
               cores_SMP = 64
            else
               print*,'rank',rank,'Unknown device type'
               stop
            end if
         case(7)
               cores_SMP = 64
         case(8)
               cores_SMP = 64
         case default
               print*,'rank',rank,'Unknown device type',
     &         ' update getDeviceProp routine for compute capability'
               cores_SMP = 128
      end select

      nSPcores = cores_SMP*nSMP

      end subroutine

      ! Enables polarisation precomputation following
      ! the available memory
      logical function allow_precomp(n) result(auth)
      implicit none
      integer,intent(in) :: n
#ifdef _OPENACC
      integer(mipk) free,total,esti
      integer ierr
      esti = 200*n*7*ti_p
      ierr = cudaMemGetinfo(free,total)
      auth = (esti.lt.(free*9)/10)
#else
      auth = .true.
#endif
      end function

      integer function get_GridDim(n,bDim,opt) result(GDim)
      integer,intent(in):: bDim,n
      integer,optional :: opt
      integer,parameter:: maxB=2**14
      if (present(opt)) then
         GDim = min(((n-1)/ishft(bDim,1))+1,maxB)
      else
         GDim = min(((n-1)*WARP_SIZE/bDim)+1,maxB)
      end if
      end function

      subroutine dev_sort_i4(n,array,stream)
      use thrust, only: thrust_sort
      implicit none
      integer,intent(in   ):: n
      integer,intent(inout):: array(n)
      integer(cuda_stream_kind),intent(in) :: stream
      integer ierr

      !ierr = cudaStreamSynchronize(stream)
!$acc host_data use_device(array)
      call thrust_sort(array,n,stream)
!$acc end host_data
      end subroutine
#endif
      !---------------------------------------
      ! memory Set routines
      !---------------------------------------
      subroutine mem_set_i4(dst,val,n,stream,offset)
      integer(4) dst(*)
      integer(4),intent(in):: val
      integer(mipk),intent(in):: n
      integer(cuda_stream_kind),intent(in):: stream
      integer(mipk),optional:: offset
      integer ierr
      integer(mipk) i,offset_
 14   format(A,I4,A,/,4X,A)

      offset_ = zero8
      if (present(offset)) offset_=offset
#ifdef _OPENACC
!$acc host_data use_device(dst)
      if (stream.eq.zero8) then
         ierr = cudaMemset(dst(offset_+1),val,n)
         ierr = ierr + cudaStreamSynchronize(stream)
      else
         ierr = cudaMemsetAsync(dst(offset_+1),val,n,stream)
      end if
!$acc end host_data
      if (ierr.ne.cudasuccess) then
         write(*,14) 'error',ierr,'detected in utilgpu_mem_set_i4'
     &              ,cudageterrorstring(ierr)
      end if
#else
      do i = offset_+1,offset_+n
         dst(i) = val
      end do
#endif
      end subroutine
      subroutine mem_set_i8(dst,val,n,stream,offset)
      integer(8) dst(*)
      integer(8),intent(in):: val
      integer(mipk),intent(in):: n
      integer(cuda_stream_kind),intent(in):: stream
      integer(mipk),optional:: offset
      integer ierr
      integer(mipk) i,offset_
 14   format(A,I4,A,/,4X,A)

      offset_ = zero8
      if (present(offset)) offset_=offset
#ifdef _OPENACC
!$acc host_data use_device(dst)
      if (stream.eq.zero8) then
         ierr = cudaMemset(dst(offset_+1),val,n)
         ierr = ierr + cudaStreamSynchronize(stream)
      else
         ierr = cudaMemsetAsync(dst(offset_+1),val,n,stream)
      end if
!$acc end host_data
      if (ierr.ne.cudasuccess) then
         write(*,14) 'error',ierr,'detected in utilgpu_mem_set_i8'
     &              ,cudageterrorstring(ierr)
      end if
#else
      do i = offset_+1,offset_+n
         dst(i) = val
      end do
#endif
      end subroutine
      subroutine mem_set_r4(dst,val,n,stream,offset)
      real(4) dst(*)
      real(4),intent(in):: val
      integer(mipk),intent(in):: n
      integer(cuda_stream_kind),intent(in):: stream
      integer(mipk),optional:: offset
      integer ierr
      integer(mipk) i,offset_
 14   format(A,I4,A,/,4X,A)

      offset_ = zero8
      if (present(offset)) offset_=offset
#ifdef _OPENACC
!$acc host_data use_device(dst)
      if (stream.eq.zero8) then
         ierr = cudaMemset(dst(offset_+1),val,n)
         ierr = ierr + cudaStreamSynchronize(stream)
      else
         ierr = cudaMemsetAsync(dst(offset_+1),val,n,stream)
      end if
!$acc end host_data
      if (ierr.ne.cudasuccess) then
         write(*,14) 'error',ierr,'detected in utilgpu_mem_set_r4'
     &              ,cudageterrorstring(ierr)
      end if
#else
      do i = offset_+1,offset_+n
         dst(i) = val
      end do
#endif
      end subroutine
      subroutine mem_set_r4_2d(dst,val,n,stream,offset)
      real(4) dst(n,1)
      real(4),intent(in):: val
      integer(mipk),intent(in):: n
      integer(cuda_stream_kind),intent(in):: stream
      integer(mipk),optional:: offset
      integer ierr
      integer(mipk) i,offset_
   14   format(A,I4,A,/,4X,A)

      offset_ = zero8
      if (present(offset)) offset_=offset
#ifdef _OPENACC
!$acc host_data use_device(dst)
      if (stream.eq.zero8) then
         ierr = cudaMemset(dst(offset_+1,1),val,n)
         ierr = ierr + cudaStreamSynchronize(stream)
      else
         ierr = cudaMemsetAsync(dst(offset_+1,1),val,n,stream)
      end if
!$acc end host_data
      if (ierr.ne.cudasuccess) then
         write(*,14) 'error',ierr,'detected in utilgpu_mem_set_r4'
     &               ,cudageterrorstring(ierr)
      end if
#else
      do i = offset_+1,offset_+n
         dst(i,1) = val
      end do
#endif
      end subroutine
      subroutine mem_set_r4_3d(dst,val,n,stream,offset)
      real(4) dst(n,1,1)
      real(4),intent(in):: val
      integer(mipk),intent(in):: n
      integer(cuda_stream_kind),intent(in):: stream
      integer(mipk),optional:: offset
      integer ierr
      integer(mipk) i,offset_
   14   format(A,I4,A,/,4X,A)

      offset_ = zero8
      if (present(offset)) offset_=offset
#ifdef _OPENACC
!$acc host_data use_device(dst)
      if (stream.eq.zero8) then
         ierr = cudaMemset(dst(offset_+1,1,1),val,n)
         ierr = ierr + cudaStreamSynchronize(stream)
      else
         ierr = cudaMemsetAsync(dst(offset_+1,1,1),val,n,stream)
      end if
!$acc end host_data
      if (ierr.ne.cudasuccess) then
         write(*,14) 'error',ierr,'detected in utilgpu_mem_set_r4'
     &               ,cudageterrorstring(ierr)
      end if
#else
      do i = offset_+1,offset_+n
         dst(i,1,1) = val
      end do
#endif
      end subroutine
      subroutine mem_set_r4_5d(dst,val,n,stream,offset)
      real(4) dst(n,1,1,1,1)
      real(4),intent(in):: val
      integer(mipk),intent(in):: n
      integer(cuda_stream_kind),intent(in):: stream
      integer(mipk),optional:: offset
      integer ierr
      integer(mipk) i,offset_
   14   format(A,I4,A,/,4X,A)

      offset_ = zero8
      if (present(offset)) offset_=offset
#ifdef _OPENACC
!$acc host_data use_device(dst)
      if (stream.eq.zero8) then
         ierr = cudaMemset(dst(offset_+1,1,1,1,1),val,n)
         ierr = ierr + cudaStreamSynchronize(stream)
      else
         ierr = cudaMemsetAsync(dst(offset_+1,1,1,1,1),val,n,stream)
      end if
!$acc end host_data
      if (ierr.ne.cudasuccess) then
         write(*,14) 'error',ierr,'detected in utilgpu_mem_set_r4'
     &               ,cudageterrorstring(ierr)
      end if
#else
      do i = offset_+1,offset_+n
         dst(i,1,1,1,1) = val
      end do
#endif
      end subroutine
      subroutine mem_set_r8(dst,val,n,stream,offset)
      real(8) dst(*)
      real(8),intent(in):: val
      integer(mipk),intent(in):: n
      integer(cuda_stream_kind),intent(in):: stream
      integer(mipk),optional:: offset
      integer ierr
      integer(mipk) i,offset_
 14   format(A,I4,A,/,4X,A)

      offset_ = zero8
      if (present(offset)) offset_=offset
#ifdef _OPENACC
!$acc host_data use_device(dst)
      if (stream.eq.zero8) then
         ierr = cudaMemset(dst(offset_+1),val,n)
         ierr = ierr + cudaStreamSynchronize(stream)
      else
         ierr = cudaMemsetAsync(dst(offset_+1),val,n,stream)
      end if
!$acc end host_data
      if (ierr.ne.cudasuccess) then
         write(*,14) 'error',ierr,'detected in utilgpu_mem_set_r8'
     &              ,cudageterrorstring(ierr)
      end if
#else
      do i = offset_+1,offset_+n
         dst(i) = val
      end do
#endif
      end subroutine
      subroutine mem_set_r8_2d(dst,val,n,stream,offset)
      real(8) dst(n,1)
      real(8),intent(in):: val
      integer(mipk),intent(in):: n
      integer(cuda_stream_kind),intent(in):: stream
      integer(mipk),optional:: offset
      integer ierr
      integer(mipk) i,offset_
   14   format(A,I4,A,/,4X,A)

      offset_ = zero8
      if (present(offset)) offset_=offset
#ifdef _OPENACC
!$acc host_data use_device(dst)
      if (stream.eq.zero8) then
         ierr = cudaMemset(dst(offset_+1,1),val,n)
         ierr = ierr + cudaStreamSynchronize(stream)
      else
         ierr = cudaMemsetAsync(dst(offset_+1,1),val,n,stream)
      end if
!$acc end host_data
      if (ierr.ne.cudasuccess) then
         write(*,14) 'error',ierr,'detected in utilgpu_mem_set_r8'
     &               ,cudageterrorstring(ierr)
      end if
#else
      do i = offset_+1,offset_+n
         dst(i,1) = val
      end do
#endif
      end subroutine
      subroutine mem_set_r8_3d(dst,val,n,stream,offset)
      real(8) dst(n,1,1)
      real(8),intent(in):: val
      integer(mipk),intent(in):: n
      integer(cuda_stream_kind),intent(in):: stream
      integer(mipk),optional:: offset
      integer ierr
      integer(mipk) i,offset_
   14   format(A,I4,A,/,4X,A)

      offset_ = zero8
      if (present(offset)) offset_=offset
#ifdef _OPENACC
!$acc host_data use_device(dst)
      if (stream.eq.zero8) then
         ierr = cudaMemset(dst(offset_+1,1,1),val,n)
         ierr = ierr + cudaStreamSynchronize(stream)
      else
         ierr = cudaMemsetAsync(dst(offset_+1,1,1),val,n,stream)
      end if
!$acc end host_data
      if (ierr.ne.cudasuccess) then
         write(*,14) 'error',ierr,'detected in utilgpu_mem_set_r8'
     &               ,cudageterrorstring(ierr)
      end if
#else
      do i = offset_+1,offset_+n
         dst(i,1,1) = val
      end do
#endif
      end subroutine
      subroutine mem_set_r8_5d(dst,val,n,stream,offset)
      real(8) dst(n,1,1,1,1)
      real(8),intent(in):: val
      integer(mipk),intent(in):: n
      integer(cuda_stream_kind),intent(in):: stream
      integer(mipk),optional:: offset
      integer ierr
      integer(mipk) i,offset_
   14   format(A,I4,A,/,4X,A)

      offset_ = zero8
      if (present(offset)) offset_=offset
#ifdef _OPENACC
!$acc host_data use_device(dst)
      if (stream.eq.zero8) then
         ierr = cudaMemset(dst(offset_+1,1,1,1,1),val,n)
         ierr = ierr + cudaStreamSynchronize(stream)
      else
         ierr = cudaMemsetAsync(dst(offset_+1,1,1,1,1),val,n,stream)
      end if
!$acc end host_data
      if (ierr.ne.cudasuccess) then
         write(*,14) 'error',ierr,'detected in utilgpu_mem_set_r4'
     &               ,cudageterrorstring(ierr)
      end if
#else
      do i = offset_+1,offset_+n
         dst(i,1,1,1,1) = val
      end do
#endif
      end subroutine

      !---------------------------------------
      ! memory Copy routines
      !---------------------------------------
      subroutine mem_move_i4(dst,src,n,stream)
      integer(4) src(*),dst(*)
      integer(mipk),intent(in):: n
      integer(cuda_stream_kind),intent(in):: stream
      integer ierr
      integer(mipk) i

#ifdef _OPENACC
!$acc host_data use_device(dst,src)
      ierr = cudaMemCpyAsync(dst,src,n,cudaMemcpyDeviceToDevice,stream)
      if(stream.eq.zero8) ierr = ierr + cudaStreamSynchronize(stream)
!$acc end host_data
      if (ierr.ne.cudasuccess) then
         write(*,*) 'error',ierr,'detected in utilgpu_mem_move_i4'
         write(*,*) cudageterrorstring(ierr)
      end if
#else
      do i = 1,n
         dst(i) = src(i)
      end do
#endif
      end subroutine

      subroutine mem_move_i8(dst,src,n,stream)
      integer(8) src(*),dst(*)
      integer(mipk),intent(in):: n
      integer(cuda_stream_kind),intent(in):: stream
      integer ierr
      integer(mipk) i

#ifdef _OPENACC
!$acc host_data use_device(dst,src)
      ierr = cudamemcpyasync(dst,src,n,cudaMemcpyDeviceToDevice,stream)
      if(stream.eq.zero8) ierr = ierr + cudaStreamSynchronize(stream)
!$acc end host_data
      if (ierr.ne.cudasuccess) then
         write(*,*) 'error',ierr,'detected in utilgpu_mem_move_i8'
         write(*,*) cudageterrorstring(ierr)
      end if
#else
      do i = 1,n
         dst(i) = src(i)
      end do
#endif
      end subroutine

      subroutine mem_move_r4(dst,src,n,stream)
      real(4) src(*),dst(*)
      integer(mipk),intent(in):: n
      integer(cuda_stream_kind),intent(in):: stream
      integer ierr
      integer(mipk) i

#ifdef _OPENACC
!$acc host_data use_device(dst,src)
      ierr = cudamemcpyasync(dst,src,n,cudaMemcpyDeviceToDevice,stream)
      if(stream.eq.zero8) ierr = ierr + cudaStreamSynchronize(stream)
!$acc end host_data
      if (ierr.ne.cudasuccess) then
         write(*,*) 'error',ierr,'detected in utilgpu_mem_move_r4'
         write(*,*) cudageterrorstring(ierr)
      end if
#else
      do i = 1,n
         dst(i) = src(i)
      end do
#endif
      end subroutine

      subroutine mem_move_r4_2d(dst,src,n,stream)
      integer(mipk),intent(in):: n
      real(4) src(n,1),dst(n,1)
      integer(cuda_stream_kind),intent(in):: stream
      integer ierr
      integer(mipk) i
   
#ifdef _OPENACC
!$acc host_data use_device(dst,src)
      ierr = cudamemcpyasync(dst,src,n,cudaMemcpyDeviceToDevice,stream)
      if(stream.eq.zero8) ierr = ierr + cudaStreamSynchronize(stream)
!$acc end host_data
      if (ierr.ne.cudasuccess) then
         write(*,*) 'error',ierr,'detected in utilgpu_mem_move_r4'
         write(*,*) cudageterrorstring(ierr)
      end if
#else
      do i = 1,n
         dst(i,1) = src(i,1)
      end do
#endif
      end subroutine

      subroutine mem_move_r8(dst,src,n,stream)
      integer(mipk),intent(in):: n
      real(8) src(*),dst(*)
      integer(cuda_stream_kind),intent(in):: stream
      integer ierr
      integer(mipk) i

#ifdef _OPENACC
!$acc host_data use_device(dst,src)
      ierr = cudamemcpyasync(dst,src,n,cudaMemcpyDeviceToDevice,stream)
      if(stream.eq.zero8) ierr = ierr + cudaStreamSynchronize(stream)
!$acc end host_data
      if (ierr.ne.cudasuccess) then
         write(*,*) 'error',ierr,'detected in utilgpu_mem_move_r8'
         write(*,*) cudageterrorstring(ierr)
      end if
#else
      do i = 1,n
         dst(i) = src(i)
      end do
#endif
      end subroutine

      subroutine mem_move_r8_2d(dst,src,n,stream)
      integer(mipk),intent(in):: n
      real(8) src(n,1),dst(n,1)
      integer(cuda_stream_kind),intent(in):: stream
      integer ierr
      integer(mipk) i

#ifdef _OPENACC
!$acc host_data use_device(dst,src)
      ierr = cudamemcpyasync(dst,src,n,cudaMemcpyDeviceToDevice,stream)
      if(stream.eq.zero8) ierr = ierr + cudaStreamSynchronize(stream)
!$acc end host_data
      if (ierr.ne.cudasuccess) then
         write(*,*) 'error',ierr,'detected in utilgpu_mem_move_r8'
         write(*,*) cudageterrorstring(ierr)
      end if
#else
      do i = 1,n
         dst(i,1) = src(i,1)
      end do
#endif
      end subroutine

      subroutine get_atom_glob_id(a_kind,a_id,glob,ntyp,natom)
      implicit none
      integer,intent(in):: natom,ntyp
      integer,intent(in):: a_kind(ntyp),a_id(natom)
      integer,intent(inout):: glob(ntyp)
      integer i
!$acc parallel loop async default(present)
      do i = 1,ntyp
         glob(i) = a_id(a_kind(i))
      end do
      end subroutine

      subroutine reduce_buffer_m(r_buff,sz_buff,redVal,queue)
      implicit none
      integer  ,intent(in ):: sz_buff, queue
      real(r_p)            :: r_buff(*)
      real(r_p),intent(out):: redVal
      integer i
!$acc parallel loop present(redVal,r_buff) async(queue)
      do i = 1,sz_buff
         redVal    = redVal + r_buff(i)
         r_buff(i) = 0
      end do
      end subroutine

      subroutine reduce_energy_virial0(energy,vxx,vxy,vxz,vyy,vyz,vzz
     &           ,queue)
      implicit none
      integer,intent(in):: queue
      real(r_p) energy
      real(r_p) vxx,vxy,vxz,vyy,vyz,vzz
      integer i
      !Make sure data are set zero before calling this routine

!$acc parallel loop async(queue) present(energy,vxx,vxy,vxz,
!$acc&    vyy,vyz,vzz,ered_buff,vred_buff)
      do i = 1,RED_BUFF_SIZE
         energy = energy + ered_buff(i)
         vxx    = vxx + real(vred_buff(i),r_p)
         vxy    = vxy + real(vred_buff(  RED_BUFF_SIZE+i),r_p)
         vxz    = vxz + real(vred_buff(2*RED_BUFF_SIZE+i),r_p)
         vyy    = vyy + real(vred_buff(3*RED_BUFF_SIZE+i),r_p)
         vyz    = vyz + real(vred_buff(4*RED_BUFF_SIZE+i),r_p)
         vzz    = vzz + real(vred_buff(5*RED_BUFF_SIZE+i),r_p)

         ered_buff(i)                 = 0
         vred_buff(i)                 = 0
         vred_buff(  RED_BUFF_SIZE+i) = 0
         vred_buff(2*RED_BUFF_SIZE+i) = 0
         vred_buff(3*RED_BUFF_SIZE+i) = 0
         vred_buff(4*RED_BUFF_SIZE+i) = 0
         vred_buff(5*RED_BUFF_SIZE+i) = 0
      end do
      end subroutine

      subroutine reduce_energy_virial1(energy,vxx,vxy,vxz,vyy,vyz,vzz
     &           ,ered_b,queue)
      implicit none
      integer,intent(in):: queue
      ener_rtyp,intent(inout):: ered_b(RED_BUFF_SIZE)
      ener_rtyp energy
      real(r_p) vxx,vxy,vxz,vyy,vyz,vzz
      integer i
      !Make sure data are set zero before calling this routine

!$acc parallel loop async(queue) present(energy,vxx,vxy,vxz,
!$acc&    vyy,vyz,vzz,ered_b,vred_buff)
      do i = 1,RED_BUFF_SIZE
         energy = energy + ered_b(i)
         vxx    = vxx + real(vred_buff(0*RED_BUFF_SIZE+i),r_p)
         vxy    = vxy + real(vred_buff(1*RED_BUFF_SIZE+i),r_p)
         vxz    = vxz + real(vred_buff(2*RED_BUFF_SIZE+i),r_p)
         vyy    = vyy + real(vred_buff(3*RED_BUFF_SIZE+i),r_p)
         vyz    = vyz + real(vred_buff(4*RED_BUFF_SIZE+i),r_p)
         vzz    = vzz + real(vred_buff(5*RED_BUFF_SIZE+i),r_p)

         ered_b   (i)                 = 0
         vred_buff(i)                 = 0
         vred_buff(  RED_BUFF_SIZE+i) = 0
         vred_buff(2*RED_BUFF_SIZE+i) = 0
         vred_buff(3*RED_BUFF_SIZE+i) = 0
         vred_buff(4*RED_BUFF_SIZE+i) = 0
         vred_buff(5*RED_BUFF_SIZE+i) = 0
      end do
      end subroutine

      subroutine reduce_energy_virial2(energy,vxx,vxy,vxz,vyy,vyz,vzz
     &           ,ered_b,vred_b,queue)
      implicit none
      integer,intent(in):: queue
      real(r_p),intent(inout):: ered_b(  RED_BUFF_SIZE)
      real(t_p),intent(inout):: vred_b(6*RED_BUFF_SIZE)
      real(r_p) energy
      real(r_p) vxx,vxy,vxz,vyy,vyz,vzz
      integer i
      !Make sure data are set zero before calling this routine

!$acc parallel loop async(queue) present(energy,vxx,vxy,vxz,
!$acc&    vyy,vyz,vzz,ered_b,vred_buff)
      do i = 1,RED_BUFF_SIZE
         energy = energy + ered_b(i)
         vxx    = vxx + real(vred_b(0*RED_BUFF_SIZE+i),r_p)
         vxy    = vxy + real(vred_b(1*RED_BUFF_SIZE+i),r_p)
         vxz    = vxz + real(vred_b(2*RED_BUFF_SIZE+i),r_p)
         vyy    = vyy + real(vred_b(3*RED_BUFF_SIZE+i),r_p)
         vyz    = vyz + real(vred_b(4*RED_BUFF_SIZE+i),r_p)
         vzz    = vzz + real(vred_b(5*RED_BUFF_SIZE+i),r_p)

         ered_b(i)                 = 0
         vred_b(i)                 = 0
         vred_b(  RED_BUFF_SIZE+i) = 0
         vred_b(2*RED_BUFF_SIZE+i) = 0
         vred_b(3*RED_BUFF_SIZE+i) = 0
         vred_b(4*RED_BUFF_SIZE+i) = 0
         vred_b(5*RED_BUFF_SIZE+i) = 0
      end do
      end subroutine

      subroutine reduce_energy_action0(energy,n,queue)
      implicit none
      integer,intent(in):: queue
      real(r_p) energy
      integer   n
      integer i
      !Make sure data are set zero before calling this routine

!$acc parallel loop async(queue) present(energy,n,
!$acc&    ered_buff,vred_buff)
      do i = 1,RED_BUFF_SIZE
         energy = energy + ered_buff(i)
         n      = n + nred_buff(i)
         ered_buff(i) = 0
         nred_buff(i) = 0
      end do
      end subroutine

      subroutine reduce_energy_action1(energy,n,ered_b,queue)
      implicit none
      integer,intent(in):: queue
      ener_rtyp ered_b(RED_BUFF_SIZE)
      ener_rtyp energy
      integer   n
      integer i
      !Make sure data are set zero before calling this routine

!$acc parallel loop async(queue) present(energy,n,
!$acc&    ered_buff,vred_buff)
      do i = 1,RED_BUFF_SIZE
         energy = energy + ered_b(i)
         n      = n      + nred_buff(i)
         ered_b (i)   = 0
         nred_buff(i) = 0
      end do
      end subroutine

      subroutine zero_mdred_buffers(queue_)
      implicit none
      integer,optional::queue_
      integer i,queue

      queue = dir_queue
      if (present(queue_)) queue=queue_

!$acc parallel loop async(queue)
!$acc&         default(present)
      do i = 1,6*RED_BUFF_SIZE
         if (i<RED_BUFF_SIZE+1) then
            ered_buff(i)=0
            ered_buf1(i)=0
            nred_buff(i)=0
         end if
         vred_buff(i)=0
         vred_buf1(i)=0
      end do
      end subroutine

      subroutine zero_evir_red_buffer(queue_)
      implicit none
      integer,optional::queue_
      integer i,queue

      queue = dir_queue
      if (present(queue_)) queue=queue_

!$acc parallel loop async(queue)
!$acc&         present(ered_buff,ered_buf1,vred_buff)
      do i = 1,6*RED_BUFF_SIZE
         if (i<RED_BUFF_SIZE+1) then
            ered_buff(i)=0
            ered_buf1(i)=0
         end if
         vred_buff(i) = 0
      end do
      end subroutine

      subroutine zero_en_red_buffer(queue_)
      implicit none
      integer,optional::queue_
      integer i,queue

      if   (present(queue_)) then
             queue = queue_
      else;  queue = dir_queue; end if

!$acc parallel loop async(queue)
!$acc&         present(ered_buff,ered_buf1,nred_buff)
      do i = 1,RED_BUFF_SIZE
         ered_buff(i)=0
         ered_buf1(i)=0
         nred_buff(i)=0
      end do

      end subroutine

      subroutine transposez_r6(src,dst,loc,sloc,siz,queue)
      implicit none
      integer  ,intent(in)   :: siz,queue,sloc,loc(sloc)
      real(t_p),intent(inout):: src(siz*6)
      real(t_p),intent(out)  :: dst(6*siz)
      integer   i,j,idxd,idxs

!$acc host_data use_device(src,dst,loc)
!$acc parallel loop collapse(2) deviceptr(src,dst,loc) async(queue)
      do i = 1,sloc; do j = 0,5
         idxd = (loc(i)-1)*6+j+1
         idxs = j*siz +i
        dst( idxd ) = dst( idxd ) + src( idxs )
        src( idxs ) = 0
      end do; end do
!$acc end host_data

      end subroutine

      subroutine transpose_az_l3(src,dst,loc,sloc,siz,queue)
      implicit none
      integer   ,parameter    :: n1=3
      integer   ,intent(in)   :: siz,queue,sloc,loc(sloc)
      integer(8),intent(inout):: src(siz*n1)
      integer(8),intent(out)  :: dst(n1*siz)
      integer   i,j,idxd,idxs
!$acc host_data use_device(src,dst,loc)
!$acc parallel loop collapse(2) deviceptr(src,dst,loc) async(queue)
      do i = 1,sloc; do j = 0,n1-1
         idxd = (loc(i)-1)*n1+j+1
         idxs = j*siz +i
        dst( idxd ) = dst( idxd ) + src( idxs )
        src( idxs ) = 0
      end do; end do
!$acc end host_data
      end subroutine
      subroutine transpose_az_d3(src,dst,loc,sloc,siz,queue)
      implicit none
      integer  ,parameter    :: n1=3
      integer  ,intent(in)   :: siz,queue,sloc,loc(sloc)
      real(8)  ,intent(inout):: src(siz*n1)
      real(8)  ,intent(out)  :: dst(n1*siz)
      integer   i,j,idxd,idxs
!$acc host_data use_device(src,dst,loc)
!$acc parallel loop collapse(2) deviceptr(src,dst,loc) async(queue)
      do i = 1,sloc; do j = 0,n1-1
         idxd = (loc(i)-1)*n1+j+1
         idxs = j*siz +i
        dst( idxd ) = dst( idxd ) + src( idxs )
        src( idxs ) = 0
      end do; end do
!$acc end host_data
      end subroutine
      subroutine transpose_az_f3(src,dst,loc,sloc,siz,queue)
      implicit none
      integer  ,parameter    :: n1=3
      integer  ,intent(in)   :: siz,queue,sloc,loc(sloc)
      real(4)  ,intent(inout):: src(siz*n1)
      real(4)  ,intent(out)  :: dst(n1*siz)
      integer   i,j,idxd,idxs
!$acc host_data use_device(src,dst,loc)
!$acc parallel loop collapse(2) deviceptr(src,dst,loc) async(queue)
      do i = 1,sloc; do j = 0,n1-1
         idxd = (loc(i)-1)*n1+j+1
         idxs = j*siz +i
        dst( idxd ) = dst( idxd ) + src( idxs )
        src( idxs ) = 0
      end do; end do
!$acc end host_data
      end subroutine

      end module utilgpu

