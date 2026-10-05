!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!
!  subroutine fft_setup: setup the arrays for the computations of
!  3d ffts with 2decomp&FFT
!
#include "tinker_macro.h"
subroutine fft_setup(nrec,rank_bis,istart1,iend1,isize1,jstart1, &
   jend1,jsize1,kstart1,kend1,ksize1,istart2,iend2,isize2,jstart2,     &
   jend2,jsize2,kstart2,kend2,ksize2,nx,ny,nz,comm_loc,ngrid1,ngrid2)

   use decomp_2d
   use decomp_2d_fft
   use domdec   ,only: ranktot
   use fft      ,only: in, out, is_fftInit
   use mpi
   use inform   ,only: verbose
   use tinMemory,only: s_cufft
   use sizes    ,only: tinkerdebug
#ifdef _OPENACC
   use mdevice  ,only: direct_work,reciprocal_work
   use decomp_2d_cufft
#endif

   implicit none
   integer :: ierror,i,j,k,nrec,rank_bis,comm_loc,ngrid1,ngrid2
   integer :: istart1(*),iend1(*),isize1(*)
   integer :: jstart1(*),jend1(*),jsize1(*)
   integer :: kstart1(*),kend1(*),ksize1(*)
   integer :: istart2(*),iend2(*),isize2(*)
   integer :: jstart2(*),jend2(*),jsize2(*)
   integer :: kstart2(*),kend2(*),ksize2(*)

   integer :: nx, ny, nz
   integer :: p_row, p_col
!complex(KIND(0.0D0)), dimension(nx,ny,nz) :: in1, out1
   real(t_p) time0,time1
#ifdef _OPENACC
#endif

   p_row = ngrid1
   p_col = ngrid2
   dcp_verbose = merge(verbose,.false.,ranktot.eq.0.and..not.is_fftInit)

! Finalize fft library if necessary
   if (is_fftInit) then
      call fft_final
      is_fftInit = .false.
   end if

   call decomp_2d_init(nx,ny,nz,p_row,p_col,comm_loc)
!
   istart1(rank_bis+1) = xstart(1)
   iend1(rank_bis+1)   = xend(1)
   isize1(rank_bis+1)  = xsize(1)
   jstart1(rank_bis+1) = xstart(2)
   jend1(rank_bis+1)   = xend(2)
   jsize1(rank_bis+1)  = xsize(2)
   kstart1(rank_bis+1) = xstart(3)
   kend1(rank_bis+1)   = xend(3)
   ksize1(rank_bis+1)  = xsize(3)
!
   istart2(rank_bis+1) = zstart(1)
   iend2(rank_bis+1)   = zend(1)
   isize2(rank_bis+1)  = zsize(1)
   jstart2(rank_bis+1) = zstart(2)
   jend2(rank_bis+1)   = zend(2)
   jsize2(rank_bis+1)  = zsize(2)
   kstart2(rank_bis+1) = zstart(3)
   kend2(rank_bis+1)   = zend(3)
   ksize2(rank_bis+1)  = zsize(3)
!
   call decomp_2d_fft_init
!
#ifdef _OPENACC
   call get_rec_dir_queue(reciprocal_work%queue,direct_work%queue)
   call init_cufft_engine
   s_cufft = decomp2d_cufftGetSize()
#endif
   is_fftInit = .true.

   if (tinkerdebug.gt.0) then
      if (rank_bis.eq.0) write(*,*) '--- fft_setup'
      do i = 0,nrec-1
         if (rank_bis.eq.i) then
13          format(' --fft_dd',I4,'in ',6I8)
            write(*,13) i, xstart(1),xend(1),xstart(2),xend(2),xstart(3),xend(3)
         end if
         call MPI_Barrier(comm_loc,k)
      end do
      do i = 0,nrec-1
         if (rank_bis.eq.i) then
14          format(9X,I4,'out',6I8)
            write(*,14) i, zstart(1),zend(1),zstart(2),zend(2),zstart(3),zend(3)
         end if
         call MPI_Barrier(comm_loc,k)
      end do
   end if
end
!
!
!  subroutine fft2d_frontmpi: compute forward 3dFFT in parallel with 2decomp&FFT
!
subroutine fft2d_frontmpi(qgridin,qgridout,n1mpimax,n2mpimax,n3mpimax)
   use decomp_2d
   use decomp_2d_fft
   use inform   ,only: deb_Path
   use timestat ,only: timer_enter,timer_exit,timer_ffts,quiet_timers
   implicit none
   include 'mpif.h'
   complex(t_p), allocatable, dimension(:,:,:) :: in, out
   real(t_p) qgridin(2,n1mpimax,n2mpimax,n3mpimax,*)
   real(t_p) qgridout(2,zsize(1),zsize(2),zsize(3))
   real(t_p) time0,time1
   integer i,j,k,n1mpimax,n2mpimax,n3mpimax
! input is X-pencil data
! output is Z-pencil data
   call timer_enter( timer_ffts )
   if (deb_Path) write(*,'(5x,a)') '>> fft2d_frontmpi'
   allocate (in(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
   allocate (out(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)))
!
! fill in array with qgridin
!
   do k = 1, xsize(3)
      do j = 1, xsize(2)
         do i = 1, xsize(1)
            in(xstart(1)+i-1,xstart(2)+j-1,xstart(3)+k-1) = cmplx(qgridin(1,i,j,k,1), &
               qgridin(2,i,j,k,1),kind=t_p)
         end do
      end do
   end do
!
! ===== 3D forward FFT =====
   time0 = mpi_wtime()
   call decomp_2d_fft_3d(in, out, DECOMP_2D_FFT_FORWARD)
   time1 = mpi_wtime()
!
! fill qgridout array with out
!
   do k = 1, zsize(3)
      do j = 1, zsize(2)
         do i = 1, zsize(1)
            qgridout(1,i,j,k) = real(out(zstart(1)+i-1,zstart(2)+j-1,zstart(3)+k-1))
            qgridout(2,i,j,k) = aimag(out(zstart(1)+i-1,zstart(2)+j-1,zstart(3)+k-1))
         end do
      end do
   end do
!
   deallocate (in)
   deallocate (out)
   decomp2d_mpi_fcall=1
   if (deb_Path) write(*,'(5x,a)') '<< fft2d_frontmpi'
   call timer_exit( timer_ffts,quiet_timers )
end
!
#ifdef _OPENACC
subroutine cufft2d_frontmpi(qgridin,qgridout,n1mpimax,n2mpimax,n3mpimax)
   use decomp_2d
   use decomp_2d_fft
   use decomp_2d_cufft
   use inform   ,only: deb_Path
   use fft
   use timestat ,only: timer_enter,timer_exit,timer_ffts,quiet_timers
   use mdevice  ,only: reciprocal_work
   implicit none
   !include 'mpif.h'
   real(t_p) qgridin(2,n1mpimax,n2mpimax,n3mpimax,*)
   real(t_p) qgridout(2,zsize(1),zsize(2),zsize(3))
   real(t_p) time0,time1
   integer i,j,k,n1mpimax,n2mpimax,n3mpimax
   logical,save::first_in=.true.


   call timer_enter( timer_ffts )
   if (deb_Path) write(*,'(5x,a)') '>> cufft2d_frontmpi'

   if (first_in) then
      first_in = .false.
      if (.not.grid_pointer_assoc) call malloc_FFtgrid_p
   end if

   if (.not.grid_pointer_assoc) then
      !
      ! fill in array with qgridin
      !
      call r2c_grid(qgridin,n1mpimax*n2mpimax*n3mpimax,in, &
         reciprocal_work)
   else
      call associate_grid_pointer(qgridin,qgridout)
   end if

   !
   ! ===== 3D forward FFT =====
   !time0 = mpi_wtime()
   call decomp_2d_cufft_3d(in, out, DECOMP_2D_FFT_FORWARD)
   !time1 = mpi_wtime()

   !
   ! fill qgridout array with out
   !
   if (.not.grid_pointer_assoc) &
      call c2r_grid(out,zsize(1)*zsize(2)*zsize(3),qgridout, &
         reciprocal_work)
   !
   if (deb_Path) write(*,'(5x,a)') '<< cufft2d_frontmpi'
   call timer_exit( timer_ffts,quiet_timers )
end
#endif
!
!  subroutine fft2d_backmpi: compute backward 3dFFT in parallel with 2decomp&FFT
!
subroutine fft2d_backmpi(qgridin,qgridout,n1mpimax,n2mpimax,n3mpimax)
   use decomp_2d
   use decomp_2d_fft
   use inform   ,only: deb_Path
   use timestat ,only: timer_enter,timer_exit,timer_ffts,quiet_timers
   implicit none
   include 'mpif.h'
   complex(t_p), allocatable, dimension(:,:,:) :: in, out
   real(t_p) qgridin(2,n1mpimax,n2mpimax,n3mpimax,*)
   real(t_p) qgridout(2,zsize(1),zsize(2),zsize(3))
   real(t_p) time0,time1
   integer i,j,k,n1mpimax,n2mpimax,n3mpimax
! input is Z-pencil data
! output is X-pencil data
   call timer_enter( timer_ffts )
   if (deb_Path) write(*,'(5x,a)') '>> fft2d_backmpi'
   allocate (in(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
   allocate (out(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)))
!
! fill out array with qgridout
!
   do k = 1, zsize(3)
      do j = 1, zsize(2)
         do i = 1, zsize(1)
            out(zstart(1)+i-1,zstart(2)+j-1,zstart(3)+k-1) = cmplx(qgridout(1,i,j,k), &
               qgridout(2,i,j,k),kind=t_p)
         end do
      end do
   end do
!
! ===== 3D backward FFT =====
   time0 = mpi_wtime()
   call decomp_2d_fft_3d(out, in, DECOMP_2D_FFT_BACKWARD)
   time1 = mpi_wtime()
!
! fill qgridin array with in
!
   do k = 1, xsize(3)
      do j = 1, xsize(2)
         do i = 1, xsize(1)
            qgridin(1,i,j,k,1) = real(in(xstart(1)+i-1,xstart(2)+j-1,xstart(3)+k-1))
            qgridin(2,i,j,k,1) = aimag(in(xstart(1)+i-1,xstart(2)+j-1,xstart(3)+k-1))
         end do
      end do
   end do
!
   deallocate (in)
   deallocate (out)
   decomp2d_mpi_fcall=1
   if (deb_Path) write(*,'(5x,a)') '<< fft2d_backmpi'
   call timer_exit( timer_ffts,quiet_timers )
end
!
#ifdef _OPENACC
subroutine cufft2d_backmpi(qgridin,qgridout,n1mpimax,n2mpimax,n3mpimax)
   use decomp_2d
   use decomp_2d_fft
   use decomp_2d_cufft
   use inform   ,only: deb_Path
   use fft
   use timestat ,only: timer_enter,timer_exit,timer_ffts,quiet_timers
   use mdevice  ,only: reciprocal_work
   implicit none
   !include 'mpif.h'
   real(t_p) qgridin(2,n1mpimax,n2mpimax,n3mpimax,*)
   real(t_p) qgridout(2,zsize(1),zsize(2),zsize(3))
   real(t_p) time0,time1
   integer, intent(in):: n1mpimax,n2mpimax,n3mpimax
   integer i,j,k,l
   integer zstart1,zstart2,zstart3,xstart1,xstart2,xstart3


   call timer_enter( timer_ffts )
   if (deb_Path) write(*,'(5x,a)') '>> cufft2d_backmpi'
   !
   ! fill out array with qgridout
   !
   if (.not.grid_pointer_assoc) &
      call r2c_grid(qgridout,zsize(1)*zsize(2)*zsize(3),out, &
         reciprocal_work)
   !
   ! ===== 3D backward FFT =====
   !
   call decomp_2d_cufft_3d(out, in, DECOMP_2D_FFT_BACKWARD)
   !
   ! fill qgridin array with in
   !
   if (.not.grid_pointer_assoc) &
      call c2r_grid(in,n1mpimax*n2mpimax*n3mpimax,qgridin, &
         reciprocal_work)

   if (deb_Path) write(*,'(5x,a)') '<< cufft2d_backmpi'
   call timer_exit( timer_ffts,quiet_timers )
end
#endif
!
subroutine malloc_FFtgrid_p()
   use decomp_2d
   use fft ,only: in,out
   implicit none
   if (associated(in)) then
      print*, "FFT grid pointer should not be associated at this point..."
   end if
   if (.not.associated(in))  &
      allocate (in(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
   if (.not.associated(out)) &
      allocate (out(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)))
!$acc enter data create(in,out)
end subroutine
!
subroutine associate_grid_pointer(qgridin,qgridout)
   use decomp_2d
   use decomp_2d_fft
   use fft ,only: in,out,in_cptr,out_cptr,n1mpimax,n2mpimax &
      ,n3mpimax
   use iso_c_binding
   implicit none
   real(t_p),target:: qgridin (*)
   real(t_p),target:: qgridout(*)

   in_cptr = c_loc(qgridin)
   out_cptr = c_loc(qgridout)

   if (associated(in)) then
      !$acc exit data detach(in,out) async
      nullify(in)
      nullify(out)
   end if

   call c_f_pointer( in_cptr, in,(/n1mpimax,n2mpimax,n3mpimax/))
   call c_f_pointer(out_cptr,out,(/zsize(1),zsize(2),zsize(3)/))

   in (xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) => in
   out(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)) => out

!$acc enter data attach(in,out) async
end subroutine
!
subroutine fft_final
   use decomp_2d
   use decomp_2d_fft
   use fft
#ifdef _OPENACC
   use decomp_2d_cufft

   call finalize_cufft_engine
#endif
   call decomp_2d_fft_finalize
   call decomp_2d_finalize
end

!     Subroutine commGridFront: Exchange FFT grid between reciprocal process
!                 before call to fft_frontmpi
!
subroutine commGridFront( qgrid,rule )
   use parallel, only: get_domain_comm, get_recip_comm
   use domdec  ,only: nproc,nrec_recep,nrec_send,rank,prec_send,prec_recep&
      &
   use mdevice ,only: reciprocal_work
   use mpi
   use inform
   use pme     ,only: GridDecomp1d
   use pme1, only: nfloor,qgridmpi,r_comm,r_wait,front_send_req,front_recv_req,&
      &front_split_snd,front_split_rcv
   use potent  ,only: use_pmecore
   use fft
   use timestat,only: timer_enter,timer_exit,timer_recreccomm&
      &,quiet_timers
   implicit none
   integer,intent(in)::rule
   real(t_p) qgrid(2,n1mpimax,n2mpimax,n3mpimax,nrec_send+1)

   integer i,ierr,tag0,tag1,commloc,proc,MPI_TYPE,MSG_sz
   integer mpi_status1(MPI_STATUS_SIZE)
   integer qIdx3
   logical decomp1d
   parameter(tag0=0,tag1=1)
!     character*3 ich

   if (nrec_send.eq.0) return
   if (.not.allocated(front_send_req)) then
      allocate(front_send_req(max(1,nproc)),front_recv_req(max(1,nproc)),&
         &front_split_snd(max(1,nproc)),front_split_rcv(max(1,nproc)))
      front_send_req = MPI_REQUEST_NULL
      front_recv_req = MPI_REQUEST_NULL
      front_split_snd = MPI_REQUEST_NULL
      front_split_rcv = MPI_REQUEST_NULL
   end if
   call timer_enter(timer_recreccomm)
   if (use_pmecore) then
      commloc  = get_recip_comm()
   else
      commloc  = get_domain_comm()
   end if
   decomp1d = GridDecomp1d.and.(2*nfloor.lt.n3mpimax)

   if (rule.eq.r_comm) then
      if (decomp1d) then
         MPI_TYPE = MPI_TPREC
         MSG_sz   = 2*n1mpimax*n2mpimax*nfloor
         qIdx3    = n3mpimax-nfloor+1
      else
         MPI_TYPE = MPI_TPREC
         MSG_sz   = 2*n1mpimax*n2mpimax*n3mpimax
      end if

      !MPI : Begin reception
!$acc host_data use_device(qgridmpi,qgrid)
      do i = 1, nrec_recep
         call MPI_IRECV(qgridmpi(1,1,1,1,i),MSG_sz,MPI_TYPE&
            &,prec_recep(i),tag0,commloc,front_recv_req(i),ierr)
         if (decomp1d)&
            &call MPI_IRECV(qgridmpi(1,1,1,qIdx3,i),MSG_sz,MPI_TYPE&
            &,prec_recep(i),tag1,commloc,front_split_rcv(i),ierr)
      end do

      !MPI : begin sending
!$acc wait(reciprocal_work%queue)
      do i = 1, nrec_send
         call MPI_ISEND(qgrid(1,1,1,1,i+1),MSG_sz,MPI_TYPE&
            &,prec_send(i),tag0,commloc,front_send_req(i),ierr)
         if (decomp1d)&
            &call MPI_ISEND(qgrid(1,1,1,qIdx3,i+1),MSG_sz,MPI_TYPE&
            &,prec_send(i),tag1,commloc,front_split_snd(i),ierr)
      end do
!$acc end host_data

   else if (rule.eq.r_wait) then
      MSG_sz   = 2*n1mpimax*n2mpimax*n3mpimax
      call MPI_Waitall(nrec_recep,front_recv_req,MPI_STATUSES_IGNORE,ierr)
      if (decomp1d) then
         call MPI_Waitall(nrec_recep,front_split_rcv,MPI_STATUSES_IGNORE,ierr)
         call MPI_Waitall(nrec_send ,front_split_snd,MPI_STATUSES_IGNORE,ierr)
      end if
      call MPI_Waitall(nrec_send ,front_send_req,MPI_STATUSES_IGNORE,ierr)

      !Reduction
      do i = 1, nrec_recep
         call aaddgpuAsync(MSG_sz&
            &,qgrid(1,1,1,1,1),qgridmpi(1,1,1,1,i)&
            &,qgrid(1,1,1,1,1))
      end do
   else
12    format("commGridFront !! Unknown configuration",I4)
      if (rank.Eq.0) then
         write(0,12) rule
         call fatal
      end if
   end if
   call timer_exit(timer_recreccomm,quiet_timers)
end subroutine
!
!     Subroutine commGridBack: Exchange FFT grid between reciprocal process
!                 after call to fft_frontmpi
!
subroutine commGridBack( qgrid,rule )
   use parallel, only: get_domain_comm, get_recip_comm
   use domdec  ,only: nproc,nrec_recep,nrec_send,rank,prec_send,prec_recep&
      &
   use mdevice ,only: reciprocal_work
   use mpi
   use pme     ,only: GridDecomp1d
   use pme1, only: nfloor,r_comm,r_wait,back_send_req,back_recv_req,&
      &back_split_snd,back_split_rcv
   use potent  ,only: use_pmecore
   use fft
   use timestat,only: timer_enter,timer_exit,timer_recreccomm&
      &,quiet_timers
   implicit none
   integer,intent(in)::rule
   real(t_p) qgrid(2,n1mpimax,n2mpimax,n3mpimax,nrec_send+1)

   integer i,j,ierr,tag0,tag1,commloc,proc,MPI_TYPE,MSG_szs,MSG_szr
   integer off,off1,qIdx3
   logical decomp1d
   parameter(tag0=0,tag1=1)

   if (nrec_send.eq.0) return
   if (.not.allocated(back_send_req)) then
      allocate(back_send_req(max(1,nproc)),back_recv_req(max(1,nproc)),&
         &back_split_snd(max(1,nproc)),back_split_rcv(max(1,nproc)))
      back_send_req = MPI_REQUEST_NULL
      back_recv_req = MPI_REQUEST_NULL
      back_split_snd = MPI_REQUEST_NULL
      back_split_rcv = MPI_REQUEST_NULL
   end if
   call timer_enter(timer_recreccomm)
   if (use_pmecore) then
      commloc  = get_recip_comm()
   else
      commloc  = get_domain_comm()
   end if
   decomp1d = GridDecomp1d.and.(2*nfloor.lt.n3mpimax)

   if (rule.eq.r_comm) then
      if (decomp1d) then
         MPI_TYPE = MPI_TPREC
         MSG_szs  = 2*n1mpimax*n2mpimax*(nfloor)
         MSG_szr  = MSG_szs
         qIdx3    = n3mpimax-nfloor+1
      else
         MPI_TYPE = MPI_TPREC
         MSG_szs  = 2*n1mpimax*n2mpimax*n3mpimax
         MSG_szr  = MSG_szs
      end if

      !MPI : Begin reception
!$acc    host_data use_device(qgrid)
      do i = 1, nrec_send
         call MPI_IRECV(qgrid(1,1,1,1,i+1),MSG_szs,MPI_TYPE,&
            &prec_send(i),tag0,commloc,back_recv_req(i),ierr)
         if (decomp1d)&
            &call MPI_IRECV(qgrid(1,1,1,qIdx3,i+1),MSG_szs,MPI_TYPE,&
            &prec_send(i),tag1,commloc,back_split_rcv(i),ierr)
      end do

      !MPI : begin sending
!$acc wait(reciprocal_work%queue)
      do i = 1, nrec_recep
         call MPI_ISEND(qgrid(1,1,1,1,1),MSG_szr,MPI_TYPE,&
            &prec_recep(i),tag0,commloc,back_send_req(i),ierr)
         if (decomp1d)&
            &call MPI_ISEND(qgrid(1,1,1,qIdx3,1),MSG_szr,MPI_TYPE,&
            &prec_recep(i),tag1,commloc,back_split_snd(i),ierr)
      end do
!$acc end host_data
   else if (rule.eq.r_wait) then
      call MPI_Waitall(nrec_send,back_recv_req,MPI_STATUSES_IGNORE,ierr)
      if (decomp1d) then
         call MPI_Waitall(nrec_send,back_split_rcv,MPI_STATUSES_IGNORE,ierr)
         call MPI_Waitall(nrec_recep,back_split_snd,MPI_STATUSES_IGNORE,ierr)
      end if
      call MPI_Waitall(nrec_recep,back_send_req,MPI_STATUSES_IGNORE,ierr)
   else
12    format("commGridBack !! Unknown configuration",I4)
      if (rank.Eq.0) then
         write(0,12) rule
         call fatal
      end if
   end if
   call timer_exit(timer_recreccomm,quiet_timers)
end subroutine

subroutine fftTrC_tmat_alloc(nrhs,npolebloc,dodiag,argmu,argef)
   use domdec     ,only: rank
   use decomp_2d
   use interfaces ,only: tmatxb_pmegpu
   implicit none
   integer,intent(in)::nrhs,npolebloc
   logical,intent(in)::dodiag
   real(t_p),dimension(3,nrhs,npolebloc),target::argmu,argef
   procedure(tmatxb_pmegpu) :: fftTrC_tmat
   tin_nrhs   = nrhs
   tin_dodiag = dodiag
   mu => argmu
   ef => argef
   decomp2d_WhileWait1 => fftTrC_tmat
end subroutine

subroutine fftTrC_tmat(nrhs,dodiag,mu,ef)
   use domdec
   use interfaces ,only: tmatxb_cp
   use mdevice, only: direct_work,reciprocal_work
   use mdevice, only: begin_direct_overlap,end_direct_overlap
   use mpole      ,only: npolebloc
   implicit none
   integer,intent(in)::nrhs
   logical,intent(in)::dodiag
   real(t_p),intent(in ),dimension(:,:,:)::mu
   real(t_p),intent(out),dimension(:,:,:)::ef


   if (direct_work%queue.ne.reciprocal_work%queue) then
#ifdef _OPENACC
      call begin_direct_overlap
#endif
      call tmatxb_cp(nrhs,.true.,mu,ef)
#ifdef _OPENACC
      call end_direct_overlap
#endif
   end if
end subroutine

subroutine fftTrC_tmat_free
   use domdec,only:rank
   use decomp_2d
   implicit none
   tin_nrhs   = 0
   tin_dodiag = .false.
   nullify(mu)
   nullify(ef)
   nullify(decomp2d_WhileWait1)
end subroutine
