!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!
!  subroutine fft_setup: setup the arrays for the computations of 
!  3d ffts with 2decomp&FFT
!
#include "tinker_precision.h"
subroutine fft_setup(nrec,rank_bis,istart1,iend1,isize1,jstart1, &
jend1,jsize1,kstart1,kend1,ksize1,istart2,iend2,isize2,jstart2,     &
jend2,jsize2,kstart2,kend2,ksize2,nx,ny,nz,comm_loc,ngrid1,ngrid2)

use decomp_2d
use decomp_2d_fft
use fft      ,only: in, out, is_fftInit
use mpi
use tinMemory,only: s_cufft
use sizes    ,only: tinkerdebug
#ifdef _OPENACC
use utilgpu  ,only: rec_queue,dir_queue
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

p_row = ngrid1 
p_col = ngrid2

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
call get_rec_dir_queue(rec_queue,dir_queue)
call init_cufft_engine
s_cufft = decomp2d_cufftGetSize()
#endif
is_fftInit = .true.

if (tinkerdebug.gt.0) then
   if (rank_bis.eq.0) write(*,*) '--- fft_setup'
   do i = 0,nrec-1
      if (rank_bis.eq.i) then
13       format(I4,'in ',6I8)
         write(*,13) i, xstart(1),xend(1),xstart(2),xend(2),xstart(3),xend(3)
      end if
      call MPI_Barrier(comm_loc,k)
   end do
   do i = 0,nrec-1
      if (rank_bis.eq.i) then
14       format(I4,'out',6I8)
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
  use utilgpu, only :rec_queue
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
     call r2c_grid(qgridin,n1mpimax*n2mpimax*n3mpimax,in,rec_queue)
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
     call c2r_grid(out,zsize(1)*zsize(2)*zsize(3),qgridout,rec_queue)
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
  use utilgpu, only:rec_queue
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
     call r2c_grid(qgridout,zsize(1)*zsize(2)*zsize(3),out,rec_queue)
  !
  ! ===== 3D backward FFT =====
  !time0 = mpi_wtime()
  call decomp_2d_cufft_3d(out, in, DECOMP_2D_FFT_BACKWARD)
  !time1 = mpi_wtime()
  !
  ! fill qgridin array with in
  !
  if (.not.grid_pointer_assoc) &
     call c2r_grid(in,n1mpimax*n2mpimax*n3mpimax,qgridin,rec_queue)

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
