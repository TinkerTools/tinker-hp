!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!
!  subroutine fft_setup: setup the arrays for the computations of 
!  3d ffts with 2decomp&FFT
!
subroutine fft_setup(nrec,rank_bis,istart1,iend1,isize1,jstart1, &
jend1,jsize1,kstart1,kend1,ksize1,istart2,iend2,isize2,jstart2,  &
jend2,jsize2,kstart2,kend2,ksize2,nx,ny,nz,comm_loc,ngrid1,ngrid2)

use decomp_2d
use decomp_2d_fft

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

complex(KIND(0.0D0)), allocatable, dimension(:,:,:) :: in, out

complex(KIND(0.0D0)), dimension(nx,ny,nz) :: in1, out1
real*8 time0,time1,mpi_wtime

p_row = ngrid1 
p_col = ngrid2
call decomp_2d_init(nx,ny,nz,p_row,p_col,comm_loc)
!
istart1(rank_bis+1) = xstart(1)
iend1(rank_bis+1) = xend(1)
isize1(rank_bis+1) = xsize(1)
jstart1(rank_bis+1) = xstart(2)
jend1(rank_bis+1) = xend(2)
jsize1(rank_bis+1) = xsize(2)
kstart1(rank_bis+1) = xstart(3)
kend1(rank_bis+1) = xend(3)
ksize1(rank_bis+1) = xsize(3)
!
istart2(rank_bis+1) = zstart(1)
iend2(rank_bis+1) = zend(1)
isize2(rank_bis+1) = zsize(1)
jstart2(rank_bis+1) = zstart(2)
jend2(rank_bis+1) = zend(2)
jsize2(rank_bis+1) = zsize(2)
kstart2(rank_bis+1) = zstart(3)
kend2(rank_bis+1) = zend(3)
ksize2(rank_bis+1) = zsize(3)
call decomp_2d_fft_init
return
end
!
!
!  subroutine fft2d_frontmpi: compute forward 3dFFT in parallel with 2decomp&FFT
!
subroutine fft2d_frontmpi(qgridin,qgridout,n1mpimax,n2mpimax,n3mpimax)
use decomp_2d
use decomp_2d_fft
implicit none
include 'mpif.h'
complex(KIND(0.0D0)), allocatable, dimension(:,:,:) :: in, out
real*8 qgridin(2,n1mpimax,n2mpimax,n3mpimax,*)
real*8 qgridout(2,zsize(1),zsize(2),zsize(3))
real*8 time0,time1
integer i,j,k,n1mpimax,n2mpimax,n3mpimax
! input is X-pencil data
! output is Z-pencil data
allocate (in(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
allocate (out(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)))
!
! fill in array with qgridin
!
do k = 1, xsize(3)
  do j = 1, xsize(2)
    do i = 1, xsize(1)
      in(xstart(1)+i-1,xstart(2)+j-1,xstart(3)+k-1) = dcmplx(qgridin(1,i,j,k,1), &
   qgridin(2,i,j,k,1))
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
return
end
!
!
!  subroutine fft2d_frontmpi: compute backward 3dFFT in parallel with 2decomp&FFT
!
subroutine fft2d_backmpi(qgridin,qgridout,n1mpimax,n2mpimax,n3mpimax)
use decomp_2d
use decomp_2d_fft
implicit none
include 'mpif.h'
complex(KIND(0.0D0)), allocatable, dimension(:,:,:) :: in, out
real*8 qgridin(2,n1mpimax,n2mpimax,n3mpimax,*)
real*8 qgridout(2,zsize(1),zsize(2),zsize(3))
real*8 time0,time1
integer i,j,k,n1mpimax,n2mpimax,n3mpimax
! input is Z-pencil data
! output is X-pencil data
allocate (in(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
allocate (out(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)))
!
! fill out array with qgridout
!
do k = 1, zsize(3)
  do j = 1, zsize(2)
    do i = 1, zsize(1)
      out(zstart(1)+i-1,zstart(2)+j-1,zstart(3)+k-1) = dcmplx(qgridout(1,i,j,k), &
   qgridout(2,i,j,k))
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
return
end
!
subroutine fft_final
use decomp_2d
use decomp_2d_fft
call decomp_2d_fft_finalize
call decomp_2d_finalize
return
end 
