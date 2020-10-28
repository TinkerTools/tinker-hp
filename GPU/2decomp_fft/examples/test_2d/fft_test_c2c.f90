!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Main test program for the FFT interface
!  - use input data from a NAG FFT library for validation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program fft_test_c2c

use decomp_2d
use decomp_2d_fft

implicit none

integer, parameter :: nx=128, ny=128, nz=128
integer, parameter :: p_row=0, p_col=0

complex(mytype), allocatable, dimension(:,:,:) :: in, out

complex(mytype), dimension(nx,ny,nz) :: in1, out1
integer :: ierror, i,j,k
real*8 time0,time1,mpi_wtime

interface
   subroutine assemble_global(ndir,local,global,nx,ny,nz)
     use decomp_2d
     integer, intent(IN) :: ndir
     integer, intent(IN) :: nx,ny,nz
     complex(mytype), dimension(:,:,:), intent(IN) :: local
     complex(mytype), dimension(nx,ny,nz), intent(OUT) :: global
   end subroutine assemble_global
end interface

call MPI_INIT(ierror)
call decomp_2d_init(nx,ny,nz,p_row,p_col)
if (nrank.eq.0) then
write(*,*) 'xstart = ',xstart
write(*,*) 'ystart = ',ystart
write(*,*) 'zstart = ',zstart
write(*,*) 'p_row = ',p_row
write(*,*) 'p_col = ',p_col
end if
call decomp_2d_fft_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! (1) Testing the complex-to-complex interface (c2c) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  input is X-pencil data
! output is Z-pencil data
allocate (in(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
allocate (out(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)))


! each processor gets its local portion of global data
do k=xstart(3),xend(3)
   do j=xstart(2),xend(2)
      do i=xstart(1),xend(1)
         in(i,j,k) = 1.0d0
      end do
   end do
end do

! ===== 3D forward FFT =====
time0 = mpi_wtime()
call decomp_2d_fft_3d(in, out, DECOMP_2D_FFT_FORWARD)
time1 = mpi_wtime()
if (nrank.eq.0) write(*,*) 'temps fft forward = ',time1 - time0

! normalisation - note FFTW doesn't normalise 
do k=zstart(3),zend(3)
   do j=zstart(2),zend(2)
      do i=zstart(1),zend(1)
         out(i,j,k) = out(i,j,k) / sqrt(real(nx*ny*nz))
      end do
   end do
end do

call assemble_global(3,out,out1,nx,ny,nz)

! ===== 3D inverse FFT =====
time0 = mpi_wtime()
call decomp_2d_fft_3d(out, in, DECOMP_2D_FFT_BACKWARD)
time1 = mpi_wtime()
if (nrank.eq.0) write(*,*) 'temps fft backward = ',time1 - time0

! normalisation - note FFTW doesn't normalise
do k=xstart(3),xend(3)
   do j=xstart(2),xend(2)
      do i=xstart(1),xend(1)
         in(i,j,k) = in(i,j,k) / sqrt(real(nx*ny*nz))
      end do
   end do
end do

call assemble_global(1,in,in1,nx,ny,nz)


call decomp_2d_fft_finalize
call decomp_2d_finalize
call MPI_FINALIZE(ierror)

end program fft_test_c2c


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Collect data from each processor and assemble into a global array
! at the master rank
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine assemble_global(ndir,local,global,nx,ny,nz)
  
  use decomp_2d
  use MPI
  
  implicit none
  
  integer, intent(IN) :: ndir  ! 1 = X-pencil; 3 = Z-pencil
  integer, intent(IN) :: nx,ny,nz
  complex(mytype), dimension(:,:,:), intent(IN) :: local
  complex(mytype), dimension(nx,ny,nz), intent(OUT) :: global
  
  complex(mytype), allocatable, dimension(:,:,:) :: rbuf
  integer, dimension(9) :: sbuf1, rbuf1
  
  integer :: ierror, i,j,k,m, i1,i2,j1,j2,k1,k2, count
  integer, dimension(MPI_STATUS_SIZE) :: status
  
  if (nrank==0) then
     ! master writes its own data to a global array
     if (ndir==3) then  ! Z-pencil 
        i1 = zstart(1)
        i2 = zend(1)
        j1 = zstart(2)
        j2 = zend(2)
        k1 = zstart(3)
        k2 = zend(3)
     else if (ndir==1) then  ! X-pencil
        i1 = xstart(1)
        i2 = xend(1)
        j1 = xstart(2)
        j2 = xend(2)
        k1 = xstart(3)
        k2 = xend(3)
     end if
     do k=k1,k2
        do j=j1,j2
           do i=i1,i2
              ! 'local' is assumbed shape array
              ! but it is OK as starting index for rank 0 always 1
              global(i,j,k)=local(i,j,k)
           end do
        end do
     end do
     ! then loop through all other ranks to collect data
     do m=1,nproc-1
        CALL MPI_RECV(rbuf1,9,MPI_INTEGER,m,m,MPI_COMM_WORLD, &
             status,ierror)
        allocate(rbuf(rbuf1(1):rbuf1(2),rbuf1(4):rbuf1(5), &
             rbuf1(7):rbuf1(8)))
        CALL MPI_RECV(rbuf,rbuf1(3)*rbuf1(6)*rbuf1(9),complex_type,m, &
             m+nproc,MPI_COMM_WORLD,status,ierror)
        do k=rbuf1(7),rbuf1(8)
           do j=rbuf1(4),rbuf1(5)
              do i=rbuf1(1),rbuf1(2)
                 global(i,j,k)=rbuf(i,j,k)
              end do
           end do
        end do
        deallocate(rbuf)
     end do
  else
     ! slaves send data to mater
     if (ndir==3) then  ! Z-pencil
        sbuf1(1) = zstart(1)
        sbuf1(2) = zend(1)
        sbuf1(3) = zsize(1)
        sbuf1(4) = zstart(2)
        sbuf1(5) = zend(2)
        sbuf1(6) = zsize(2)
        sbuf1(7) = zstart(3)
        sbuf1(8) = zend(3)
        sbuf1(9) = zsize(3)
        count = zsize(1)*zsize(2)*zsize(3)
     else if (ndir==1) then  ! X-pencil
        sbuf1(1) = xstart(1)
        sbuf1(2) = xend(1)
        sbuf1(3) = xsize(1)
        sbuf1(4) = xstart(2)
        sbuf1(5) = xend(2)
        sbuf1(6) = xsize(2)
        sbuf1(7) = xstart(3)
        sbuf1(8) = xend(3)
        sbuf1(9) = xsize(3)
        count = xsize(1)*xsize(2)*xsize(3)
     end if
     ! send partition information
     CALL MPI_SEND(sbuf1,9,MPI_INTEGER,0,nrank,MPI_COMM_WORLD,ierror)
     ! send data array
     CALL MPI_SEND(local,count,complex_type,0, &
          nrank+nproc,MPI_COMM_WORLD,ierror)
  end if
  
  return
end subroutine assemble_global


