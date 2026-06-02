!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This program computes a distributed 3D FFT using
!   -  2DECOMP&FFT's decomposition
!   -  Global Arrays toolkit for communication
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program ga_test

  use decomp_2d
  use MPI

  implicit none

  ! GA header files
  ! require preprocessing, but not a problem for 2DECOMP applications
#include "mafdecls.fh"
#include "global.fh"

  include "fftw3.f"

  integer, parameter :: nx=16, ny=16, nz=16
  integer, parameter :: p_row=2, p_col=2
  integer, parameter :: NFFT=20   ! number of independent FFTs

  integer :: i,j,k, m, nmax, ierror
  real(mytype) :: tmp1, tmp2
  logical :: status
  double precision :: t1, t2

  integer :: ga1, ga2, ga3

  ! FFTW plans for the 1D forward/backward transforms
  integer*8, save :: x_plan_f, x_plan_b
  integer*8, save :: y_plan_f, y_plan_b
  integer*8, save :: z_plan_f, z_plan_b

  ! dummy array used for planning
  complex(mytype), allocatable, dimension(:,:,:) :: buf1
  complex(mytype), allocatable, dimension(:,:) :: buf2

  ! input/output of the FFT
  complex(mytype), allocatable, dimension(:,:,:) :: in, out, wk2


  ! initialisation
  call MPI_INIT(ierror)
  call ga_initialize()
  call decomp_2d_init(nx,ny,nz,p_row,p_col)

  ! Create global arrays that map to 2DECOMP's pencils
  call get_global_array(ga1, 1, ga_complex_type)
  call get_global_array(ga2, 2, ga_complex_type)
  call get_global_array(ga3, 3, ga_complex_type)
  
  ! ===== planning 1D FFTs in X =====
  allocate(buf1(xsize(1),xsize(2),xsize(3)))

#ifdef DOUBLE_PREC
  call dfftw_plan_many_dft(x_plan_f, 1, xsize(1), &
         xsize(2)*xsize(3), buf1, xsize(1), 1, &
         xsize(1), buf1, xsize(1), 1, xsize(1), &
         FFTW_FORWARD,  FFTW_MEASURE)
  call dfftw_plan_many_dft(x_plan_b, 1, xsize(1), &
         xsize(2)*xsize(3), buf1, xsize(1), 1, &
         xsize(1), buf1, xsize(1), 1, xsize(1), &
         FFTW_BACKWARD, FFTW_MEASURE)
#else
  call sfftw_plan_many_dft(x_plan_f, 1, xsize(1), &
         xsize(2)*xsize(3), buf1, xsize(1), 1, &
         xsize(1), buf1, xsize(1), 1, xsize(1), &
         FFTW_FORWARD,  FFTW_MEASURE)
  call sfftw_plan_many_dft(x_plan_b, 1, xsize(1), &
         xsize(2)*xsize(3), buf1, xsize(1), 1, &
         xsize(1), buf1, xsize(1), 1, xsize(1), &
         FFTW_BACKWARD, FFTW_MEASURE)
#endif

  deallocate(buf1)

  ! ===== planning 1D FFTs in Y =====
  allocate(buf2(ysize(1),ysize(2)))
  
#ifdef DOUBLE_PREC
  call dfftw_plan_many_dft(y_plan_f, 1, ysize(2), ysize(1), &
         buf2, ysize(2), ysize(1), 1, buf2, ysize(2), &
         ysize(1), 1, FFTW_FORWARD,   FFTW_MEASURE)
  call dfftw_plan_many_dft(y_plan_b, 1, ysize(2), ysize(1), &
         buf2, ysize(2), ysize(1), 1, buf2, ysize(2), &
         ysize(1), 1, FFTW_BACKWARD,  FFTW_MEASURE)
#else
  call sfftw_plan_many_dft(y_plan_f, 1, ysize(2), ysize(1), &
         buf2, ysize(2), ysize(1), 1, buf2, ysize(2), &
         ysize(1), 1, FFTW_FORWARD,   FFTW_MEASURE)
  call sfftw_plan_many_dft(y_plan_b, 1, ysize(2), ysize(1), &
         buf2, ysize(2), ysize(1), 1, buf2, ysize(2), &
         ysize(1), 1, FFTW_BACKWARD,  FFTW_MEASURE)
#endif
  
  deallocate(buf2)
  
  ! ===== planning 1D FFTs in Z =====
  allocate(buf1(zsize(1),zsize(2),zsize(3)))
  
#ifdef DOUBLE_PREC
  call dfftw_plan_many_dft(z_plan_f, 1, zsize(3), &
         zsize(1)*zsize(2), buf1, zsize(3), &
         zsize(1)*zsize(2), 1, buf1, zsize(3), &
         zsize(1)*zsize(2), 1, FFTW_FORWARD,  FFTW_MEASURE)
  call dfftw_plan_many_dft(z_plan_b, 1, zsize(3), &
         zsize(1)*zsize(2), buf1, zsize(3), &
         zsize(1)*zsize(2), 1, buf1, zsize(3), &
         zsize(1)*zsize(2), 1, FFTW_BACKWARD, FFTW_MEASURE)
#else
  call sfftw_plan_many_dft(z_plan_f, 1, zsize(3), &
         zsize(1)*zsize(2), buf1, zsize(3), &
         zsize(1)*zsize(2), 1, buf1, zsize(3), &
         zsize(1)*zsize(2), 1, FFTW_FORWARD,  FFTW_MEASURE)
  call sfftw_plan_many_dft(z_plan_b, 1, zsize(3), &
         zsize(1)*zsize(2), buf1, zsize(3), &
         zsize(1)*zsize(2), 1, buf1, zsize(3), &
         zsize(1)*zsize(2), 1, FFTW_BACKWARD, FFTW_MEASURE)
#endif
  
  deallocate(buf1)


  allocate( in(xsize(1),xsize(2),xsize(3)))   ! x-pencil input
  allocate(out(zsize(1),zsize(2),zsize(3)))   ! z-pencil output
  allocate(wk2(ysize(1),ysize(2),ysize(3)))   ! y-pencil intermediate

  t1 = MPI_WTIME() ! start timer

  do m=1,NFFT

     do k=1,xsize(3)
        do j=1,xsize(2)
           do i=1,xsize(1)
              tmp1 = real(xstart(1)+i-1, mytype) / real(nx, mytype) &
                   * real(xstart(2)+j-1, mytype) / real(ny, mytype) &
                   * real(xstart(3)+k-1, mytype) / real(nz, mytype) &
                   * real(m, mytype) / real(NFFT, mytype)
              in(i,j,k) = cmplx(tmp1, 0._mytype, mytype)
           end do
        end do
     end do

     ! 1D FFT in X
#ifdef DOUBLE_PREC
     call dfftw_execute_dft(x_plan_f, in, in)
#else
     call sfftw_execute_dft(x_plan_f, in, in)
#endif

     ! ===== Swap X --> Y =====
     !call transpose_x_to_y(in,wk2)
     call nga_put(ga1, xstart, xend, in, xsize)
     call ga_copy(ga1, ga2)
     call nga_get(ga2, ystart, yend, wk2, ysize)

     ! ===== 1D FFTs in Y =====
     do k=1,ysize(3)
#ifdef DOUBLE_PREC
        call dfftw_execute_dft(y_plan_f, wk2(:,:,k), wk2(:,:,k))
#else
        call sfftw_execute_dft(y_plan_f, wk2(:,:,k), wk2(:,:,k))
#endif
     end do

     ! ===== Swap Y --> Z =====
     !call transpose_y_to_z(wk2,out)
     call nga_put(ga2, ystart, yend, wk2, ysize)
     call ga_copy(ga2, ga3)
     call nga_get(ga3, zstart, zend, out, zsize)

     ! ===== 1D FFTs in Z =====
#ifdef DOUBLE_PREC
     call dfftw_execute_dft(z_plan_f, out, out)
#else
     call sfftw_execute_dft(z_plan_f, out, out)
#endif

     if (nrank==0) write(*,*) 'TEST ', m, out(1:2,1:2,1:2)

  end do ! NFFT

  t2 = MPI_WTIME() - t1
  call MPI_ALLREDUCE(t2,t1,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
       MPI_COMM_WORLD,ierror)
  t1 = t1 / real(nproc, mytype)

  if (nrank==0) then
     write(*,*) 'Average Forward FFT Time(sec): ', t1
  end if

  ! clean up
#ifdef DOUBLE_PREC
  call dfftw_destroy_plan(x_plan_f)
  call dfftw_destroy_plan(x_plan_b)
  call dfftw_destroy_plan(y_plan_f)
  call dfftw_destroy_plan(y_plan_b)
  call dfftw_destroy_plan(z_plan_f)
  call dfftw_destroy_plan(z_plan_b)
#else
  call sfftw_destroy_plan(x_plan_f)
  call sfftw_destroy_plan(x_plan_b)
  call sfftw_destroy_plan(y_plan_f)
  call sfftw_destroy_plan(y_plan_b)
  call sfftw_destroy_plan(z_plan_f)
  call sfftw_destroy_plan(z_plan_b)
#endif

  status = ga_destroy(ga1)
  status = ga_destroy(ga2)
  status = ga_destroy(ga3)

  call decomp_2d_finalize 
  call ga_terminate()
  call MPI_FINALIZE(ierror)
  deallocate(in,out,wk2)
  

end program ga_test
