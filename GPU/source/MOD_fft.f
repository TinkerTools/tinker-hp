c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #####################################################################
c     ##                                                                 ##
c     ##  module fft  --  values and options for Fast Fourier transform  ##
c     ##                                                                 ##
c     #####################################################################
c
c
c     maxtable   maximum size of the FFT table intermediate array
c     maxprime   maximum number of prime factors of FFT dimension
c     n1mpimax   maximum number of points along the divided first dimension of the fft grid 
c     n2mpimax   maximum number of points along the divided second dimension of the fft grid 
c     n3mpimax   maximum number of points along the divided third dimension of the fft grid 
c     ngrid1     first dimension of the 2d proc grid for the 2decomp library (2d pencils decomposition)
c     ngrid2     second dimension of the 2d proc grid for the 2decomp library (2d pencils decomposition)
c     istart1    list of beginning indexes of the first axis of the grid among the procs 
c     jstart1    list of beginning indexes of the second axis of the grid among the procs
c     kstart1    list of beginning indexes of the third axis of the grid among the procs
c     iend1    list of ending indexes of the first axis of the grid among the procs
c     jend1    list of ending indexes of the second axis of the grid among the procs
c     kend1    list of ending indexes of the third axis of the grid among the procs 
c     isize1   list of size of the first axis of the grid among the procs
c     jsize1   list of size of the second axis of the grid among the procs
c     ksize1   list of size of the third axis of the grid among the procs
c     istart2    list of beginning indexes of the first axis of the grid among the procs for the transposed grid
c     jstart2    list of beginning indexes of the second axis of the grid among the procs for the transposed grid
c     kstart2    list of beginning indexes of the third axis of the grid among the procs for the transposed grid
c     iend2    list of ending indexes of the first axis of the grid among the procs for the transposed grid
c     jend2    list of ending indexes of the second axis of the grid among the procs for the transposed grid
c     kend2    list of ending indexes of the third axis of the grid among the procs for the transposed grid
c     isize2   list of size of the first axis of the transposed grid among the procs
c     jsize2   list of size of the second axis of the transposed grid among the procs
c     ksize2   list of size of the third axis of the transposed grid among the procs
c     input    is X-pencil data
c     output   is Z-pencil data
c     in       temporary complex grid to be send to transform with FFT mpi
c     out      temporary complex grid to be send to transform with FFT mpi
c     in_cptr,out_cptr  C generic pointer used for casting
c     grid_pointer_assoc  decide whether or not in & out should be allocated or associated
c
#include "tinker_precision.h"
      module fft
      use iso_c_binding,only: c_ptr
      use sizes        ,only: maxfft
      implicit none
      integer maxtable
      integer maxprime
      integer n1mpimax,n2mpimax,n3mpimax
      parameter (maxtable=4*maxfft)
      parameter (maxprime=15)
      integer ngrid1,ngrid2
      logical,parameter:: grid_pointer_assoc=.true.
      logical:: is_fftInit=.false.
      !DIR$ ATTRIBUTES ALIGN:64:: istart1, iend1
      integer, allocatable, target :: istart1(:), iend1(:)
      !DIR$ ATTRIBUTES ALIGN:64:: jstart1, jend1
      integer, allocatable, target :: jstart1(:), jend1(:)
      !DIR$ ATTRIBUTES ALIGN:64:: kstart1, kend1
      integer, allocatable, target :: kstart1(:), kend1(:)
      !DIR$ ATTRIBUTES ALIGN:64:: isize1 ,jsize1,ksize1
      integer, allocatable, target :: isize1(:),jsize1(:),ksize1(:)
      !DIR$ ATTRIBUTES ALIGN:64:: istart2, iend2
      integer, allocatable, target :: istart2(:), iend2(:)
      !DIR$ ATTRIBUTES ALIGN:64:: jstart2, jend2
      integer, allocatable, target :: jstart2(:), jend2(:)
      !DIR$ ATTRIBUTES ALIGN:64:: kstart2, kend2
      integer, allocatable, target :: kstart2(:), kend2(:)
      !DIR$ ATTRIBUTES ALIGN:64:: isize2 ,jsize2,ksize2
      integer, allocatable, target :: isize2(:),jsize2(:),ksize2(:)
      complex(t_p), pointer, contiguous :: in(:,:,:), out(:,:,:)
      type(c_ptr) in_cptr,out_cptr
      save
!$acc declare create(n1mpimax,n2mpimax,n3mpimax,
!$acc&   istart2,iend2,jstart2,jend2,kstart2,kend2,
!$acc&   istart1,iend1,jstart1,jend1,kstart1,kend1,
!$acc&   isize2,jsize2,ksize2)

      interface
        subroutine associate_grid_pointer(qgridin,qgridout)
        real(t_p),target:: qgridin (*)
        real(t_p),target:: qgridout(*)
        end subroutine
      end interface

      contains

c
c     Convert real grid to complex grid for fft mpi
c
      subroutine r2c_grid(real_g,size_grid,cmplx_g,rec_queue)
      implicit none
      complex(t_p),intent(out):: cmplx_g(*)
      real   (t_p),intent(in) :: real_g(*)
      integer     ,intent(in) :: size_grid,rec_queue
      integer i

!$acc parallel loop async(rec_queue) present(cmplx_g,real_g)
      do i = 1, size_grid
         cmplx_g(i) = cmplx(real_g(2*(i-1)+1),real_g(2*i),kind=t_p)
      end do
      end subroutine
c
c     Convert complex grid to real grid for fft mpi
c
      subroutine c2r_grid(cmplx_g,size_grid,real_g,rec_queue)
      implicit none
      complex(t_p),intent(in) :: cmplx_g(*)
      integer     ,intent(in) :: size_grid,rec_queue
      real   (t_p),intent(out):: real_g(*)
      integer i

!$acc parallel loop async(rec_queue) present(cmplx_g,real_g)
      do i = 1, 2*size_grid
         if (btest(i,0)) then
            real_g(i) = real(cmplx_g((i-1)/2+1),kind=t_p)
         else
            real_g(i) = aimag(cmplx_g((i-1)/2+1))
         end if
      end do
      end subroutine

      subroutine free_FFTgrid_p
      implicit none

      if (.not.grid_pointer_assoc.and.associated(in)) then
         !$acc exit data delete(in,out) async
         deallocate(in)
         deallocate(out)
      end if
      end subroutine

      end
