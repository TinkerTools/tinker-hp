!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #####################################################################
!     ##                                                                 ##
!     ##  module fft  --  values and options for Fast Fourier transform  ##
!     ##                                                                 ##
!     #####################################################################
!
!
!     maxtable   maximum size of the FFT table intermediate array
!     maxprime   maximum number of prime factors of FFT dimension
!     n1mpimax   maximum number of points along the divided first dimension of the fft grid
!     n2mpimax   maximum number of points along the divided second dimension of the fft grid
!     n3mpimax   maximum number of points along the divided third dimension of the fft grid
!     ngrid1     first dimension of the 2d proc grid for the 2decomp library (2d pencils decomposition)
!     ngrid2     second dimension of the 2d proc grid for the 2decomp library (2d pencils decomposition)
!     istart1    list of beginning indexes of the first axis of the grid among the procs
!     jstart1    list of beginning indexes of the second axis of the grid among the procs
!     kstart1    list of beginning indexes of the third axis of the grid among the procs
!     iend1    list of ending indexes of the first axis of the grid among the procs
!     jend1    list of ending indexes of the second axis of the grid among the procs
!     kend1    list of ending indexes of the third axis of the grid among the procs
!     isize1   list of size of the first axis of the grid among the procs
!     jsize1   list of size of the second axis of the grid among the procs
!     ksize1   list of size of the third axis of the grid among the procs
!     istart2    list of beginning indexes of the first axis of the grid among the procs for the transposed grid
!     jstart2    list of beginning indexes of the second axis of the grid among the procs for the transposed grid
!     kstart2    list of beginning indexes of the third axis of the grid among the procs for the transposed grid
!     iend2    list of ending indexes of the first axis of the grid among the procs for the transposed grid
!     jend2    list of ending indexes of the second axis of the grid among the procs for the transposed grid
!     kend2    list of ending indexes of the third axis of the grid among the procs for the transposed grid
!     isize2   list of size of the first axis of the transposed grid among the procs
!     jsize2   list of size of the second axis of the transposed grid among the procs
!     ksize2   list of size of the third axis of the transposed grid among the procs
!     input    is X-pencil data
!     output   is Z-pencil data
!     in       temporary complex grid to be send to transform with FFT mpi
!     out      temporary complex grid to be send to transform with FFT mpi
!     in_cptr,out_cptr  C generic pointer used for casting
!     grid_pointer_assoc  decide whether or not in & out should be allocated or associated
!
#include "tinker_macro.h"
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
!$acc declare create(n1mpimax,n2mpimax,n3mpimax, &
!$acc   istart2,iend2,jstart2,jend2,kstart2,kend2, &
!$acc   istart1,iend1,jstart1,jend1,kstart1,kend1, &
!$acc   isize2,jsize2,ksize2)

   interface
      subroutine associate_grid_pointer(qgridin,qgridout)
         real(t_p),target:: qgridin (*)
         real(t_p),target:: qgridout(*)
      end subroutine
   end interface

contains

!
!     Convert real grid to complex grid for fft mpi
!
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
!
!     Convert complex grid to real grid for fft mpi
!
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
