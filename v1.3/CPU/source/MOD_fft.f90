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
module fft
   use sizes
   implicit none
   integer :: n1mpimax,n2mpimax,n3mpimax
   integer :: ngrid1 !<maximum number of points along the divided first dimension of the fft grid
   integer :: ngrid2 !<maximum number of points along the divided second dimension of the fft grid
   integer, allocatable :: istart1(:) !<list of beginning indexes of the first axis of the grid among the procs
   integer, allocatable :: iend1(:) !<list of ending indexes of the first axis of the grid among the procs
   integer, allocatable :: jstart1(:) !<list of beginning indexes of the second axis of the grid among the procs
   integer, allocatable :: jend1(:) !<list of ending indexes of the second axis of the grid among the procs
   integer, allocatable :: kstart1(:) !<list of beginning indexes of the third axis of the grid among the procs
   integer, allocatable :: kend1(:) !<list of ending indexes of the third axis of the grid among the procs
   integer, allocatable :: isize1(:) !<list of size of the first axis of the grid among the procs
   integer, allocatable :: jsize1(:) !<list of size of the second axis of the grid among the procs
   integer, allocatable :: ksize1(:) !<list of size of the third axis of the grid among the procs
   integer, allocatable :: istart2(:) !<list of beginning indexes of the first axis of the grid among the procs for the transposed grid

   integer, allocatable :: iend2(:) !<list of ending indexes of the first axis of the grid among the procs for the transposed grid
   integer, allocatable :: jstart2(:) !<list of beginning indexes of the second axis of the grid among the procs for the transposed grid

   integer, allocatable :: jend2(:) !<list of ending indexes of the second axis of the grid among the procs for the transposed grid
   integer, allocatable :: kstart2(:) !<list of beginning indexes of the third axis of the grid among the procs for the transposed grid

   integer, allocatable :: kend2(:) !<list of ending indexes of the third axis of the grid among the procs for the transposed grid
   integer, allocatable :: isize2(:) !<list of size of the first axis of the transposed grid among the procs
   integer, allocatable :: jsize2(:) !<list of size of the second axis of the transposed grid among the procs
   integer, allocatable :: ksize2(:) !<list of size of the third axis of the transposed grid among the procs
   save
end
