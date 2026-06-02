!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ####################################################################
!     ##                                                                ##
!     ##  module chunks  --  values for PME grid spatial decomposition  ##
!     ##                                                                ##
!     ####################################################################
!
module chunks
   implicit none
   integer :: nlpts !<PME grid points to the left of center point
   integer :: nrpts !<PME grid points to the right of center point
   integer :: grdoff !<offset for index into B-spline coefficients
   save
end
