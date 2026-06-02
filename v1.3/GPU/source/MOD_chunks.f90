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
!
!     nlpts      PME grid points to the left of center point
!     nrpts      PME grid points to the right of center point
!     grdoff     offset for index into B-spline coefficients
!
!
module chunks
   implicit none
   integer nlpts,nrpts,grdoff
!$acc declare create(nlpts,nrpts,grdoff)
   save
end
