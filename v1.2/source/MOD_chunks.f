c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ####################################################################
c     ##                                                                ##
c     ##  module chunks  --  values for PME grid spatial decomposition  ##
c     ##                                                                ##
c     ####################################################################
c
c
c     nlpts      PME grid points to the left of center point
c     nrpts      PME grid points to the right of center point
c     grdoff     offset for index into B-spline coefficients
c
c
      module chunks
      implicit none
      integer nlpts,nrpts,grdoff
      save
      end
