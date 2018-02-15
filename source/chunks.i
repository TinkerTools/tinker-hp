c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################
c     ##                                                           ##
c     ##  chunks.i  --  values for PME grid spatial decomposition  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     nlpts      PME grid points to the left of center point
c     nrpts      PME grid points to the right of center point
c     grdoff     offset for index into B-spline coefficients
c
c
      integer nlpts,nrpts,grdoff
      common /chunks/ nlpts,nrpts,grdoff
