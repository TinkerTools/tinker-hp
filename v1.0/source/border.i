c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###########################################################
c     ##                                                       ##
c     ##  border.i  --  bond orders for a conjugated pisystem  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     pbpl   pi-bond orders for bonds in "planar" pisystem
c     pnpl   pi-bond orders for bonds in "nonplanar" pisystem
c
c
      real*8 pbpl,pnpl
      common /border/ pbpl(maxbnd),pnpl(maxbnd)
