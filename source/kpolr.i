c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #############################################################
c     ##                                                         ##
c     ##  kpolr.i  --  forcefield parameters for polarizability  ##
c     ##                                                         ##
c     #############################################################
c
c
c     polr   dipole polarizability parameters for each atom type
c     athl   Thole polarizability damping value for each atom type
c     pgrp   connected types in polarization group of each atom type
c
c
      integer pgrp
      real*8 polr,athl
      common /kpolr/ polr(maxtyp),athl(maxtyp),pgrp(maxval,maxtyp)
