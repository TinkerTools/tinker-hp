c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ####################################################################
c     ##                                                                ##
c     ##  module ewald  --  parameters and options for Ewald summation  ##
c     ##                                                                ##
c     ####################################################################
c
c
c     aewald     Ewald convergence coefficient value (Ang-1)
c     ewald_tiny tiniest Ewald convergence coefficient value (Ang-1)
c     aeewald    Ewald convergence coefficient value (Ang-1) for electrostatics
c     apewald    Ewald convergence coefficient value (Ang-1) for polarization
c     adwald     Ewald convergence coefficient value (Ang-1) for dispersion
c     boundary   Ewald boundary condition; none, tinfoil or vacuum
c
c
#include "tinker_macro.h"
      module ewald
      implicit none
      real(t_p) aewald,adewald,aeewald,apewald,ewald_tiny
      character*7 boundary
      parameter( ewald_tiny=1d-6 )
      end
