c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #####################################################################
c     ##                                                                 ##
c     ##  module moldyn  --  velocity and acceleration on MD trajectory  ##
c     ##                                                                 ##
c     #####################################################################
c
c
c     v       current velocity of each atom along the x,y,z-axes
c     a       current acceleration of each atom along x,y,z-axes
c     aalt    alternate acceleration of each atom along x,y,z-axes
c     aalt2   alternate 2 acceleration of each atom along x,y,z-axes
c     dshort  short range (bonded) forces time step for respa and respa1 integrator
c     nalt    number of inner short range time steps for respa integrator and intermediate for respa1
c     dinter    intermediate range forces time step for respa1 integrator
c     nalt2    number of inner short range time steps for respa1 integrator
c
c
#include "tinker_precision.h"
      module moldyn
      implicit none
      real(r_p),allocatable :: v(:,:),a(:,:),aalt(:,:),aalt2(:,:)
      real(r_p) dshort,dinter
      integer nalt,nalt2
      end
