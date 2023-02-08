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
c     velx    current velocity of each atom along the x-axes
c     vely    current velocity of each atom along the y-axes
c     velz    current velocity of each atom along the z-axes
c     a       current acceleration of each atom along x,y,z-axes
c     accx    current acceleration of each atom along x-axes
c     accy    current acceleration of each atom along y-axes
c     accz    current acceleration of each atom along z-axes
c     aalt    alternate acceleration of each atom along x,y,z-axes
c     aalt2   alternate 2 acceleration of each atom along x,y,z-axes
c     dshort  short range (bonded) forces time step for respa and respa1 integrator
c     nalt    number of inner short range time steps for respa integrator and intermediate for respa1
c     dinter    intermediate range forces time step for respa1 integrator
c     nalt2    number of inner short range time steps for respa1 integrator
c     stepfast index of inner fast timestep eval (Bonded terms) for multi-timestep integrator
c     stepint  index of inner intermediate timestep for short range
c     step_c   index of outer timestep
c
c
#include "tinker_macro.h"
      module moldyn
      implicit none
      real(r_p),allocatable,target :: v(:,:),a(:,:),aalt(:,:),aalt2(:,:)
      real(r_p),pointer :: velx(:),vely(:),velz(:)
      real(r_p),pointer :: a_x(:),a_y(:),a_z(:)
      real(r_p),pointer :: aalt_x(:),aalt_y(:),aalt_z(:)
      real(r_p) dshort,dinter
      integer nalt,nalt2
      integer stepfast,stepint,step_c
      end
