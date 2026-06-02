!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #####################################################################
!     ##                                                                 ##
!     ##  module moldyn  --  velocity and acceleration on MD trajectory  ##
!     ##                                                                 ##
!     #####################################################################
!
!
!     v       current velocity of each atom along the x,y,z-axes
!     velx    current velocity of each atom along the x-axes
!     vely    current velocity of each atom along the y-axes
!     velz    current velocity of each atom along the z-axes
!     a       current acceleration of each atom along x,y,z-axes
!     accx    current acceleration of each atom along x-axes
!     accy    current acceleration of each atom along y-axes
!     accz    current acceleration of each atom along z-axes
!     aalt    alternate acceleration of each atom along x,y,z-axes
!     aalt2   alternate 2 acceleration of each atom along x,y,z-axes
!     dshort  short range (bonded) forces time step for respa and respa1 integrator
!     nalt    number of inner short range time steps for respa integrator and intermediate for respa1
!     dinter    intermediate range forces time step for respa1 integrator
!     nalt2    number of inner short range time steps for respa1 integrator
!     stepfast index of inner fast timestep eval (Bonded terms) for multi-timestep integrator
!     stepint  index of inner intermediate timestep for short range
!     step_c   index of outer timestep
!
!
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
