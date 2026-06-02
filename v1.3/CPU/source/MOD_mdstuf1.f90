!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###################################################################
!     ##                                                               ##
!     ##  module mdstuf1  --  control of molecular dynamics trajectory  ##
!     ##                                                               ##
!     ###################################################################
!
!
module mdstuf1
   implicit none
   real*8 :: etot !<holds the system total energy at current timestep
   real*8 :: epot !<holds the system potential energy at current timestep
   real*8 :: eksum !<holds the system kinetic energy at current timestep
   real*8 :: ealt !<holds the alt energy at current timestep (mts integrators)
   real*8 :: ealt2 !<holds the alt2 energy at current timestep (mts integrators)
   real*8 :: eml !<holds the system ml energy at current timestep
   real*8 :: ealtml !<holds the system alt ml energy at current timestep (mts integrators)
   real*8 :: temp !<current temperature
   real*8 :: pres !<current pressure
   real*8 :: ekin(3,3) !<kinetic energy tensor
   real*8 :: stress(3,3) !<stress tensor
   real*8 :: viralt(3,3) !<virial tensor
   real*8 :: viralt2(3,3)!<alt virial tensor (mts integrators)
   real*8,allocatable::derivs(:,:) !<stores forces computed by gradient routines
   save
end module mdstuf1
