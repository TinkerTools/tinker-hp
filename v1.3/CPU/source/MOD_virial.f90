!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###############################################################
!     ##                                                           ##
!     ##  module virial  --  components of internal virial tensor  ##
!     ##                                                           ##
!     ###############################################################
!
!
!
module virial
   implicit none
   real*8 :: dedv !<derivative of energy with respect to volume
   real*8 :: vir(3,3) !<total internal virial Cartesian tensor components
   real*8 :: virsave(3,3) !<stored internal virial Cartesian tensor components
   logical :: virnum !<stored internal virial Cartesian tensor components computed with finite differences
   logical :: kin_instant !<compute instantaneous kinetic energy
   save
end
