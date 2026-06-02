!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas !<cutoff distance for short range  real space Ewald dispersion at Austin
!
!     ###################################################################
!     ##                                                               ##
!     ##  module cutoff  --  cutoff distances for energy interactions  ##
!     ##                                                               ##
!     ###################################################################
!
!
!
module cutoff
   implicit none
   real*8 :: vdwcut !<cutoff distance for van der Waals interactions
   real*8 :: chgcut !<cutoff distance for charge-charge interactions
   real*8 :: mpolecut !<cutoff distance for atomic multipole interactions
   real*8 :: vdwtaper !<distance at which van der Waals switching begins
   real*8 :: shortheal !<healing length for switching short range van der Waals and mpole terms
   real*8 :: chgtaper !<distance at which charge-charge switching begins
   real*8 :: mpoletaper !<distance at which atomic multipole switching begins
   real*8 :: disptaper !<distance at which dispersion switching begins
   real*8 :: ctrntaper !<distance at which charge transfer switching begins
   real*8 :: reptaper !<distance at which repulsion switching begins
   real*8 :: ewaldcut !<cutoff distance for direct space Ewald summation
   real*8 :: ewaldshortcut !<cutoff distance for short range direct space Ewald summation
   real*8 :: vdwshortcut !<cutoff distance for short range van der Waals interactions
   real*8 :: mpoleshortcut !<cutoff distance for short range atomic multipole interactions
   real*8 :: chgshortcut !<cutoff distance for short range charge-charge interactions
   real*8 :: ctrncut !<cutoff distance for charge transfer interactions
   real*8 :: dispcut !<cutoff distance for dispersions interactions
   real*8 :: repcut !<cutoff distance for Pauli repulsions interactions
   real*8 :: ctrnshortcut !<cutoff distance for short range charge transfer interactions
   real*8 :: dispshortcut !<cutoff distance for short range dispersions interactions
   real*8 :: repshortcut !<cutoff distance for short range Pauli repulsions interactions
   real*8 :: dewaldcut !<cutoff distance for real space Ewald dispersion
   real*8 :: dewaldshortcut !<cutoff distance for short range  real space Ewald dispersion
   real*8 :: ddcut !<additional cutoff for domain decomposition
   logical :: use_ewald !<logical flag governing use of Ewald summation
   logical :: use_list !<logical flag governing use of any neighbor lists
   logical :: use_vlist !<logical flag governing use of vdw neighbor list
   logical :: use_mlist !<logical flag governing use of multipole neighbor list
   logical :: use_clist !<logical flag governing use of charge neighbor list
   logical :: use_dlist !<logical flag governing use of dispersion neighbor list
   logical :: use_shortmlist !<logical flag governing use of short range mulitpole neighbor list
   logical :: use_shortclist !<logical flag governing use of short range charge neighbor list
   logical :: use_shortvlist !<logical flag governing use of short range van der Waals neighbor list
   logical :: use_shortdlist !<logical flag governing use of short range dispersion neighbor list
   save
end
