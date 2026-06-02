!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###################################################################
!     ##                                                               ##
!     ##  module cutoff  --  cutoff distances for energy interactions  ##
!     ##                                                               ##
!     ###################################################################
!
!
!     vdwcut      cutoff distance for van der Waals interactions
!     vdwshortcut  cutoff distance for short range direct space Ewald summation
!     chgcut      cutoff distance for charge-charge interactions
!     chgshortcut   cutoff distance for short range charge-charge interactions
!     mpolecut    cutoff distance for atomic multipole interactions
!     mpoleshortcut    cutoff distance for short range atomic multipole interactions
!     repcut    cutoff distance for Pauli repulsions interactions
!     repshortcut    cutoff distance for short range Pauli repulsions interactions
!     dispcut    cutoff distance for dispersions interactions
!     dispshortcut    cutoff distance for short range dispersions interactions
!     ctrncut    cutoff distance for charge transfer interactions
!     ctrnshortcut    cutoff distance for short range charge transfer interactions
!     vdwtaper    distance at which van der Waals switching begins
!     shortheal   healing length for switching short range van der Waals and mpole terms
!     chgtaper    distance at which charge-charge switching begins
!     mpoletaper  distance at which atomic multipole switching begins
!     ctrntaper    distance at which charge transfer switching begins
!     disptaper    distance at which dispersion switching begins
!     reptaper    distance at which repulsion switching begins
!     ewaldcut    cutoff distance for direct space Ewald summation
!     ewaldshortcut  cutoff distance for short range direct space Ewald summation
!     dewaldcut   cutoff distance for real space Ewald dispersion
!     dewaldshortcut   cutoff distance for short range  real space Ewald dispersion
!     use_ewald   logical flag governing use of Ewald summation
!     use_list    logical flag governing use of any neighbor lists
!     use_vlist   logical flag governing use of vdw neighbor list
!     use_mlist   logical flag governing use of multipole neighbor list
!     use_dlist   logical flag governing use of dispersion neighbor list
!     use_shortmlist   logical flag governing use of short range mulitpole neighbor list
!     use_shortclist   logical flag governing use of short range charge neighbor list
!     use_shortvlist   logical flag governing use of short range vdw neighbor list
!     use_shortdlist   logical flag governing use of short range dispersion neighbor list
!
!
#include "tinker_macro.h"
module cutoff
   implicit none
   real(t_p) vdwcut,chgcut
   real(t_p) mpolecut
   real(t_p) vdwtaper,shortheal,chgtaper
   real(t_p) mpoletaper,disptaper,ctrntaper,reptaper
   real(t_p) ewaldcut,ewaldshortcut
   real(t_p) vdwshortcut,mpoleshortcut,chgshortcut
   real(t_p) ctrncut,dispcut,repcut
   real(t_p) ctrnshortcut,dispshortcut,repshortcut
   real(t_p) dewaldcut,dewaldshortcut
   real(t_p) ddcut
   logical use_ewald,use_lights
   logical use_list,use_vlist
   logical use_mlist,use_clist
   logical use_dlist
   logical use_shortmlist,use_shortclist,use_shortvlist
   logical use_shortdlist
end
