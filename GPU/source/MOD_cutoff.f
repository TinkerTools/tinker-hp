c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###################################################################
c     ##                                                               ##
c     ##  module cutoff  --  cutoff distances for energy interactions  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     vdwcut      cutoff distance for van der Waals interactions
c     vdwshortcut  cutoff distance for short range direct space Ewald summation
c     chgcut      cutoff distance for charge-charge interactions
c     chgshortcut   cutoff distance for short range charge-charge interactions
c     mpolecut    cutoff distance for atomic multipole interactions
c     mpoleshortcut    cutoff distance for short range atomic multipole interactions
c     vdwtaper    distance at which van der Waals switching begins
c     shortheal   healing length for switching short range van der Waals and mpole terms
c     chgtaper    distance at which charge-charge switching begins
c     mpoletaper  distance at which atomic multipole switching begins
c     ewaldcut    cutoff distance for direct space Ewald summation
c     ewaldshortcut  cutoff distance for short range direct space Ewald summation
c     use_ewald   logical flag governing use of Ewald summation
c     use_list    logical flag governing use of any neighbor lists
c     use_vlist   logical flag governing use of vdw neighbor list
c     use_mlist   logical flag governing use of multipole neighbor list
c     use_shortmlist   logical flag governing use of short range mulitpole neighbor list
c     use_shortclist   logical flag governing use of short range charge neighbor list
c     use_shortvlist   logical flag governing use of short range vdw neighbor list
c
c
#include "tinker_precision.h"
      module cutoff
      implicit none
      real(t_p) vdwcut,chgcut
      real(t_p) mpolecut
      real(t_p) vdwtaper,shortheal,chgtaper
      real(t_p) mpoletaper
      real(t_p) ewaldcut,ewaldshortcut
      real(t_p) vdwshortcut,mpoleshortcut,chgshortcut
      real(t_p) ddcut
      logical use_ewald,use_lights
      logical use_list,use_vlist
      logical use_mlist,use_clist
      logical use_shortmlist,use_shortclist,use_shortvlist
      save
      end
