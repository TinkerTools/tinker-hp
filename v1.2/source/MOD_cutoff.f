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
c     repcut    cutoff distance for Pauli repulsions interactions
c     repshortcut    cutoff distance for short range Pauli repulsions interactions
c     dispcut    cutoff distance for dispersions interactions
c     dispshortcut    cutoff distance for short range dispersions interactions
c     ctrncut    cutoff distance for charge transfer interactions
c     ctrnshortcut    cutoff distance for short range charge transfer interactions
c     vdwtaper    distance at which van der Waals switching begins
c     shortheal   healing length for switching short range van der Waals and mpole terms
c     chgtaper    distance at which charge-charge switching begins
c     mpoletaper  distance at which atomic multipole switching begins
c     ctrntaper    distance at which charge transfer switching begins
c     disptaper    distance at which dispersion switching begins
c     reptaper    distance at which repulsion switching begins
c     ewaldcut    cutoff distance for direct space Ewald summation
c     ewaldshortcut  cutoff distance for short range direct space Ewald summation
c     dewaldcut   cutoff distance for real space Ewald dispersion
c     dewaldshortcut   cutoff distance for short range  real space Ewald dispersion
c     use_ewald   logical flag governing use of Ewald summation
c     use_list    logical flag governing use of any neighbor lists
c     use_vlist   logical flag governing use of vdw neighbor list
c     use_mlist   logical flag governing use of multipole neighbor list
c     use_dlist   logical flag governing use of dispersion neighbor list
c     use_shortmlist   logical flag governing use of short range mulitpole neighbor list
c     use_shortclist   logical flag governing use of short range charge neighbor list
c     use_shortvlist   logical flag governing use of short range vdw neighbor list
c     use_shortdlist   logical flag governing use of short range dispersion neighbor list
c
c
      module cutoff
      implicit none
      real*8 vdwcut,chgcut
      real*8 mpolecut
      real*8 vdwtaper,shortheal,chgtaper
      real*8 mpoletaper,disptaper,ctrntaper,reptaper
      real*8 ewaldcut,ewaldshortcut
      real*8 vdwshortcut,mpoleshortcut,chgshortcut
      real*8 ctrncut,dispcut,repcut
      real*8 ctrnshortcut,dispshortcut,repshortcut
      real*8 dewaldcut,dewaldshortcut
      real*8 ddcut
      logical use_ewald,use_lights
      logical use_list,use_vlist
      logical use_mlist,use_clist
      logical use_dlist
      logical use_shortmlist,use_shortclist,use_shortvlist
      logical use_shortdlist
      save
      end
