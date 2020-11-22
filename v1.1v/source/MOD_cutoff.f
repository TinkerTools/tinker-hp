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
c     chgcut      cutoff distance for charge-charge interactions
c     mpolecut    cutoff distance for atomic multipole interactions
c     vdwtaper    distance at which van der Waals switching begins
c     chgtaper    distance at which charge-charge switching begins
c     mpoletaper  distance at which atomic multipole switching begins
c     ewaldcut    cutoff distance for direct space Ewald summation
c     use_ewald   logical flag governing use of Ewald summation
c     use_list    logical flag governing use of any neighbor lists
c     use_vlist   logical flag governing use of vdw neighbor list
c     use_mlist   logical flag governing use of multipole neighbor list
c
c
      module cutoff
      implicit none
      real*8 vdwcut,chgcut
      real*8 mpolecut
      real*8 vdwtaper,chgtaper
      real*8 mpoletaper
      real*8 ewaldcut
      real*8 ddcut
      logical use_ewald,use_lights
      logical use_list,use_vlist
      logical use_mlist,use_clist
      save
      end
