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
c     disptaper  distance at which dispersion switching begins
c     ewaldcut    cutoff distance for direct space Ewald summation
c     ddcut       additional cutoff distance for parallelism
c     repcut      cutoff distance for repulsion interactions
c     dispcut     cutoff distance for dispersion interactions
c     ctransfercut  cutoff distance for charge transfer interactions
c     use_ewald   logical flag governing use of Ewald summation
c     use_list    logical flag governing use of any neighbor lists
c     use_vlist   logical flag governing use of vdw neighbor list
c     use_mlist   logical flag governing use of multipole neighbor list
c     use_replist logical flag governing use of repulsion neighbor list
c     use_displist logical flag governing use of dispersion neighbor list
c     use_ctransferlist logical flag governing use of charge transfer neighbor list
c
c
      module cutoff
      implicit none
      real*8 vdwcut,chgcut
      real*8 mpolecut
      real*8 vdwtaper,chgtaper
      real*8 mpoletaper,disptaper
      real*8 ewaldcut
      real*8 ddcut
      real*8 repcut,dispcut,ctransfercut
      real*8 mpolectcut
      logical use_ewald,use_lights
      logical use_list,use_vlist
      logical use_mlist,use_clist
      logical use_replist,use_displist,use_ctransferlist
      save
      end
