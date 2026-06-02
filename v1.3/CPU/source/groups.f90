!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     "groups" tests a set of atoms to see if all are members of a
!     single atom group or a pair of atom groups; if so, then the
!     correct intra- or intergroup weight is assigned
!
!     note the default group-based interaction weight is 1.0; only
!     interactions involving two or fewer groups can be scaled
!
!
!> @brief 
!> tests a set of atoms to see if all are members of a
!> single atom group or a pair of atom groups; if so, then the
!> correct intra- or intergroup weight is assigned
!> @param[in] weigh: weight for the group
!> @param[in] ia: index of the first atom
!> @param[in] ib: index of the second atom
!> @param[in] ic: index of the third atom
!> @param[in] id: index of the fourth atom
!> @param[in] ie: index of the fifth atom
!> @param[in] ig: index of the sixth atom
subroutine groups (weigh,ia,ib,ic,id,ie,ig)
   use group
   implicit none
   integer ia,ib,ic
   integer id,ie,ig
   integer iga,igb,igc
   integer igd,ige,igg
   integer nset
   integer gmax,gmin
   real*8 weigh
!
!
!     determine the number of atoms in the set to be compared
!
   nset = 0
   weigh = 1.0d0
   if (ig .ne. 0) then
      nset = 6
   else if (ie .ne. 0) then
      nset = 5
   else if (id .ne. 0) then
      nset = 4
   else if (ic .ne. 0) then
      nset = 3
   else if (ib .ne. 0) then
      nset = 2
   else if (ia .ne. 0) then
      nset = 1
   end if
!
!     check group membership for a set containing one atom
!
   if (nset .eq. 1) then
      iga = grplist(ia)
      weigh = wgrp(iga,iga)
!
!     check group membership for a set containing two atoms
!
   else if (nset .eq. 2) then
      iga = grplist(ia)
      igb = grplist(ib)
      weigh = wgrp(iga,igb)
!
!     check group membership for a set containing three atoms
!
   else if (nset .eq. 3) then
      iga = grplist(ia)
      igb = grplist(ib)
      igc = grplist(ic)
      if (iga.eq.igb .or. igb.eq.igc) then
         weigh = wgrp(iga,igc)
      else if (iga .eq. igc) then
         weigh = wgrp(iga,igb)
      end if
!
!     check group membership for a set containing four atoms
!
   else if (nset .eq. 4) then
      iga = grplist(ia)
      igb = grplist(ib)
      igc = grplist(ic)
      igd = grplist(id)
      gmin = min(iga,igb,igc,igd)
      gmax = max(iga,igb,igc,igd)
      if ((iga.eq.gmin .or. iga.eq.gmax) .and.&
      &(igb.eq.gmin .or. igb.eq.gmax) .and.&
      &(igc.eq.gmin .or. igc.eq.gmax) .and.&
      &(igd.eq.gmin .or. igd.eq.gmax))  weigh = wgrp(gmin,gmax)
!
!     check group membership for a set containing five atoms
!
   else if (nset .eq. 5) then
      iga = grplist(ia)
      igb = grplist(ib)
      igc = grplist(ic)
      igd = grplist(id)
      ige = grplist(ie)
      gmin = min(iga,igb,igc,igd,ige)
      gmax = max(iga,igb,igc,igd,ige)
      if ((iga.eq.gmin .or. iga.eq.gmax) .and.&
      &(igb.eq.gmin .or. igb.eq.gmax) .and.&
      &(igc.eq.gmin .or. igc.eq.gmax) .and.&
      &(igd.eq.gmin .or. igd.eq.gmax) .and.&
      &(ige.eq.gmin .or. ige.eq.gmax))  weigh = wgrp(gmin,gmax)
!
!     check group membership for a set containing six atoms
!
   else if (nset .eq. 6) then
      iga = grplist(ia)
      igb = grplist(ib)
      igc = grplist(ic)
      igd = grplist(id)
      ige = grplist(ie)
      igg = grplist(ig)
      gmin = min(iga,igb,igc,igd,ige,igg)
      gmax = max(iga,igb,igc,igd,ige,igg)
      if ((iga.eq.gmin .or. iga.eq.gmax) .and.&
      &(igb.eq.gmin .or. igb.eq.gmax) .and.&
      &(igc.eq.gmin .or. igc.eq.gmax) .and.&
      &(igd.eq.gmin .or. igd.eq.gmax) .and.&
      &(ige.eq.gmin .or. ige.eq.gmax) .and.&
      &(igg.eq.gmin .or. igg.eq.gmax))  weigh = wgrp(gmin,gmax)
   end if
!
   return
end
