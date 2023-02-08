c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     "groups" tests a set of atoms to see if all are members of a
c     single atom group or a pair of atom groups; if so, then the
c     correct intra- or intergroup weight is assigned
c
c     note the default group-based interaction weight is 1.0; only
c     interactions involving two or fewer groups can be scaled
c
#ifndef GROUPS_INC_F
#define GROUPS_INC_F
#include "tinker_macro.h"

      M_subroutine
     &             groups1_inl (weigh,ia,ngrp,grplist,wgrp)
!$acc routine seq
      use sizes ,only: maxgrp
      implicit none
      integer  ,intent(in) :: ia,ngrp
      integer  ,intent(in) :: grplist(*)
      real(t_p),intent(in) :: wgrp(ngrp+1,*)
      real(t_p),intent(out):: weigh
      integer iga
      iga   = grplist(ia)+1
      weigh = wgrp(iga,iga)
      end subroutine

      M_subroutine
     &             groups2_inl (weigh,ia,ib,ngrp,grplist,wgrp)
!$acc routine seq
      use sizes ,only: maxgrp
      implicit none
      integer  ,intent(in) :: ia,ib,ngrp
      integer  ,intent(in) :: grplist(*)
      real(t_p),intent(in) :: wgrp(ngrp+1,*)
      real(t_p),intent(out):: weigh
      weigh = wgrp(grplist(ia)+1,grplist(ib)+1)
      end subroutine

      M_subroutine
     &             groups3_inl (weigh,ia,ib,ic,ngrp,grplist,wgrp)
!$acc routine seq
      use sizes ,only: maxgrp
      implicit none
      integer  ,intent(in) :: ia,ib,ic,ngrp
      integer  ,intent(in) :: grplist(*)
      real(t_p),intent(in) :: wgrp(ngrp+1,*)
      real(t_p),intent(out):: weigh
      integer iga,igb,igc
      weigh = 1.0
      iga   = grplist(ia)
      igb   = grplist(ib)
      igc   = grplist(ic)
      if (iga.eq.igb .or. igb.eq.igc) then
         weigh = wgrp(iga+1,igc+1)
      else if (iga .eq. igc) then
         weigh = wgrp(iga+1,igb+1)
      end if
      end subroutine

      M_subroutine
     &             groups4_inl (weigh,ia,ib,ic,id,ngrp,grplist,wgrp)
!$acc routine seq
      use sizes ,only: maxgrp
      implicit none
      integer  ,intent(in) :: ia,ib,ic,id,ngrp
      integer  ,intent(in) :: grplist(*)
      real(t_p),intent(in) :: wgrp(ngrp+1,*)
      real(t_p),intent(out):: weigh
      integer iga,igb,igc,igd,gmin,gmax
      iga  = grplist(ia)
      igb  = grplist(ib)
      igc  = grplist(ic)
      igd  = grplist(id)
      gmin = min(iga,igb,igc,igd)
      gmax = max(iga,igb,igc,igd)
      weigh= 1.0
      if ((iga.eq.gmin .or. iga.eq.gmax) .and.
     &    (igb.eq.gmin .or. igb.eq.gmax) .and.
     &    (igc.eq.gmin .or. igc.eq.gmax) .and.
     &    (igd.eq.gmin .or. igd.eq.gmax)) weigh = wgrp(gmin+1,gmax+1)
      end subroutine

      subroutine groups5_inl (weigh,ia,ib,ic,id,ie,ngrp,grplist,wgrp)
!$acc routine seq
      use sizes ,only: maxgrp
      implicit none
      integer  ,intent(in ):: ia,ib,ic,id,ie,ngrp
      integer  ,intent(in ):: grplist(*)
      real(t_p),intent(in ):: wgrp(ngrp+1,*)
      real(t_p),intent(out):: weigh

      integer iga,igb,igc,igd,ige
      integer gmax,gmin
c
c     set containing five atoms
c
      iga = grplist(ia)
      igb = grplist(ib)
      igc = grplist(ic)
      igd = grplist(id)
      ige = grplist(ie)
      gmin = min(iga,igb,igc,igd,ige)
      gmax = max(iga,igb,igc,igd,ige)
      if ((iga.eq.gmin .or. iga.eq.gmax) .and.
     &    (igb.eq.gmin .or. igb.eq.gmax) .and.
     &    (igc.eq.gmin .or. igc.eq.gmax) .and.
     &    (igd.eq.gmin .or. igd.eq.gmax) .and.
     &    (ige.eq.gmin .or. ige.eq.gmax)) weigh = wgrp(gmin+1,gmax+1)
c
      end
      M_subroutine
     &           groups6_inl (weigh,ia,ib,ic,id,ie,ig,ngrp,grplist,wgrp)
!$acc routine seq
      use sizes ,only: maxgrp
      implicit none
      integer  ,intent(in ):: ia,ib,ic,id,ie,ig,ngrp
      integer  ,intent(in ):: grplist(*)
      real(t_p),intent(in ):: wgrp(ngrp+1,*)
      real(t_p),intent(out):: weigh

      integer nset
      integer iga,igb,igc
      integer igd,ige,igg
      integer gmax,gmin
c
c     set containing six atoms
c
      iga = grplist(ia)
      igb = grplist(ib)
      igc = grplist(ic)
      igd = grplist(id)
      ige = grplist(ie)
      igg = grplist(ig)
      gmin = min(iga,igb,igc,igd,ige,igg)
      gmax = max(iga,igb,igc,igd,ige,igg)
      if ((iga.eq.gmin .or. iga.eq.gmax) .and.
     &    (igb.eq.gmin .or. igb.eq.gmax) .and.
     &    (igc.eq.gmin .or. igc.eq.gmax) .and.
     &    (igd.eq.gmin .or. igd.eq.gmax) .and.
     &    (ige.eq.gmin .or. ige.eq.gmax) .and.
     &    (igg.eq.gmin .or. igg.eq.gmax)) weigh = wgrp(gmin+1,gmax+1)
c
      end

      subroutine groups_inl (weigh,ia,ib,ic,id,ie,ig,ngrp,grplist,wgrp)
!$acc routine seq
      use sizes ,only: maxgrp
      implicit none
      integer,intent(in):: ia,ib,ic,id,ie,ig,ngrp
      integer,intent(in):: grplist(*)
      real(t_p),intent(in) :: wgrp(ngrp+1,*)
      real(t_p),intent(out):: weigh

      integer nset
      integer iga,igb,igc
      integer igd,ige,igg
      integer gmax,gmin
c
c
c     determine the number of atoms in the set to be compared
c
      nset = 0
      weigh = 1.0
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
c
c     check group membership for a set containing one atom
c
      if (nset .eq. 1) then
         iga = grplist(ia)+1
         weigh = wgrp(iga,iga)
c
c     check group membership for a set containing two atoms
c
      else if (nset .eq. 2) then
         iga = grplist(ia)+1
         igb = grplist(ib)+1
         weigh = wgrp(iga,igb)
c
c     check group membership for a set containing three atoms
c
      else if (nset .eq. 3) then
         iga = grplist(ia)
         igb = grplist(ib)
         igc = grplist(ic)
         if (iga.eq.igb .or. igb.eq.igc) then
            weigh = wgrp(iga+1,igc+1)
         else if (iga .eq. igc) then
            weigh = wgrp(iga+1,igb+1)
         end if
c
c     check group membership for a set containing four atoms
c
      else if (nset .eq. 4) then
         iga = grplist(ia)
         igb = grplist(ib)
         igc = grplist(ic)
         igd = grplist(id)
         gmin = min(iga,igb,igc,igd)
         gmax = max(iga,igb,igc,igd)
         if ((iga.eq.gmin .or. iga.eq.gmax) .and.
     &       (igb.eq.gmin .or. igb.eq.gmax) .and.
     &       (igc.eq.gmin .or. igc.eq.gmax) .and.
     &       (igd.eq.gmin .or. igd.eq.gmax)) weigh = wgrp(gmin+1,gmax+1)
c
c     check group membership for a set containing five atoms
c
      else if (nset .eq. 5) then
         iga = grplist(ia)
         igb = grplist(ib)
         igc = grplist(ic)
         igd = grplist(id)
         ige = grplist(ie)
         gmin = min(iga,igb,igc,igd,ige)
         gmax = max(iga,igb,igc,igd,ige)
         if ((iga.eq.gmin .or. iga.eq.gmax) .and.
     &       (igb.eq.gmin .or. igb.eq.gmax) .and.
     &       (igc.eq.gmin .or. igc.eq.gmax) .and.
     &       (igd.eq.gmin .or. igd.eq.gmax) .and.
     &       (ige.eq.gmin .or. ige.eq.gmax)) weigh = wgrp(gmin+1,gmax+1)
c
c     check group membership for a set containing six atoms
c
      else if (nset .eq. 6) then
         iga = grplist(ia)
         igb = grplist(ib)
         igc = grplist(ic)
         igd = grplist(id)
         ige = grplist(ie)
         igg = grplist(ig)
         gmin = min(iga,igb,igc,igd,ige,igg)
         gmax = max(iga,igb,igc,igd,ige,igg)
         if ((iga.eq.gmin .or. iga.eq.gmax) .and.
     &       (igb.eq.gmin .or. igb.eq.gmax) .and.
     &       (igc.eq.gmin .or. igc.eq.gmax) .and.
     &       (igd.eq.gmin .or. igd.eq.gmax) .and.
     &       (ige.eq.gmin .or. ige.eq.gmax) .and.
     &       (igg.eq.gmin .or. igg.eq.gmax)) weigh = wgrp(gmin+1,gmax+1)
      end if
c
      end
#endif
