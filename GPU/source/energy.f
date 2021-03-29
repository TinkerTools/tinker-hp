c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #############################################################
c     ##                                                         ##
c     ##  function energy  --  evaluates energy terms and total  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "energy" calls the subroutines to calculate the potential
c     energy terms and sums up to form the total energy
c
c
#include "tinker_precision.h"

      module energy_inl
      contains
#include "convert.f.inc"
      end module

      function energy ()
      use action
      use energy_inl
      use sizes
      use domdec,only: rank
      use energi
      use inform,only: deb_Path,abort
      use iounit
      use potent
      use tinheader
      use utilgpu,only: inf_r
      use vdwpot
      implicit none
      real(r_p) energy
      logical tinker_isnan_m
      logical ::first_in=.true.

      if (deb_Path) write(*,*) "energy"

      if (first_in) then
         call create_action_data_ondevice
         first_in=.false.
      end if
c
c     zero out each of the potential energy components
c
!$acc data present(eb,eba,eub,eopb,et,ept,eat,ett,ebt,ea
!$acc&      ,eaa,eopd,eid,eit,ec,ev,em,ep,eg,ex,esum
!$acc&      ,ev_r,ec_r,em_r,ep_r,eb_r
!$acc&      ,nev,nec,nem,nep,nev_,nec_,nem_,nep_)

!$acc serial async
      eb   = 0.0_re_p
      eb_r = 0
      ea   = 0.0_re_p
      eba  = 0.0_re_p
      eub  = 0.0_re_p
      eaa  = 0.0_re_p
      eopb = 0.0_re_p
      eopd = 0.0_re_p
      eid  = 0.0_re_p
      eit  = 0.0_re_p
      et   = 0.0_re_p
      ept  = 0.0_re_p
      ebt  = 0.0_re_p
      eat  = 0.0_re_p
      ett  = 0.0_re_p
      ev   = 0.0_re_p
      ev_r = 0
      ec   = 0.0_re_p
      ec_r = 0
      em   = 0.0_re_p
      em_r = 0
      ep   = 0.0_re_p
      ep_r = 0
      eg   = 0.0_re_p
      ex   = 0.0_re_p
      nec  = 0
      nev  = 0
      nem  = 0
      nep  = 0
      nev_ = 0.0
      nem_ = 0.0
      nep_ = 0.0
      nec_ = 0.0
!$acc end serial
c
c     call the local geometry energy component routines
c
      if (use_bond)    call ebond
      if (use_angle)   call eangle
      if (use_strbnd)  call estrbnd
      if (use_urey)    call eurey
      if (use_angang)  call eangang
      if (use_opbend)  call eopbend
      if (use_opdist)  call eopdist
      if (use_improp)  call eimprop
      if (use_imptor)  call eimptor
      if (use_tors)    call etors
      if (use_pitors)  call epitors
      if (use_strtor)  call estrtor
      if (use_angtor)  call eangtor
      if (use_tortor)  call etortor
c
c     call the van der Waals energy component routines
c
      if (use_vdw) then
         if (vdwtyp .eq. 'LENNARD-JONES')  call elj3gpu
         if (vdwtyp .eq. 'BUFFERED-14-7')  call ehal3gpu
      end if
c
c     call the electrostatic energy component routines
c
      if (use_charge) call echarge3gpu
      if (use_mpole)  call empole3gpu
      if (use_polar)  call epolar3gpu
c
c     call any miscellaneous energy component routines
c
      if (use_geom)   call egeom
      if (use_extra)  call extra
c
c     sum up to give the total potential energy
c
!$acc serial async
      if (ev_r.ne.0) ev = ev + enr2en(ev_r)
      if (eb_r.ne.0) eb = eb + enr2en(eb_r)
      esum = eb + ea + eba + eub + eaa + eopb + eopd + eid + eit
     &          + et + ept + ebt + eat + ett  + ev   + ec  + em
     &          + ep + eg + ex + eg
!$acc end serial
!$acc update host(esum) async
!$acc wait
      energy = esum
c
c     check for an illegal value for the total energy
c
      if (tinker_isnan_m(esum).or.esum.eq.inf_r) then
         write (iout,10) esum
   10    format (/,' ENERGY  --  Illegal Value for the Total',
     &             ' Potential Energy',F16.6)
         call info_energy(0)
         call fatal
      end if
!$acc end data
      end
