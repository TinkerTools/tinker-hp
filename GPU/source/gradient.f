c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine gradient  --  find energy & gradient components  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "gradient" calls subroutines to calculate the potential energy
c     and first derivatives with respect to Cartesian coordinates
c
c
#include "tinker_precision.h"
      subroutine gradient (energy,derivs)
      use atoms
      use atmlst
      use atmtyp
      use bitor
      use couple
      use deriv
      use domdec
      use energi
      use inter
      use iounit
      use inform    ,only:abort
      use ktrtor
      use mpi
      use potent
      use pitors
      use strbnd
      use tinheader ,only:ti_p,re_p
      use timestat
      use tors
      use tortor
      use utilgpu   ,only:def_queue,inf
      use sizes
      use vdwpot
      use virial
      implicit none
      integer i,j
      real(r_p) energy
      real(r_p) derivs(3,nbloc)
      real(r_p) time0,time1,time2
      logical tinker_isnan_m
      logical,save::first_in_gradient=.true.

      if (first_in_gradient) then
         first_in_gradient=.false.
      end if

!     Energy scalars used in gradient kernels
!$acc data present(derivs,energy)
!$acc&     present(eb,eba,eub,eopb,et,ept,ett,ebt,ea
!$acc&           ,eaa,eopd,eid,eit,ec,ecrec,ev,em,emrec,ep,eprec
!$acc&           ,eg,ex,esum,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz
!$acc&           ,ensmd,einter,ePaMD,eDaMD,eW1aMD,vir,desum)

      call timer_enter( timer_fmanage )
!$acc serial async
      eb       = 0.0_re_p  ! ebond
      eba      = 0.0_re_p  ! estrbnd
      eub      = 0.0_re_p  ! eurey
      eopb     = 0.0_re_p  ! eopbend
      et       = 0.0_re_p  ! etors
      ept      = 0.0_re_p  ! epitors
      ett      = 0.0_re_p  ! etortor 
      ebt      = 0.0_re_p  ! estrtor
      ea       = 0.0_re_p  ! eangle
      ev       = 0.0_re_p  ! ehal1
      em       = 0.0_re_p  ! empole
      emrec    = 0.0_re_p
      ep       = 0.0_re_p  ! epolar
      eprec    = 0.0_re_p
      eaa      = 0.0_re_p  ! eangang 
      eopd     = 0.0_re_p  ! eopdist
      eid      = 0.0_re_p  ! eimprop
      eit      = 0.0_re_p  ! eimptor
      ec       = 0.0_re_p  ! echarge
      ecrec    = 0.0_re_p
      eg       = 0.0_re_p  ! egeom
      ensmd    = 0.0_re_p  ! esmd
      ex       = 0.0_re_p  ! extra
      einter   = 0.0_re_p  ! intermolecular energy
      ePaMD    = 0.0_re_p
      eDaMD    = 0.0_re_p
      eW1aMD   = 0.0_re_p
c
c     zero out the virial
c
      g_vxx    = 0.0_re_p
      g_vxy    = 0.0_re_p
      g_vxz    = 0.0_re_p
      g_vyy    = 0.0_re_p
      g_vyz    = 0.0_re_p
      g_vzz    = 0.0_re_p
      vir(1,1) = 0.0_re_p
      vir(2,1) = 0.0_re_p
      vir(3,1) = 0.0_re_p
      vir(1,2) = 0.0_re_p
      vir(2,2) = 0.0_re_p
      vir(3,2) = 0.0_re_p
      vir(1,3) = 0.0_re_p
      vir(2,3) = 0.0_re_p
      vir(3,3) = 0.0_re_p
!$acc end serial
c
c     zero out each of the first derivative components
c
      call resetForces_p
      call timer_exit ( timer_fmanage,quiet_timers )
!
!     ----------------------------------------
!     ******      Bonded Section        ******
!     ----------------------------------------
!
c
c     call the local geometry energy and gradient routines
c
      if (use_bond) call timer_enter( timer_bonded )
      if (use_bond)     call ebond1gpu
      if (use_strbnd)   call estrbnd1gpu
      if (use_urey)     call eurey1gpu
      if (use_angang)   call eangang1     !no acc
      if (use_opbend)   call eopbend1gpu
      if (use_opdist)   call eopdist1     !no acc
      if (use_improp)   call eimprop1gpu
      if (use_imptor)   call eimptor1     !no acc
      if (use_tors)     call etors1gpu
      if (use_pitors)   call epitors1gpu
      if (use_strtor)   call estrtor1gpu
      if (use_tortor)   call etortor1gpu
      if (use_angle)    call eangle1gpu
      ! miscellaneous energy
      if (use_geom)     call egeom1gpu
      if (use_extra)    call extra1
      if (use_bond) call timer_exit ( timer_bonded )

!
!     ----------------------------------------
!     ******     Non Bonded Section     ******
!     ----------------------------------------
!
      if (use_vdw.or.use_mpole.or.use_charge)
     &   call timer_enter(timer_nonbonded)
c
c     call the van der Waals energy and gradient routines
c
      if (use_vdw) then
         call timer_enter( timer_vdw )
         if (vdwtyp .eq. 'LENNARD-JONES')  call elj1gpu
         if (vdwtyp .eq. 'BUFFERED-14-7')  call ehal1gpu
         call timer_exit ( timer_vdw )
      end if
c
c     call the electrostatic energy and gradient routines
c
      if (use_mpole.or.use_charge) call timer_enter( timer_elec )
      if (use_charge)   call echarge1gpu
      if (use_mpole)    call empole1gpu
      if (use_mpole.or.use_charge) call timer_exit ( timer_elec )

      if (use_polar)  then
         call timer_enter( timer_polar )
                        call epolar1gpu
         call timer_exit ( timer_polar )
      endif

      if (use_smd_velconst.or.use_smd_forconst) call esmd1

      if (use_vdw.or.use_mpole.or.use_charge)
     &   call timer_exit(timer_nonbonded)
 


!     ----------------------------------------
!     ******   Compute Partial Forces   ******
!     ----------------------------------------

      call timer_enter( timer_fmanage )
      call addForces_p

!$acc parallel loop collapse(2) async
      do i = 1, nloc
         do j = 1, 3
            derivs(j,i) = derivs(j,i) + desum(j,i)
         end do
      end do

!$acc serial async
c
c     sum up the virial component
c
      vir(1,1) = vir(1,1) + g_vxx
      vir(2,1) = vir(2,1) + g_vxy
      vir(3,1) = vir(3,1) + g_vxz
      vir(1,2) = vir(1,2) + g_vxy
      vir(2,2) = vir(2,2) + g_vyy
      vir(3,2) = vir(3,2) + g_vyz
      vir(1,3) = vir(1,3) + g_vxz
      vir(2,3) = vir(2,3) + g_vyz
      vir(3,3) = vir(3,3) + g_vzz
c
c     sum up to get the total energy and first derivatives
c
      esum = eit + eopd + eopb + eaa + eub + eba + ea + eb + em  + ep
     &      + ec + ev   + et   + ept + ebt + ett + eg + ex + eid + ensmd
      energy = esum
!$acc end serial
      call timer_exit ( timer_fmanage,quiet_timers )

!$acc end data
c
c     check for an illegal value for the total energy
c
!$acc update host(esum) async
      if (tinkerdebug.gt.0) then
!$acc wait
         if (tinker_isnan_m(esum).or.esum.eq.3*huge(0.0_re_p)) then
            write  (0,10) esum
   10       format (/,' GRADIENT  --  Illegal Value for the Total',
     &                ' Potential Energy',F16.6)
            abort=.true.
            call info_energy(0)
         end if
      end if
      end
