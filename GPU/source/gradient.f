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
#include "tinker_macro.h"

      module gradient_inl
      contains
#include "convert.f.inc"
      end module

      subroutine gradient (energy,derivs)
      use ani
      use atoms
      use atmlst
      use atmtyp
      use bitor
      use cell
      use colvars
      use couple
      use deriv
      use domdec
      use energi
      use gradient_inl
      use group
      use inter
      use iounit
      use inform    ,only:abort,minmaxone,deb_Force
      use ktrtor
      use mpi
      use pitors
      use plumed
      use potent
      use strbnd
      use tinheader ,only:ti_p,re_p
      use timestat
      use tors
      use tortor
      use utilgpu   ,only:rec_queue,def_queue,inf_r,zero_mdred_buffers
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
      real(t_p), allocatable, save :: wgrp_save(:,:)

      if (first_in_gradient) then
         first_in_gradient=.false.
      end if

      if(use_ml_embedding .and. use_mlpot) then
        call save_wgrp
        call set_embedding_weights()
      endif

      call timer_enter( timer_fmanage )
      if (calc_e.or.use_virial) then    ! Energy and virial computation
!     Energy scalars used in gradient kernels
!$acc host_data use_device(energy
!$acc&          ,eb,eba,eub,eopb,et,ept,ett,eat,ebt,ea
!$acc&          ,eaa,eopd,eid,eit,ec,ecrec,ev,em,emrec,ep,eprec
!$acc&          ,edsp,edsprec,er,ect,emlpot
!$acc&          ,eg,ex,esum,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz
!$acc&          ,ev_r,ec_r,em_r,ep_r,eb_r
!$acc&          ,ensmd,einter,ePaMD,eDaMD,eW1aMD)
!$acc serial deviceptr(eb,eba,eub,eopb,et,ept,ett,eat,ebt,ea
!$acc&      ,eaa,eopd,eid,eit,ec,ecrec,ev,em,emrec,ep,eprec
!$acc&      ,edsp,edsprec,er,ect,emlpot
!$acc&      ,eg,ex,esum,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz
!$acc&      ,ev_r,ec_r,em_r,ep_r,eb_r
!$acc&      ,ensmd,einter,ePaMD,eDaMD,eW1aMD) async(rec_queue)
      eb       = 0.0_re_p  ! ebond
      eb_r     = 0
      eba      = 0.0_re_p  ! estrbnd
      eub      = 0.0_re_p  ! eurey
      eopb     = 0.0_re_p  ! eopbend
      et       = 0.0_re_p  ! etors
      ept      = 0.0_re_p  ! epitors
      ett      = 0.0_re_p  ! etortor 
      eat      = 0.0_re_p  ! eangtor
      ebt      = 0.0_re_p  ! estrtor
      ea       = 0.0_re_p  ! eangle
      emlpot   = 0.0_re_p  ! emlpot
      ev       = 0.0_re_p  ! ehal1
      ev_r     = 0
      edsp     = 0         ! edisp
      edsprec  = 0         ! edisp rec
      er       = 0         ! erepel
      ect      = 0         ! echgtrn
      em       = 0.0_re_p  ! empole
      emrec    = 0.0_re_p  ! empole (reciprocal)
      em_r     = 0
      ep       = 0.0_re_p  ! epolar
      eprec    = 0.0_re_p  ! epolar (reciprocal)
      ep_r     = 0
      eaa      = 0.0_re_p  ! eangang 
      eopd     = 0.0_re_p  ! eopdist
      eid      = 0.0_re_p  ! eimprop
      eit      = 0.0_re_p  ! eimptor
      ec       = 0.0_re_p  ! echarge
      ec_r     = 0
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
!$acc end host_data
      call zero_mdred_buffers(rec_queue)
      end if  ! Energy and virial computation
 
      !zero out each of the first derivative components
      if (deb_Force) call zero_forces

      if (use_lambdadyn.and.
     &   (use_vdw.or.use_mpole.or.use_charge.or.use_polar)) then ! Zero Halmitonian derivative
!$acc serial async present(delambdae,delambdav)
         delambdae = 0.0
         delambdav = 0.0
!$acc end serial
      end if
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
#ifdef _OPENACC
      if (use_bond.and.fuse_bonded) then
         call eliaison1cu
         goto 26
      end if
#endif
      if (use_bond)     call ebond1gpu
      if (use_urey)     call eurey1gpu
      if (use_opbend)   call eopbend1gpu
      if (use_strbnd)   call estrbnd1gpu
      if (use_angle)    call eangle1gpu
      if (use_tors)     call etors1gpu
      if (use_pitors)   call epitors1gpu
      if (use_tortor)   call etortor1gpu
      if (use_improp)   call eimprop1gpu
      if (use_imptor)   call eimptor1gpu
      if (use_angtor)   call eangtor1

 26   continue
      if (use_angang)   call eangang1     !no acc
      if (use_opdist)   call eopdist1     !no acc
      if (use_strtor)   call estrtor1gpu
      ! miscellaneous energy
      if (use_geom)     call egeom1gpu
      if (use_extra)    call extra1
      if (use_bond) call timer_exit ( timer_bonded )
      
      if (use_mlpot) then
         call timer_enter(timer_b1)
         call ml_potential(.true.)
         call timer_exit(timer_b1)
      end if

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
      if (use_repuls)   call erepel1
      if (use_disp)     call edisp1
c
c     alter partial charges and multipoles for charge flux
c
      if (use_chgflx)  call alterchg1
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
      if (use_chgtrn)   call echgtrn1gpu

      if (use_smd_velconst.or.use_smd_forconst) call esmd1

      if (use_group) call switch_group(.TRUE.)

      if (use_vdw.or.use_mpole.or.use_charge)
     &   call timer_exit(timer_nonbonded)
 


!     ----------------------------------------
!     ******   Add Partial Forces   ******
!     ----------------------------------------

      call timer_enter( timer_fmanage )
      if (ftot_l) then
         call add_forces( de_tot )
      else
         call add_forces1( derivs )
      end if

      ! Apply colvars
      if (use_colvars) call colvars_run(derivs)

      if (calc_e.or.use_virial) then     ! Energy and virial computation
!$acc host_data use_device(energy
!$acc&         ,eb,eba,eub,eopb,et,ept,ett,eat,ebt,ea
!$acc&         ,eaa,eopd,eid,eit,ec,ecrec,ev,em,emrec,ep,eprec
!$acc&         ,edsp,edsprec,er,ect,emlpot
!$acc&         ,eg,ex,esum,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz
!$acc&         ,ev_r,ec_r,em_r,ep_r,eb_r
!$acc&         ,ensmd,einter,ePaMD,eDaMD,eW1aMD)
!$acc serial deviceptr(energy,eb,eba,eub,eopb,et,ept,ett,eat,ebt,ea
!$acc&      ,eaa,eopd,eid,eit,ec,ecrec,ev,em,emrec,ep,eprec
!$acc&      ,edsp,edsprec,er,ect,emlpot
!$acc&      ,eg,ex,esum,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz
!$acc&      ,ev_r,ec_r,em_r,ep_r,eb_r
!$acc&      ,ensmd,einter,ePaMD,eDaMD,eW1aMD) async(rec_queue)
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
      if (ev_r.ne.0) ev = ev + enr2en(ev_r)
      if (eb_r.ne.0) eb = eb + enr2en(eb_r)
      esum = eit + eopd + eopb + eaa + eub + eba + ea + eb + em  + ep
     &      + ec + ev   + et   + ept + ebt + eat +ett + eg + ex + eid
     &      + ensmd + emlpot + enr2en(edsp + er + ect)
      energy = esum
!$acc end serial
!$acc end host_data
      end if        ! Energy and virial computation
      call timer_exit ( timer_fmanage,quiet_timers )

      ! Apply plumed bias
      if (lplumed) call eplumed(energy,derivs)

      if(use_ml_embedding .and. use_mlpot) then
        call load_wgrp()
      endif

c
c     check for an illegal value for the total energy
c
      if (tinkerdebug.gt.0.and.calc_e) then
!$acc update host(esum) async(rec_queue)
!$acc wait
         if (tinker_isnan_m(esum).or.esum.eq.inf_r) then
            write  (0,10) esum
   10       format (/,' GRADIENT  --  Illegal Value for the Total',
     &                ' Potential Energy',F16.6)
            abort=.true.
            call info_energy(0)
         end if
      end if
      end
