c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ####################################################################
c     ##                                                                ##
c     ##  module subderiv                                               ##
c     ##  --  Routines associated to derivatives component managment    ##
c     ##                                                                ##
c     ####################################################################
c
#include "tinker_macro.h"

      submodule(deriv) subderiv
      use atoms    ,only:x,y,z,n
      use argue    ,only:arg
      use colvars  ,only:use_colvars,ncvatoms,decv,decv_tot
      use domdec
      use dcdmod   ,only:dcdio
      use inform   ,only:deb_Path,deb_Force,abort,n_fwriten
     &             ,dint1,dint2,mmo=>minmaxone
      use,intrinsic:: iso_c_binding ,only: c_ptr,c_f_pointer,c_loc
      use mpi
      use mdstuf   ,only:integrate
      use mdstuf1  ,only:epot
      use moldyn   ,only:stepfast,stepint,nalt,step_c
      use neigh
      use output   ,only:archive,new_restart,f_mdsave
      use potent
      use sizes
      use timestat
      use tinheader,only:ti_p,re_p
      use tinMemory,only:prmem_request,prmem_requestm,mipk
     &             ,prmem_int8_req1,prmemGetAllocSize,extra_alloc
      use utilgpu  ,only:mem_move,mem_set,rec_stream,rec_queue
     &             ,Tinker_shellEnv
      use utilcomm
      implicit none
      logical abortall
      integer inte(2)
      integer:: n_adjust=0
      mdyn_rtyp,parameter::zero_md=0
      real(r_p),parameter::zero_rp=0
      logical n_debf
      logical:: mem_alloc_deriv_fcall=.true.
      real(r_p) f_ulim

      contains

#include "convert.f.inc"

      module subroutine ConfigPotential
      implicit none
      integer watchPot
      watchPot = 0

      PotentialAmoeba = use_bond.and.use_angle.and.use_strbnd
     &      .and.use_urey.and..not.use_angang.and.use_opbend
     &      .and..not.use_opdist.and..not.use_improp.and..not.use_imptor
     &      .and.use_tors.and.use_pitors.and..not.use_strtor
     &      .and.use_tortor.and..not.use_angtor.and.use_vdw
     &      .and..not.use_chgtrn.and..not.use_disp.and..not.use_repuls
     &      .and..not.use_charge.and.use_mpole.and.use_polar
     &      .and..not.use_extra
      PotentialAmoeba18 = use_bond.and.use_angle.and.use_strbnd
     &      .and.use_urey.and..not.use_angang.and.use_opbend
     &      .and..not.use_opdist.and..not.use_improp.and..not.use_imptor
     &      .and.use_tors.and.use_pitors.and.use_strtor
     &      .and..not.use_angtor.and.use_tortor.and.use_vdw
     &      .and..not.use_chgtrn.and..not.use_disp.and..not.use_repuls
     &      .and..not.use_charge.and.use_mpole.and.use_polar
     &      .and..not.use_extra
      PotentialAmoeba181= use_bond.and.use_angle.and.use_strbnd
     &      .and.use_urey.and..not.use_angang.and.use_opbend
     &      .and..not.use_opdist.and..not.use_improp.and..not.use_imptor
     &      .and.use_tors.and.use_pitors.and.use_strtor
     &      .and.use_angtor.and..not.use_tortor.and.use_vdw
     &      .and..not.use_chgtrn.and..not.use_disp.and..not.use_repuls
     &      .and..not.use_charge.and.use_mpole.and.use_polar
     &      .and..not.use_extra
      PotentialAmoeba182= use_bond.and.use_angle.and.use_strbnd
     &      .and.use_urey.and..not.use_angang.and.use_opbend
     &      .and..not.use_opdist.and..not.use_improp.and..not.use_imptor
     &      .and.use_tors.and.use_pitors.and.use_strtor
     &      .and.use_angtor.and.use_tortor.and.use_vdw
     &      .and..not.use_chgtrn.and..not.use_disp.and..not.use_repuls
     &      .and..not.use_charge.and.use_mpole.and.use_polar
     &      .and..not.use_extra
      PotentialWaterAmoeba = use_bond.and.use_angle.and..not.use_strbnd
     &      .and.use_urey.and..not.use_angang.and..not.use_opbend
     &      .and..not.use_opdist.and..not.use_improp.and..not.use_imptor
     &      .and..not.use_tors.and..not.use_pitors.and..not.use_strtor
     &      .and..not.use_tortor.and..not.use_angtor.and.use_vdw
     &      .and..not.use_chgtrn.and..not.use_disp.and..not.use_repuls
     &      .and.use_mpole.and.use_polar.and..not.use_charge
     &      .and..not.use_extra
      PotentialWaterCharmm = use_bond.and.use_angle.and..not.use_strbnd
     &      .and.use_urey.and..not.use_angang.and..not.use_opbend
     &      .and..not.use_opdist.and..not.use_improp.and..not.use_imptor
     &      .and..not.use_tors.and..not.use_pitors.and..not.use_strtor
     &      .and..not.use_tortor.and..not.use_angtor.and.use_vdw
     &      .and..not.use_chgtrn.and..not.use_disp.and..not.use_repuls
     &      .and..not.use_mpole.and..not.use_polar.and.use_charge
     &      .and..not.use_extra
      PotentialCharmm = use_bond.and.use_angle.and..not.use_strbnd
     &      .and.use_urey.and..not.use_angang.and..not.use_opbend
     &      .and..not.use_opdist.and.use_improp.and..not.use_imptor
     &      .and.use_tors.and..not.use_pitors.and..not.use_strtor
     &      .and..not.use_tortor.and..not.use_angtor.and.use_vdw
     &      .and..not.use_chgtrn.and..not.use_disp.and..not.use_repuls
     &      .and..not.use_mpole.and..not.use_polar.and.use_charge
     &      .and..not.use_extra

      if (PotentialAmoeba) then
         if (deb_Path) print*,'Using Amoeba Potential'
         PotentialAll=.false.
         watchPot = watchPot + 1
      end if
      if (PotentialAmoeba18) then
         if (deb_Path) print*,'Using Amoeba18 Potential'
         PotentialAll=.false.
         watchPot = watchPot + 1
      end if
      if (PotentialAmoeba181) then
         if (deb_Path) print*,'Using Amoeba18 (1) Potential'
         PotentialAll=.false.
         watchPot = watchPot + 1
      end if
      if (PotentialAmoeba182) then
         if (deb_Path) print*,'Using Amoeba18 (2) Potential'
         PotentialAll=.false.
         watchPot = watchPot + 1
      end if
      if (PotentialWaterAmoeba) then
         if (deb_Path) print*,'Using Water Amoeba Potential'
         PotentialAll=.false.
         watchPot = watchPot + 1
      end if
      if (PotentialWaterCharmm) then
         PotentialAll=.false.
         watchPot = watchPot + 1
      end if
      if (PotentialCharmm) then
         PotentialAll=.false.
         if (deb_Path) print*,'Using Charmm Potential'
         watchPot = watchPot + 1
      end if
      if (PotentialAll) then
         if (deb_Path) print*,'Using Default Potential'
         watchPot = watchPot + 1
      end if

      if (watchPot.gt.1) then
11       format ( " WARNING !!! ConfigPotential Routine ",/,
     &   " Found two of more differents configurations on running",
     &   " potentials ",/,
     &   " Switched to default potential by reserving space for all",
     &   " forces ")
         print 11

         PotentialAll         = .true.
         PotentialAmoeba      = .false.
         PotentialAmoeba18    = .false.
         PotentialCharmm      = .false.
         PotentialWaterAmoeba = .false.
         PotentialWaterCharmm = .false.
      end if

      block
      logical OK,respa1_l,smd_l,gamd_l

      fdebs_l    = .false.
      tdes_l     = .false.
      ftot_l     = .false.

      respa1_l   = (index(integrate,'RESPA1').gt.0)
      smd_l      = use_smd_velconst.or.use_smd_forconst
      gamd_l     = use_gamd.or.use_amd_ene.or.use_amd_dih.or.
     &             use_amd_wat1

      ftot_l = (integrate.eq.'RESPA'.or.respa1_l
     &      .or.integrate.eq.'BAOAB'.or.integrate.eq.'BAOABRESPA'
     &      .or.integrate.eq.'VERLET'.or.integrate.eq.'BBK'
     &      .or.integrate.eq.'BEEMAN')

      fdebs_l    = merge(.true.,.false.,ftot_l.and..not.gamd_l
     &                   .and..not.deb_Force)

      OK         = .not.fdebs_l
      n_debf     = .not.deb_Force

      dr_nb0     = 0
      dr_nbnbr   = 0
      dr_nb1     = 1  !( derivs )

      dr_stride  = 0
      dr_stride3 = 0

      dr_nb0 = dr_nb0 + merge(1,0,use_bond.and.OK)
      dr_nb0 = dr_nb0 + merge(1,0,use_angle.and.OK)
      dr_nb0 = dr_nb0 + merge(1,0,use_urey.and.OK)
      dr_nb0 = dr_nb0 + merge(1,0,use_strbnd.and.OK)
      dr_nb0 = dr_nb0 + merge(1,0,use_opbend.and.OK)
      dr_nb0 = dr_nb0 + merge(1,0,use_opdist.and.OK)
      dr_nb0 = dr_nb0 + merge(1,0,use_angang.and.OK)
      dr_nb0 = dr_nb0 + merge(1,0,use_improp.and.OK)
      dr_nb0 = dr_nb0 + merge(1,0,use_imptor.and.OK)
      dr_nb0 = dr_nb0 + merge(1,0,use_tors.and.OK)
      dr_nb0 = dr_nb0 + merge(1,0,use_pitors.and.OK)
      dr_nb0 = dr_nb0 + merge(1,0,use_strtor.and.OK)
      dr_nb0 = dr_nb0 + merge(1,0,use_angtor.and.OK)
      dr_nb0 = dr_nb0 + merge(1,0,use_tortor.and.OK)
      dr_nb0 = dr_nb0 + merge(1,0,use_geom.and.OK)
      dr_nb0 = dr_nb0 + merge(1,0,(use_mlpot.or.use_ml_embedding)
     &                        .and.OK)
      dr_nb0 = dr_nb0 + merge(1,0,fdebs_l)
      dr_nbb = dr_nb0

      dr_nb0 = dr_nb0 + 1  !( de_ws0 )
      dr_nb0 = dr_nb0 + merge(2,0,gamd_l)

      dr_nbnbr = dr_nbnbr + merge(1,0,use_charge)  !( decrec )
      dr_nbnbr = dr_nbnbr + merge(1,0,use_mpole )  !( demrec )
      dr_nbnbr = dr_nbnbr + merge(1,0,use_polar )  !( deprec )
      dr_nbnbr = dr_nbnbr + merge(1,0,use_disp  )  !( dedsprec )

      dr_nb1 = dr_nb1 + merge(1,0,use_vdw   )
      dr_nb1 = dr_nb1 + merge(1,0,use_disp  )
      dr_nb1 = dr_nb1 + merge(1,0,use_repuls)
      dr_nb1 = dr_nb1 + merge(1,0,use_charge)
      dr_nb1 = dr_nb1 + merge(1,0,use_chgtrn)
      dr_nb1 = dr_nb1 + merge(1,0,use_mpole )
      dr_nb1 = dr_nb1 + merge(1,0,use_polar )
      dr_nb1 = dr_nb1 + merge(1,0,smd_l     )
      dr_nbnb= dr_nb1 - 1  !( Remove de_tot )
      dr_nb1 = dr_nb1 + merge(1,0,respa1_l  ) !( desave )  
      dr_nb1 = dr_nb1 + 2  !( de_ws1,de_ws2 )

      if (tinkerdebug.gt.0.and.rank.eq.0) then
 12      format(
     &        /,'  --- Forces buffers ---',
     &        /,' buff0 ',I4,' buff1 ',I4,
     &        /,' bonded',I4,' non bonded ',I4,' nrec',I4
     &        /,' ftot  ',L4,' fdebs      ',L4,' tdes',L4)
         write(*,12) dr_nb0,dr_nb1,dr_nbb,dr_nbnb,dr_nbnbr
     &              ,ftot_l,fdebs_l,tdes_l
      end if

      !set f_ulim
      call Tinker_shellEnv("DEBUG_FMAX",f_ulim,-1.0_re_p)
      if (tinkerdebug.gt.0.and.rank.eq.0.and.f_ulim.gt.0.0) then
  13     format(" Set a limit of ",F8.2," on forces L_infty norm")
         write(*,13) f_ulim
      end if
      end block

      end subroutine

      module subroutine mem_alloc_deriv(opt)
      implicit none
      integer,optional:: opt
      integer(mipk) i,j
      integer opt_
      integer(mipk) siz0,of0,siz1,of1
      integer(mipk) ofr,sizr
      logical OK,save_ex_alloc

      if (present(opt)) then
         opt_ = opt
      else
         opt_ = cBond + cNBond
      end if

      OK = .not.fdebs_l

      if (nbloc<=dr_stride) goto 20
      if (deb_Path) write(*,*) "mem_alloc_deriv",opt_,nbloc

      dr_stride   = prmemGetAllocSize(nbloc)
      dr_stride3  = 3*dr_stride

      siz0       = dr_nb0*int(dr_stride3,mipk)
      siz1       = dr_nb1*int(dr_stride3,mipk)
      of0        = 0
      of1        = 0

      save_ex_alloc = extra_alloc
      extra_alloc   = .false.
      call prmem_requestm(de_buff0,siz0,async=.false.)
      call prmem_requestm(de_buff1,siz1,async=.false.)
      extra_alloc   = save_ex_alloc

      call mem_set( de_buff0,zero_rp,siz0,rec_stream )
      call mem_set( de_buff1,zero_md,siz1,rec_stream )

      dr_obb = of0
      if (use_bond) then
         deb(1:3,1:dr_stride) => de_buff0(of0+1:of0+dr_stride3)
         of0 = of0 + merge(dr_stride3,0,OK)
      end if
      if (use_angle) then
         dea(1:3,1:dr_stride) => de_buff0(of0+1:of0+dr_stride3)
         of0 = of0 + merge(dr_stride3,0,OK)
      end if
      if (use_urey) then
         deub(1:3,1:dr_stride) => de_buff0(of0+1:of0+dr_stride3)
         of0 = of0 + merge(dr_stride3,0,OK)
      end if
      if (use_strbnd) then
         deba(1:3,1:dr_stride) => de_buff0(of0+1:of0+dr_stride3)
         of0 = of0 + merge(dr_stride3,0,OK)
      end if
      if (use_opbend) then
         deopb(1:3,1:dr_stride) => de_buff0(of0+1:of0+dr_stride3)
         of0 = of0 + merge(dr_stride3,0,OK)
      end if
      if (use_opdist) then
         deopd(1:3,1:dr_stride) => de_buff0(of0+1:of0+dr_stride3)
         of0 = of0 + merge(dr_stride3,0,OK)
      end if
      if (use_angang) then
         deaa(1:3,1:dr_stride) => de_buff0(of0+1:of0+dr_stride3)
         of0 = of0 + merge(dr_stride3,0,OK)
      end if
      if (use_improp) then
         deid(1:3,1:dr_stride) => de_buff0(of0+1:of0+dr_stride3)
         of0 = of0 + merge(dr_stride3,0,OK)
      end if
      if (use_imptor) then
         deit(1:3,1:dr_stride) => de_buff0(of0+1:of0+dr_stride3)
         of0 = of0 + merge(dr_stride3,0,OK)
      end if
      if (use_tors) then
         det(1:3,1:dr_stride) => de_buff0(of0+1:of0+dr_stride3)
         of0 = of0 + merge(dr_stride3,0,OK)
      end if
      if (use_pitors) then
         dept(1:3,1:dr_stride) => de_buff0(of0+1:of0+dr_stride3)
         of0 = of0 + merge(dr_stride3,0,OK)
      end if
      if (use_tortor) then
         dett(1:3,1:dr_stride) => de_buff0(of0+1:of0+dr_stride3)
         of0 = of0 + merge(dr_stride3,0,OK)
      end if
      if (use_strtor) then
         debt(1:3,1:dr_stride) => de_buff0(of0+1:of0+dr_stride3)
         of0 = of0 + merge(dr_stride3,0,OK)
      end if
      if (use_angtor) then
         deat(1:3,1:dr_stride) => de_buff0(of0+1:of0+dr_stride3)
         of0 = of0 + merge(dr_stride3,0,OK)
      end if
      if (use_mlpot.or.use_ml_embedding) then
         dmlpot(1:3,1:dr_stride) => de_buff0(of0+1:of0+dr_stride3)
         of0 = of0 + merge(dr_stride3,0,OK)
      end if
      if (use_geom) then
         deg(1:3,1:dr_stride) => de_buff0(of0+1:of0+dr_stride3)
         of0 = of0 + merge(dr_stride3,0,OK)
      end if
      if (fdebs_l) of0 = of0 + dr_stride3

      de_ws0 (1:3,1:dr_stride) => de_buff0(of0+1:of0+dr_stride3)
      of0 = of0 + dr_stride3

      if (use_gamd.or.use_amd_ene.or.use_amd_dih.or.use_amd_wat1) then
         deamdD (1:3,1:dr_stride) => de_buff0(of0+1:of0+dr_stride3)
         of0 = of0 + dr_stride3
         deW1aMD(1:3,1:dr_stride) => de_buff0(of0+1:of0+dr_stride3)
         of0 = of0 + dr_stride3
      else
          deamdD(1:3,1:1) => de_buff0(of0-dr_stride3+1:)
         deW1aMD(1:3,1:1) => de_buff0(of0-dr_stride3+1:)
      end if

      ! ------------------------------------------------------------------------
      !         de_buff1 (mdyn_rtyp)

      de_tot (1:3,1:dr_stride) => de_buff1(of1+1:of1+dr_stride3)
      derivx(1:dr_stride3) =>de_buff1(of1            +1:of1+dr_stride3)
      derivy(1:dr_stride) => de_buff1(of1+1*dr_stride+1:of1+2*dr_stride)
      derivz(1:dr_stride) => de_buff1(of1+2*dr_stride+1:of1+3*dr_stride)
      of1 = of1 + dr_stride3

      dr_obnb = of1
      if (use_vdw) then
         dev(1:3,1:dr_stride) => de_buff1(of1+1:of1+dr_stride3)
         devx(1:dr_stride)=> de_buff1(of1+0*dr_stride+1:of1+1*dr_stride)
         devy(1:dr_stride)=> de_buff1(of1+1*dr_stride+1:of1+2*dr_stride)
         devz(1:dr_stride)=> de_buff1(of1+2*dr_stride+1:of1+3*dr_stride)
         of1 = of1 + dr_stride3
      end if
      if (use_disp) then
         dedsp(1:3,1:dr_stride) => de_buff1(of1+1:of1+dr_stride3)
         of1 = of1 + dr_stride3
      end if
      if (use_repuls) then
         der(1:3,1:dr_stride) => de_buff1(of1+1:of1+dr_stride3)
         of1 = of1 + dr_stride3
      end if
      if (use_charge) then
         dec(1:3,1:dr_stride) => de_buff1(of1+1:of1+dr_stride3)
         decx(1:dr_stride)=> de_buff1(of1+0*dr_stride+1:of1+1*dr_stride)
         decy(1:dr_stride)=> de_buff1(of1+1*dr_stride+1:of1+2*dr_stride)
         decz(1:dr_stride)=> de_buff1(of1+2*dr_stride+1:of1+3*dr_stride)
         of1 = of1 + dr_stride3
      end if
      if (use_chgtrn) then
         dect(1:3,1:dr_stride) => de_buff1(of1+1:of1+dr_stride3)
         of1 = of1 + dr_stride3
      end if
      if (use_mpole ) then
         dem(1:3,1:dr_stride) => de_buff1(of1+1:of1+dr_stride3)
         demx(1:dr_stride)=> de_buff1(of1+0*dr_stride+1:of1+1*dr_stride)
         demy(1:dr_stride)=> de_buff1(of1+1*dr_stride+1:of1+2*dr_stride)
         demz(1:dr_stride)=> de_buff1(of1+2*dr_stride+1:of1+3*dr_stride)
         of1 = of1 + dr_stride3
      end if
      if (use_polar ) then
         dep(1:3,1:dr_stride) => de_buff1(of1+1:of1+dr_stride3)
         depx(1:dr_stride)=> de_buff1(of1+0*dr_stride+1:of1+1*dr_stride)
         depy(1:dr_stride)=> de_buff1(of1+1*dr_stride+1:of1+2*dr_stride)
         depz(1:dr_stride)=> de_buff1(of1+2*dr_stride+1:of1+3*dr_stride)
         of1 = of1 + dr_stride3
      end if
      if (use_smd_velconst.or.use_smd_forconst) then
         desmd  (1:3,1:dr_stride) => de_buff1(of1+1:of1+dr_stride3)
         of1 = of1 + dr_stride3
      end if
      if (integrate.eq.'RESPA1'.or.integrate.eq.'BAOABRESPA1') then
         desave (1:3,1:dr_stride) => de_buff1(of1+1:of1+dr_stride3)
         of1 = of1 + dr_stride3
      end if

      dr_obws = of1
      de_ws1 (1:3,1:dr_stride) => de_buff1(of1+1:of1+dr_stride3)
      d_x  (1:dr_stride) => de_buff1(of1+0*dr_stride+1:of1+1*dr_stride)
      d_y  (1:dr_stride) => de_buff1(of1+1*dr_stride+1:of1+2*dr_stride)
      d_z  (1:dr_stride) => de_buff1(of1+2*dr_stride+1:of1+3*dr_stride)
      of1 = of1 + dr_stride3
      de_ws2 (1:3,1:dr_stride) => de_buff1(of1+1:of1+dr_stride3)
      de1x (1:dr_stride) => de_buff1(of1+0*dr_stride+1:of1+1*dr_stride)
      de1y (1:dr_stride) => de_buff1(of1+1*dr_stride+1:of1+2*dr_stride)
      de1z (1:dr_stride) => de_buff1(of1+2*dr_stride+1:of1+3*dr_stride)
      of1 = of1 + dr_stride3

      if (of0.ne.siz0.or.of1.ne.siz1) then
 11      format('ERROR! mem_alloc_deriv',
     &        /,' end offset unequal to size',4I10)
         print 11, of0,siz0,of1,siz1
         call fatal
      end if

 20   continue
      if (nlocrec2<=dr_strider.or..not.btest(opt_,idNBond)) goto 30
      if (deb_Path) write(*,*)
     &   "mem_alloc_deriv_rec",opt_,nlocrec2

      dr_strider  = prmemGetAllocSize(nlocrec2)
      dr_strider3 = 3*dr_strider

      sizr        = dr_nbnbr*int(dr_strider3,mipk)
      ofr         = 0
      dr_obnbr    = ofr

      save_ex_alloc = extra_alloc
      extra_alloc   = .false.
      call prmem_requestm(de_buffr,sizr,async=.false.)
      extra_alloc   = save_ex_alloc
      call mem_set( de_buffr,zero_rp,sizr,rec_stream )

      if (use_charge) then
         decrec(1:3,1:dr_strider) => de_buffr(ofr+1:ofr+dr_strider3)
         ofr = ofr + dr_strider3
      end if
      if (use_mpole) then
         demrec(1:3,1:dr_strider) => de_buffr(ofr+1:ofr+dr_strider3)
         ofr = ofr + dr_strider3
      end if
      if (use_polar) then
         deprec(1:3,1:dr_strider) => de_buffr(ofr+1:ofr+dr_strider3)
         ofr = ofr + dr_strider3
      end if
      if (use_disp ) then
       dedsprec(1:3,1:dr_strider) => de_buffr(ofr+1:ofr+dr_strider3)
         ofr = ofr + dr_strider3
      end if

      if (ofr.ne.sizr) then
 21      format('ERROR! mem_alloc_deriv_rec',
     &        /,' end offset unequal to size',2I10)
         print 21, ofr,sizr
         call fatal
      end if
      if(use_mlpot) then
!$acc update host(dmlpot)
      end if

 30   continue
      if (mem_alloc_deriv_fcall) call zero_forces_host
      mem_alloc_deriv_fcall=.false.
      end subroutine

      module subroutine mem_free_deriv
      implicit none

      if (associated(deb))    nullify (deb)
      if (associated(dea))    nullify (dea)
      if (associated(deba))   nullify (deba)
      if (associated(deub))   nullify (deub)
      if (associated(deaa))   nullify (deaa)
      if (associated(deopb))  nullify (deopb)
      if (associated(deopd))  nullify (deopd)
      if (associated(deid))   nullify (deid)
      if (associated(det))    nullify (det)
      if (associated(dept))   nullify (dept)
      if (associated(deit))   nullify (deit)
      if (associated(deat))   nullify (deat)
      if (associated(debt))   nullify (debt)
      if (associated(dett))   nullify (dett)
      if (associated(deg))    nullify (deg)
      if (associated(dmlpot)) nullify (dmlpot)
      if (associated(dev))    nullify (dev)
      if (associated(dec))    nullify (dec)
      if (associated(dem))    nullify (dem)
      if (associated(dep))    nullify (dep)
      if (associated(dect))   nullify (dect)
      if (associated(der))    nullify (der)
      if (associated(dedsp))  nullify (dedsp)
      if (associated(dedsprec)) nullify (dedsprec)
      if (associated(decrec)) nullify (decrec)
      if (associated(demrec)) nullify (demrec)
      if (associated(deprec)) nullify (deprec)
      if (associated(desave)) nullify (desave)
      if (associated(desmd))  nullify (desmd)

      if (allocated(de_buff0)) deallocate(de_buff0)
      if (allocated(de_buff1)) deallocate(de_buff1)
      if (allocated(de_buffr)) deallocate(de_buffr)
      end subroutine

      subroutine updateForcesHost
      implicit none

      print*, "Empty Routine subDeriv_updateForcesHost"
      end subroutine

      module subroutine zero_forces
      integer(mipk) siz0,siz1,sizr

      siz0 = dr_nb0  *int(dr_stride3 ,mipk)
      siz1 = dr_nb1  *int(dr_stride3 ,mipk)
      sizr = dr_nbnbr*int(dr_strider3,mipk)
      call mem_set( de_buff0,zero_rp,siz0,rec_stream )
      call mem_set( de_buff1,zero_md,siz1,rec_stream )
      call mem_set( de_buffr,zero_rp,sizr,rec_stream )

      end subroutine

      module subroutine zero_forces_host
      integer(mipk) siz0,siz1,sizr,i

      siz0 = dr_nb0  *int(dr_stride3 ,mipk)
      siz1 = dr_nb1  *int(dr_stride3 ,mipk)
      sizr = dr_nbnbr*int(dr_strider3,mipk)
      do i = 1,siz0; de_buff0(i) = 0; end do
      do i = 1,siz1; de_buff1(i) = 0; end do
      do i = 1,sizr; de_buffr(i) = 0; end do

      end subroutine

      module subroutine check_nzero
      integer(mipk) siz0,siz1,sizr,offset,i
      integer n0,n1,nr
      offset = dr_stride3
      siz0   = dr_nb0  *int(dr_stride3 ,mipk)
      siz1   = dr_nb1  *int(dr_stride3 ,mipk)
      sizr = dr_nbnbr*int(dr_strider3,mipk)
      n0=0;n1=0;nr=0;
!$acc parallel loop async default(present)
      do i = 1,siz0
         if (de_buff0(i).ne.0) n0=n0+1
      end do
!$acc parallel loop async default(present)
      do i = 1,siz1
         if (de_buff1(i).ne.0) n1=n1+1
      end do
!$acc parallel loop async default(present)
      do i = 1,sizr
         if (de_buffr(i).ne.0) nr=nr+1
      end do
!$acc wait
      print '(A,3I10)', 'check_zero_f',n0,n1,nr
      end subroutine

      subroutine reset_forces_d
      integer(mipk) siz0,siz1,sizr,offset

      offset = dr_stride3
      siz0   = dr_nb0  *int(dr_stride3 ,mipk)
      siz1   = dr_nb1  *int(dr_stride3 ,mipk)
      call mem_set( de_buff0,zero_rp,siz0,rec_stream )
      call mem_set( de_buff1,zero_md,siz1-offset,rec_stream,offset )

      end subroutine

      module subroutine zero_forces_rec
      integer(mipk) sizr

      sizr = dr_nbnbr*int(dr_strider3,mipk)
      call mem_set( de_buffr,zero_rp,sizr,rec_stream )

      end subroutine

      module subroutine get_ftot(derivs,nbloc_)
      integer   nbloc_
      real(r_p) derivs(3,nbloc_)
      integer i,j

      if (tdes_l) then
!$acc parallel loop async default(present)
         do i = 1,nbloc_
            derivs(1,i) = mdr2md(derivx(i))
            derivs(2,i) = mdr2md(derivy(i))
            derivs(3,i) = mdr2md(derivz(i))
         end do
      else
!$acc parallel loop collapse(2) async default(present)
         do i = 1,nbloc_; do j = 1,3
            derivs(j,i) = mdr2md(de_tot(j,i))
         end do; end do
      end if
      end subroutine

      ! Perform on mdyn_rtyp
      module subroutine add_forces(deadd)
      implicit none
      mdyn_rtyp deadd(*)
      integer i,j,id1,id2,idx
      integer(mipk) stride
      real(r_p) tot_
      mdyn_rtyp totr

      stride = dr_stride3


      !save the values of the short range real space polarizable forces
      if (shortnonbonded_l.and.stepint.eq.nalt)
     &   call mem_move(desave,dep,stride,rec_stream)

      if (tdes_l) then  ! (n,3) buffer format

!$acc host_data use_device(deadd,de_buff0,de_buff1)
      if (bonded_l.and..not.nonbonded_l) then
         if (.not.fdebs_l) then ! Check merge and transposition
            if(deb_Path) write(*,*) "   add_forces bonded"
!$acc parallel loop async(rec_queue) deviceptr(de_buff0,deadd)
            do i = 1,3*nbloc
               tot_ = 0
!$acc loop seq
               do j = 0,dr_nbb-1
                  idx  = dr_obb + i + j*stride
                  tot_ = tot_ + de_buff0(idx)
                  if(n_debf) de_buff0(idx) = 0
               end do
               id1 = mod(i-1,3)
               id2 = (i-1)/3
               idx = id1*dr_stride + id2+1 
               deadd(idx) = deadd(idx) + rp2mdr(tot_)
            end do
         end if
      else if (nonbonded_l.and..not.bonded_l) then
         if(deb_Path) write(*,*) "   add_forces non bonded"
!$acc parallel loop async(rec_queue) deviceptr(de_buff1,deadd)
         do i = 1,3*nbloc
            totr = 0
!$acc loop seq
            do j = 0,dr_nbnb-1
               idx  = dr_obnb + i + j*stride
               totr = totr + de_buff1(idx)
               if(n_debf) de_buff1(idx) = 0
            end do
            id1 = mod(i-1,3)
            id2 = (i-1)/3
            idx = id1*dr_stride + id2+1 
            deadd(idx) = deadd(idx) + totr
         end do
      else
         if(deb_Path) write(*,*) "   add_forces"
!$acc parallel loop async(rec_queue)
!$acc&         deviceptr(de_buff0,de_buff1,deadd)
         do i = 1,3*nbloc
            tot_ = 0
            totr = 0
!$acc loop seq
            do j = 0,dr_nbb-1
               idx  = dr_obb + i + j*stride
               tot_ = tot_ + de_buff0(idx)
               if(n_debf) de_buff0(idx) = 0
            end do
!$acc loop seq
            do j = 0,dr_nbnb-1
               idx  = dr_obnb + i + j*stride
               totr = totr + de_buff1(idx)
               if(n_debf) de_buff1(idx) = 0
            end do
            id1 = mod(i-1,3)
            id2 = (i-1)/3
            idx = id1*dr_stride + id2+1 
            deadd(idx) = deadd(idx) + totr + rp2mdr(tot_)
         end do
      end if
!$acc end host_data

      else  ! (n,3) buffer format

!$acc host_data use_device(deadd,de_buff0,de_buff1)
      if (bonded_l.and..not.nonbonded_l) then
         if(deb_Path) write(*,*) "   add_forces bonded"
!$acc parallel loop async(rec_queue) deviceptr(de_buff0,deadd)
         do i = 1,3*nbloc
            tot_ = 0
!$acc loop seq
            do j = 0,dr_nbb-1
               idx  = dr_obb + i + j*stride
               tot_ = tot_ + de_buff0(idx)
               if(n_debf) de_buff0(idx) = 0
            end do
            deadd(i) = deadd(i) + rp2mdr(tot_)
         end do
      else if (nonbonded_l.and..not.bonded_l) then
         if(deb_Path) write(*,*) "   add_forces non bonded"
!$acc parallel loop async(rec_queue) deviceptr(de_buff1,deadd)
         do i = 1,3*nbloc
            totr = 0
!$acc loop seq
            do j = 0,dr_nbnb-1
               idx  = dr_obnb + i + j*stride
               totr = totr + de_buff1(idx)
               if(n_debf) de_buff1(idx) = 0
            end do
            deadd(i) = deadd(i) + totr
         end do
      else
         if(deb_Path) write(*,*) "   add_forces"
!$acc parallel loop async(rec_queue)
!$acc&         deviceptr(de_buff0,de_buff1,deadd)
         do i = 1,3*nbloc
            tot_ = 0
            totr = 0
!$acc loop seq
            do j = 0,dr_nbb-1
               idx  = dr_obb + i + j*stride
               tot_ = tot_ + de_buff0(idx)
               if(n_debf) de_buff0(idx) = 0
            end do
!$acc loop seq
            do j = 0,dr_nbnb-1
               idx  = dr_obnb + i + j*stride
               totr = totr + de_buff1(idx)
               if(n_debf) de_buff1(idx) = 0
            end do
            deadd(i) = deadd(i) + totr + rp2mdr(tot_)
         end do
      end if
!$acc end host_data

      end if  ! (n,3) buffer format
      end subroutine

      ! Perform on real(r_p)
      module subroutine add_forces1(deadd)
      real(r_p) deadd(*)
      integer i,j
      integer(mipk) idx,stride
      real(r_p) tot_
      mdyn_rtyp totr

      stride = dr_stride3

      !save the values of the short range real space polarizable forces
      if (shortnonbonded_l.and.stepint.eq.nalt)
     &   call mem_move(desave,dep,stride,rec_stream)

!$acc host_data use_device(de_buff0,de_buff1,deadd)
      if (bonded_l.and..not.nonbonded_l) then
!$acc parallel loop async(rec_queue) deviceptr(de_buff0,deadd)
         do i = 1,3*nbloc
            tot_ = 0
!$acc loop seq
            do j = 0,dr_nbb-1
               idx  = dr_obb + i + j*stride
               tot_ = tot_ + de_buff0(idx)
               if(n_debf) de_buff0(idx) = 0
            end do
            deadd(i) = (tot_)
         end do
      else if (nonbonded_l.and..not.bonded_l) then
!$acc parallel loop async(rec_queue) deviceptr(de_buff1,deadd)
         do i = 1,3*nbloc
            totr = 0
!$acc loop seq
            do j = 0,dr_nbnb-1
               idx  = dr_obnb + i + j*stride
               totr = totr + de_buff1(idx)
               if(n_debf) de_buff1(idx) = 0
            end do
            deadd(i) = mdr2md(totr)
         end do
      else
!$acc parallel loop async(rec_queue) deviceptr(de_buff0,de_buff1,deadd)
         do i = 1,3*nbloc
            tot_ = 0
            totr = 0
!$acc loop seq
            do j = 0,dr_nbb-1
               idx  = dr_obb + i + j*stride
               tot_ = tot_ + de_buff0(idx)
               if(n_debf) de_buff0(idx) = 0
            end do
!$acc loop seq
            do j = 0,dr_nbnb-1
               idx  = dr_obnb + i + j*stride
               totr = totr + de_buff1(idx)
               if(n_debf) de_buff1(idx) = 0
            end do
            deadd(i) = mdr2md(totr) + tot_
         end do
      end if
!$acc end host_data
      end subroutine

      ! Perform on mdyn_rtyp
      module subroutine add_forces_rec(deadd)
      mdyn_rtyp deadd(3*nlocrec2)
      integer i,j,id1,id2,ilocrec,stride
      integer(mipk) idx
      real(r_p) tot_
      mdyn_rtyp totr

      if(deb_Path) write(*,*) "   add_forces_rec"

      stride = dr_strider3
      
!$acc host_data use_device(deadd,de_buffr)
      if (tdes_l) then

!$acc parallel loop async(rec_queue) deviceptr(deadd,de_buffr)
      do i = 1,3*nlocrec2
         tot_ = 0
!$acc loop seq
         do j = 0,dr_nbnbr-1
            idx  = dr_obnbr + i + j*stride
            tot_ = tot_ + de_buffr(idx)
            if(n_debf) de_buffr(idx) = 0
         end do
         id1 = mod(i-1,3)
         id2 = (i-1)/3
         idx = id1*dr_stride + id2+1 
         deadd(idx) = deadd(idx) + rp2mdr(tot_)
      end do

      else

!$acc parallel loop async(rec_queue) deviceptr(deadd,de_buffr)
      do i = 1,3*nlocrec2
         tot_ = 0
!$acc loop seq
         do j = 0,dr_nbnbr-1
            idx  = dr_obnbr + i + j*stride
            tot_ = tot_ + de_buffr(idx)
            if(n_debf) de_buffr(idx) = 0
         end do
         deadd(i) = deadd(i) + rp2mdr(tot_)
      end do

      end if
!$acc end host_data
      end subroutine

      ! Perform on real(r_p)
      module subroutine add_forces_rec1(deadd)
      real(r_p) deadd(3*nlocrec2)
      integer i,j,ilocrec,stride
      integer(mipk) idx
      real(r_p) tot_
      mdyn_rtyp totr

      if(deb_Path) write(*,*) "   add_forces_rec1"

      stride = dr_strider3
!$acc host_data use_device(deadd,de_buffr)
!$acc parallel loop async(rec_queue) deviceptr(deadd,de_buffr)
      do i = 1,3*nlocrec2
         tot_ = 0
!$acc loop seq
         do j = 0,dr_nbnbr-1
            idx  = dr_obnbr + i + j*stride
            tot_ = tot_ + de_buffr(idx)
            if(n_debf) de_buffr(idx) = 0
         end do
         deadd(i) = deadd(i) + (tot_)
      end do
!$acc end host_data
      end subroutine

      module subroutine sum_forces_rec1(desum)
      real(r_p) desum(3*nlocrec2)
      integer i,j,ilocrec,stride
      integer(mipk) idx
      real(r_p) tot_
      mdyn_rtyp totr

      if(deb_Path) write(*,*) "   sum_add_forces_rec1"

      stride = dr_strider3
!$acc parallel loop async(rec_queue) default(present)
      do i = 1,3*nlocrec2
         tot_ = 0
!$acc loop seq
         do j = 0,dr_nbnbr-1
            idx  = dr_obnbr + i + j*stride
            tot_ = tot_ + de_buffr(idx)
            if(n_debf) de_buffr(idx) = 0
         end do
         desum(i) = (tot_)
      end do
      end subroutine

      module subroutine add_forces_rec_1d(deadd)
      mdyn_rtyp deadd(3*nbloc)
      integer i,j,id1,id2,ilocrec,i1,i2,stride
      integer(mipk) idx
      real(r_p) tot_
      mdyn_rtyp totr

      if(deb_Path) write(*,*) "    add_forces_rec_1d"

      stride = dr_strider3

      if (tdes_l) then

!$acc parallel loop async(rec_queue) default(present)
      do i = 0,3*nbloc-1
         tot_    = 0
         i1      = mod(i,3)
         i2      = i/3 + 1
         ilocrec = locrec(glob(i2)) -1
         if (ilocrec.eq.-1) cycle
!$acc loop seq
         do j = 0,dr_nbnbr-1
            idx  = dr_obnbr + i1+1+ilocrec*3 + j*stride
            tot_ = tot_ + de_buffr(idx)
            if(n_debf) de_buffr(idx) = 0
         end do
         idx        = i2 + i1*dr_stride  
         deadd(idx) = deadd(idx) + rp2mdr(tot_)
      end do

      else

!$acc parallel loop async(rec_queue) default(present)
      do i = 1,3*nbloc
         tot_    = 0
         i1      = mod(i-1,3) +1
         ilocrec = locrec(glob((i-1)/3+1)) -1
         if (ilocrec.eq.-1) cycle
!$acc loop seq
         do j = 0,dr_nbnbr-1
            idx  = dr_obnbr + ilocrec*3+i1 + j*stride
            tot_ = tot_ + de_buffr(idx)
            if(n_debf) de_buffr(idx) = 0
         end do
         deadd(i) = deadd(i) + rp2mdr(tot_)
      end do

      end if
      end subroutine

      module subroutine add_forces_rec_1d1(deadd,derec)
      mdyn_rtyp deadd(3*nbloc)
      real(r_p) derec(3*nlocrec2)
      integer i,j,id1,id2,ilocrec,i1,i2,stride
      integer(mipk) idx,idx1
      real(r_p) tot_
      mdyn_rtyp totr

      if(deb_Path) write(*,*) "    add_forces_rec_1d"

      stride = dr_strider3

      if (tdes_l) then

!$acc parallel loop async(rec_queue) default(present)
      do i = 0,3*nbloc-1
         i1      = mod(i,3)
         i2      = i/3 + 1
         ilocrec = locrec(glob(i2)) -1
         if (ilocrec.eq.-1) cycle
         idx1    = i1+1 + ilocrec*3
         idx        = i2 + i1*dr_stride  
         deadd(idx) = deadd(idx) + rp2mdr(derec(idx1))
      end do

      else

!$acc parallel loop async(rec_queue) default(present)
      do i = 1,3*nbloc
         i1      = mod(i-1,3) +1
         ilocrec = locrec(glob((i-1)/3+1)) -1
         if (ilocrec.eq.-1) cycle
         idx1    = ilocrec*3 + i1
         deadd(i) = deadd(i) + rp2mdr(derec(idx1))
      end do

      end if
      end subroutine

      module subroutine add_forces_rec1_1d(deadd)
      real(r_p) deadd(3*nbloc)
      integer i,j,ilocrec,i1,i2,stride
      integer(mipk) idx
      real(r_p) tot_
      mdyn_rtyp totr

      if(deb_Path) write(*,*) "    add_forces_rec1_1d1"

      stride = dr_strider3
!$acc parallel loop async(rec_queue) default(present)
      do i = 1,3*nbloc
         tot_    = 0
         i1      = mod(i-1,3) +1
         ilocrec = locrec(glob((i-1)/3+1)) -1
         if (ilocrec.eq.-1) cycle
!$acc loop seq
         do j = 0,dr_nbnbr-1
            idx  = dr_obnbr + ilocrec*3+i1 + j*stride
            tot_ = tot_ + de_buffr(idx)
            de_buffr(idx) = 0
         end do
         deadd(i) = deadd(i) + (tot_)
      end do
      end subroutine

      module subroutine add_forces_rec1_1d1(deadd,derec)
      real(r_p) deadd(3*nbloc), derec(3*nlocrec2)
      integer i,j,ilocrec,i1,i2,stride
      integer(mipk) idx1
      real(r_p) tot_
      mdyn_rtyp totr

      if(deb_Path) write(*,*) "    add_forces_rec1_1d1"

      stride = dr_strider3
!$acc parallel loop async(rec_queue) default(present)
      do i = 1,3*nbloc
         i1       = mod(i-1,3) +1
         ilocrec  = locrec(glob((i-1)/3+1)) -1
         if (ilocrec.eq.-1) cycle
         idx1     = ilocrec*3 + i1
         deadd(i) = deadd(i) + derec(idx1)
      end do
      end subroutine

      module subroutine remove_desave(derivs)
      real(r_p) derivs(3,nbloc)
      integer i,j
      if (ftot_l) then
         if (tdes_l) then
!$acc parallel loop
!$acc&         default(present) async
            do i = 1,nbloc;
               derivx(i) = derivx(i) - desave(1,i)
               derivy(i) = derivy(i) - desave(2,i)
               derivz(i) = derivz(i) - desave(3,i)
            end do;
         else
!$acc host_data use_device(de_tot,desave)
!$acc parallel loop collapse(2) 
!$acc&         deviceptr(de_tot,desave) async
            do i = 1,nbloc; do j = 1,3
               de_tot(j,i) = de_tot(j,i) - desave(j,i)
            end do; end do
!$acc end host_data
         end if

      else

!$acc host_data use_device(derivs,desave)
!$acc parallel loop collapse(2) 
!$acc&         deviceptr(derivs,desave) async
      do i = 1,nbloc; do j = 1,3
         derivs(j,i) = derivs(j,i) - mdr2md(desave(j,i))
      end do; end do
!$acc end host_data

      end if
      end subroutine

      module subroutine resetForcesAMD
      implicit none
      integer i,j
      if (deb_Path) print*,'resetForcesAMD'
!$acc parallel loop collapse(2) default(present) async
      do i = 1, nbloc; do j = 1, 3
         deamdD(j,i) = 0.0_re_p
        deW1aMD(j,i) = 0.0_re_p
      end do; end do
      end subroutine

      ! attach commforces mpi_buffer
      subroutine load_mpi_buffer( opt )
      implicit none
      integer,intent(in):: opt
      integer i,j

      if      (opt.eq.0) then
         if (nproc.eq.1.or..not.tdes_l) then
            de_mpi(1:3,1:dr_stride) => de_buff1(1:dr_stride3)
         else
            de_mpi(1:3,1:dr_stride) => 
     &      de_buff1(dr_obws+1:dr_obws+dr_stride3)
!$acc parallel loop collapse(2) async default(present)
            do i = 1,dr_stride; do j = 1,3
               de_mpi(j,i) = derivx(i+(j-1)*dr_stride)
            end do; end do
         end if
      else if (opt.eq.1) then
         if (nproc.ne.1.and.tdes_l) then
!$acc parallel loop async default(present)
            do i=1,dr_stride
                 derivx(i) = de_mpi(1,i)
                 derivy(i) = de_mpi(2,i)
                 derivz(i) = de_mpi(3,i)
               de_mpi(1,i) = 0
               de_mpi(2,i) = 0
               de_mpi(3,i) = 0
            end do
         end if
         nullify(de_mpi)
      end if

      end subroutine

      ! communicate forces between processes
      module subroutine comm_forces_dd(derivs,opt)
      implicit none
      mdyn_rtyp derivs(3,*)
      integer,optional:: opt
      integer i,j,k,tag,ierr,siz,opt_
      integer nsend,nrecv,psend,precv,pbeg,dlen
      integer commloc,reqsend(nproc),reqrecv(nproc)
      integer status(MPI_STATUS_SIZE)
      parameter(tag=0)

      integer  ,pointer:: p_send(:),p_recv(:),p_beg(:),d_len(:)
      mdyn_rtyp,pointer:: buffer(:,:,:)
c
      if (nproc.eq.1.or.((use_pmecore).and.(rank.gt.ndir-1))) return
c
      opt_ = cNBond
      if (present(opt)) opt_=opt
c
      commloc = COMM_TINKER
      if      (opt_.eq.cBond) then
         nsend     =  nneig_send
         nrecv     =  nneig_recep
         p_send(1:nproc) => pneig_send (:)
         p_recv(1:nproc) => pneig_recep(:)
         p_beg (1:nproc) => bufbeg(:)
         d_len (1:nproc) => domlen(:)
      else if (opt_.eq.cSNBond) then
         nsend     =  nbigshort_send
         nrecv     =  nbigshort_recep
         p_send(1:nproc) => pbigshort_send (:)
         p_recv(1:nproc) => pbigshort_recep(:)
         p_beg (1:nproc) => bufbeg(:)
         d_len (1:nproc) => domlen(:)
      else if (opt_.eq.cNBond) then
         nsend     =  nbig_send
         nrecv     =  nbig_recep
         p_send(1:nproc) => pbig_send (:)
         p_recv(1:nproc) => pbig_recep(:)
         p_beg (1:nproc) => bufbeg(:)
         d_len (1:nproc) => domlen(:)
      else
 16      format(' ERROR! routine comm_forces_dd' 
     &       ,/,' --- Available options ',3I3
     &       ,/,' ---   UNKNOWN OPTION  ',I3   )
         write(0,16) cBond,cSNBond,cNBond,opt_
         call fatal
      end if
c
      if (nsend.eq.0.and.nrecv.eq.0) return
c
      if (deb_Path) write(*,*) '   >> comm_forces_dd',opt_
      call timer_enter( timer_dirbondfcomm )
c
      ! Create exchange buffer
      siz = 3*max(1,nloc)*nsend
      call prmem_requestm(buffMpi_p1,siz,async=.false.)
      call c_f_pointer(c_loc(buffMpi_p1),buffer,[3,max(1,nloc),nsend])
      !buffer(1:3,1:max(1,nloc),1:nsend) => buffMpi_p1(1:siz)
c
      !MPI : Start reception in buffer
!$acc host_data use_device(buffer,derivs)
      do i = 1, nsend
         psend = p_send(i)
         dlen  = 3*nloc
         call MPI_IRECV(buffer(1,1,i),dlen,MPI_MDTYP,psend,tag
     $                 ,COMM_TINKER,reqrecv(i),ierr)
      end do
c
      !MPI : Start sending
!$acc wait(rec_queue)
      do i = 1, nrecv
         precv = p_recv(i)
         pbeg  = p_beg ( precv+1 )
         dlen  = d_len ( precv+1 )*3
         call MPI_ISEND(derivs(1,pbeg),dlen,MPI_MDTYP,precv,tag
     $                 ,COMM_TINKER,reqsend(i),ierr)
      end do
!$acc end host_data
c
      ! Wait for Exchanges to end
      do i = 1,nsend; call MPI_WAIT(reqrecv(i),status,ierr); end do;
      do i = 1,nrecv; call MPI_WAIT(reqsend(i),status,ierr); end do;
c
c     MPI : move in global arrays
c
!$acc parallel loop gang vector collapse(2) async
!$acc&         present(derivs,buffer)
      do j = 1, nloc; do k = 1, 3; do i = 1, nsend;
         derivs(k,j) = derivs(k,j) + buffer(k,j,i)
      end do; end do; end do
c
      nullify(buffer)

      if (deb_Path) write(*,*) '   << comm_forces_dd',opt_
 20   continue
      call timer_exit( timer_dirbondfcomm,quiet_timers )
      end subroutine comm_forces_dd

      ! communicate forces between processes
      module subroutine comm_forces_dd1(derivs,opt)
      implicit none
      real(r_p) derivs(3,*)
      integer,optional:: opt
      integer i,j,k,tag,ierr,siz,opt_
      integer nsend,nrecv,psend,precv,pbeg,dlen
      integer commloc,reqsend(nproc),reqrecv(nproc)
      integer status(MPI_STATUS_SIZE)
      parameter(tag=0)

      integer  ,pointer:: p_send(:),p_recv(:),p_beg(:),d_len(:)
      real(r_p),pointer:: buffer(:,:,:)
c
      if (nproc.eq.1.or.((use_pmecore).and.(rank.gt.ndir-1))) return
c
      opt_ = cNBond
      if (present(opt)) opt_=opt
c
      commloc = COMM_TINKER
      if      (opt_.eq.cBond) then
         nsend     =  nneig_send
         nrecv     =  nneig_recep
         p_send(1:nproc) => pneig_send (:)
         p_recv(1:nproc) => pneig_recep(:)
         p_beg (1:nproc) => bufbeg(:)
         d_len (1:nproc) => domlen(:)
      else if (opt_.eq.cSNBond) then
         nsend     =  nbigshort_send
         nrecv     =  nbigshort_recep
         p_send(1:nproc) => pbigshort_send (:)
         p_recv(1:nproc) => pbigshort_recep(:)
         p_beg (1:nproc) => bufbeg(:)
         d_len (1:nproc) => domlen(:)
      else if (opt_.eq.cNBond) then
         nsend     =  nbig_send
         nrecv     =  nbig_recep
         p_send(1:nproc) => pbig_send (:)
         p_recv(1:nproc) => pbig_recep(:)
         p_beg (1:nproc) => bufbeg(:)
         d_len (1:nproc) => domlen(:)
      else
 16      format(' ERROR! routine comm_forces_dd1' 
     &       ,/,' --- Available options ',3I3
     &       ,/,' ---   UNKNOWN OPTION  ',I3   )
         write(0,16) cBond,cSNBond,cNBond,opt_
         call fatal
      end if
c
      if (nsend.eq.0.and.nrecv.eq.0) return
c
      if (deb_Path) write(*,*) '   >> comm_forces_dd1',opt_
      call timer_enter( timer_dirbondfcomm )
c
      ! Create exchange buffer
      siz = 3*max(1,nloc)*nsend
      call prmem_requestm(buffMpi_p1,siz,async=.false.)
      buffer(1:3,1:max(1,nloc),1:nsend) => buffMpi_p1(1:siz)
c
      !MPI : Start reception in buffer
!$acc host_data use_device(buffer,derivs)
      do i = 1, nsend
         psend = p_send(i)
         dlen  = 3*nloc
         call MPI_IRECV(buffer(1,1,i),dlen,MPI_RPREC,psend,tag
     $                 ,COMM_TINKER,reqrecv(i),ierr)
      end do
c
      !MPI : Start sending
!$acc wait(rec_queue)
      do i = 1, nrecv
         precv = p_recv(i)
         pbeg  = p_beg ( precv+1 )
         dlen  = d_len ( precv+1 )*3
         call MPI_ISEND(derivs(1,pbeg),dlen,MPI_RPREC,precv,tag
     $                 ,COMM_TINKER,reqsend(i),ierr)
      end do
!$acc end host_data
c
      ! Wait for Exchanges to end
      do i = 1,nsend; call MPI_WAIT(reqrecv(i),status,ierr); end do;
      do i = 1,nrecv; call MPI_WAIT(reqsend(i),status,ierr); end do;
c
c     MPI : move in global arrays
c
!$acc parallel loop gang vector collapse(2) async(rec_queue)
!$acc&         present(derivs,buffer)
      do j = 1, nloc; do k = 1, 3; do i = 1, nsend;
         derivs(k,j) = derivs(k,j) + buffer(k,j,i)
      end do; end do; end do
c
      nullify(buffer)

      if (deb_Path) write(*,*) '   << comm_forces_dd1',opt_
 20   continue
      call timer_exit( timer_dirbondfcomm,quiet_timers )
      end subroutine comm_forces_dd1

      ! Communicate reciprocal forces between processes
      module subroutine comm_forces_rec(de_rec,opt)
      implicit none
      real(r_p) de_rec(3,nlocrec2)
      integer,optional:: opt
      integer i,j,k,tag,ierr,siz,opt_,commloc
      integer nsend,nrecv,psend,precv,pbeg,dlen
      integer reqsend(nproc),reqrecv(nproc)
      integer status(MPI_STATUS_SIZE)
      parameter(tag=0)

      integer  ,pointer:: p_send(:),p_recv(:),p_beg(:),d_len(:)
      real(r_p),pointer:: buffer(:,:,:)
c
      if (nproc.eq.1) return
c
      nsend     =  nrec_send
      nrecv     =  nrec_recep
      if (nsend.eq.0.and.nrecv.eq.0) return
c
      if (deb_Path) write(*,*) '   >> comm_forces_rec'
      call timer_enter( timer_dirreccomm )
c
      p_send(1:nproc) => prec_send (:)
      p_recv(1:nproc) => prec_recep(:)
      p_beg (1:nproc) => bufbegrec(:)
      d_len (1:nproc) => domlenrec(:)
      commloc = merge(comm_rec,COMM_TINKER,use_pmecore)
c
      ! Create exchange buffer
      siz = max(size(buffMpi_p1), 3*max(1,nlocrec)*nsend)
      call prmem_requestm(buffMpi_p1,siz,async=.false.)
      buffer(1:3,1:max(1,nlocrec),1:nsend) => buffMpi_p1(1:siz)
c
      !MPI : Start reception in buffer
!$acc host_data use_device(buffer,de_rec)
      do i = 1, nsend
         psend = p_send(i)
         dlen  = 3*nlocrec
         call MPI_IRECV(buffer(1,1,i),dlen,MPI_RPREC,psend,tag
     $                 ,commloc,reqrecv(i),ierr)
      end do
c
      !MPI : Start sending
!$acc wait(rec_queue)
      do i = 1, nrecv
         precv = p_recv(i)
         pbeg  = p_beg ( precv+1 )
         dlen  = d_len ( precv+1 )*3
         call MPI_ISEND(de_rec(1,pbeg),dlen,MPI_RPREC,precv,tag
     $                 ,commloc,reqsend(i),ierr)
      end do
!$acc end host_data
c
      ! Wait for Exchanges to end
      do i = 1,nsend; call MPI_WAIT(reqrecv(i),status,ierr); end do;
      do i = 1,nrecv; call MPI_WAIT(reqsend(i),status,ierr); end do;
c
c     MPI : move in global arrays
c
!$acc parallel loop gang vector collapse(2) async
!$acc&         present(de_rec,buffer)
      do j = 1,nlocrec; do k = 1,3; do i = 1,nsend;
         de_rec(k,j) = de_rec(k,j) + buffer(k,j,i)
      end do; end do; end do
c
      nullify(buffer)

      if (deb_Path) write(*,*) '   << comm_forces_rec'
      call timer_exit( timer_dirreccomm,quiet_timers )
      end subroutine comm_forces_rec

      !
      ! Communicate reciprocal forces to direct ones
      module subroutine comm_forces_recdir(derivs,de_rec,opt)
      implicit none
      mdyn_rtyp,intent(inout):: derivs(dr_stride3)
      real(r_p),intent(in)::    de_rec(3,nlocrec2)
      integer  ,optional:: opt

      integer reqsend(nproc),reqrec(nproc)
      integer i,iglob,j,k,iloc,jloc,jlocrec,jglob
     &       ,iproc,ibufbeg,idlen,sz1,sz2,opt_,
     &       rankloc,tag,ierr,status(MPI_STATUS_SIZE)
      real(r_p),pointer:: buffer(:,:,:),buffers(:,:)
      parameter(tag=0)

      call timer_enter( timer_dirreccomm )
      if (deb_Path) write(*,*) '   << comm_forces_recdir'

      ! Reserve buffers space
      rankloc =   merge(rank_bis,rank,use_pmecore)
      sz1     = 3*max(1,nloc)*nrecdir_send1
      sz2     = 3*max(1,nblocrecdir)
      call prmem_requestm(buffMpi_p1,max(size(buffMpi_p1),sz1))
      call prmem_requestm(buffMpi_p2,max(size(buffMpi_p2),sz2))
      buffer(1:3,1:max(1,nloc),1:nrecdir_send1) => buffMpi_p1(1:sz1)
      buffers(1:3,1:max(1,nblocrecdir))         => buffMpi_p2(1:sz2)
 
      !MPI : begin reception in buffer
!$acc host_data use_device(buffer)
      do i = 1, nrecdir_send1
         if (precdir_send1(i).ne.rank) then
           !tag = nproc*rank + precdir_send1(i) + 1
           call MPI_IRECV(buffer(1,1,i),3*nloc,MPI_RPREC
     &         ,precdir_send1(i),tag,COMM_TINKER,reqrec(i),ierr)
         end if
      end do
!$acc end host_data
 
      !Move in buffers
      do i = 1, nrecdir_recep1
        if (precdir_recep1(i).ne.rank) then
          ibufbeg = bufbeg(precdir_recep1(i)+1)-1 
          idlen   = domlen(precdir_recep1(i)+1)
!$acc parallel loop collapse(2) async(rec_queue)
!$acc&    present(repartrec,locrec,glob,de_rec,buffers)
          do j = 1, idlen
            do k = 1,3
              jloc = ibufbeg+j
              jglob = glob(jloc)
              if (repartrec(jglob).eq.rankloc) then
                jlocrec = locrec(jglob)
                buffers(k,jloc) = de_rec(k,jlocrec)
              else
                buffers(k,jloc) = 0
              end if
            end do
          end do
        end if
      end do
c
      ! MPI : Start Sending
!$acc wait(rec_queue)
!$acc host_data use_device(buffers)
      do i = 1, nrecdir_recep1
        if (precdir_recep1(i).ne.rank) then
          !tag = nproc*precdir_recep1(i) + rank + 1
          call MPI_ISEND(buffers(1,bufbeg(precdir_recep1(i)+1))
     &        ,3*domlen(precdir_recep1(i)+1),MPI_RPREC
     &        ,precdir_recep1(i),tag,COMM_TINKER,reqsend(i),ierr)
        end if
      end do
!$acc end host_data
c
      do i = 1, nrecdir_send1
        if (precdir_send1(i).ne.rank)
     &  call MPI_WAIT(reqrec(i),status,ierr)
      end do
c
c     MPI : move in output
c
      if (tdes_l) then

      do i = 1, nrecdir_send1
         if (precdir_send1(i).ne.rank) then
!$acc parallel loop collapse(2) async default(present)
            do k = 0,2; do j = 1,nloc;
               derivs(j+k*dr_stride) =
     &         derivs(j+k*dr_stride) + rp2mdr(buffer(k+1,j,i))
            end do; end do
         else
!$acc parallel loop collapse(2) async default(present)
            do k = 0,2; do j = 1,nlocrec;
               iglob = globrec(j)
               iloc  = loc(iglob)
               if (repart(iglob).eq.rank) then
                  derivs(iloc+k*dr_stride) =
     &            derivs(iloc+k*dr_stride) + rp2mdr(de_rec(k+1,j))
               end if
            end do; end do
         end if
      end do

      else

      do i = 1, nrecdir_send1
         if (precdir_send1(i).ne.rank) then
!$acc parallel loop collapse(2) async default(present)
            do j = 1,nloc; do k = 1,3
               derivs(k+(j-1)*3) = derivs(k+(j-1)*3) 
     &                           + rp2mdr(buffer(k,j,i))
            end do; end do
         else
!$acc parallel loop collapse(2) async default(present)
            do j = 1, nlocrec; do k = 1,3
               iglob = globrec(j)
               iloc  = loc(iglob)
               if (repart(iglob).eq.rank) then
                  derivs(k+(iloc-1)*3) = derivs(k+(iloc-1)*3)
     &                                 + rp2mdr(de_rec(k,j))
               end if
            end do; end do
         end if
      end do

      end if

      do i = 1, nrecdir_recep1
        if (precdir_recep1(i).ne.rank) then
        call MPI_WAIT(reqsend(i),status,ierr)
        end if
      end do

      nullify(buffer)
      nullify(buffers)

      if (deb_Path) write(*,*) '   << comm_forces_recdir'
      call timer_exit( timer_dirreccomm,quiet_timers )
      end subroutine comm_forces_recdir

      module subroutine comm_forces_recdir1(derivs,de_rec,opt)
      implicit none
      real(r_p),intent(inout):: derivs(3,nbloc)
      real(r_p),intent(in)::    de_rec(3,nlocrec2)
      integer  ,optional:: opt

      integer reqsend(nproc),reqrec(nproc)
      integer i,iglob,j,k,iloc,jloc,jlocrec,jglob
     &       ,iproc,ibufbeg,idlen,sz1,sz2,opt_,
     &       rankloc,tag,ierr,status(MPI_STATUS_SIZE)
      real(r_p),pointer:: buffer(:,:,:),buffers(:,:)
      parameter(tag=0)

      call timer_enter( timer_dirreccomm )
      if (deb_Path) write(*,*) '   << comm_forces_recdir1'

      ! Reserve buffers space
      rankloc =   merge(rank_bis,rank,use_pmecore)
      sz1     = 3*max(1,nloc)*nrecdir_send1
      sz2     = 3*max(1,nblocrecdir)
      call prmem_requestm(buffMpi_p1,max(size(buffMpi_p1),sz1))
      call prmem_requestm(buffMpi_p2,max(size(buffMpi_p2),sz2))
      buffer(1:3,1:max(1,nloc),1:nrecdir_send1) => buffMpi_p1(1:sz1)
      buffers(1:3,1:max(1,nblocrecdir))         => buffMpi_p2(1:sz2)
 
      !MPI : begin reception in buffer
!$acc host_data use_device(buffer)
      do i = 1, nrecdir_send1
         if (precdir_send1(i).ne.rank) then
           !tag = nproc*rank + precdir_send1(i) + 1
           call MPI_IRECV(buffer(1,1,i),3*nloc,MPI_RPREC
     &         ,precdir_send1(i),tag,COMM_TINKER,reqrec(i),ierr)
         end if
      end do
!$acc end host_data
 
      !Move in buffers
      do i = 1, nrecdir_recep1
        if (precdir_recep1(i).ne.rank) then
          ibufbeg = bufbeg(precdir_recep1(i)+1)-1 
          idlen   = domlen(precdir_recep1(i)+1)
!$acc parallel loop collapse(2) async(rec_queue)
!$acc&    present(repartrec,locrec,glob,de_rec,buffers)
          do j = 1, idlen
            do k = 1,3
              jloc = ibufbeg+j
              jglob = glob(jloc)
              if (repartrec(jglob).eq.rankloc) then
                jlocrec = locrec(jglob)
                buffers(k,jloc) = de_rec(k,jlocrec)
              else
                buffers(k,jloc) = 0
              end if
            end do
          end do
        end if
      end do
c
      ! MPI : Start Sending
!$acc wait(rec_queue)
!$acc host_data use_device(buffers)
      do i = 1, nrecdir_recep1
        if (precdir_recep1(i).ne.rank) then
          !tag = nproc*precdir_recep1(i) + rank + 1
          call MPI_ISEND(buffers(1,bufbeg(precdir_recep1(i)+1))
     &        ,3*domlen(precdir_recep1(i)+1),MPI_RPREC
     &        ,precdir_recep1(i),tag,COMM_TINKER,reqsend(i),ierr)
        end if
      end do
!$acc end host_data
c
      do i = 1, nrecdir_send1
        if (precdir_send1(i).ne.rank)
     &  call MPI_WAIT(reqrec(i),status,ierr)
      end do
c
c     MPI : move in output
c
      do i = 1, nrecdir_send1
         if (precdir_send1(i).ne.rank) then
!$acc parallel loop collapse(2) async default(present)
            do j = 1,nloc; do k = 1,3
               derivs(k,j) = derivs(k,j) + buffer(k,j,i)
            end do; end do
         else
!$acc parallel loop collapse(2) async default(present)
            do j = 1, nlocrec; do k = 1,3
               iglob = globrec(j)
               iloc  = loc(iglob)
               if (repart(iglob).eq.rank) then
                  derivs(k,iloc) = derivs(k,iloc) + de_rec(k,j)
               end if
            end do; end do
         end if
      end do

      do i = 1, nrecdir_recep1
        if (precdir_recep1(i).ne.rank) then
        call MPI_WAIT(reqsend(i),status,ierr)
        end if
      end do

      nullify(buffer)
      nullify(buffers)

      if (deb_Path) write(*,*) '   << comm_forces_recdir1'
      call timer_exit( timer_dirreccomm,quiet_timers )
      end subroutine comm_forces_recdir1

      module subroutine comm_forces(derivs,opt)
      implicit  none
      real(r_p) derivs(3,*)
      integer  ,optional:: opt
      integer   opt_

      opt_ = cNBond
      if (present(opt)) opt_=opt
      call timer_enter( timer_fcomm )

      if (opt_.eq.cNBond) call commforcesrec(derivs,opt_)

      if (ftot_l) then
         call load_mpi_buffer(0)
         call comm_forces_dd(de_mpi,opt_)
         call load_mpi_buffer(1)
      else
         call comm_forces_dd1(derivs,opt_)
      end if

      call timer_exit( timer_fcomm,quiet_timers )
      end subroutine

      subroutine info_Forces_one( array,name )
      implicit none
      character(*),intent(in) :: name
      real(r_p),intent(inout):: array(:,:)
      integer i,j
      real(r_p) mini,maxi
      real(8) norm_l1,temp

      mini = huge(mini)
      maxi = tiny(maxi)
      norm_l1 = 0.0d0
!$acc data present(array)

      call commforce_one(array)

!$acc wait
!$acc parallel loop
      do i = 1,nloc
         do j = 1,3
            mini = min( mini,array(j,i) )
            maxi = max( maxi,array(j,i) )
            norm_l1 = norm_l1 + abs(array(j,i))
         end do
      end do

!$acc end data
      if (rank.eq.0) then
         call MPI_REDUCE(MPI_IN_PLACE,mini,1,MPI_RPREC,
     &                   MPI_MIN,0,MPI_COMM_WORLD,i)
         call MPI_REDUCE(MPI_IN_PLACE,maxi,1,MPI_RPREC,
     &                   MPI_MAX,0,MPI_COMM_WORLD,i)
         call MPI_REDUCE(MPI_IN_PLACE,norm_l1,1,MPI_RPREC,
     &                   MPI_SUM,0,MPI_COMM_WORLD,i)
      else
         call MPI_REDUCE(mini,mini,1,MPI_RPREC,
     &                   MPI_MIN,0,MPI_COMM_WORLD,i)
         call MPI_REDUCE(maxi,maxi,1,MPI_RPREC,
     &                   MPI_MAX,0,MPI_COMM_WORLD,i)
         call MPI_REDUCE(norm_l1,norm_l1,1,MPI_RPREC,
     &                   MPI_SUM,0,MPI_COMM_WORLD,i)
      end if

 30   format(a7,a3,3F20.8)
      if (rank.eq.0) print 30,name,'=>',mini,maxi,norm_l1
      end subroutine

      module subroutine prtEForces(des,etot)
      use atoms  ,only: n
      use domdec
      use argue
      implicit none
      integer i, freeunit, iunit
      character(30) Ffile
      real(r_p) des(:,:)
      real(r_p) etot

      if (nproc.ne.1) then
         print '(10x,A)', "WARNING prtForces routine !!!"
         print '(10x,A)', "Not Operational in parallel execution " 
         return
      end if
      print*, 'prtForces'

!$acc wait
!$acc update host(des,etot)
      iunit = freeunit()
#if TINKER_SINGLE_PREC
      Ffile = trim(arg(1))//"_fs.txt"
#elif TINKER_MIXED_PREC
#  ifdef USE_DETERMINISTIC_REDUCTION
      Ffile = trim(arg(1))//"_ff.txt"
#  else
      Ffile = trim(arg(1))//"_fm.txt"
#  endif
#else
      Ffile = trim(arg(1))//"_fd.txt"
#endif
      call version(Ffile,'new')
      open( unit=iunit,file=Ffile,status='new' )
 12   format(3F18.10)
 13   format(F30.10)
      write (iunit,'(I0)') n
      write (iunit,13 ) etot

      do i = 1,n
         write(iunit,12) des(1:3,i)
      end do

      close(iunit)
      call fatal
      end subroutine


      ! Debug routine on Forces
      module subroutine info_forces(rule)
      !use atoms
      implicit none
      integer,intent(in) :: rule
      integer,parameter:: nf=28  !number of forces
      integer i,j,sze,ids
      real(8) mmx(3*nf),maxi
      logical tinker_isnan_m, save_arc,save_dcdio,focus_nbond,abortf
      real(r_p) dt
      integer,save:: doin = 1

      enum,bind(C)
      enumerator ::commBonded
      enumerator ::commShortNonBonded
      enumerator ::commNonBonded
      end enum

      mmx    = 0
      sze    = 3*nloc
      focus_nbond = merge(.true.,.false.,f_ulim.lt.700.0)

      if (deb_Path) write(*,*) 'info_forces', rule

      !TODO Add desave and desmd to debug functions
      !TODO Reduce communication amount in this function
!$acc wait
      if (btest(rule,idBond)) then

      if (.not.fdebs_l) then
         if(use_bond)   call comm_forces_dd1(deb  ,cBond)
         if(use_angle)  call comm_forces_dd1(dea  ,cBond)
         if(use_strbnd) call comm_forces_dd1(deba ,cBond)
         if(use_urey)   call comm_forces_dd1(deub ,cBond)
         if(use_angtor) call comm_forces_dd1(deat ,cBond)
         if(use_improp) call comm_forces_dd1(deid ,cBond)
         if(use_imptor) call comm_forces_dd1(deit ,cBond)
         if(use_tors)   call comm_forces_dd1(det  ,cBond)
         if(use_pitors) call comm_forces_dd1(dept ,cBond)
         if(use_tortor) call comm_forces_dd1(dett ,cBond)
         if(use_opbend) call comm_forces_dd1(deopb,cBond)
         if(use_strtor) call comm_forces_dd1(debt ,cBond)
         if(use_opdist) call comm_forces_dd1(deopd,cBond)
         if(use_angang) call comm_forces_dd1(deaa ,cBond)
         if(use_geom)   call comm_forces_dd1(deg  ,cBond)
         if(use_extra)  call comm_forces_dd1(dex  ,cBond)
      else
         !call commforce_one(de_buff0,commBonded)
      end if

      end if

      if (btest(rule,idSNBond)) then

      if (use_vdw) then
         call comm_forces_dd(dev ,cSNBond)
         call minmaxone1(mmx(16),mmx(nf+16),mmx(2*nf+16),dev
     &                  ,sze,'dev');
      end if
      if (use_charge) then
         call comm_forces_dd(dec ,cSNBond)
         call minmaxone1(mmx(17),mmx(nf+17),mmx(2*nf+17),dec
     &                  ,sze,'dec');
      end if
      if (use_mpole) then
         call comm_forces_dd(dem ,cSNBond)
         call minmaxone1(mmx(19),mmx(nf+19),mmx(2*nf+19),dem
     &                  ,sze,'dem');
      end if
      if (use_polar) then
         call comm_forces_dd(dep ,cSNBond)
         call minmaxone1(mmx(21),mmx(nf+21),mmx(2*nf+21),dep
     &                  ,sze,'dep');
      end if

      end if

      if (btest(rule,idNBond)) then

      if (use_vdw) then
         call comm_forces_dd(dev,cNBond)
         call minmaxone1(mmx(16),mmx(nf+16),mmx(2*nf+16),dev
     &                  ,sze,'dev')
      end if

      if(use_charge) then
         call comm_forces_dd(dec,cNBond)
         call minmaxone1(mmx(17),mmx(nf+17),mmx(2*nf+17),dec
     &                  ,sze,'dec');
         call commforcesrec1(dec,decrec)
         call minmaxone1(mmx(18),mmx(nf+18),mmx(2*nf+18),dec
     &                  ,sze,'decsum');
      end if

      if(use_mpole) then
         call comm_forces_dd(dem,cNBond)
         call minmaxone1(mmx(19),mmx(nf+19),mmx(2*nf+19),dem
     &                  ,sze,'dem');
         call commforcesrec1(dem,demrec)
         call minmaxone1(mmx(20),mmx(nf+20),mmx(2*nf+20),dem
     &                  ,sze,'demsum');
      end if

      if(use_polar) then
         call comm_forces_dd(dep,cNBond)
         call minmaxone1(mmx(21),mmx(nf+21),mmx(2*nf+21),dep
     &                  ,sze,'dep');
         call commforcesrec1(dep,deprec)
         call minmaxone1(mmx(22),mmx(nf+22),mmx(2*nf+22),dep
     &                  ,sze,'depsum');
      end if

      if (use_chgtrn) then
         call comm_forces_dd(dect,cNBond)
         call minmaxone1(mmx(28),mmx(nf+28),mmx(2*nf+28),dect
     &                  ,sze,'dect');
      end if

      end if

      if (btest(rule,idBond)) then

      if (.not.fdebs_l) then
      if(use_bond)
     &call minmaxone(mmx(01),mmx(nf+01),mmx(2*nf+01),deb
     &     ,sze,'deb');
      if(use_angle)
     &call minmaxone(mmx(02),mmx(nf+02),mmx(2*nf+02),dea
     &     ,sze,'dea');
      if(use_strbnd)
     &call minmaxone(mmx(03),mmx(nf+03),mmx(2*nf+03),deba
     &     ,sze,'deba');
      if(use_urey)
     &call minmaxone(mmx(04),mmx(nf+04),mmx(2*nf+04),deub
     &     ,sze,'deub');
      if(use_angang)
     &call minmaxone(mmx(05),mmx(nf+05),mmx(2*nf+05),deaa
     &     ,sze,'deaa');
      if(use_improp)
     &call minmaxone(mmx(06),mmx(nf+06),mmx(2*nf+06),deid
     &     ,sze,'deid');
      if(use_imptor)
     &call minmaxone(mmx(07),mmx(nf+07),mmx(2*nf+07),deit
     &     ,sze,'deit');
      if(use_tors)
     &call minmaxone(mmx(08),mmx(nf+08),mmx(2*nf+08),det
     &     ,sze,'det');
      if(use_pitors)
     &call minmaxone(mmx(09),mmx(nf+09),mmx(2*nf+09),dept
     &     ,sze,'dept');
      if(use_strtor)
     &call minmaxone(mmx(10),mmx(nf+10),mmx(2*nf+10),debt
     &     ,sze,'debt');
      if(use_tortor)
     &call minmaxone(mmx(11),mmx(nf+11),mmx(2*nf+11),dett
     &     ,sze,'dett');
      if(use_angtor)
     &call minmaxone(mmx(25),mmx(nf+25),mmx(2*nf+25),deat
     &     ,sze,'deat');
      if(use_opbend)
     &call minmaxone(mmx(12),mmx(nf+12),mmx(2*nf+12),deopb
     &     ,sze,'deopb');
      if(use_opdist)
     &call minmaxone(mmx(13),mmx(nf+13),mmx(2*nf+13),deopd
     &     ,sze,'deopd');
      if(use_geom)
     &call minmaxone(mmx(14),mmx(nf+14),mmx(2*nf+14),deg
     &     ,sze);
      if(use_extra)
     &call minmaxone(mmx(15),mmx(nf+15),mmx(2*nf+15),dex
     &     ,sze,'dex');
      if(use_mlpot.or.use_ml_embedding)
     &call minmaxone(mmx(28),mmx(nf+28),mmx(2*nf+28),dmlpot
     &     ,sze,'dmlpot');

      else
      call minmaxone1(mmx(01),mmx(nf+01),mmx(2*nf+01),de_tot
     &     ,sze,'de_tot');
      end if

      end if

      if (btest(rule,idSNBond)) then

      if (abort) then
         dint1 = minval(inte); dint2=maxval(inte)
         call searchpair(nshortvlst,shortvlst,maxvlst
     &                  ,dint1,dint2)
      end if

      else if (btest(rule,idNBond)) then

c     if (abort.and.vlst_enable) then
c        dint1 = minval(inte); dint2=maxval(inte)
c        call searchpair(nshortvlst,shortvlst,maxvlst,
c    &        dint1,dint2)
c     end if

      end if

c     if (allocated(deamdD)) then
c        call minmaxone(mmx(23),mmx(nf+23),mmx(2*nf+23),deamdD ,3*nloc);
c        call minmaxone(mmx(24),mmx(nf+24),mmx(2*nf+24),deW1aMD,3*nloc);
c     end if

      if (use_colvars.and.ncvatoms.gt.0) then
!$acc data copyin(decv,decv_tot)
      call minmaxone(mmx(26),mmx(nf+26),mmx(2*nf+26),decv_tot
     &              ,3*ncvatoms,'decolv');
      call minmaxone(mmx(27),mmx(nf+27),mmx(2*nf+27),decv
     &              ,3*ncvatoms,'decolv');
!$acc end data
      end if

      if (nproc.gt.1) then
      if (rank.eq.0) then
         call MPI_REDUCE(MPI_IN_PLACE,mmx,nf,MPI_RPREC,
     &                   MPI_MIN,0,MPI_COMM_WORLD,i)
         call MPI_REDUCE(MPI_IN_PLACE,mmx(nf+1),nf,MPI_RPREC,
     &                   MPI_MAX,0,MPI_COMM_WORLD,i)
         call MPI_REDUCE(MPI_IN_PLACE,mmx(2*nf+1),nf,MPI_RPREC,
     &                   MPI_SUM,0,MPI_COMM_WORLD,i)
      else
         call MPI_REDUCE(mmx,mmx,nf,MPI_RPREC,
     &                   MPI_MIN,0,MPI_COMM_WORLD,i)
         call MPI_REDUCE(mmx(nf+1),mmx(nf+1),nf,MPI_RPREC,
     &                   MPI_MAX,0,MPI_COMM_WORLD,i)
         call MPI_REDUCE(mmx(2*nf+1),mmx(2*nf+1),nf,MPI_RPREC,
     &                   MPI_SUM,0,MPI_COMM_WORLD,i)
      end if
      end if

      maxi=0
      ids = merge(2,1,fdebs_l)
      do i = ids,26; maxi = max(mmx(nf+i),maxi); end do
      if (focus_nbond) maxi = maxval(mmx(nf+16:nf+21))

      if ((maxi.gt.f_ulim).and.f_ulim.gt.0)
     &   then
         n_adjust = n_adjust + 1
 11      format('info_forces: above frc uplimit',F8.2,' cnt',I10)
         if(rank.eq.0) write(*,11) f_ulim,n_adjust

         if (maxi.gt.f_ulim+7) then
         read(arg(3),*) dt
         dt = 1d-3*dt
         new_restart = .true.
         f_mdsave    = .true.
         save_dcdio  = dcdio
         save_arc    = archive
         dcdio       = .false.
         archive     = .false.
         n_fwriten   = n_fwriten + 1
         call mdsave(step_c,dt,epot)
         new_restart = .false.
         f_mdsave    = .false.
         dcdio       = save_dcdio
         archive     = save_arc
         end if
      end if

      abortf = merge(.true.,.false.,maxi.gt.10000)
      if (abortf) then
 24   format(' info_forces ! Abnormal forces detected ',/)
         write(0,24)
         abort=.true.
      end if

      if (rank.eq.0) then
 30      format(a10,3F20.8)
 40      format(80('='))
         print 40
         if (fdebs_l) then
         if(mmx(2*nf+01)/=0.0) print 30,'de_tot =>',
     &      mmx(01),mmx(01+nf),mmx(01+nf*2)
         else
         if(mmx(2*nf+01)/=0.0) print 30,'deb    =>',
     &      mmx(01),mmx(01+nf),mmx(01+nf*2)
         end if
         if(mmx(2*nf+02)/=0.0) print 30,'dea    =>',
     &      mmx(02),mmx(02+nf),mmx(02+nf*2)
         if(mmx(2*nf+03)/=0.0) print 30,'deba   =>',
     &      mmx(03),mmx(03+nf),mmx(03+nf*2)
         if(mmx(2*nf+04)/=0.0) print 30,'deub   =>',
     &      mmx(04),mmx(04+nf),mmx(04+nf*2)
         if(mmx(2*nf+05)/=0.0) print 30,'deaa   =>',
     &      mmx(05),mmx(05+nf),mmx(05+nf*2)
         if(mmx(2*nf+06)/=0.0) print 30,'deid   =>',
     &      mmx(06),mmx(06+nf),mmx(06+nf*2)
         if(mmx(2*nf+07)/=0.0) print 30,'deit   =>',
     &      mmx(07),mmx(07+nf),mmx(07+nf*2)
         if(mmx(2*nf+08)/=0.0) print 30,'det    =>',
     &      mmx(08),mmx(08+nf),mmx(08+nf*2)
         if(mmx(2*nf+09)/=0.0) print 30,'dept   =>',
     &      mmx(09),mmx(09+nf),mmx(09+nf*2)
         if(mmx(2*nf+10)/=0.0) print 30,'debt   =>',
     &      mmx(10),mmx(10+nf),mmx(10+nf*2)
         if(mmx(2*nf+11)/=0.0) print 30,'dett   =>',
     &      mmx(11),mmx(11+nf),mmx(11+nf*2)
         if(mmx(2*nf+25)/=0.0) print 30,'deat   =>',
     &      mmx(25),mmx(25+nf),mmx(25+nf*2)
         if(mmx(2*nf+12)/=0.0) print 30,'deopb  =>',
     &      mmx(12),mmx(12+nf),mmx(12+nf*2)
         if(mmx(2*nf+13)/=0.0) print 30,'deopd  =>',
     &      mmx(13),mmx(13+nf),mmx(13+nf*2)
         if(mmx(2*nf+14)/=0.0) print 30,'deg    =>',
     &      mmx(14),mmx(14+nf),mmx(14+nf*2)
         if(mmx(2*nf+15)/=0.0) print 30,'dex    =>',
     &      mmx(15),mmx(15+nf),mmx(15+nf*2)
         if(mmx(2*nf+16)/=0.0) print 30,'dev    =>',
     &      mmx(16),mmx(16+nf),mmx(16+nf*2)
         if(mmx(2*nf+17)/=0.0) print 30,'dec    =>',
     &      mmx(17),mmx(17+nf),mmx(17+nf*2)
         if(mmx(2*nf+18)/=0.0) print 30,'decsum =>',
     &      mmx(18),mmx(18+nf),mmx(18+nf*2)
         if(mmx(2*nf+19)/=0.0) print 30,'dem    =>',
     &      mmx(19),mmx(19+nf),mmx(19+nf*2)
         if(mmx(2*nf+20)/=0.0) print 30,'demsum =>',
     &      mmx(20),mmx(20+nf),mmx(20+nf*2)
         if(mmx(2*nf+21)/=0.0) print 30,'dep    =>',
     &      mmx(21),mmx(21+nf),mmx(21+nf*2)
         if(mmx(2*nf+22)/=0.0) print 30,'depsum =>',
     &      mmx(22),mmx(22+nf),mmx(22+nf*2)
         if(mmx(2*nf+28)/=0.0) print 30,'dect   =>',
     &      mmx(28),mmx(28+nf),mmx(28+nf*2)
         if(mmx(2*nf+23)/=0.0) print 30,'deamdD =>',
     &      mmx(23),mmx(23+nf),mmx(23+nf*2)
         if(mmx(2*nf+24)/=0.0) print 30,'deW1aMD=>',
     &      mmx(24),mmx(24+nf),mmx(24+nf*2)
         if(mmx(2*nf+26)/=0.0) print 30,'declvi =>',
     &      mmx(26),mmx(26+nf),mmx(26+nf*2)
         if(mmx(2*nf+27)/=0.0) print 30,'declvo =>',
     &      mmx(27),mmx(27+nf),mmx(27+nf*2)
         if(mmx(2*nf+28)/=0.0) print 30,'dmlpot =>',
     &      mmx(28),mmx(28+nf),mmx(28+nf*2)
      end if

      call reset_forces_d
      if (btest(rule,idNBond)) then
         call zero_forces_rec
      end if

      doin = doin +1
      end subroutine

      subroutine minmaxone1( mi,ma,on,vector,sz,name )
      implicit none
      integer sz
      integer,parameter::lgli=10
      real(8) mi,ma,on
      mdyn_rtyp vector(*)
      character(*),optional,intent(in)::name
      integer i,j,i1,i2,iglob,cap,cap1,cap2
      integer gli(lgli,nproc)
      real(8) val

      abortall = .false.
      if (present(name)) then
         cap = 1
         gli = 0
         if (name.eq.'devi') then
!$acc parallel loop default(present) copy(abort,cap,gli)
            do i = 1, sz/3
               cap2  = 0
               iglob = glob(i)
!$acc loop seq
               do j  = 1,3
               val   = mdr2md(vector((i-1)*3+j))
               if (abs(val).gt.90.0_re_p) then
                  print*,j,iglob,rank,real(x(iglob),4),real(y(iglob),4)
     &                  ,real(z(iglob),4),real(val,4)
!$acc atomic write
                  abort=.true.
                  cap2 = cap2 + 1
               end if
               end do
               if (cap2.gt.0) then
!$acc atomic capture
               cap1 = cap
               cap  = cap + 1
!$acc end atomic
               if (cap1.le.lgli) gli(cap1,rank+1) = iglob
               end if
            end do
            do i = 1, nproc
               if (rank.eq.i-1) abortall=abort
               call MPI_BCAST(abortall,1,MPI_LOGICAL,i-1,COMM_TINKER,i1)
               if (abortall) then
                  abort = .true.
                  call MPI_AllGather(MPI_IN_PLACE,lgli,MPI_DATATYPE_NULL
     $                              ,gli,lgli,MPI_INT,COMM_TINKER,i1)
                  exit
               end if
            end do

            if (abort) then
            
            cap  = 0
            cap1 = 0
            do i1 = 1, nproc; do i2 = 1,lgli
               if ( gli(i2,i1).ne.0 ) then
                  if (cap.eq.0) then
                     cap = 1
                     inte(cap) = gli(i2,i1)
                  else if (cap.gt.0.and.abs(gli(i2,i1)-inte(1)).lt.5)
     &                 then
                     cap = cap + 1
                     inte(2) = gli(i2,i1)
                  else
                     cap = cap + 1
                  end if
               end if
            end do; end do
            
            if (cap.ne.2.and.rank.eq.0) then
               print*,' more than one interactions found '
     &               ,rank,cap
            do i = 1,nproc
               do j = 1,lgli
                  if (gli(j,i).ne.0) write(*,'(I10,$)') gli(j,i)
               end do
               print*
            end do
            !print*, 'interaction', inte,rank
            end if
            
            end if
         end if

      end if

!$acc parallel loop present(vector(1:sz))
      do i = 1, sz
         val = mdr2md(vector(i))
         mi = min( mi,val )
         ma = max( ma,val )
         on = on + abs(val)
      end do
      end subroutine

      subroutine minmaxone( mi,ma,on,vector,sz,name )
      implicit none
      integer sz
      real(8) mi,ma,on
      real(8) vector(*)
      character(*),optional,intent(in)::name
      integer i,j
      real(8) val

      abortall = .false.

!$acc parallel loop present(vector(1:sz))
      do i = 1, sz
         val = (vector(i))
         mi = min( mi,val )
         ma = max( ma,val )
         on = on + abs(val)
      end do
      end subroutine
      end submodule
