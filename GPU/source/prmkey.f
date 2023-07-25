c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine prmkey  --  interpret force field keywords  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "prmkey" parses a text string to extract keywords related to
c     force field potential energy functional forms and constants
c
c
#include "tinker_precision.h"
      subroutine prmkey (text)
      use angpot
      use bndpot
      use bound
      use chgpot
      use ctrpot
      use cutoff
      use divcon
      use fields
      use interfaces,only:sub_config,itrf_legacy
     &              , itrf_adapted
      use mplpot
      use neigh     ,only: buffMatb
      use polpot
      use potent
      use reppot
      use torpot
      use random_mod
      use precompute_pole,only: polar_precomp
      use urypot
      use USampling ,only: US_enable
      use utilgpu
      use uprior
      use vdwpot
      implicit none
      integer next,ivalue
      character*4 value
      character*25 keyword
      character*240 text
      character*240 record
      character*240 string
c
c
c     parse the line to extract any possible keyword
c
      record = text
      next = 1
      call upcase (record)
      call gettext (record,keyword,next)
      string = record(next:240)
c
c     select the individual force field potential terms
c
      if (keyword(1:9) .eq. 'BONDTERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_bond = .true.
         if (value .eq. 'NONE')  use_bond = .false.
      else if (keyword(1:10) .eq. 'ANGLETERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_angle = .true.
         if (value .eq. 'NONE')  use_angle = .false.
      else if (keyword(1:11) .eq. 'STRBNDTERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_strbnd = .true.
         if (value .eq. 'NONE')  use_strbnd = .false.
      else if (keyword(1:9) .eq. 'UREYTERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_urey = .true.
         if (value .eq. 'NONE')  use_urey = .false.
      else if (keyword(1:11) .eq. 'ANGANGTERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_angang = .true.
         if (value .eq. 'NONE')  use_angang = .false.
      else if (keyword(1:11) .eq. 'OPBENDTERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_opbend = .true.
         if (value .eq. 'NONE')  use_opbend = .false.
      else if (keyword(1:11) .eq. 'OPDISTTERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_opdist = .true.
         if (value .eq. 'NONE')  use_opdist = .false.
      else if (keyword(1:11) .eq. 'IMPROPTERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_improp = .true.
         if (value .eq. 'NONE')  use_improp = .false.
      else if (keyword(1:12) .eq. 'IMPTORSTERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_imptor = .true.
         if (value .eq. 'NONE')  use_imptor = .false.
      else if (keyword(1:12) .eq. 'TORSIONTERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_tors = .true.
         if (value .eq. 'NONE')  use_tors = .false.
      else if (keyword(1:11) .eq. 'PITORSTERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_pitors = .true.
         if (value .eq. 'NONE')  use_pitors = .false.
      else if (keyword(1:11) .eq. 'STRTORTERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_strtor = .true.
         if (value .eq. 'NONE')  use_strtor = .false.
      else if (keyword(1:11) .eq. 'TORTORTERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_tortor = .true.
         if (value .eq. 'NONE')  use_tortor = .false.
      else if (keyword(1:8) .eq. 'VDWTERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_vdw = .true.
         use_vlist = .true.
         if (value .eq. 'NONE')  then
           use_vdw = .false.
           use_vlist = .false.
         end if
      else if (keyword(1:8) .eq. 'DISPTERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_disp = .true.
         use_dlist = .true.
         if (value .eq. 'NONE')  then
           use_disp = .false.
           use_dlist = .false.
         end if
      else if (keyword(1:14) .eq. 'REPULSIONTERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_repuls = .true.
         if (value .eq. 'NONE')  use_repuls = .false.
      else if (keyword(1:15) .eq. 'DISPERSIONTERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_disp = .true.
         if (value .eq. 'NONE')  use_disp = .false.
      else if (keyword(1:11) .eq. 'CHGTRNTERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_chgtrn = .true.
         if (value .eq. 'NONE')  use_chgtrn = .false.
      else if (keyword(1:11) .eq. 'CHARGETERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_charge = .true.
         use_clist = .true.
         if (value .eq. 'NONE')  then
           use_charge = .false.
           use_clist = .false.
         end if
      else if (keyword(1:10) .eq. 'MPOLETERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_mpole = .true.
         use_mlist = .true.
         if (value .eq. 'NONE')  then
           use_mpole = .false.
           use_mlist = .false.
         end if
!$acc update device(use_mpole)
      else if (keyword(1:13) .eq. 'POLARIZETERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_polar = .true.
         if (value .eq. 'NONE')  use_polar = .false.
!$acc update device(use_polar)
      else if (keyword(1:12) .eq. 'SOLVATETERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_solv = .true.
         if (value .eq. 'NONE')  use_solv = .false.
      else if (keyword(1:13) .eq. 'RESTRAINTERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_geom = .true.
         if (value .eq. 'NONE')  use_geom = .false.
      else if (keyword(1:10) .eq. 'EXTRATERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_extra = .true.
         if (value .eq. 'NONE')  use_extra = .false.
      else if (keyword(1:5) .eq. 'EMTP ') then
         use_emtp = .true.
         if (value .eq. 'NONE')  use_emtp = .false.
      else if (keyword(1:8) .eq. 'POLYMER ') then
         use_polymer = .true.
      end if
c
c     select the name of the force field parameter set
c
      if (keyword(1:11) .eq. 'FORCEFIELD ') then
         call getword (record,forcefield,next)
c
c     set control parameters for bond stretching potentials
c
      else if (keyword(1:9) .eq. 'BONDTYPE ') then
         call getword (record,bndtyp,next)
         if (bndtyp.eq.'MORSE') then
            bndtyp_i = BND_MORSE
         else if (bndtyp.eq.'MORSE4') then
            bndtyp_i = BND_MORSE4
         else if (bndtyp.eq.'HARMONIC') then
            bndtyp_i = BND_HARMONIC
         else
            print*, 'WARNING;:prmkey.f !!!'
     &            , '  Unrecognize BONDTYPE in keyfile !!!', bndtyp
            bndtyp_i = BND_NO_TYPE
         end if
      else if (keyword(1:9) .eq. 'BONDUNIT ') then
         read (string,*,err=10,end=10)  bndunit
      else if (keyword(1:11) .eq. 'BOND-CUBIC ') then
         read (string,*,err=10,end=10)  cbnd
      else if (keyword(1:13) .eq. 'BOND-QUARTIC ') then
         read (string,*,err=10,end=10)  qbnd
c
c     set control parameters for bond angle bending potentials
c
      else if (keyword(1:10) .eq. 'ANGLEUNIT ') then
         read (string,*,err=10,end=10)  angunit
      else if (keyword(1:12) .eq. 'ANGLE-CUBIC ') then
         read (string,*,err=10,end=10)  cang
      else if (keyword(1:14) .eq. 'ANGLE-QUARTIC ') then
         read (string,*,err=10,end=10)  qang
      else if (keyword(1:13) .eq. 'ANGLE-PENTIC ') then
         read (string,*,err=10,end=10)  pang
      else if (keyword(1:13) .eq. 'ANGLE-SEXTIC ') then
         read (string,*,err=10,end=10)  sang
c
c     set control parameters for stretch-bend potentials
c
      else if (keyword(1:11) .eq. 'STRBNDUNIT ') then
         read (string,*,err=10,end=10)  stbnunit
c
c     set control parameters for Urey-Bradley potentials
c
      else if (keyword(1:9) .eq. 'UREYUNIT ') then
         read (string,*,err=10,end=10)  ureyunit
      else if (keyword(1:11) .eq. 'UREY-CUBIC ') then
         read (string,*,err=10,end=10)  cury
      else if (keyword(1:13) .eq. 'UREY-QUARTIC ') then
         read (string,*,err=10,end=10)  qury
c
c     set control parameters for out-of-plane bend potentials
c
      else if (keyword(1:11) .eq. 'OPBENDTYPE ') then
         call getword (record,opbtyp,next)
         if( opbtyp.eq.'ALLINGER') opbtypInt = OPB_ALLINGER
      else if (keyword(1:11) .eq. 'OPBENDUNIT ') then
         read (string,*,err=10,end=10)  opbunit
      else if (keyword(1:13) .eq. 'OPBEND-CUBIC ') then
         read (string,*,err=10,end=10)  copb
      else if (keyword(1:15) .eq. 'OPBEND-QUARTIC ') then
         read (string,*,err=10,end=10)  qopb
      else if (keyword(1:14) .eq. 'OPBEND-PENTIC ') then
         read (string,*,err=10,end=10)  popb
      else if (keyword(1:14) .eq. 'OPBEND-SEXTIC ') then
         read (string,*,err=10,end=10)  sopb
c
c     set control parameters for out-of-plane distance potentials
c
      else if (keyword(1:11) .eq. 'OPDISTUNIT ') then
         read (string,*,err=10,end=10)  opdunit
      else if (keyword(1:13) .eq. 'OPDIST-CUBIC ') then
         read (string,*,err=10,end=10)  copd
      else if (keyword(1:15) .eq. 'OPDIST-QUARTIC ') then
         read (string,*,err=10,end=10)  qopd
      else if (keyword(1:14) .eq. 'OPDIST-PENTIC ') then
         read (string,*,err=10,end=10)  popd
      else if (keyword(1:14) .eq. 'OPDIST-SEXTIC ') then
         read (string,*,err=10,end=10)  sopd
c
c     set control parameters for other local geometry potentials
c
      else if (keyword(1:11) .eq. 'ANGANGUNIT ') then
         read (string,*,err=10,end=10)  aaunit
      else if (keyword(1:11) .eq. 'IMPROPUNIT ') then
         read (string,*,err=10,end=10)  idihunit
      else if (keyword(1:11) .eq. 'IMPTORUNIT ') then
         read (string,*,err=10,end=10)  itorunit
      else if (keyword(1:12) .eq. 'TORSIONUNIT ') then
         read (string,*,err=10,end=10)  torsunit
      else if (keyword(1:11) .eq. 'PITORSUNIT ') then
         read (string,*,err=10,end=10)  ptorunit
      else if (keyword(1:11) .eq. 'STRTORUNIT ') then
         read (string,*,err=10,end=10)  storunit
      else if (keyword(1:11) .eq. 'ANGTORUNIT ') then
         read (string,*,err=10,end=10)  atorunit
      else if (keyword(1:11) .eq. 'TORTORUNIT ') then
         read (string,*,err=10,end=10)  ttorunit
c
c     set control parameters for van der Waals potentials
c
      else if (keyword(1:9) .eq. 'VDWINDEX ') then
         call getword (record,vdwindex,next)
      else if (keyword(1:8) .eq. 'VDWTYPE ') then
         call getword (record,vdwtyp,next)
      else if (keyword(1:11) .eq. 'RADIUSTYPE ') then
         call getword (record,radtyp,next)
      else if (keyword(1:11) .eq. 'RADIUSSIZE ') then
         call getword (record,radsiz,next)
      else if (keyword(1:11) .eq. 'RADIUSRULE ') then
         call getword (record,radrule,next)
         call set_vdwradepsrule(radrule,radrule_i)
      else if (keyword(1:12) .eq. 'EPSILONRULE ') then
         call getword (record,epsrule,next)
         call set_vdwradepsrule(epsrule,epsrule_i)
      else if (keyword(1:14) .eq. 'GAUSSTYPE ') then
         call getword (record,gausstyp,next)
      else if (keyword(1:10) .eq. 'A-EXPTERM ') then
         read (string,*,err=10,end=10)  abuck
      else if (keyword(1:10) .eq. 'B-EXPTERM ') then
         read (string,*,err=10,end=10)  bbuck
      else if (keyword(1:10) .eq. 'C-EXPTERM ') then
         read (string,*,err=10,end=10)  cbuck
      else if (keyword(1:14) .eq. 'GAMMA-HALGREN ') then
         read (string,*,err=10,end=10)  ghal
      else if (keyword(1:14) .eq. 'DELTA-HALGREN ') then
         read (string,*,err=10,end=10)  dhal
      else if (keyword(1:13) .eq. 'VDW-12-SCALE ') then
         read (string,*,err=10,end=10)  v2scale
         if (v2scale .gt. 1.0_ti_p)  v2scale = 1.0_ti_p / v2scale
!$acc update device(v2scale)
      else if (keyword(1:13) .eq. 'VDW-13-SCALE ') then
         read (string,*,err=10,end=10)  v3scale
         if (v3scale .gt. 1.0_ti_p)  v3scale = 1.0_ti_p / v3scale
!$acc update device(v3scale)
      else if (keyword(1:13) .eq. 'VDW-14-SCALE ') then
         read (string,*,err=10,end=10)  v4scale
         if (v4scale .gt. 1.0_ti_p)  v4scale = 1.0_ti_p / v4scale
!$acc update device(v4scale)
      else if (keyword(1:13) .eq. 'VDW-15-SCALE ') then
         read (string,*,err=10,end=10)  v5scale
         if (v5scale .gt. 1.0_ti_p)  v5scale = 1.0_ti_p / v5scale
!$acc update device(v5scale)
      else if (keyword(1:15) .eq. 'VDW-CORRECTION ') then
         use_vcorr = .true.
c
c     set control parameters for Pauli repulsion potential
c
      else if (keyword(1:13) .eq. 'REP-12-SCALE ') then
         read (string,*,err=10,end=10)  r2scale
         if (r2scale .gt. 1.0d0)  r2scale = 1.0_ti_p / r2scale
      else if (keyword(1:13) .eq. 'REP-13-SCALE ') then
         read (string,*,err=10,end=10)  r3scale
         if (r3scale .gt. 1.0d0)  r3scale = 1.0_ti_p / r3scale
      else if (keyword(1:13) .eq. 'REP-14-SCALE ') then
         read (string,*,err=10,end=10)  r4scale
         if (r4scale .gt. 1.0d0)  r4scale = 1.0_ti_p / r4scale
      else if (keyword(1:13) .eq. 'REP-15-SCALE ') then
         read (string,*,err=10,end=10)  r5scale
         if (r5scale .gt. 1.0d0)  r5scale = 1.0_ti_p / r5scale
c
c     set control parameters for charge-charge potentials
c
      else if (keyword(1:9) .eq. 'ELECTRIC ') then
         read (string,*,err=10,end=10)  electric
      else if (keyword(1:11) .eq. 'DIELECTRIC ') then
         read (string,*,err=10,end=10)  dielec
      else if (keyword(1:11) .eq. 'CHG-BUFFER ') then
         read (string,*,err=10,end=10)  ebuffer
      else if (keyword(1:13) .eq. 'CHG-12-SCALE ') then
         read (string,*,err=10,end=10)  c2scale
         if (c2scale .gt. 1.0_ti_p)  c2scale = 1.0_ti_p / c2scale
      else if (keyword(1:13) .eq. 'CHG-13-SCALE ') then
         read (string,*,err=10,end=10)  c3scale
         if (c3scale .gt. 1.0_ti_p)  c3scale = 1.0_ti_p / c3scale
      else if (keyword(1:13) .eq. 'CHG-14-SCALE ') then
         read (string,*,err=10,end=10)  c4scale
         if (c4scale .gt. 1.0_ti_p)  c4scale = 1.0_ti_p / c4scale
      else if (keyword(1:13) .eq. 'CHG-15-SCALE ') then
         read (string,*,err=10,end=10)  c5scale
         if (c5scale .gt. 1.0_ti_p)  c5scale = 1.0_ti_p / c5scale
      else if (keyword(1:16) .eq. 'NEIGHBOR-GROUPS ') then
         neutnbr = .true.
      else if (keyword(1:15) .eq. 'NEUTRAL-GROUPS ') then
         neutcut = .true.
c
c     set control parameters for atomic multipole potentials
c
      else if (keyword(1:12) .eq. 'PENETRATION ') then
         call getword (record,pentyp,next)
         if (pentyp.eq.'GORDON2') then
            pentyp_i = PT_GORDON2
         else if (pentyp.eq.'GORDON1') then
            pentyp_i = PT_GORDON1
         else
            pentyp_i = PT_NONE
         end if
      else if (keyword(1:15) .eq. 'MPOLE-12-SCALE ') then
         read (string,*,err=10,end=10)  m2scale
         if (m2scale .gt. 1.0_ti_p)  m2scale = 1.0_ti_p / m2scale
!$acc update device(m2scale)
      else if (keyword(1:15) .eq. 'MPOLE-13-SCALE ') then
         read (string,*,err=10,end=10)  m3scale
         if (m3scale .gt. 1.0_ti_p)  m3scale = 1.0_ti_p / m3scale
!$acc update device(m3scale)
      else if (keyword(1:15) .eq. 'MPOLE-14-SCALE ') then
         read (string,*,err=10,end=10)  m4scale
         if (m4scale .gt. 1.0_ti_p)  m4scale = 1.0_ti_p / m4scale
!$acc update device(m4scale)
      else if (keyword(1:15) .eq. 'MPOLE-15-SCALE ') then
         read (string,*,err=10,end=10)  m5scale
         if (m5scale .gt. 1.0_ti_p)  m5scale = 1.0_ti_p / m5scale
!$acc update device(m5scale)
c
c     set control parameters for polarization potentials
c
      else if (keyword(1:13) .eq. 'POLARIZATION ') then
         call getword (record,poltyp,next)
      else if (keyword(1:10) .eq. 'POLAR-FF ') then
         read (string,*,err=10,end=10)  polff
      else if (keyword(1:10) .eq. 'POLAR-PRT ') then
         read (string,*,err=10,end=10)  polprt
      else if (keyword(1:10) .eq. 'POLAR-ALG ') then
         read (string,*,err=10,end=10)  polalg
         polalgshort = polalg
      else if (keyword(1:15) .eq. 'POLAR-ALGSHORT ') then
         read (string,*,err=10,end=10)  polalgshort
      else if (keyword(1:9) .eq. 'TCGORDER ') then
         read (string,*,err=10,end=10)  tcgorder
         tcgordershort = tcgorder
      else if (keyword(1:14) .eq. 'TCGORDERSHORT ') then
         read (string,*,err=10,end=10)  tcgordershort
      else if (keyword(1:8) .eq. 'TCGPREC ') then
         read (string,*,err=10,end=10)  value
         if (value.eq.'1') tcgprec = .true.
         if (value.eq.'0') tcgprec = .false.
         tcgprecshort = tcgprec
      else if (keyword(1:13) .eq. 'TCGPRECSHORT ') then
         read (string,*,err=10,end=10)  value
         if (value.eq.'1') tcgprecshort = .true.
         if (value.eq.'0') tcgprecshort = .false.
      else if (keyword(1:9) .eq. 'TCGGUESS ') then
         read (string,*,err=10,end=10)  value
         if (value.eq.'1') tcgguess = .true.
         if (value.eq.'0') tcgguess = .false.
         tcgguessshort = tcgguess
      else if (keyword(1:14) .eq. 'TCGGUESSSHORT ') then
         read (string,*,err=10,end=10)  value
         if (value.eq.'1') tcgguessshort = .true.
         if (value.eq.'0') tcgguessshort = .false.
      else if (keyword(1:8) .eq. 'TCGPEEK ') then
         read (string,*,err=10,end=10) value 
         if (value.eq.'1') tcgpeek = .true.
         if (value.eq.'0') tcgpeek = .false.
         tcgpeekshort = tcgpeek
      else if (keyword(1:13) .eq. 'TCGPEEKSHORT ') then
         read (string,*,err=10,end=10) value 
         if (value.eq.'0') tcgpeekshort = .false.
         if (value.eq.'1') tcgpeekshort = .true.
      else if (keyword(1:15) .eq. 'TCGOMEGAFITFREQ') then
         read (string,*,err=10,end=10)  omegafitfreq
      else if (keyword(1:12) .eq. 'TCGOMEGAFIT ') then
        tcgomegafit = .true.
      else if (keyword(1:8) .eq. 'TCGOMEGA') then
         read (string,*,err=10,end=10)  tcgomega
         tcgomegashort = tcgomega
      else if (keyword(1:13) .eq. 'TCGOMEGASHORT') then
         read (string,*,err=10,end=10)  tcgomegashort
      else if (keyword(1:12) .eq. 'CLUSTER-ALG ') then
         read (string,*,err=10,end=10)  clsttype
      else if (keyword(1:11) .eq. 'BLOCK-SIZE ') then
         read (string,*,err=10,end=10)  natprblk
      else if (keyword(1:11) .eq. 'PRECOMPUTE ') then
         read (string,*,err=10,end=10)  precomp
      else if (keyword(1:12) .eq. 'NO-COM-DIIS ') then
         read (string,*,err=10,end=10)  nocomdiis
      else if (keyword(1:12) .eq. 'POLAR-GUESS ') then
         read (string,*,err=10,end=10)  polgsf
      else if (keyword(1:11) .eq. 'POLAR-ITER ') then
         read (string,*,err=10,end=10)  politer
      else if (keyword(1:10) .eq. 'POLAR-EPS ') then
         read (string,*,err=10,end=10)  poleps
      else if (keyword(1:11) .eq. 'D-EQUALS-P ') then
         dpequal = .true.
      else if (keyword(1:15) .eq. 'POLAR-12-SCALE ') then
         read (string,*,err=10,end=10)  p2scale
         if (p2scale .gt. 1.0_ti_p)  p2scale = 1.0_ti_p / p2scale
!$acc update device(p2scale)
      else if (keyword(1:15) .eq. 'POLAR-13-SCALE ') then
         read (string,*,err=10,end=10)  p3scale
         if (p3scale .gt. 1.0_ti_p)  p3scale = 1.0_ti_p / p3scale
!$acc update device(p3scale)
      else if (keyword(1:15) .eq. 'POLAR-14-SCALE ') then
         read (string,*,err=10,end=10)  p4scale
         if (p4scale .gt. 1.0_ti_p)  p4scale = 1.0_ti_p / p4scale
!$acc update device(p4scale)
      else if (keyword(1:15) .eq. 'POLAR-15-SCALE ') then
         read (string,*,err=10,end=10)  p5scale
         if (p5scale .gt. 1.0_ti_p)  p5scale = 1.0_ti_p / p5scale
!$acc update device(p5scale)
      else if (keyword(1:15) .eq. 'POLAR-12-INTRA ') then
         read (string,*,err=10,end=10)  p2iscale
         if (p2iscale .gt. 1.0_ti_p)  p2iscale = 1.0_ti_p / p2iscale
      else if (keyword(1:15) .eq. 'POLAR-13-INTRA ') then
         read (string,*,err=10,end=10)  p3iscale
         if (p3iscale .gt. 1.0_ti_p)  p3iscale = 1.0_ti_p / p3iscale
      else if (keyword(1:15) .eq. 'POLAR-14-INTRA ') then
         read (string,*,err=10,end=10)  p4iscale
         if (p4iscale .gt. 1.0_ti_p)  p4iscale = 1.0_ti_p / p4iscale
            p41scale = p4iscale
      else if (keyword(1:15) .eq. 'POLAR-15-INTRA ') then
         read (string,*,err=10,end=10)  p5iscale
         if (p5iscale .gt. 1.0_ti_p)  p5iscale = 1.0_ti_p / p5iscale
      else if (keyword(1:16) .eq. 'DIRECT-11-SCALE ') then
         read (string,*,err=10,end=10)  d1scale
         if (d1scale .gt. 1.0_ti_p)  d1scale = 1.0_ti_p / d1scale
!$acc update device(d1scale)
      else if (keyword(1:16) .eq. 'DIRECT-12-SCALE ') then
         read (string,*,err=10,end=10)  d2scale
         if (d2scale .gt. 1.0_ti_p)  d2scale = 1.0_ti_p / d2scale
!$acc update device(d2scale)
      else if (keyword(1:16) .eq. 'DIRECT-13-SCALE ') then
         read (string,*,err=10,end=10)  d3scale
         if (d3scale .gt. 1.0_ti_p)  d3scale = 1.0_ti_p / d3scale
!$acc update device(d3scale)
      else if (keyword(1:16) .eq. 'DIRECT-14-SCALE ') then
         read (string,*,err=10,end=10)  d4scale
         if (d4scale .gt. 1.0_ti_p)  d4scale = 1.0_ti_p / d4scale
!$acc update device(d4scale)
      else if (keyword(1:16) .eq. 'MUTUAL-11-SCALE ') then
         read (string,*,err=10,end=10)  u1scale
         if (u1scale .gt. 1.0_ti_p)  u1scale = 1.0_ti_p / u1scale
!$acc update device(u1scale)
      else if (keyword(1:16) .eq. 'MUTUAL-12-SCALE ') then
         read (string,*,err=10,end=10)  u2scale
         if (u2scale .gt. 1.0_ti_p)  u2scale = 1.0_ti_p / u2scale
!$acc update device(u2scale)
      else if (keyword(1:16) .eq. 'MUTUAL-13-SCALE ') then
         read (string,*,err=10,end=10)  u3scale
         if (u3scale .gt. 1.0_ti_p)  u3scale = 1.0_ti_p / u3scale
!$acc update device(u3scale)
      else if (keyword(1:16) .eq. 'MUTUAL-14-SCALE ') then
         read (string,*,err=10,end=10)  u4scale
         if (u4scale .gt. 1.0_ti_p)  u4scale = 1.0_ti_p / u4scale
!$acc update device(u4scale)
      else if (keyword(1:16) .eq. 'INDUCE-12-SCALE ') then
         read (string,*,err=10,end=10)  w2scale
         if (w2scale .gt. 1.0_ti_p)  w2scale = 1.0_ti_p / w2scale
      else if (keyword(1:16) .eq. 'INDUCE-13-SCALE ') then
         read (string,*,err=10,end=10)  w3scale
         if (w3scale .gt. 1.0_ti_p)  w3scale = 1.0_ti_p / w3scale
      else if (keyword(1:16) .eq. 'INDUCE-14-SCALE ') then
         read (string,*,err=10,end=10)  w4scale
         if (w4scale .gt. 1.0_ti_p)  w4scale = 1.0_ti_p / w4scale
      else if (keyword(1:16) .eq. 'INDUCE-15-SCALE ') then
         read (string,*,err=10,end=10)  w5scale
         if (w5scale .gt. 1.0_ti_p)  w5scale = 1.0_ti_p / w5scale
c
c     set control parameters for charge transfer potentials
c
      else if (keyword(1:15) .eq. 'CHARGETRANSFER ') then
         call getword (record,ctrntyp,next)
         ctrntyp_ID=merge(CHGT_SEPARATE,CHGT_COMBINED
     &                  ,(ctrntyp.eq.'SEPARATE'))
      else if (keyword(1:10) .eq. 'LAMBDADYN ') then
         use_lambdadyn = .true.
      else if (keyword(1:5) .eq. 'OSRW ') then
         use_lambdadyn = .true.
         use_OSRW      = .true.

      else if (keyword(1:14) .eq. 'POLAR-PREDICT ') then
         call getword (record,polpred,next)
         call upcase (polpred)
         if (polpred .eq. 'NONE') use_pred = .false.
         if (polpred .eq. 'ENAB') use_pred = .true.
      !GPU Tools
      else if (keyword(1:10) .eq. 'GPU-GANGS ')       then
         read (string,*,err=10,end=10)   gpu_gangs
!$acc update device(gpu_gangs)
         mod_gangs = .false.
      else if (keyword(1:10) .eq. 'MOD-GANGS ')       then
         read (string,*,err=10,end=10)   mod_gangs
      else if (keyword(1:11) .eq. 'NDIR-BLOCK ')       then
         read (string,*,err=10,end=10)   Ndir_block
         Ndir_block = max(1,Ndir_block)
      else if (keyword(1:16) .eq. 'NDIR-ASYNC-BLOCK ')       then
         read (string,*,err=10,end=10)   Ndir_async_block
         Ndir_async_block = max(0,min(3*Ndir_async_block/4,
     &                                  Ndir_async_block))
      else if (keyword(1:12) .eq. 'GPU-WORKERS ')     then
         read (string,*,err=10,end=10)   gpu_workers
!$acc update device(gpu_workers)
      else if (keyword(1:11) .eq. 'GPU-VECTOR ')      then
         read (string,*,err=10,end=10)   gpu_vector
!$acc update device(gpu_vector)
#ifdef _OPENACC
      else if (keyword(1:21) .eq. 'HOST-RANDOM-PLATFORM ') then
         host_rand_platform = .true.
#endif
      else if (keyword(1:18) .eq. 'POLAR-PRECOMP ') then
         read (string,*,err=10,end=10) ivalue
         polar_precomp = btest(ivalue,0)
      else if (keyword(1:9) .eq. 'RUN-MODE ') then
         if      (index(string,'LEGACY').gt.0) then
            sub_config = itrf_legacy
         else if (index(string,'ADAPTED').gt.0) then
            sub_config = itrf_adapted
         end if
      else if (keyword(1:15) .eq. 'CONFIG-ROUTINE ') then
         read (string,'(Z32)',err=10,end=10) sub_config
      else if (keyword(1:11) .eq. 'U-SAMPLING ') then
         US_enable = .true.
      else if (keyword(1:15) .eq. 'MATLIST-BUFFER ') then
         read(string,*,err=10,end=10) buffMatb
      end if
c
c     jump directly to the end if any error was detected
c
   10 continue
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine potoff  --  turn off all potential functions  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "potoff" clears the forcefield definition by turning off
c     the use of each of the potential energy functions
c
c
      subroutine potoff
      use potent
      implicit none
c
c
c     turn off the use of each of the potential energy functions
c
      use_bond = .false.
      use_angle = .false.
      use_strbnd = .false.
      use_urey = .false.
      use_angang = .false.
      use_opbend = .false.
      use_opdist = .false.
      use_improp = .false.
      use_imptor = .false.
      use_tors = .false.
      use_pitors = .false.
      use_strtor = .false.
      use_tortor = .false.
      use_vdw = .false.
      use_charge = .false.
      use_mpole = .false.
      use_polar = .false.
      use_solv = .false.
      use_geom = .false.
      use_extra = .false.
      use_emtp = .false.
      use_disp = .false.
      use_repuls = .false.
      use_chgtrn = .false.
!$acc update device (use_mpole,use_polar)
      return
      end
