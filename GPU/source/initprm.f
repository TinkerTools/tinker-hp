c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine initprm  --  initialize force field parameters  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "initprm" completely initializes a force field by setting all
c     parameters to zero and using defaults for control values
c
c
#include "tinker_precision.h"
      subroutine initprm
      use sizes
      use angpot
      use bndpot
      use atoms
      use chgpot
      use divcon
      use fields
      use kanang
      use kangs
      use katoms
      use kbonds
      use kchrge
      use khbond
      use kiprop
      use kitors
      use kmulti
      use kopbnd
      use kopdst
      use kitors
      use kmulti
      use kpitor
      use kpolr
      use kstbnd
      use ksttor
      use ktorsn
      use ktrtor
      use kurybr
      use kvdws
      use kvdwpr
      use math
      use merck
      use mplpot
      use polpot
      use torpot
      use units
      use urypot
      use vdwpot
      implicit none
      integer i,j,k
      character*3 blank3
      character*8 blank8
      character*12 blank12
      character*16 blank16
      character*20 blank20
      character*24 blank24
c
c
c     define blank character strings of various lengths
c
      blank3 = '   '
      blank8 = '        '
      blank12 = '            '
      blank16 = '                '
      blank20 = '                    '
      blank24 = '                        '
c
c     initialize strings of parameter atom types and classes
c
      do i = 1, maxnvp
         kvpr(i) = blank8
      end do
      do i = 1, maxnhb
         khb(i) = blank8
      end do
      do i = 1, maxnb
         kb(i) = blank8
      end do
      do i = 1, maxnb5
         kb5(i) = blank8
      end do
      do i = 1, maxnb4
         kb4(i) = blank8
      end do
      do i = 1, maxnb3
         kb3(i) = blank8
      end do
      do i = 1, maxnel
         kel(i) = blank12
      end do
      do i = 1, maxna
         ka(i) = blank12
      end do
      do i = 1, maxna5
         ka5(i) = blank12
      end do
      do i = 1, maxna4
         ka4(i) = blank12
      end do
      do i = 1, maxna3
         ka3(i) = blank12
      end do
      do i = 1, maxnaf
         kaf(i) = blank12
      end do
      do i = 1, maxnsb
         ksb(i) = -1
      end do
      do i = 1, maxnu
         ku(i)  = -1
      end do
      do i = 1, maxnopb
         kopb(i) = -1
      end do
      do i = 1, maxnopd
         kopd(i) = blank16
      end do
      do i = 1, maxndi
         kdi(i) = -1
      end do
      do i = 1, maxnti
         kti(i) = -1
      end do
      do i = 1, maxnt
         kt(i) = blank16
      end do
      do i = 1, maxnt5
         kt5(i) = blank16
      end do
      do i = 1, maxnt4
         kt4(i) = blank16
      end do
      do i = 1, maxnpt
         kpt(i) = -1
      end do
      do i = 1, maxnbt
         kbt(i) = -1
      end do
      do i = 1, maxntt
         ktt(i) = -1
      end do
c      do i = 1, maxnd
c         kd(i) = blank8
c      end do
c      do i = 1, maxnd5
c         kd5(i) = blank8
c      end do
c      do i = 1, maxnd4
c         kd4(i) = blank8
c      end do
c      do i = 1, maxnd3
c         kd3(i) = blank8
c      end do
      do i = 1, maxnmp
         kmp(i) = blank12
      end do
c      do i = 1, maxnpi
c         kpi(i) = blank8
c      end do
c      do i = 1, maxnpi5
c         kpi5(i) = blank8
c      end do
c      do i = 1, maxnpi4
c         kpi4(i) = blank8
c      end do
c
c     initialize values of some force field parameters
c
      forcefield = blank20
      do i = 1, maxtyp
         symbol(i) = blank3
         atmcls(i) = 0
         atmnum(i) = 0
         weight(i) = 0.0_re_p
         ligand(i) = 0
         describe(i) = blank24
         rad(i)    = 0.0_ti_p
         eps(i)    = 0.0_ti_p
         rad4(i)   = 0.0_ti_p
         eps4(i)   = 0.0_ti_p
         reduct(i) = 0.0_ti_p
         chg(i)    = 0.0_ti_p
         polr(i)   = 0.0_ti_p
         athl(i)   = 0.0_ti_p
         do j = 1, maxvalue
            pgrp(j,i) = 0
         end do
         sibfacp(1,i) = 0.0_ti_p
         sibfacp(2,i) = 0.0_ti_p
         sibfacp(3,i) = 0.0_ti_p
      end do
      do i = 1, maxclass
         do j = 1, 2
            stbn(j,i) = 0.0_ti_p
         end do
         do j = 1, 3
            anan(j,i) = 0.0_ti_p
         end do
c         electron(i) = 0.0_ti_p
c         ionize(i)   = 0.0_ti_p
c         repulse(i)  = 0.0_ti_p
      end do
      do i = 1, maxbio
         biotyp(i) = 0
      end do
c
c     initialize values of some MMFF-specific parameters
c
      do i = 1, 100
         do j = 1, 100
            mmff_kb(j,i)  = 1000.0_ti_p
            mmff_kb1(j,i) = 1000.0_ti_p
            mmff_b0(j,i)  = 1000.0_ti_p
            mmff_b1(j,i)  = 1000.0_ti_p
            bci(j,i)   = 1000.0_ti_p
            bci_1(j,i) = 1000.0_ti_p
            do k = 1, 100
               stbn_abc  (k,j,i) = 1000.0_ti_p
               stbn_cba  (k,j,i) = 1000.0_ti_p
               stbn_abc1 (k,j,i) = 1000.0_ti_p
               stbn_cba1 (k,j,i) = 1000.0_ti_p
               stbn_abc2 (k,j,i) = 1000.0_ti_p
               stbn_cba2 (k,j,i) = 1000.0_ti_p
               stbn_abc3 (k,j,i) = 1000.0_ti_p
               stbn_cba3 (k,j,i) = 1000.0_ti_p
               stbn_abc4 (k,j,i) = 1000.0_ti_p
               stbn_cba4 (k,j,i) = 1000.0_ti_p
               stbn_abc5 (k,j,i) = 1000.0_ti_p
               stbn_cba5 (k,j,i) = 1000.0_ti_p
               stbn_abc6 (k,j,i) = 1000.0_ti_p
               stbn_cba6 (k,j,i) = 1000.0_ti_p
               stbn_abc7 (k,j,i) = 1000.0_ti_p
               stbn_cba7 (k,j,i) = 1000.0_ti_p
               stbn_abc8 (k,j,i) = 1000.0_ti_p
               stbn_cba8 (k,j,i) = 1000.0_ti_p
               stbn_abc9 (k,j,i) = 1000.0_ti_p
               stbn_cba9 (k,j,i) = 1000.0_ti_p
               stbn_abc10(k,j,i) = 1000.0_ti_p
               stbn_cba10(k,j,i) = 1000.0_ti_p
               stbn_abc11(k,j,i) = 1000.0_ti_p
               stbn_cba11(k,j,i) = 1000.0_ti_p
            end do
         end do
      end do
      do i = 0, 100
         do j = 1, 100
            do k = 0, 100
               mmff_ka  (k,j,i) = 1000.0_ti_p
               mmff_ka1 (k,j,i) = 1000.0_ti_p
               mmff_ka2 (k,j,i) = 1000.0_ti_p
               mmff_ka3 (k,j,i) = 1000.0_ti_p
               mmff_ka4 (k,j,i) = 1000.0_ti_p
               mmff_ka5 (k,j,i) = 1000.0_ti_p
               mmff_ka6 (k,j,i) = 1000.0_ti_p
               mmff_ka7 (k,j,i) = 1000.0_ti_p
               mmff_ka8 (k,j,i) = 1000.0_ti_p
               mmff_ang0(k,j,i) = 1000.0_ti_p
               mmff_ang1(k,j,i) = 1000.0_ti_p
               mmff_ang2(k,j,i) = 1000.0_ti_p
               mmff_ang3(k,j,i) = 1000.0_ti_p
               mmff_ang4(k,j,i) = 1000.0_ti_p
               mmff_ang5(k,j,i) = 1000.0_ti_p
               mmff_ang6(k,j,i) = 1000.0_ti_p
               mmff_ang7(k,j,i) = 1000.0_ti_p
               mmff_ang8(k,j,i) = 1000.0_ti_p
            end do
         end do
      end do
      do i = 1, maxnt
         kt  (i) = blank16
         kt_1(i) = blank16
         kt_2(i) = blank16
         t1  (1,i) = 1000.0_ti_p
         t1  (2,i) = 1000.0_ti_p
         t2  (1,i) = 1000.0_ti_p
         t2  (2,i) = 1000.0_ti_p
         t3  (1,i) = 1000.0_ti_p
         t3  (2,i) = 1000.0_ti_p
         t1_1(1,i) = 1000.0_ti_p
         t1_1(2,i) = 1000.0_ti_p
         t2_1(1,i) = 1000.0_ti_p
         t2_1(2,i) = 1000.0_ti_p
         t3_1(1,i) = 1000.0_ti_p
         t3_1(2,i) = 1000.0_ti_p
         t1_2(1,i) = 1000.0_ti_p
         t1_2(2,i) = 1000.0_ti_p
         t2_2(1,i) = 1000.0_ti_p
         t2_2(2,i) = 1000.0_ti_p
         t3_2(1,i) = 1000.0_ti_p
         t3_2(2,i) = 1000.0_ti_p
      end do
      do i = 1, 5
         do j = 1, 500
            eqclass(j,i) = 1000
         end do
      end do
      do i = 1, 6
         do j = 1, maxtyp
            mmffarom(j,i) = 0
            mmffaromc(j,i) = 0
            mmffaroma(j,i) = 0
         end do
      end do
c
c     set default control parameters for local geometry terms
c
      bndtyp = 'HARMONIC'
      bndtyp_i = BND_HARMONIC
      bndunit = 1.0_ti_p
      cbnd = 0.0_ti_p
      qbnd = 0.0_ti_p
      angunit = 1.0_ti_p / radian**2
      cang = 0.0_ti_p
      qang = 0.0_ti_p
      pang = 0.0_ti_p
      sang = 0.0_ti_p
      stbnunit = 1.0_ti_p / radian
      ureyunit = 1.0_ti_p
      cury = 0.0_ti_p
      qury = 0.0_ti_p
      aaunit = 1.0_ti_p / radian**2
      opbtyp = 'W-D-C'
      opbunit = 1.0_ti_p / radian**2
      copb = 0.0_ti_p
      qopb = 0.0_ti_p
      popb = 0.0_ti_p
      sopb = 0.0_ti_p
      opdunit = 1.0_ti_p
      copd = 0.0_ti_p
      qopd = 0.0_ti_p
      popd = 0.0_ti_p
      sopd = 0.0_ti_p
      idihunit = 1.0_ti_p
      itorunit = 1.0_ti_p
      torsunit = 1.0_ti_p
      ptorunit = 1.0_ti_p
      storunit = 1.0_ti_p
      ttorunit = 1.0_ti_p
c
c     set default control parameters for van der Waals terms
c
      vdwindex = 'CLASS'
      vdwtyp = 'LENNARD-JONES'
      radrule = 'ARITHMETIC'
      radtyp = 'R-MIN'
      radsiz = 'RADIUS'
      epsrule = 'GEOMETRIC'
      gausstyp = 'NONE'
      ngauss = 0
      abuck = 0.0_ti_p
      bbuck = 0.0_ti_p
      cbuck = 0.0_ti_p
      ghal = 0.12_ti_p
      dhal = 0.07_ti_p
      v2scale = 0.0_ti_p
      v3scale = 0.0_ti_p
      v4scale = 1.0_ti_p
      v5scale = 1.0_ti_p
      use_vcorr = .false.
!$acc update device(ghal,dhal)
!$acc update device(v2scale,v3scale,v4scale,v5scale)
c
c     set default control parameters for charge-charge terms
c
      electric = coulomb
      dielec   = 1.0_ti_p
      ebuffer  = 0.0_ti_p
      c2scale  = 0.0_ti_p
      c3scale  = 0.0_ti_p
      c4scale  = 1.0_ti_p
      c5scale  = 1.0_ti_p
      neutnbr  = .false.
      neutcut  = .false.
c
c     set default control parameters for polarizable multipoles
c
      m2scale = 0.0_ti_p
      m3scale = 0.0_ti_p
      m4scale = 1.0_ti_p
      m5scale = 1.0_ti_p
      p2scale = 0.0_ti_p
      p3scale = 0.0_ti_p
      p4scale = 1.0_ti_p
      p5scale = 1.0_ti_p
      p41scale = 0.5_ti_p
!$acc update device(m2scale,m3scale,m4scale,m5scale,p2scale,
!$acc& p3scale,p4scale,p41scale,p5scale)
c
c     set default control parameters for induced dipoles
c
      poltyp = 'MUTUAL'
      politer = 500
      polprt = 0
      polalg = 1   !TODO 1.2 Set this to 5 when finished
      polalgshort = 5
      polgsf = 0
      polff  = 0
      tcgprec = .true.
      tcgguess = .false.
      tcgpeek = .true.
      tcgprecshort = .true.
      tcgguessshort = .false.
      tcgpeekshort = .false.
      tcgorder = 2
      tcgordershort = 1
      tcgomega = 1.0d0
      tcgomegafit = .false.
      omegafitfreq = 1000
      poleps = 0.00001_re_p
      d1scale = 0.0_ti_p
      d2scale = 1.0_ti_p
      d3scale = 1.0_ti_p
      d4scale = 1.0_ti_p
      u1scale = 1.0_ti_p
      u2scale = 1.0_ti_p
      u3scale = 1.0_ti_p
      u4scale = 1.0_ti_p
!$acc update device(d1scale,d2scale,d3scale,d4scale,
!$acc& u1scale,u2scale,u3scale,u4scale)
c
c     set default divide and conquer ji/diis parameters
c
      clsttype = 2
      precomp = 0
      nocomdiis = 0
      natprblk = 60
      return
      end
