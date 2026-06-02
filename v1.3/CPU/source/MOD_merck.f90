!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ##################################################################
!     ##                                                              ##
!     ##  module merck  --  parameters specific for MMFF force field  ##
!     ##                                                              ##
!     ##################################################################
!
module merck
   use sizes
   implicit none
   integer bt_1(500,2) !<atom pairs having MMFF Bond Type 1
   integer nlignes !<number of atom pairs having MMFF Bond Type 1
   integer eqclass(500,5) !<table of atom class equivalencies used to find default parameters if explicit values are missing (see J. Comput. Chem., 17, 490-519, '95, Table IV)
   integer :: crd(100) !<number of attached neighbors
   integer :: val(100) !<valency value                   |  see T. A. Halgren,
   integer :: pilp(100) !<if 0, no lone pair              |  J. Comput. Chem., if 1, one or more lone pair(s)  |  17, 616-645 (1995)
   integer :: mltb(100) !<multibond indicator             |
   integer :: arom(100) !<aromaticity indicator           |
   integer :: lin(100) !<linearity indicator             |
   integer :: sbmb(100) !<single- vs multiple-bond flag   |
   integer :: mmffarom(maxtyp,6) !<aromatic rings parameters
   integer :: mmffaromc(maxtyp,6) !<cationic aromatic rings parameters
   integer :: mmffaroma(maxtyp,6) !<anionic aromatic rings parameters
   real*8 :: mmff_kb(100,100),mmff_kb1(100,100)
   real*8 :: mmff_b0(100,100),mmff_b1(100,100)
   real*8 :: rad0(100),r0ref(100,100)
   real*8 :: kbref(100,100),paulel(100)
   real*8 :: mmff_ka(0:100,100,0:100),mmff_ka1(0:100,100,0:100)
   real*8 :: mmff_ka2(0:100,100,0:100)
   real*8 :: mmff_ka3(0:100,100,0:100),mmff_ka4(0:100,100,0:100)
   real*8 :: mmff_ka5(0:100,100,0:100)
   real*8 :: mmff_ka6(0:100,100,0:100),mmff_ka7(0:100,100,0:100)
   real*8 :: mmff_ka8(0:100,100,0:100)
   real*8 :: mmff_ang0(0:100,100,0:100),mmff_ang1(0:100,100,0:100)
   real*8 :: mmff_ang2(0:100,100,0:100)
   real*8 :: mmff_ang3(0:100,100,0:100),mmff_ang4(0:100,100,0:100)
   real*8 :: mmff_ang5(0:100,100,0:100)
   real*8 :: mmff_ang6(0:100,100,0:100),mmff_ang7(0:100,100,0:100)
   real*8 :: mmff_ang8(0:100,100,0:100)
   real*8 :: stbn_abc(100,100,100),stbn_cba(100,100,100)
   real*8 :: stbn_abc1(100,100,100),stbn_cba1(100,100,100)
   real*8 :: stbn_abc2(100,100,100),stbn_cba2(100,100,100)
   real*8 :: stbn_abc3(100,100,100),stbn_cba3(100,100,100)
   real*8 :: stbn_abc4(100,100,100),stbn_cba4(100,100,100)
   real*8 :: stbn_abc5(100,100,100),stbn_cba5(100,100,100)
   real*8 :: stbn_abc6(100,100,100),stbn_cba6(100,100,100)
   real*8 :: stbn_abc7(100,100,100),stbn_cba7(100,100,100)
   real*8 :: stbn_abc8(100,100,100),stbn_cba8(100,100,100)
   real*8 :: stbn_abc9(100,100,100),stbn_cba9(100,100,100)
   real*8 :: stbn_abc10(100,100,100),stbn_cba10(100,100,100)
   real*8 :: stbn_abc11(100,100,100),stbn_cba11(100,100,100)
   real*8 :: defstbn_abc(0:4,0:4,0:4),defstbn_cba(0:4,0:4,0:4)
   real*8 :: t1_1(2,0:2000),t2_1(2,0:2000),t3_1(2,0:2000)
   real*8 :: t1_2(2,0:2000),t2_2(2,0:2000),t3_2(2,0:2000)
   character*16 :: kt_1(0:2000),kt_2(0:2000)
   real*8 :: g(maxclass),alph(maxclass),nn(maxclass)
   character*1 :: da(maxclass)
   real*8 :: bci(100,100),bci_1(100,100)
   real*8 :: pbci(maxclass),fcadj(maxclass)
   save
end
