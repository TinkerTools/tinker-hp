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
!
!     bt_1       atom pairs having MMFF Bond Type 1
!     nlignes    number of atom pairs having MMFF Bond Type 1
!     eqclass    table of atom class equivalencies used to find
!                default parameters if explicit values are missing
!                (see J. Comput. Chem., 17, 490-519, '95, Table IV)
!     crd        number of attached neighbors    |
!     val        valency value                   |  see T. A. Halgren,
!     pilp       if 0, no lone pair              |  J. Comput. Chem.,
!                if 1, one or more lone pair(s)  |  17, 616-645 (1995)
!     mltb       multibond indicator             |
!     arom       aromaticity indicator           |
!     lin        linearity indicator             |
!     sbmb       single- vs multiple-bond flag   |
!     mmffarom   aromatic rings parameters
!     mmffaromc  cationic aromatic rings parameters
!     mmffaroma  anionic aromatic rings parameters
!
!
#include "tinker_macro.h"
module merck
   use sizes
   implicit none
   integer bt_1(500,2)
   integer nlignes
   integer eqclass(500,5)
   integer crd(100),val(100),pilp(100),mltb(100)
   integer arom(100),lin(100),sbmb(100)
   integer mmffarom(maxtyp,6)
   integer mmffaromc(maxtyp,6)
   integer mmffaroma(maxtyp,6)
!      common /merck1/ bt_1(500,2),nlignes,eqclass(500,5),crd(100),
!     &                val(100),pilp(100),mltb(100),arom(100),lin(100),
!     &                sbmb(100),mmffarom(maxtyp,6),
!     &                mmffaromc(maxtyp,6),mmffaroma(maxtyp,6)
!
!
!     mmff_kb   bond force constant for pairs of atom classes
!     mmff_kb1  bond force constant for class pairs with Bond Type 1
!     mmff_b0   bond length value for pairs of atom classes
!     mmff_b1   bond length value for class pairs with Bond Type 1
!     rad0      covalent atomic radius for empirical bond rules
!     paulel    Pauling electronegativities for empirical bond rules
!     r0ref     reference bond length for empirical bond rules
!     kbref     reference force constant for empirical bond rules
!
!
   real(t_p) mmff_kb(100,100),mmff_kb1(100,100)
   real(t_p) mmff_b0(100,100),mmff_b1(100,100)
   real(t_p) rad0(100),r0ref(100,100)
   real(t_p) kbref(100,100),paulel(100)
!      common /merck2/ mmff_kb(100,100),mmff_kb1(100,100),
!     &                mmff_b0(100,100),mmff_b1(100,100),rad0(100),
!     &                r0ref(100,100),kbref(100,100),paulel(100)
!
!
!     mmff_ka     angle force constant for triples of atom classes
!     mmff_ka1    angle force constant with one bond of Type 1
!     mmff_ka2    angle force constant with both bonds of Type 1
!     mmff_ka3    angle force constant for 3-membered ring
!     mmff_ka4    angle force constant for 4-membered ring
!     mmff_ka5    angle force constant for 3-ring and one Bond Type 1
!     mmff_ka6    angle force constant for 3-ring and both Bond Type 1
!     mmff_ka7    angle force constant for 4-ring and one Bond Type 1
!     mmff_ka8    angle force constant for 4-ring and both Bond Type 1
!     mmff_ang0   ideal bond angle for triples of atom classes
!     mmff_ang1   ideal bond angle with one bond of Type 1
!     mmff_ang2   ideal bond angle with both bonds of Type 1
!     mmff_ang3   ideal bond angle for 3-membered ring
!     mmff_ang4   ideal bond angle for 4-membered ring
!     mmff_ang5   ideal bond angle for 3-ring and one Bond Type 1
!     mmff_ang6   ideal bond angle for 3-ring and both Bond Type 1
!     mmff_ang7   ideal bond angle for 4-ring and one Bond Type 1
!     mmff_ang8   ideal bond angle for 4-ring and both Bond Type 1
!
!
   real(t_p) mmff_ka (0:100,100,0:100)
   real(t_p) mmff_ka1(0:100,100,0:100)
   real(t_p) mmff_ka2(0:100,100,0:100)
   real(t_p) mmff_ka3(0:100,100,0:100)
   real(t_p) mmff_ka4(0:100,100,0:100)
   real(t_p) mmff_ka5(0:100,100,0:100)
   real(t_p) mmff_ka6(0:100,100,0:100)
   real(t_p) mmff_ka7(0:100,100,0:100)
   real(t_p) mmff_ka8(0:100,100,0:100)
   real(t_p) mmff_ang0(0:100,100,0:100)
   real(t_p) mmff_ang1(0:100,100,0:100)
   real(t_p) mmff_ang2(0:100,100,0:100)
   real(t_p) mmff_ang3(0:100,100,0:100)
   real(t_p) mmff_ang4(0:100,100,0:100)
   real(t_p) mmff_ang5(0:100,100,0:100)
   real(t_p) mmff_ang6(0:100,100,0:100)
   real(t_p) mmff_ang7(0:100,100,0:100)
   real(t_p) mmff_ang8(0:100,100,0:100)
!      common /merck3/ mmff_ka(0:100,100,0:100),
!     &                mmff_ka1(0:100,100,0:100),
!     &                mmff_ka2(0:100,100,0:100),
!     &                mmff_ka3(0:100,100,0:100),
!     &                mmff_ka4(0:100,100,0:100),
!     &                mmff_ka5(0:100,100,0:100),
!     &                mmff_ka6(0:100,100,0:100),
!     &                mmff_ka7(0:100,100,0:100),
!     &                mmff_ka8(0:100,100,0:100),
!     &                mmff_ang0(0:100,100,0:100),
!     &                mmff_ang1(0:100,100,0:100),
!     &                mmff_ang2(0:100,100,0:100),
!     &                mmff_ang3(0:100,100,0:100),
!     &                mmff_ang4(0:100,100,0:100),
!     &                mmff_ang5(0:100,100,0:100),
!     &                mmff_ang6(0:100,100,0:100),
!     &                mmff_ang7(0:100,100,0:100),
!     &                mmff_ang8(0:100,100,0:100)
!
!
!     Stretch-Bend Type 0
!     stbn_abc     stretch-bend parameters for A-B-C atom classes
!     stbn_cba     stretch-bend parameters for C-B-A atom classes
!     Stretch-Bend Type 1  (A-B is Bond Type 1)
!     stbn_abc1    stretch-bend parameters for A-B-C atom classes
!     stbn_cba1    stretch-bend parameters for C-B-A atom classes
!     Stretch-Bend Type 2  (B-C is Bond Type 1)
!     stbn_abc2    stretch-bend parameters for A-B-C atom classes
!     stbn_cba2    stretch-bend parameters for C-B-A atom classes
!     Stretch-Bend Type = 3  (A-B and B-C are Bond Type 1)
!     stbn_abc3    stretch-bend parameters for A-B-C atom classes
!     stbn_cba3    stretch-bend parameters for C-B-A atom classes
!     Stretch-Bend Type 4  (both Bond Types 0, 4-membered ring)
!     stbn_abc4    stretch-bend parameters for A-B-C atom classes
!     stbn_cba4    stretch-bend parameters for C-B-A atom classes
!     Stretch-Bend Type 5  (both Bond Types 0, 3-membered ring)
!     stbn_abc5    stretch-bend parameters for A-B-C atom classes
!     stbn_cba5    stretch-bend parameters for C-B-A atom classes
!     Stretch-Bend Type 6  (A-B is Bond Type 1, 3-membered ring)
!     stbn_abc6    stretch-bend parameters for A-B-C atom classes
!     stbn_cba6    stretch-bend parameters for C-B-A atom classes
!     Stretch-Bend Type 7  (B-C is Bond Type 1, 3-membered ring)
!     stbn_abc7    stretch-bend parameters for A-B-C atom classes
!     stbn_cba7    stretch-bend parameters for C-B-A atom classes
!     Stretch-Bend Type 8  (both Bond Types 1, 3-membered ring)
!     stbn_abc8    stretch-bend parameters for A-B-C atom classes
!     stbn_cba8    stretch-bend parameters for C-B-A atom classes
!     Stretch-Bend Type 9  (A-B is Bond Type 1, 4-membered ring)
!     stbn_abc9    stretch-bend parameters for A-B-C atom classes
!     stbn_cba9    stretch-bend parameters for C-B-A atom classes
!     Stretch-Bend Type 10  (B-C is Bond Type 1, 4-membered ring)
!     stbn_abc10   stretch-bend parameters for A-B-C atom classes
!     stbn_cba10   stretch-bend parameters for C-B-A atom classes
!     Stretch-Bend Type 11  (both Bond Types 1, 4-membered ring)
!     stbn_abc11   stretch-bend parameters for A-B-C atom classes
!     stbn_cba11   stretch-bend parameters for C-B-A atom classes
!     defstbn_abc  default stretch-bend parameters for A-B-C classes
!     defstbn_cba  default stretch-bend parameters for C-B-A classes
!
!
   real(t_p) stbn_abc(100,100,100),stbn_cba(100,100,100)
   real(t_p) stbn_abc1(100,100,100),stbn_cba1(100,100,100)
   real(t_p) stbn_abc2(100,100,100),stbn_cba2(100,100,100)
   real(t_p) stbn_abc3(100,100,100),stbn_cba3(100,100,100)
   real(t_p) stbn_abc4(100,100,100),stbn_cba4(100,100,100)
   real(t_p) stbn_abc5(100,100,100),stbn_cba5(100,100,100)
   real(t_p) stbn_abc6(100,100,100),stbn_cba6(100,100,100)
   real(t_p) stbn_abc7(100,100,100),stbn_cba7(100,100,100)
   real(t_p) stbn_abc8(100,100,100),stbn_cba8(100,100,100)
   real(t_p) stbn_abc9(100,100,100),stbn_cba9(100,100,100)
   real(t_p) stbn_abc10(100,100,100),stbn_cba10(100,100,100)
   real(t_p) stbn_abc11(100,100,100),stbn_cba11(100,100,100)
   real(t_p),dimension(0:4,0:4,0:4)::defstbn_cba,defstbn_abc
!      common /merck4/ stbn_abc(100,100,100),stbn_cba(100,100,100),
!     &                stbn_abc1(100,100,100),stbn_cba1(100,100,100),
!     &                stbn_abc2(100,100,100),stbn_cba2(100,100,100),
!     &                stbn_abc3(100,100,100),stbn_cba3(100,100,100),
!     &                stbn_abc4(100,100,100),stbn_cba4(100,100,100),
!     &                stbn_abc5(100,100,100),stbn_cba5(100,100,100),
!     &                stbn_abc6(100,100,100),stbn_cba6(100,100,100),
!     &                stbn_abc7(100,100,100),stbn_cba7(100,100,100),
!     &                stbn_abc8(100,100,100),stbn_cba8(100,100,100),
!     &                stbn_abc9(100,100,100),stbn_cba9(100,100,100),
!     &                stbn_abc10(100,100,100),stbn_cba10(100,100,100),
!     &                stbn_abc11(100,100,100),stbn_cba11(100,100,100),
!     &                defstbn_abc(0:4,0:4,0:4),defstbn_cba(0:4,0:4,0:4)
!
!
!     t1_1     torsional parameters for 1-fold, MMFF Torsion Type 1
!     t1_2     torsional parameters for 1-fold, MMFF Torsion Type 2
!     t2_1     torsional parameters for 2-fold, MMFF Torsion Type 1
!     t2_2     torsional parameters for 2-fold, MMFF Torsion Type 2
!     t3_1     torsional parameters for 3-fold, MMFF Torsion Type 1
!     t3_2     torsional parameters for 3-fold, MMFF Torsion Type 2
!     kt_1     string of classes for torsions, MMFF Torsion Type 1
!     kt_2     string of classes for torsions, MMFF Torsion Type 2
!
!
   real(t_p) t1_1(2,0:2000),t2_1(2,0:2000),t3_1(2,0:2000)
   real(t_p) t1_2(2,0:2000),t2_2(2,0:2000),t3_2(2,0:2000)
   character*16 kt_1(0:2000),kt_2(0:2000)
!      common /merck5/ t1_1(2,0:2000),t2_1(2,0:2000),t3_1(2,0:2000),
!     &                t1_2(2,0:2000),t2_2(2,0:2000),t3_2(2,0:2000),
!     &                kt_1(0:2000),kt_2(0:2000)
!
!
!     g        scale factors for calculation of MMFF eps
!     alph     atomic polarizabilities for calculation of MMFF eps
!     nn       effective number of valence electrons for MMFF eps
!     da       donor/acceptor atom classes
!
!
   real(t_p) g(maxclass),alph(maxclass),nn(maxclass)
   character*1 da(maxclass)
!      common /merck6/ g(maxclass),alph(maxclass),nn(maxclass),
!     &                da(maxclass)
!
!
!     bci      bond charge increments for building atom charges
!     bci_1    bond charge increments for MMFF Bond Type 1
!     pbci     partial BCI for building missing BCI's
!     fcadj    formal charge adjustment factor
!
!
   real(t_p) bci(100,100),bci_1(100,100)
   real(t_p) pbci(maxclass),fcadj(maxclass)
!      common /merck7/ bci(100,100),bci_1(100,100),pbci(maxclass),
!     &                fcadj(maxclass)
   save
end
