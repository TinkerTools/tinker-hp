!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #####################################################################
!     ##                                                                 ##
!     ##  module ktorsn  --  forcefield parameters for torsional angles  ##
!     ##                                                                 ##
!     #####################################################################
!
!
!     maxnt    maximum number of torsional angle parameter entries
!     maxnt5   maximum number of 5-membered ring torsion entries
!     maxnt4   maximum number of 4-membered ring torsion entries
!
!     t1       torsional parameters for standard 1-fold rotation
!     t2       torsional parameters for standard 2-fold rotation
!     t3       torsional parameters for standard 3-fold rotation
!     t4       torsional parameters for standard 4-fold rotation
!     t5       torsional parameters for standard 5-fold rotation
!     t6       torsional parameters for standard 6-fold rotation
!     t15      torsional parameters for 1-fold rotation in 5-ring
!     t25      torsional parameters for 2-fold rotation in 5-ring
!     t35      torsional parameters for 3-fold rotation in 5-ring
!     t45      torsional parameters for 4-fold rotation in 5-ring
!     t55      torsional parameters for 5-fold rotation in 5-ring
!     t65      torsional parameters for 6-fold rotation in 5-ring
!     t14      torsional parameters for 1-fold rotation in 4-ring
!     t24      torsional parameters for 2-fold rotation in 4-ring
!     t34      torsional parameters for 3-fold rotation in 4-ring
!     t44      torsional parameters for 4-fold rotation in 4-ring
!     t54      torsional parameters for 5-fold rotation in 4-ring
!     t64      torsional parameters for 6-fold rotation in 4-ring
!     kt       string of atom classes for torsional angles
!     kt5      string of atom classes for 5-ring torsions
!     kt4      string of atom classes for 4-ring torsions
!
!
#include "tinker_macro.h"
module ktorsn
   implicit none
   integer maxnt,maxnt5,maxnt4
   parameter (maxnt=2000)
   parameter (maxnt5=500)
   parameter (maxnt4=500)
   real(t_p) t1(2,maxnt),t2(2,maxnt)
   real(t_p) t3(2,maxnt),t4(2,maxnt)
   real(t_p) t5(2,maxnt),t6(2,maxnt)
   real(t_p) t15(2,maxnt5),t25(2,maxnt5)
   real(t_p) t35(2,maxnt5),t45(2,maxnt5)
   real(t_p) t55(2,maxnt5),t65(2,maxnt5)
   real(t_p) t14(2,maxnt4),t24(2,maxnt4)
   real(t_p) t34(2,maxnt4),t44(2,maxnt4)
   real(t_p) t54(2,maxnt4),t64(2,maxnt4)
   character*16 kt(maxnt),kt5(maxnt5),kt4(maxnt4)
   save
end
