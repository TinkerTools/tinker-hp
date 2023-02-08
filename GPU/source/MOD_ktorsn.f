c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #####################################################################
c     ##                                                                 ##
c     ##  module ktorsn  --  forcefield parameters for torsional angles  ##
c     ##                                                                 ##
c     #####################################################################
c
c
c     maxnt    maximum number of torsional angle parameter entries
c     maxnt5   maximum number of 5-membered ring torsion entries
c     maxnt4   maximum number of 4-membered ring torsion entries
c
c     t1       torsional parameters for standard 1-fold rotation
c     t2       torsional parameters for standard 2-fold rotation
c     t3       torsional parameters for standard 3-fold rotation
c     t4       torsional parameters for standard 4-fold rotation
c     t5       torsional parameters for standard 5-fold rotation
c     t6       torsional parameters for standard 6-fold rotation
c     t15      torsional parameters for 1-fold rotation in 5-ring
c     t25      torsional parameters for 2-fold rotation in 5-ring
c     t35      torsional parameters for 3-fold rotation in 5-ring
c     t45      torsional parameters for 4-fold rotation in 5-ring
c     t55      torsional parameters for 5-fold rotation in 5-ring
c     t65      torsional parameters for 6-fold rotation in 5-ring
c     t14      torsional parameters for 1-fold rotation in 4-ring
c     t24      torsional parameters for 2-fold rotation in 4-ring
c     t34      torsional parameters for 3-fold rotation in 4-ring
c     t44      torsional parameters for 4-fold rotation in 4-ring
c     t54      torsional parameters for 5-fold rotation in 4-ring
c     t64      torsional parameters for 6-fold rotation in 4-ring
c     kt       string of atom classes for torsional angles
c     kt5      string of atom classes for 5-ring torsions
c     kt4      string of atom classes for 4-ring torsions
c
c
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
