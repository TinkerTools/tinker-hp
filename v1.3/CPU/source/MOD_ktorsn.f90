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
module ktorsn
   implicit none
   integer :: maxnt !<maximum number of torsional angle parameter entries
   integer :: maxnt5 !<maximum number of 5-membered ring torsion entries
   integer :: maxnt4 !<maximum number of 4-membered ring torsion entries
   parameter (maxnt=2000)
   parameter (maxnt5=500)
   parameter (maxnt4=500)
   real*8 :: t1(2,maxnt) !<torsional parameters for standard 1-fold rotation
   real*8 :: t2(2,maxnt) !<torsional parameters for standard 2-fold rotation
   real*8 :: t3(2,maxnt) !<torsional parameters for standard 3-fold rotation
   real*8 :: t4(2,maxnt) !<torsional parameters for standard 4-fold rotation
   real*8 :: t5(2,maxnt) !<torsional parameters for standard 5-fold rotation
   real*8 :: t6(2,maxnt) !<torsional parameters for standard 6-fold rotation
   real*8 :: t15(2,maxnt5) !<torsional parameters for 1-fold rotation in 5-ring
   real*8 :: t25(2,maxnt5) !<torsional parameters for 2-fold rotation in 5-ring
   real*8 :: t35(2,maxnt5) !<torsional parameters for 3-fold rotation in 5-ring
   real*8 :: t45(2,maxnt5) !<torsional parameters for 4-fold rotation in 5-ring
   real*8 :: t55(2,maxnt5) !<torsional parameters for 5-fold rotation in 5-ring
   real*8 :: t65(2,maxnt5) !<torsional parameters for 6-fold rotation in 5-ring
   real*8 :: t14(2,maxnt4) !<torsional parameters for 1-fold rotation in 4-ring
   real*8 :: t24(2,maxnt4) !<torsional parameters for 2-fold rotation in 4-ring
   real*8 :: t34(2,maxnt4) !<torsional parameters for 3-fold rotation in 4-ring
   real*8 :: t44(2,maxnt4) !<torsional parameters for 4-fold rotation in 4-ring
   real*8 :: t54(2,maxnt4) !<torsional parameters for 5-fold rotation in 4-ring
   real*8 :: t64(2,maxnt4) !<torsional parameters for 6-fold rotation in 4-ring
   character*16 :: kt(maxnt) !<string of atom classes for torsional angles
   character*16 :: kt5(maxnt5) !<string of atom classes for 5-ring torsions
   character*16 :: kt4(maxnt4) !<string of atom classes for 4-ring torsions
   save
end
