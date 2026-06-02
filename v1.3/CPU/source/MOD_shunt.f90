!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ####################################################################
!     ##                                                                ##
!     ##  module shunt  --  polynomial switching function coefficients  ##
!     ##                                                                ##
!     ####################################################################
!
!
module shunt
   implicit none
   real*8 :: off !<distance at which the potential energy goes to zero
   real*8 :: off2 !<square of distance at which the potential goes to zero
   real*8 :: cut !<distance at which switching of the potential begins
   real*8 :: cut2 !<square of distance at which the switching begins
   real*8 :: c0 !<zeroth order coefficient of multiplicative switch
   real*8 :: c1 !<first order coefficient of multiplicative switch
   real*8 :: c2 !<second order coefficient of multiplicative switch
   real*8 :: c3 !<third order coefficient of multiplicative switch
   real*8 :: c4 !<fourth order coefficient of multiplicative switch
   real*8 :: c5 !<fifth order coefficient of multiplicative switch
   real*8 :: f0 !<zeroth order coefficient of additive switch function
   real*8 :: f1 !<first order coefficient of additive switch function
   real*8 :: f2 !<second order coefficient of additive switch function
   real*8 :: f3 !<third order coefficient of additive switch function
   real*8 :: f4 !<fourth order coefficient of additive switch function
   real*8 :: f5 !<fifth order coefficient of additive switch function
   real*8 :: f6 !<sixth order coefficient of additive switch function
   real*8 :: f7 !<seventh order coefficient of additive switch function
   save
end
