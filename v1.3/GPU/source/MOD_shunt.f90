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
!     off    distance at which the potential energy goes to zero
!     off2   square of distance at which the potential goes to zero
!     cut    distance at which switching of the potential begins
!     cut2   square of distance at which the switching begins
!     c0     zeroth order coefficient of multiplicative switch
!     c1     first order coefficient of multiplicative switch
!     c2     second order coefficient of multiplicative switch
!     c3     third order coefficient of multiplicative switch
!     c4     fourth order coefficient of multiplicative switch
!     c5     fifth order coefficient of multiplicative switch
!     f0     zeroth order coefficient of additive switch function
!     f1     first order coefficient of additive switch function
!     f2     second order coefficient of additive switch function
!     f3     third order coefficient of additive switch function
!     f4     fourth order coefficient of additive switch function
!     f5     fifth order coefficient of additive switch function
!     f6     sixth order coefficient of additive switch function
!     f7     seventh order coefficient of additive switch function
!
!
#include "tinker_macro.h"
module shunt
   implicit none
   real(t_p) off,off2,cut,cut2
   real(t_p) c0,c1,c2,c3,c4,c5
   real(t_p) f0,f1,f2,f3,f4,f5,f6,f7
   save
end
