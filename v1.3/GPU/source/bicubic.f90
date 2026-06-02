!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ################################################################
!     ##                                                            ##
!     ##  subroutine bcuint  --  bicubic interpolation of function  ##
!     ##                                                            ##
!     ################################################################
!
!
!     "bcuint" performs a bicubic interpolation of the function
!     value on a 2D spline grid
!
!
#include "tinker_precision.h"
subroutine bcuint (y,y1,y2,y12,x1l,x1u,x2l,x2u,x1,x2,ansy)
!$acc routine seq
   implicit none
   integer,parameter::ti_p=t_p
   integer i
   real(t_p) x1,x1l,x1u
   real(t_p) x2,x2l,x2u
   real(t_p) t,u,ansy
   real(t_p) y(4),y12(4)
   real(t_p) y1(4),y2(4)
   real(t_p) c(4,4)
!$acc routine(bcucof) seq
!
!
!     get coefficients, then perform bicubic interpolation
!
   call bcucof (y,y1,y2,y12,x1u-x1l,x2u-x2l,c)
   t = (x1-x1l) / (x1u-x1l)
   u = (x2-x2l) / (x2u-x2l)
   ansy = 0.0_ti_p
   do i = 4, 1, -1
      ansy = t*ansy + ((c(i,4)*u+c(i,3))*u+c(i,2))*u + c(i,1)
   end do
   return
end
!
!
!     #################################################################
!     ##                                                             ##
!     ##  subroutine bcuint1  --  bicubic interpolation of gradient  ##
!     ##                                                             ##
!     #################################################################
!
!
!     "bcuint1" performs a bicubic interpolation of the function
!     value and gradient along the directions of a 2D spline grid
!
!     literature reference:
!
!     W. H. Press, S. A. Teukolsky, W. T. Vetterling and B. P.
!     Flannery, Numerical Recipes (Fortran), 2nd Ed., Cambridge
!     University Press, 1992, Section 3.6
!
!
subroutine bcuint1 (y,y1,y2,y12,x1l,x1u,x2l,x2u,&
   &x1,x2,ansy,ansy1,ansy2)
!$acc routine seq
   implicit none
   integer,parameter::ti_p=t_p
   integer i
   real(t_p) x1,x1l,x1u
   real(t_p) x2,x2l,x2u
   real(t_p) t,u,ansy
   real(t_p) ansy1,ansy2
   real(t_p) y(4),y12(4)
   real(t_p) y1(4),y2(4)
   real(t_p) c(4,4)
!$acc routine(bcucof) seq
!
!
!     get coefficients, then perform bicubic interpolation
!
   call bcucof (y,y1,y2,y12,x1u-x1l,x2u-x2l,c)
   t = (x1-x1l) / (x1u-x1l)
   u = (x2-x2l) / (x2u-x2l)
   ansy = 0.0_ti_p
   ansy1 = 0.0_ti_p
   ansy2 = 0.0_ti_p
   do i = 4, 1, -1
      ansy = t*ansy + ((c(i,4)*u+c(i,3))*u+c(i,2))*u + c(i,1)
      ansy1 = u*ansy1 +(3.0_ti_p*c(4,i)*t+2.0_ti_p*c(3,i))*t + c(2,i)
      ansy2 = t*ansy2 +(3.0_ti_p*c(i,4)*u+2.0_ti_p*c(i,3))*u + c(i,2)
   end do
   ansy1 = ansy1 / (x1u-x1l)
   ansy2 = ansy2 / (x2u-x2l)
   return
end
!
!
!     ################################################################
!     ##                                                            ##
!     ##  subroutine bcuint2  --  bicubic interpolation of Hessian  ##
!     ##                                                            ##
!     ################################################################
!
!
!     "bcuint2" performs a bicubic interpolation of the function value,
!     gradient and Hessian along the directions of a 2D spline grid
!
!
subroutine bcuint2 (y,y1,y2,y12,x1l,x1u,x2l,x2u,x1,x2,ansy,&
   &ansy1,ansy2,ansy12,ansy11,ansy22)
   implicit none
   integer,parameter::ti_p=t_p
   integer i
   real(t_p) zero,two
   real(t_p) x1,x1l,x1u,x2,x2l,x2u
   real(t_p) ansy,ansy1,ansy2
   real(t_p) ansy11,ansy22,ansy12
   real(t_p) y(4),y1(4),y2(4),y12(4)
   real(t_p) t,u,c(4,4)
   parameter(zero=0.0_ti_p,two=2.0_ti_p)
!
!
!     get coefficients, then perform bicubic interpolation
!
   call bcucof (y,y1,y2,y12,x1u-x1l,x2u-x2l,c)
   t = (x1-x1l) / (x1u-x1l)
   u = (x2-x2l) / (x2u-x2l)
   ansy = zero
   ansy1 = zero
   ansy2 = zero
   ansy11 = zero
   ansy22 = zero
   do i = 4, 1, -1
      ansy = t*ansy + ((c(i,4)*u+c(i,3))*u+c(i,2))*u + c(i,1)
      ansy1 = u*ansy1 + (3.0_ti_p*c(4,i)*t+two*c(3,i))*t + c(2,i)
      ansy2 = t*ansy2 + (3.0_ti_p*c(i,4)*u+two*c(i,3))*u + c(i,2)
      ansy11 = u*ansy11 + 6.0_ti_p*c(4,i)*t + two*c(3,i)
      ansy22 = t*ansy22 + 6.0_ti_p*c(i,4)*u + two*c(i,3)
   end do
   ansy12 = 3.0_ti_p*t*t*((3.0_ti_p*c(4,4)*u+two*c(4,3))*u+c(4,2))&
      &+ two*t*((3.0_ti_p*c(3,4)*u+two*c(3,3))*u+c(3,2))&
      &+ (3.0_ti_p*c(2,4)*u+two*c(2,3))*u + c(2,2)
   ansy1 = ansy1 / (x1u-x1l)
   ansy2 = ansy2 / (x2u-x2l)
   ansy11 = ansy11 / ((x1u-x1l)*(x1u-x1l))
   ansy22 = ansy22 / ((x2u-x2l)*(x2u-x2l))
   ansy12 = ansy12 / ((x1u-x1l)*(x2u-x2l))
   return
end
!
!
!     #################################################################
!     ##                                                             ##
!     ##  subroutine bcucof  --  bicubic interpolation coefficients  ##
!     ##                                                             ##
!     #################################################################
!
!
!     "bcucof" determines the coefficient matrix needed for bicubic
!     interpolation of a function, gradients and cross derivatives
!
!     literature reference:
!
!     W. H. Press, S. A. Teukolsky, W. T. Vetterling and B. P.
!     Flannery, Numerical Recipes (Fortran), 2nd Ed., Cambridge
!     University Press, 1992, Section 3.6
!
!
subroutine bcucof (y,y1,y2,y12,d1,d2,c)
!$acc routine seq
   implicit none
   integer,parameter::ti_p=t_p
   integer i,j,k
   real(t_p) xx,d1,d2,d1d2
   real(t_p) y(4),y12(4)
   real(t_p) y1(4),y2(4)
   real(t_p) x(16),cl(16)
   real(t_p) c(4,4)
   real(t_p) wt(16,16)
   save wt
   data wt / 1.0d0, 0.0d0,-3.0d0, 2.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,&
      &-3.0d0, 0.0d0, 9.0d0,-6.0d0, 2.0d0, 0.0d0,-6.0d0, 4.0d0,&
      &0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,&
      &3.0d0, 0.0d0,-9.0d0, 6.0d0,-2.0d0, 0.0d0, 6.0d0,-4.0d0,&
      &0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,&
      &0.0d0, 0.0d0, 9.0d0,-6.0d0, 0.0d0, 0.0d0,-6.0d0, 4.0d0,&
      &0.0d0, 0.0d0, 3.0d0,-2.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,&
      &0.0d0, 0.0d0,-9.0d0, 6.0d0, 0.0d0, 0.0d0, 6.0d0,-4.0d0,&
      &0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0,-3.0d0, 2.0d0,&
      &-2.0d0, 0.0d0, 6.0d0,-4.0d0, 1.0d0, 0.0d0,-3.0d0, 2.0d0,&
      &0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,&
      &-1.0d0, 0.0d0, 3.0d0,-2.0d0, 1.0d0, 0.0d0,-3.0d0, 2.0d0,&
      &0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,&
      &0.0d0, 0.0d0,-3.0d0, 2.0d0, 0.0d0, 0.0d0, 3.0d0,-2.0d0,&
      &0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 3.0d0,-2.0d0,&
      &0.0d0, 0.0d0,-6.0d0, 4.0d0, 0.0d0, 0.0d0, 3.0d0,-2.0d0,&
      &0.0d0, 1.0d0,-2.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,&
      &0.0d0,-3.0d0, 6.0d0,-3.0d0, 0.0d0, 2.0d0,-4.0d0, 2.0d0,&
      &0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,&
      &0.0d0, 3.0d0,-6.0d0, 3.0d0, 0.0d0,-2.0d0, 4.0d0,-2.0d0,&
      &0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,&
      &0.0d0, 0.0d0,-3.0d0, 3.0d0, 0.0d0, 0.0d0, 2.0d0,-2.0d0,&
      &0.0d0, 0.0d0,-1.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,&
      &0.0d0, 0.0d0, 3.0d0,-3.0d0, 0.0d0, 0.0d0,-2.0d0, 2.0d0,&
      &0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0,-2.0d0, 1.0d0,&
      &0.0d0,-2.0d0, 4.0d0,-2.0d0, 0.0d0, 1.0d0,-2.0d0, 1.0d0,&
      &0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,&
      &0.0d0,-1.0d0, 2.0d0,-1.0d0, 0.0d0, 1.0d0,-2.0d0, 1.0d0,&
      &0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,&
      &0.0d0, 0.0d0, 1.0d0,-1.0d0, 0.0d0, 0.0d0,-1.0d0, 1.0d0,&
      &0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,-1.0d0, 1.0d0,&
      &0.0d0, 0.0d0, 2.0d0,-2.0d0, 0.0d0, 0.0d0,-1.0d0, 1.0d0 /
!
!
!     pack a temporary vector of corner values
!
   d1d2 = d1 * d2
   do i = 1, 4
      x(i) = y(i)
      x(i+4) = y1(i) * d1
      x(i+8) = y2(i) * d2
      x(i+12) = y12(i) * d1d2
   end do
!
!     matrix multiply by the stored weight table
!
   do i = 1, 16
      xx = 0.0_ti_p
      do k = 1, 16
         xx = xx + wt(i,k)*x(k)
      end do
      cl(i) = xx
   end do
!
!     unpack the result into the coefficient table
!
   j = 0
   do i = 1, 4
      do k = 1, 4
         j = j + 1
         c(i,k) = cl(j)
      end do
   end do
   return
end
