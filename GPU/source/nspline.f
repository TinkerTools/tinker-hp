c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##################################################################
c     ##                                                              ##
c     ##   subroutine nspline  --  nonperiodic natural cubic spline   ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "nspline" computes coefficients for an nonperiodic cubic spline
c     with natural boundary conditions where the first and last second
c     derivatives are already known
c
c
#include "tinker_precision.h"
      subroutine nspline (n,x0,y0,s1,s2,h,g,dy,dla,dmu)
      implicit none
      integer,parameter::ti_p=t_p
      integer i,n
      real(t_p) t,y21,y2n
      real(t_p) x0(0:*)
      real(t_p) y0(0:*)
      real(t_p) s1(0:*)
      real(t_p) s2(0:*)
      real(t_p) h(0:*)
      real(t_p) g(0:*)
      real(t_p) dy(0:*)
      real(t_p) dla(0:*)
      real(t_p) dmu(0:*)
c
c
c     set first and last second deriviatives to zero
c
      y21 = 0.0_ti_p
      y2n = 0.0_ti_p
c
c     find the intervals to be used
c
      do i = 0, n-1
         h(i) = x0(i+1) - x0(i)
         dy(i) = (y0(i+1)-y0(i)) / h(i)
      end do
c
c     calculate the spline coeffcients
c
      do i = 1, n-1
         dla(i) = h(i) / (h(i)+h(i-1))
         dmu(i) = 1.0_ti_p - dla(i)
         g(i) = 3.0_ti_p * (dla(i)*dy(i-1)+dmu(i)*dy(i))
      end do
c
c     set the initial value via natural boundary condition
c
      dla(n) = 1.0_ti_p
      dla(0) = 0.0_ti_p
      dmu(n) = 0.0_ti_p
      dmu(0) = 1.0_ti_p
      g(0) = 3.0_ti_p*dy(0) - 0.5_ti_p*h(0)*y21
      g(n) = 3.0_ti_p*dy(n-1) + 0.5_ti_p*h(n-1)*y2n
c
c     solve the triagonal system of linear equations
c
      dmu(0) = 0.5_ti_p * dmu(0)
      g(0) = 0.5_ti_p * g(0)
      do i = 1, n
         t = 2.0_ti_p - dmu(i-1)*dla(i)
         dmu(i) = dmu(i) / t
         g(i) = (g(i)-g(i-1)*dla(i)) / t
      end do
      do i = n-1, 0, -1
         g(i) = g(i) - dmu(i)*g(i+1)
      end do
c
c     get the first derivative at each grid point
c
      do i = 0, n
         s1(i) = g(i)
      end do
c
c     get the second derivative at each grid point
c
      s2(0) = y21
      s2(n) = y2n
      do i = 1, n-1
         s2(i) = 6.0_ti_p*(y0(i+1)-y0(i))/(h(i)*h(i))
     &              - 4.0_ti_p*s1(i)/h(i) - 2.0_ti_p*s1(i+1)/h(i)
      end do
      return
      end
