!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ##################################################################
!     ##                                                              ##
!     ##  subroutine cspline  --  periodic interpolating cube spline  ##
!     ##                                                              ##
!     ##################################################################
!
!
!     "cspline" computes the coefficients for a periodic interpolating
!     cubic spline
!
!     literature reference:
!
!     G. Engeln-Mullges and F. Uhlig, Numerical Algorithms with Fortran,
!     Springer Verlag, 1996, Section 10.1.2  [see routine "isplpe"]
!
!
#include "tinker_precision.h"
subroutine cspline (n,xn,fn,b,c,d,h,du,dm,rc,rs)
   use tinheader
   use iounit
   implicit none
   integer i,n,iflag
   real(t_p) eps,average
   real(t_p) temp1,temp2
   real(t_p) xn(0:*)
   real(t_p) fn(0:*)
   real(t_p) b(0:*)
   real(t_p) c(0:*)
   real(t_p) d(0:*)
   real(t_p) h(0:*)
   real(t_p) du(0:*)
   real(t_p) dm(0:*)
   real(t_p) rc(0:*)
   real(t_p) rs(0:*)
!
!
!     check the periodicity of fn, and for subsequent call
!
   eps = 0.000001_ti_p
   if (abs(fn(n)-fn(0)) .gt. eps) then
      write (iout,10)  fn(0),fn(n)
10    format (' CSPLINE  --  Warning, Non-Periodic Input',&
         &' Values',2f12.5)
   end if
   average = 0.5_ti_p * (fn(0) + fn(n))
   fn(0) = average
   fn(n) = average
!
!     get auxiliary variables and matrix elements on first call
!
   do i = 0, n-1
      h(i) = xn(i+1) - xn(i)
   end do
   h(n) = h(0)
   do i = 1, n-1
      du(i) = h(i)
   end do
   du(n) = h(0)
   do i = 1, n
      dm(i) = 2.0_ti_p * (h(i-1)+h(i))
   end do
!
!     compute the right hand side
!
   temp1 = (fn(1)-fn(0)) / h(0)
   do i = 1, n-1, 1
      temp2 = (fn(i+1)-fn(i)) / h(i)
      rs(i)  = 3.0_ti_p * (temp2-temp1)
      temp1 = temp2
   end do
   rs(n) = 3.0_ti_p * ((fn(1)-fn(0))/h(0)-temp1)
!
!     solve the linear system with factorization
!
   call cytsy (n,dm,du,rc,rs,c,iflag)
   if (iflag .ne. 1)  return
!
!     compute remaining spline coefficients
!
   c(0) = c(n)
   do i = 0, n-1
      b(i) = (fn(i+1)-fn(i))/h(i) -&
         &h(i)/3.0_ti_p*(c(i+1)+2.0_ti_p*c(i))
      d(i) = (c(i+1)-c(i)) / (3.0_ti_p*h(i))
   end do
   b(n) = (fn(1)-fn(n))/h(n) - h(n)/3.0_ti_p*(c(1)+2.0_ti_p*c(n))
   return
end
!
!
!     #############################################################
!     ##                                                         ##
!     ##  subroutine cytsy  --  solve cyclic tridiagonal system  ##
!     ##                                                         ##
!     #############################################################
!
!
!     "cytsy" solves a system of linear equations for a cyclically
!     tridiagonal, symmetric, positive definite matrix
!
!     literature reference:
!
!     G. Engeln-Mullges and F. Uhlig, Numerical Algorithms with Fortran,
!     Springer Verlag, 1996, Section 4.11.2
!
!
subroutine cytsy (n,dm,du,cr,rs,x,iflag)
   implicit none
   integer n,iflag
   real(t_p) dm(0:*)
   real(t_p) du(0:*)
   real(t_p) cr(0:*)
   real(t_p) rs(0:*)
   real(t_p) x(0:*)
!
!
!     factorization of the input matrix
!
   iflag = -2
   if (n .lt. 3)  return
   call cytsyp (n,dm,du,cr,iflag)
!
!     update and back substitute as necessary
!
   if (iflag .eq. 1)  call cytsys (n,dm,du,cr,rs,x)
   return
end
!
!
!     #################################################################
!     ##                                                             ##
!     ##  subroutine cytsyp  --  tridiagonal Cholesky factorization  ##
!     ##                                                             ##
!     #################################################################
!
!
!     "cytsyp" finds the Cholesky factors of a cyclically tridiagonal
!     symmetric, positive definite matrix given by two vectors
!
!     literature reference:
!
!     G. Engeln-Mullges and F. Uhlig, Numerical Algorithms with Fortran,
!     Springer Verlag, 1996, Section 4.11.2
!
!
subroutine cytsyp (n,dm,du,cr,iflag)
   implicit none
   integer,parameter::ti_p=t_p
   integer i,n,iflag
   real(t_p) eps,row,d
   real(t_p) temp1,temp2
   real(t_p) dm(0:*)
   real(t_p) du(0:*)
   real(t_p) cr(0:*)
!
!
!     set error bound and test for condition n greater than 2
!
   eps = 0.00000001_ti_p
   iflag = -2
   if (n .lt. 3)  return
!
!     checking to see if matrix is positive definite
!
   row = abs(dm(1)) + abs(du(1)) + abs(du(n))
   if (row .eq. 0.0_ti_p) then
      iflag = 0
      return
   end if
   d = 1.0_ti_p / row
   if (dm(1) .lt. 0.0_ti_p) then
      iflag = -1
      return
   else if (abs(dm(1))*d .le. eps) then
      iflag = 0
      return
   end if
!
!     factoring a while checking for a positive definite and strong
!     nonsingular matrix a
!
   temp1 = du(1)
   du(1) = du(1) / dm(1)
   cr(1) = du(n) / dm(1)
   do i = 2, n-1
      row = abs(dm(i)) + abs(du(i)) + abs(temp1)
      if (row .eq. 0.0_ti_p) then
         iflag = 0
         return
      end if
      d = 1.0_ti_p / row
      dm(i) = dm(i) - temp1*du(i-1)
      if (dm(i) .lt. 0.0_ti_p) then
         iflag = -1
         return
      else if (abs(dm(i))*d .le. eps) then
         iflag = 0
         return
      end if
      if (i .lt. (n-1)) then
         cr(i) = -temp1 * cr(i-1) / dm(i)
         temp1 = du(i)
         du(i) = du(i) / dm(i)
      else
         temp2 = du(i)
         du(i) = (du(i) - temp1*cr(i-1)) / dm(i)
      end if
   end do
   row = abs(du(n)) + abs(dm(n)) + abs(temp2)
   if (row .eq. 0.0_ti_p) then
      iflag = 0
      return
   end if
   d = 1.0_ti_p / row
   dm(n) = dm(n) - dm(n-1)*du(n-1)*du(n-1)
   temp1 = 0.0_ti_p
   do i = 1, n-2
      temp1 = temp1 + dm(i)*cr(i)*cr(i)
   end do
   dm(n) = dm(n) - temp1
   if (dm(n) .lt. 0) then
      iflag = -1
      return
   else if (abs(dm(n))*d .le. eps) then
      iflag = 0
      return
   end if
   iflag = 1
   return
end
!
!
!     ################################################################
!     ##                                                            ##
!     ##  subroutine cytsys  --  tridiagonal solution from factors  ##
!     ##                                                            ##
!     ################################################################
!
!
!     "cytsys" solves a cyclically tridiagonal linear system
!     given the Cholesky factors
!
!     literature reference:
!
!     G. Engeln-Mullges and F. Uhlig, Numerical Algorithms with Fortran,
!     Springer Verlag, 1996, Section 4.11.2
!
!
subroutine cytsys (n,dm,du,cr,rs,x)
   implicit none
   integer i,n
   real(t_p) sum,temp
   real(t_p) dm(0:*)
   real(t_p) du(0:*)
   real(t_p) cr(0:*)
   real(t_p) rs(0:*)
   real(t_p) x(0:*)
!
!
!     updating phase
!
   temp = rs(1)
   rs(1) = temp / dm(1)
   sum = cr(1) * temp
   do i = 2, n-1
      temp = rs(i) - du(i-1)*temp
      rs(i) = temp / dm(i)
      if (i .ne. (n-1))  sum = sum + cr(i)*temp
   end do
   temp = rs(n) - du(n-1)*temp
   temp = temp - sum
   rs(n) = temp / dm(n)
!
!     back substitution phase
!
   x(n) = rs(n)
   x(n-1) = rs(n-1) - du(n-1)*x(n)
   do i = n-2, 1, -1
      x(i) = rs(i) - du(i)*x(i+1) - cr(i)*x(n)
   end do
   return
end

!-----------------------------------------------------------------------

subroutine cubic_spline (n, x, y, b, c, d)
   integer n
   real(t_p) x(n), y(n), b(n), c(n), d(n)
!     (adapted from https://www.netlib.org/fmm/spline.f)
!  the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed
!  for a cubic interpolating spline
!
!    s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
!
!    for  x(i) .le. x .le. x(i+1)
!
!  input..
!
!    n = the number of data points or knots (n.ge.2)
!    x = the abscissas of the knots in strictly increasing order
!    y = the ordinates of the knots
!
!  output..
!
!    b, c, d  = arrays of spline coefficients as defined above.
!
!  using  p  to denote differentiation,
!
!    y(i) = s(x(i))
!    b(i) = sp(x(i))
!    c(i) = spp(x(i))/2
!    d(i) = sppp(x(i))/6  (derivative from the right)
!
!  the accompanying function subprogram  seval  can be used
!  to evaluate the spline.
!
!
   integer nm1, ib, i
   real(t_p) t
!
   nm1 = n-1
   if ( n .lt. 2 ) return
   if ( n .lt. 3 ) then
      b(1) = (y(2)-y(1))/(x(2)-x(1))
      c(1) = 0.
      d(1) = 0.
      b(2) = b(1)
      c(2) = 0.
      d(2) = 0.
      return
   endif
!
!  set up tridiagonal system
!
!  b = diagonal, d = offdiagonal, c = right hand side.
!
   d(1) = x(2) - x(1)
   c(2) = (y(2) - y(1))/d(1)
   do i = 2, nm1
      d(i) = x(i+1) - x(i)
      b(i) = 2.*(d(i-1) + d(i))
      c(i+1) = (y(i+1) - y(i))/d(i)
      c(i) = c(i+1) - c(i)
   enddo
!
!  end conditions.  third derivatives at  x(1)  and  x(n)
!  obtained from divided differences
!
   b(1) = -d(1)
   b(n) = -d(n-1)
   c(1) = 0.
   c(n) = 0.
   if ( n /= 3 ) then
      c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
      c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
      c(1) = c(1)*d(1)**2/(x(4)-x(1))
      c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
   endif
!
!  forward elimination
!
   do i = 2, n
      t = d(i-1)/b(i-1)
      b(i) = b(i) - t*d(i-1)
      c(i) = c(i) - t*c(i-1)
   enddo
!
!  back substitution
!
   c(n) = c(n)/b(n)
   do ib = 1, nm1
      i = n-ib
      c(i) = (c(i) - d(i)*c(i+1))/b(i)
   enddo
!
!  c(i) is now the sigma(i) of the text
!
!  compute polynomial coefficients
!
   b(n) = (y(n) - y(nm1))/d(nm1) + d(nm1)*(c(nm1) + 2.*c(n))
   do i = 1, nm1
      b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.*c(i))
      d(i) = (c(i+1) - c(i))/d(i)
      c(i) = 3.*c(i)
   enddo
   c(n) = 3.*c(n)
   d(n) = d(n-1)

end subroutine cubic_spline

function seval(n, u, x, y, b, c, d)
   integer n
   real(t_p) seval
   real(t_p)  u, x(n), y(n), b(n), c(n), d(n)
!     (adapted from https://www.netlib.org/fmm/seval.f)
!
!  this subroutine evaluates the cubic spline function
!
!    seval = y(i) + b(i)*(u-x(i)) + c(i)*(u-x(i))**2 + d(i)*(u-x(i))**3
!
!    where  x(i) .lt. u .lt. x(i+1), using horner's rule
!
!  if  u .lt. x(1) then  i = 1  is used.
!  if  u .ge. x(n) then  i = n  is used.
!
!  input..
!
!    n = the number of data points
!    u = the abscissa at which the spline is to be evaluated
!    x,y = the arrays of data abscissas and ordinates
!    b,c,d = arrays of spline coefficients computed by spline
!
!  if  u  is not in the same interval as the previous call, then a
!  binary search is performed to determine the proper interval.
!
   integer j, k
   real(t_p) dx
   logical found
   integer, save :: i=1
   if ( i .ge. n ) i = 1
   found = u >= x(i) .and. u<=x(i+1)

   if ( .not. found ) then
      ! binary search
      i = 1
      j = n+1
      do
         k = (i+j)/2
         if ( u .lt. x(k) ) j = k
         if ( u .ge. x(k) ) i = k
         if ( j <= i+1 ) exit
      enddo
   endif
!
!  evaluate spline
!
   dx = u - x(i)
   seval = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
   return
end
