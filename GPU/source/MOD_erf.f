c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  module erf  --  data need to cumpute erf and erfc function ##
c     ##                                                             ##
c     #################################################################
c
c     literature reference:
c
c     W. J. Cody, "Rational Chebyshev Approximations for the Error
c     Function", Mathematics of Computation, 631-638, 1969
c
c     adapted from an original program written by W. J. Cody,
c     Mathematics and Computer Science Division, Argonne National
c     Laboratory, Argonne, IL 60439
c
c     machine-dependent constants:
c
c     xtiny   argument below which erf(x) may be represented by
c             2*x/sqrt(pi) and above which x*x won't underflow;
c             a conservative value is the largest machine number
c             X such that 1.0 + X = 1.0 to machine precision
c
c     xbig    largest argument acceptable for erfc; solution to
c             the equation:  W(x) * (1-0.5/x**2) = XMIN, where
c             W(x) = exp(-x*x)/[x*sqrt(pi)]
c
c
#include "tinker_precision.h"
      module erf_mod
      use tinheader ,only:ti_p
      implicit none
      private
      real(t_p) erf_xtiny,erf_xbig,erf_sqrpi
      real(t_p) erf_a(5),erf_b(4)
      real(t_p) erf_c(9),erf_d(8)
      real(t_p) erf_p(6),erf_q(5)
      real(t_p) a_inv(4),b_inv(4)
      real(t_p) c_inv(4),d_inv(2)
c
c     mathematical and machine-dependent constants
c
      parameter(
     &   erf_sqrpi = 5.6418958354775628695d-1 ,
#ifdef SINGLE
     &   erf_xtiny = 5.96d-8 ,
#elif defined(MIXED)
     &   erf_xtiny = 5.96d-8 ,
#else
     &   erf_xtiny = 1.11d-16 ,
#endif
     &   erf_xbig  = 26.543_ti_p ,
c
c     coefficients for approximation to erf in first interval
c
     &   erf_a=(/ 3.16112374387056560d0,  1.13864154151050156d2,
     &            3.77485237685302021d2,  3.20937758913846947d3,
     &            1.85777706184603153d-1 /) ,
     &   erf_b=(/ 2.36012909523441209d1,  2.44024637934444173d2,
     &            1.28261652607737228d3,  2.84423683343917062d3 /),
c
c     coefficients for approximation to erfc in second interval
c
     &   erf_c=(/ 5.64188496988670089d-1, 8.88314979438837594d0,
     &            6.61191906371416295d1 , 2.98635138197400131d2,
     &            8.81952221241769090d2 , 1.71204761263407058d3,
     &            2.05107837782607147d3 , 1.23033935479799725d3,
     &            2.15311535474403846d-8 /),
     &   erf_d=(/ 1.57449261107098347d1,  1.17693950891312499d2,
     &            5.37181101862009858d2,  1.62138957456669019d3,
     &            3.29079923573345963d3,  4.36261909014324716d3,
     &            3.43936767414372164d3,  1.23033935480374942d3 /),
c
c     coefficients for approximation to erfc in third interval
c
     &   erf_p=(/ 3.05326634961232344d-1, 3.60344899949804439d-1,
     &            1.25781726111229246d-1, 1.60837851487422766d-2,
     &            6.58749161529837803d-4, 1.63153871373020978d-2 /),
     &   erf_q=(/ 2.56852019228982242d0 , 1.87295284992346047d0,
     &            5.27905102951428412d-1, 6.05183413124413191d-2,
     &            2.33520497626869185d-3 /) )
c
c     coefficients for approximation to erfinv in central range
c
      parameter(
     &   a_inv= (/  0.886226899_ti_p, -1.645349621_ti_p,
     &              0.914624893_ti_p, -0.140543331_ti_p /),
     &   b_inv= (/ -2.118377725_ti_p,  1.442710462_ti_p,
     &             -0.329097515_ti_p,  0.012229801_ti_p /),
c
c     coefficients for approximation to erfinv near endpoints
c
     &   c_inv= (/ -1.970840454_ti_p, -1.624906493_ti_p,
     &              3.429567803_ti_p,  1.641345311_ti_p /),
     &   d_inv= (/  3.543889200_ti_p,  1.637067800_ti_p /) )

!$acc declare copyin(erf_xtiny,erf_xbig,erf_sqrpi,erf_a,erf_b,erf_c,
!$acc&        erf_d,erf_p,erf_q,a_inv,b_inv,c_inv,d_inv)

      public :: tinker_erf,tinker_erfc,vderfc_acc,erfcore_acc,erfinv
      contains
c
c     ##################################################################
c     ##                                                              ##
c     ##  function erf      --  evaluate the standard error function  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "erf_acc" computes a numerical approximation to the value of
c     the error function via a Chebyshev approximation on a gpu
c
c
         function tinker_erf (x)
!$acc routine seq
         implicit none
         integer mode
         real(t_p) tinker_erf,x
         real(t_p) result
c
c        compute the error function via Chebyshev fitting
c
         mode = 0
         call erfcore_acc (x,result,mode)
         tinker_erf = result
         end function
c
c
c     ####################################################################
c     ##                                                                ##
c     ##  function erfc      --  evaluate complementary error function  ##
c     ##                                                                ##
c     ####################################################################
c
c
c     "erfc_acc" computes a numerical approximation to the value of the
c     complementary error function via a Chebyshev approximation a gpu
c
c
         function tinker_erfc (x)
!$acc routine seq
         implicit none
         integer mode
         real(t_p) tinker_erfc,x
         real(t_p) result
c
c     get the complementary error function via Chebyshev fitting
c
         mode = 1
         call erfcore_acc (x,result,mode)
         tinker_erfc = result
         end function
c
c
c     ########################################################################
c     ##                                                                    ##
c     ##  subroutine vderfc_acc  --  evaluate complementary error function  ##
c     ##                                                                    ##
c     ########################################################################
c
c
c     "vderfc_acc" computes a numerical approximation to the value of the
c     complementary error function via a Chebyshev approximation a gpu
c
c
         subroutine vderfc_acc (n,nmax,x,res)
!$acc routine vector
         implicit none
         integer,intent(in):: n,nmax
         real(t_p) x(nmax),res(nmax)
         integer i
c
c        get the complementary error function via Chebyshev fitting
c
!$acc loop vector
         do i=1,n
            call erfcore_acc (x(i),res(i),1)
         end do
         end subroutine
c
c
c     #####################################################################
c     ##                                                                 ##
c     ##  subroutine erfcore_acc  --  erf and erfc via Chebyshev approx  ##
c     ##                                                                 ##
c     #####################################################################
c
c
c     "erfcore_acc" evaluates erf(x) or erfc(x) for a real argument x;
c     when called with mode set to 0 it returns erf, a mode of 1
c     returns erfc; uses rational functions that approximate erf(x)
c     and erfc(x) to at least 18 significant decimal digits
c
c
         subroutine erfcore_acc (arg,result,mode)
!$acc routine seq
         implicit none
         integer,intent(in)::mode
         integer i
         real(t_p),intent(in)::arg
         real(t_p),intent(out)::result
         real(t_p) x,y,ysq
         real(t_p) del
         real(t_p) xnum,xden
c
c
c        store the argument and its absolute value
c
         x = arg
         y = abs(x)
c
c        evaluate error function for |x| less than 0.46875
c
         if (y .le. 0.46875_ti_p) then
            ysq = 0.0_ti_p
            if (y .gt. erf_xtiny)  ysq = y * y
            xnum = erf_a(5) * ysq
            xden = ysq
            do i = 1, 3
               xnum = (xnum + erf_a(i)) * ysq
               xden = (xden + erf_b(i)) * ysq
            end do
            result = x * (xnum + erf_a(4)) / (xden + erf_b(4))
            if (mode .ne. 0)  result = 1.0_ti_p - result
c
c        get complementary error function for 0.46875 <= |x| <= 4.0
c
         else if (y .le. 4.0_ti_p) then
            xnum = erf_c(9) * y
            xden = y
            do i = 1, 7
               xnum = (xnum + erf_c(i)) * y
               xden = (xden + erf_d(i)) * y
            end do
            result = (xnum + erf_c(8)) / (xden + erf_d(8))
            ysq = aint(16.0_ti_p*y) / 16.0_ti_p
            del = (y-ysq) * (y+ysq)
c           result = exp(-ysq*ysq) * exp(-del) * result
            result = exp(-ysq*ysq-del) * result
            if (mode .eq. 0) then
               result = 1.0_ti_p - result
               if (x .lt. 0.0_ti_p)  result = -result
            else
               if (x .lt. 0.0_ti_p)  result = 2.0_ti_p - result
            end if
c
c        get complementary error function for |x| greater than 4.0
c
         else
            result = 0.0_ti_p
            if (y .lt. erf_xbig) then
               ysq = 1.0_ti_p / (y * y)
               xnum = erf_p(6) * ysq
               xden = ysq
               do i = 1, 4
                  xnum = (xnum + erf_p(i)) * ysq
                  xden = (xden + erf_q(i)) * ysq
               end do
               result = ysq * (xnum + erf_p(5)) / (xden + erf_q(5))
               result = (erf_sqrpi - result) / y
               ysq = aint(16.0_ti_p*y) / 16.0_ti_p
               del = (y-ysq) * (y+ysq)
               result = exp(-ysq*ysq-del) * result
            end if
            if (mode .eq. 0) then
               result = 1.0_ti_p - result
               if (x .lt. 0.0_ti_p)  result = -result
            else
               if (x .lt. 0.0_ti_p)  result = 2.0_ti_p - result
            end if
         end if
         end subroutine
c
c
c     ################################################################
c     ##                                                            ##
c     ##  function erfinv  --  evaluate the error function inverse  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "erfinv" evaluates the inverse of the error function for
c     an argument in the range (-1,1) using a rational function
c     approximation followed by cycles of Newton-Raphson correction
c
c     adapted from the pseudocode for the Matlab function of the
c     same name; Matlab, version 4.2c, March 1995
c
c
         function erfinv (x)
!$acc routine seq
         use math
         implicit none
         real(t_p):: erfinv
         real(t_p),intent(in ):: x
         real(t_p):: y,z
c
c
c        get an initial estimate for the inverse error function
c
         if (abs(x) .le. 0.7_ti_p) then
            y = x * x
            z = x * (((a_inv(4)*y+a_inv(3))*y+a_inv(2))*y+a_inv(1))
     &            / ((((b_inv(4)*y+b_inv(3))*y+b_inv(2))*y+b_inv(1))
     &              *y+1.0_ti_p)
         else if (x.gt.0.7_ti_p .and. x.lt.1.0_ti_p) then
            y = sqrt(-log((1.0_ti_p-x)/2.0_ti_p))
            z = (((c_inv(4)*y+c_inv(3))*y+c_inv(2))*y+c_inv(1)) 
     &        / ((d_inv(2)*y+d_inv(1))*y+1.0_ti_p)
         else if (x.lt.-0.7_ti_p .and. x.gt.-1.0_ti_p) then
            y = sqrt(-log((1.0_ti_p+x)/2.0_ti_p))
            z = -(((c_inv(4)*y+c_inv(3))*y+c_inv(2))*y+c_inv(1)) 
     &      / ((d_inv(2)*y+d_inv(1))*y+1.0_ti_p)
         else
            print*,' ERFINV  --  Illegal Argument to Inverse',
     &             ' Error Function',x
            stop
         end if
c
c        use two steps of Newton-Raphson correction to increase accuracy
c
         z = z - (tinker_erf(z) - x) / (2.0_ti_p/sqrtpi * exp(-z*z))
         z = z - (tinker_erf(z) - x) / (2.0_ti_p/sqrtpi * exp(-z*z))
         erfinv = z
         end function

      end
