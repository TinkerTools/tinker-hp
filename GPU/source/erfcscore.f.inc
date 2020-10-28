         subroutine erfcore_inl (arg,result,mode)
!$acc routine seq
         implicit none
         integer,intent(in)::mode
         real*4 ,intent(in)::arg
         real*4 ,intent(out)::result
         integer i
         real*4 x,y,ysq
         real*4 del
         real*4 xnum,xden
c
c        store the argument and its absolute value
c
         x = arg
         y = abs(x)
c
c        evaluate error function for |x| less than 0.46875
c
         if (y .le. 0.46875) then
            ysq = 0.0
            if (y .gt. erf_xtiny)  ysq = y * y
            xnum = erf_a(5) * ysq
            xden = ysq
            do i = 1, 3
               xnum = (xnum + erf_a(i)) * ysq
               xden = (xden + erf_b(i)) * ysq
            end do
            result = x * (xnum + erf_a(4)) / (xden + erf_b(4))
            if (mode .ne. 0)  result = 1.0 - result
c
c        get complementary error function for 0.46875 <= |x| <= 4.0
c
         else if (y .le. 4.0) then
            xnum = erf_c(9) * y
            xden = y
            do i = 1, 7
               xnum = (xnum + erf_c(i)) * y
               xden = (xden + erf_d(i)) * y
            end do
            result = (xnum + erf_c(8)) / (xden + erf_d(8))
            ysq = aint(16.0*y) / 16.0
            del = (y-ysq) * (y+ysq)
c           result = exp(-ysq*ysq) * exp(-del) * result
            result = exp(-ysq*ysq-del) * result
            if (mode .eq. 0) then
               result = 1.0 - result
               if (x .lt. 0.0)  result = -result
            else
               if (x .lt. 0.0)  result = 2.0 - result
            end if
c
c        get complementary error function for |x| greater than 4.0
c
         else
            result = 0.0
            if (y .lt. erf_xbig) then
               ysq = 1.0 / (y * y)
               xnum = erf_p(6) * ysq
               xden = ysq
               do i = 1, 4
                  xnum = (xnum + erf_p(i)) * ysq
                  xden = (xden + erf_q(i)) * ysq
               end do
               result = ysq * (xnum + erf_p(5)) / (xden + erf_q(5))
               result = (erf_sqrpi - result) / y
               ysq = aint(16.0*y) / 16.0
               del = (y-ysq) * (y+ysq)
               result = exp(-ysq*ysq-del) * result
            end if
            if (mode .eq. 0) then
               result = 1.0 - result
               if (x .lt. 0.0)  result = -result
            else
               if (x .lt. 0.0)  result = 2.0 - result
            end if
         end if
         end subroutine
