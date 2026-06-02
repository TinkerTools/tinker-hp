!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!> @brief 
!> computes prime factor decomposition of an integer
!> @param[in]: integer of interest
!> @param[d]: prime factors
!> @param[i]: number of prime factors
subroutine prime(n,d,i)
   implicit none
   integer n,i,div,next,rest
   integer d(*)
   i = 1
   div = 2
   next = 3
   rest = n
   do while (rest.ne.1)
      do while (mod(rest,div).eq.0)
         d(i) = div
         i = i + 1
         rest = rest / div
      end do
      div = next
      next = next + 2
   end do
   return
end

