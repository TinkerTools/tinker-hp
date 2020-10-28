c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine search  --  perform unidimensional line search  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "search" is a unidimensional line search based upon parabolic
c     extrapolation and cubic interpolation using both function and
c     gradient values
c
c     variables used by the routine :
c
c     f       function value at the best line search point
c     x       current values of variables during line search
c     g       gradient at the current point during line search
c     p       initial search vector, unchanged by this routine
c     s       scaled search vector at current line search point
c     angle   angle between search and negative gradient vector
c
c     parameters used by the routine :
c
c     stpmin   minimum step length in current line search direction
c     stpmax   maximum step length in current line search direction
c     cappa    stringency of line search (0=tight < cappa < 1=loose)
c     slpmax   projected gradient above which stepsize is reduced
c     angmax   maximum angle between search direction and -gradient
c     intmax   maximum number of interpolations during line search
c
c     status codes upon return :
c
c     Success     normal termination after satisfying "cappa" test
c     ScaleStep   normal termination after a step size rescaling
c     ReSearch    normal termination after a reinterpolation
c     WideAngle   large angle between search direction and -gradient
c     BadIntpln   unsatisfied "cappa" test after two searches
c     IntplnErr   function value increase or serious gradient error
c
c
#include "tinker_precision.h"
      subroutine search (n,f,g,x,p,f_move,angle,ncalls,
     &                          fgvalue,status)
      use domdec
      use linmin
      use math
      use mpi
      use sizes     ,only:tinkerdebug
      use tinheader ,only:ti_p,re_p
      implicit none
      integer i,n,iglob,j
      integer ncalls
      integer intpln
      integer ierr
      real(r_p) fgvalue
      real(r_p) f
      real(r_p) f_move
      real(r_p) s_norm,g_norm
      real(r_p) cosang,angle
      real(r_p) step,parab
      real(r_p) cube,cubstp
      real(r_p) sss,ttt
      real(r_p) f_0,f_1
      real(r_p) f_a,f_b,f_c
      real(r_p) sg_0,sg_1
      real(r_p) sg_a,sg_b,sg_c
      real(r_p) x(*)
      real(r_p) g(*)
      real(r_p) p(*)
      real(r_p), allocatable :: x_0(:)
      real(r_p), allocatable :: s(:)
      logical restart
      character*9 status
      character*9 blank
      external fgvalue
c
c
c     use default parameters for the line search if needed
c
      blank = '         '
      if (stpmin .eq. 0.0_re_p)  stpmin = 1.0d-16
      if (stpmax .eq. 0.0_re_p)  stpmax = 2.0_re_p
      if (cappa .eq. 0.0_re_p)  cappa = 0.1_re_p
      if (slpmax .eq. 0.0_re_p)  slpmax = 10000.0_re_p
      if (angmax .eq. 0.0_re_p)  angmax = 180.0_re_p
      if (intmax .eq. 0)  intmax = 5

      if (rank.eq.0.and.tinkerdebug) write(*,'(x,a)') "***search"
c
c     perform dynamic allocation of some local arrays
c
      allocate (x_0(3*n))
      allocate (s(3*n))
!$acc enter data create(x_0,s) async
c
c     copy the search direction into a new vector
c
!$acc parallel loop collapse(2) async present(glob,p,s)
      do i = 1, nloc
        do j = 1, 3
          iglob = glob(i)
          s(3*(iglob-1)+j) = p(3*(iglob-1)+j)
        end do
      end do
c
c     compute the length of gradient and search direction
c
      g_norm = 0.0_re_p
      s_norm = 0.0_re_p
!$acc parallel loop collapse(2) async present(g,s)
      do i = 1, nloc
        do j = 1, 3
          iglob = glob(i)
          g_norm = g_norm + g(3*(iglob-1)+j)*g(3*(iglob-1)+j)
          s_norm = s_norm + s(3*(iglob-1)+j)*s(3*(iglob-1)+j)
        end do
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,g_norm,1,MPI_RPREC,
     $     MPI_SUM,COMM_TINKER,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,s_norm,1,MPI_RPREC,
     $     MPI_SUM,COMM_TINKER,ierr)
      g_norm = sqrt(g_norm)
      s_norm = sqrt(s_norm)
c
c     store initial function, then normalize the
c     search vector and find projected gradient
c
      f_0 = f
      sg_0 = 0.0_re_p
!$acc parallel loop collapse(2) async present(x_0,x,s,g)
      do i = 1, nloc
        do j = 1, 3
          iglob = glob(i)
          x_0(3*(iglob-1)+j) = x(3*(iglob-1)+j)
          s(3*(iglob-1)+j) = s(3*(iglob-1)+j) / s_norm
          sg_0 = sg_0 + s(3*(iglob-1)+j)*g(3*(iglob-1)+j)
        end do
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,sg_0,1,MPI_RPREC,
     $     MPI_SUM,COMM_TINKER,ierr)
c
c     check the angle between the search direction
c     and the negative gradient vector
c
      cosang = -sg_0 / g_norm
      cosang = min(1.0_re_p,max(-1.0_re_p,cosang))
      angle = radian * acos(cosang)
      if (angle .gt. angmax) then
         status = 'WideAngle'
!$acc exit data delete(x_0,s) async
         deallocate (x_0)
         deallocate (s)
         return
      end if
c
c     set the initial stepsize to the length of the passed
c     search vector, or based on previous function decrease
c
      step = 2.0_re_p * abs(f_move/sg_0)
      step = min(step,s_norm)
      if (step .gt. stpmax)  step = stpmax
      if (step .lt. stpmin)  step = stpmin
c
c     beginning of the parabolic extrapolation procedure
c
   10 continue
      restart = .true.
      intpln = 0
      f_b = f_0
      sg_b = sg_0
c
c     replace last point by latest and take another step
c
   20 continue
      f_a = f_b
      sg_a = sg_b
!$acc parallel loop async collapse(2) present(x,s)
      do i = 1, nloc
        do j = 1, 3
          iglob = glob(i)
          x(3*(iglob-1)+j) = x(3*(iglob-1)+j) + step*s(3*(iglob-1)+j)
        end do
      end do
c
c     send s,x_0,xx vector among the neighbors, the dd is going to change
c
      call sendvecmin(s)
      call sendvecmin(x_0)
      call sendvecmin(x)
c
c     get new function and projected gradient following a step
c
      ncalls = ncalls + 1
      f_b = fgvalue(x,g)
c
      sg_b = 0.0_re_p
!$acc parallel loop collapse(2) async present(s,g)
      do i = 1, nloc
        do j = 1, 3
          iglob = glob(i)
          sg_b = sg_b + s(3*(iglob-1)+j)*g(3*(iglob-1)+j)
        end do
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,sg_b,1,MPI_RPREC,
     $     MPI_SUM,COMM_TINKER,ierr)
c
c     scale stepsize if initial gradient change is too large
c
      if (abs(sg_b/sg_a).ge.slpmax .and. restart) then
!$acc parallel loop collapse(2) async present(x,x_0)
         do i = 1, nloc
           do j = 1, 3
             iglob = glob(i)
             x(3*(iglob-1)+j) = x_0(3*(iglob-1)+j)
           end do
         end do
         step = step / 10.0_re_p
         status = 'ScaleStep'
         goto 10
      end if
      restart = .false.
c
c     return if the gradient is small and function decreases
c
      if (abs(sg_b/sg_0).le.cappa .and. f_b.lt.f_a) then
         f = f_b
         if (status .eq. blank)  status = ' Success '
!$acc exit data delete(x_0,s) async
         deallocate (x_0)
         deallocate (s)
         return
      end if
c
c     interpolate if gradient changes sign or function increases
c
      if (sg_b*sg_a.lt.0.0_re_p .or. f_b.gt.f_a)  goto 30
c
c     if the finite difference curvature is negative double the step;
c     or if  step < parabolic estimate < 4*step  use this estimate,
c     otherwise truncate to step or 4*step, respectively
c
      step = 2.0_re_p * step
      if (sg_b .gt. sg_a) then
         parab = (f_a-f_b) / (sg_b-sg_a)
         if (parab .gt. 2.0_re_p*step)  parab = 2.0_re_p * step
         if (parab .lt. 0.5_re_p*step)  parab = 0.5_re_p * step
         step = parab
      end if
      if (step .gt. stpmax)  step = stpmax
      goto 20
c
c     beginning of the cubic interpolation procedure
c
   30 continue
      intpln = intpln + 1
      sss = 3.0_re_p*(f_b-f_a)/step - sg_a - sg_b
      ttt = sss*sss - sg_a*sg_b
      if (ttt .lt. 0.0_re_p) then
         f = f_b
         status = 'IntplnErr'
!$acc exit data delete(x_0,s) async
         deallocate (x_0)
         deallocate (s)
         return
      end if
      ttt = sqrt(ttt)
      cube = step * (sg_b+ttt+sss)/(sg_b-sg_a+2.0_re_p*ttt)
      if (cube.lt.0.0_re_p .or. cube.gt.step) then
         f = f_b
         status = 'IntplnErr'
!$acc exit data delete(x_0,s) async
         deallocate (x_0)
         deallocate (s)
         return
      end if
!$acc parallel loop collapse(2) async present(x,s)
      do i = 1, nloc
        do j = 1, 3
          iglob = glob(i)
          x(3*(iglob-1)+j) = x(3*(iglob-1)+j) - cube*s(3*(iglob-1)+j)
        end do
      end do
      call sendvecmin(x)
c
c     get new function and gradient, then test for termination
c
      ncalls = ncalls + 1
      f_c = fgvalue (x,g)
      sg_c = 0.0_re_p
!$acc parallel loop collapse(2) async present(glob,s,g)
      do i = 1, nloc
        do j = 1, 3
          iglob = glob(i)
          sg_c = sg_c + s(3*(iglob-1)+j)*g(3*(iglob-1)+j)
        end do
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,sg_c,1,MPI_RPREC,
     $     MPI_SUM,COMM_TINKER,ierr)
      if (abs(sg_c/sg_0) .le. cappa) then
         f = f_c
         if (status .eq. blank)  status = ' Success '
!$acc exit data delete(x_0,s) async
         deallocate (x_0)
         deallocate (s)
         return
      end if
c
c     get the next pair of bracketing points by replacing one
c     of the current brackets with the interpolated point
c
      if (f_c.le.f_a .or. f_c.le.f_b) then
         cubstp = min(abs(cube),abs(step-cube))
         if (cubstp.ge.stpmin .and. intpln.lt.intmax) then
c
c     if the current brackets have slopes of opposite sign,
c     then substitute the interpolated point for the bracket
c     point with slope of same sign as the interpolated point
c
            if (sg_a*sg_b .lt. 0.0_re_p) then
               if (sg_a*sg_c .lt. 0.0_re_p) then
                  f_b = f_c
                  sg_b = sg_c
                  step = step - cube
               else
                  f_a = f_c
                  sg_a = sg_c
                  step = cube
!$acc parallel loop collapse(2) async present(x,s,glob)
                  do i = 1, nloc
                     do j = 1, 3
                       iglob = glob(i)
                       x(3*(iglob-1)+j) = x(3*(iglob-1)+j) + 
     $                               cube*s(3*(iglob-1)+j)
                     end do
                  end do
               end if
c
c     if current brackets have slope of same sign, then replace
c     the far bracket if the interpolated point has a slope of
c     the opposite sign or a lower function value than the near
c     bracket, otherwise replace the far bracket point
c
            else
               if (sg_a*sg_c.lt.0.0_re_p .or. f_a.le.f_c) then
                  f_b  = f_c
                  sg_b = sg_c
                  step = step - cube
               else
                  f_a  = f_c
                  sg_a = sg_c
                  step = cube
!$acc parallel loop collapse(2) async present(x,s,glob)
                  do i = 1, nloc
                    do j = 1, 3
                      iglob = glob(i)
                      x(3*(iglob-1)+j) = x(3*(iglob-1)+j) + 
     $                              cube*s(3*(iglob-1)+j)
                    end do
                  end do
               end if
            end if
            goto 30
         end if
      end if
c
c     interpolation has failed, reset to best current point
c
      f_1 = min(f_a,f_b,f_c)
      if (f_1 .eq. f_a) then
         sg_1 = sg_a
!$acc parallel loop collapse(2) async present(x,s,glob)
         do i = 1, nloc
           do j = 1, 3
             iglob = glob(i)
             x(3*(iglob-1)+j) = x(3*(iglob-1)+j) + 
     $          (cube-step)*s(3*(iglob-1)+j)
           end do
         end do
      else if (f_1 .eq. f_b) then
         sg_1 = sg_b
!$acc parallel loop collapse(2) async present(x,s,glob)
         do i = 1, nloc
           do j = 1, 3
             iglob = glob(i)
             x(3*(iglob-1)+j) = x(3*(iglob-1)+j) + cube*s(3*(iglob-1)+j)
           end do
         end do
      else if (f_1 .eq. f_c) then
         sg_1 = sg_c
      end if
c
c     try to restart from best point with smaller stepsize
c
      if (f_1 .gt. f_0) then
         ncalls = ncalls + 1
         f = fgvalue (x,g)
         status = 'IntplnErr'
!$acc exit data delete(x_0,s) async
         deallocate (x_0)
         deallocate (s)
         return
      end if
      f_0  = f_1
      sg_0 = sg_1
      if (sg_1 .gt. 0.0_re_p) then
!$acc parallel loop collapse(2) async present(s,glob)
         do i = 1, nloc
           do j = 1, 3
             iglob = glob(i)
             s(3*(iglob-1)+j) = -s(3*(iglob-1)+j)
           end do
         end do
         sg_0 = -sg_1
      end if
      step = max(cube,step-cube) / 10.0_re_p
      if (step .lt. stpmin)  step = stpmin
c
c     if already restarted once, then return with best point
c
      if (status .eq. ' ReSearch') then
         ncalls = ncalls + 1
         f = fgvalue (x,g)
         status = 'BadIntpln'
!$acc exit data delete(x_0,s) async
         deallocate (x_0)
         deallocate (s)
         return
      else
         status = ' ReSearch'
         goto 10
      end if
      end
