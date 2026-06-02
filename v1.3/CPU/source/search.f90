!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #################################################################
!     ##                                                             ##
!     ##  subroutine search  --  perform unidimensional line search  ##
!     ##                                                             ##
!     #################################################################
!
!
!     "search" is a unidimensional line search based upon parabolic
!     extrapolation and cubic interpolation using both function and
!     gradient values
!
!     variables used by the routine :
!
!     f       function value at the best line search point
!     x       current values of variables during line search
!     g       gradient at the current point during line search
!     p       initial search vector, unchanged by this routine
!     s       scaled search vector at current line search point
!     angle   angle between search and negative gradient vector
!
!     parameters used by the routine :
!
!     stpmin   minimum step length in current line search direction
!     stpmax   maximum step length in current line search direction
!     cappa    stringency of line search (0=tight < cappa < 1=loose)
!     slpmax   projected gradient above which stepsize is reduced
!     angmax   maximum angle between search direction and -gradient
!     intmax   maximum number of interpolations during line search
!
!     status codes upon return :
!
!     Success     normal termination after satisfying "cappa" test
!     ScaleStep   normal termination after a step size rescaling
!     ReSearch    normal termination after a reinterpolation
!     WideAngle   large angle between search direction and -gradient
!     BadIntpln   unsatisfied "cappa" test after two searches
!     IntplnErr   function value increase or serious gradient error
!
!
!> @brief 
!> is a unidimensional line search based upon parabolic
!> extrapolation and cubic interpolation using both function and
!> gradient values
!> @param no params
!> @params[in]     f:       function value at the best line search point
!> @params[in]     x:       current values of variables during line search
!> @params[in]     g:       gradient at the current point during line search
!> @params[in]     p:       initial search vector, unchanged by this routine
!> @params[in]     s:       scaled search vector at current line search point
!> @params[in]     angle:   angle between search and negative gradient vector
!> @params[in]     ncalls:   number of calls to the routine
!> @params[in]     fgvalue:  current value of the energy and gradients
!> @params[in]     status:  status of the line search
subroutine search (n,f,g,x,p,f_move,angle,ncalls,&
&fgvalue,status)
   use domdec
   use linmin
   use math
   use mpi
   implicit none
   integer i,n,iglob,j
   integer ncalls
   integer intpln
   integer ierr
   real*8 fgvalue
   real*8 f,f_move
   real*8 s_norm,g_norm
   real*8 cosang,angle
   real*8 step,parab
   real*8 cube,cubstp
   real*8 sss,ttt
   real*8 f_0,f_1
   real*8 f_a,f_b,f_c
   real*8 sg_0,sg_1
   real*8 sg_a,sg_b,sg_c
   real*8 x(*)
   real*8 g(*)
   real*8 p(*)
   real*8, allocatable :: x_0(:)
   real*8, allocatable :: s(:)
   logical restart
   character*9 status
   character*9 blank
   external fgvalue
!
!
!     use default parameters for the line search if needed
!
   blank = '         '
   if (stpmin .eq. 0.0d0)  stpmin = 1.0d-16
   if (stpmax .eq. 0.0d0)  stpmax = 2.0d0
   if (cappa .eq. 0.0d0)  cappa = 0.1d0
   if (slpmax .eq. 0.0d0)  slpmax = 10000.0d0
   if (angmax .eq. 0.0d0)  angmax = 180.0d0
   if (intmax .eq. 0)  intmax = 5
!
!     perform dynamic allocation of some local arrays
!
   allocate (x_0(3*n))
   allocate (s(3*n))
!
!     copy the search direction into a new vector
!
   do i = 1, nloc
      iglob = glob(i)
      do j = 1, 3
         s(3*(iglob-1)+j) = p(3*(iglob-1)+j)
      end do
   end do
!
!     compute the length of gradient and search direction
!
   g_norm = 0.0d0
   s_norm = 0.0d0
   do i = 1, nloc
      iglob = glob(i)
      do j = 1, 3
         g_norm = g_norm + g(3*(iglob-1)+j)*g(3*(iglob-1)+j)
         s_norm = s_norm + s(3*(iglob-1)+j)*s(3*(iglob-1)+j)
      end do
   end do
   call MPI_ALLREDUCE(MPI_IN_PLACE,g_norm,1,MPI_REAL8,&
   &MPI_SUM,COMM_TINKER,ierr)
   call MPI_ALLREDUCE(MPI_IN_PLACE,s_norm,1,MPI_REAL8,&
   &MPI_SUM,COMM_TINKER,ierr)
   g_norm = sqrt(g_norm)
   s_norm = sqrt(s_norm)
!
!     store initial function, then normalize the
!     search vector and find projected gradient
!
   f_0 = f
   sg_0 = 0.0d0
   do i = 1, nloc
      iglob = glob(i)
      do j = 1, 3
         x_0(3*(iglob-1)+j) = x(3*(iglob-1)+j)
         s(3*(iglob-1)+j) = s(3*(iglob-1)+j) / s_norm
         sg_0 = sg_0 + s(3*(iglob-1)+j)*g(3*(iglob-1)+j)
      end do
   end do
   call MPI_ALLREDUCE(MPI_IN_PLACE,sg_0,1,MPI_REAL8,&
   &MPI_SUM,COMM_TINKER,ierr)
!
!     check the angle between the search direction
!     and the negative gradient vector
!
   cosang = -sg_0 / g_norm
   cosang = min(1.0d0,max(-1.0d0,cosang))
   angle = radian * acos(cosang)
   if (angle .gt. angmax) then
      status = 'WideAngle'
      deallocate (x_0)
      deallocate (s)
      return
   end if
!
!     set the initial stepsize to the length of the passed
!     search vector, or based on previous function decrease
!
   step = 2.0d0 * abs(f_move/sg_0)
   step = min(step,s_norm)
   if (step .gt. stpmax)  step = stpmax
   if (step .lt. stpmin)  step = stpmin
!
!     beginning of the parabolic extrapolation procedure
!
10 continue
   restart = .true.
   intpln = 0
   f_b = f_0
   sg_b = sg_0
!
!     replace last point by latest and take another step
!
20 continue
   f_a = f_b
   sg_a = sg_b
   do i = 1, nloc
      iglob = glob(i)
      do j = 1, 3
         x(3*(iglob-1)+j) = x(3*(iglob-1)+j) + step*s(3*(iglob-1)+j)
      end do
   end do
!
!     send s,x_0,x,g and p vector among the neighbors, the dd is going to change
!
   call sendvecmin(s)
   call sendvecmin(x_0)
   call sendvecmin(x)
   call sendvecmin(g)
   call sendvecmin(p)
!
!     get new function and projected gradient following a step
!
   ncalls = ncalls + 1
   f_b = fgvalue(x,g)
!
   sg_b = 0.0d0
   do i = 1, nloc
      iglob = glob(i)
      do j = 1, 3
         sg_b = sg_b + s(3*(iglob-1)+j)*g(3*(iglob-1)+j)
      end do
   end do
   call MPI_ALLREDUCE(MPI_IN_PLACE,sg_b,1,MPI_REAL8,&
   &MPI_SUM,COMM_TINKER,ierr)
!
!     scale stepsize if initial gradient change is too large
!
   if (abs(sg_b/sg_a).ge.slpmax .and. restart) then
      do i = 1, nloc
         iglob = glob(i)
         do j = 1, 3
            x(3*(iglob-1)+j) = x_0(3*(iglob-1)+j)
         end do
      end do
      step = step / 10.0d0
      status = 'ScaleStep'
      goto 10
   end if
   restart = .false.
!
!     return if the gradient is small and function decreases
!
   if (abs(sg_b/sg_0).le.cappa .and. f_b.lt.f_a) then
      f = f_b
      if (status .eq. blank)  status = ' Success '
      deallocate (x_0)
      deallocate (s)
      return
   end if
!
!     interpolate if gradient changes sign or function increases
!
   if (sg_b*sg_a.lt.0.0d0 .or. f_b.gt.f_a)  goto 30
!
!     if the finite difference curvature is negative double the step;
!     or if  step < parabolic estimate < 4*step  use this estimate,
!     otherwise truncate to step or 4*step, respectively
!
   step = 2.0d0 * step
   if (sg_b .gt. sg_a) then
      parab = (f_a-f_b) / (sg_b-sg_a)
      if (parab .gt. 2.0d0*step)  parab = 2.0d0 * step
      if (parab .lt. 0.5d0*step)  parab = 0.5d0 * step
      step = parab
   end if
   if (step .gt. stpmax)  step = stpmax
   goto 20
!
!     beginning of the cubic interpolation procedure
!
30 continue
   intpln = intpln + 1
   sss = 3.0d0*(f_b-f_a)/step - sg_a - sg_b
   ttt = sss*sss - sg_a*sg_b
   if (ttt .lt. 0.0d0) then
      f = f_b
      status = 'IntplnErr'
      deallocate (x_0)
      deallocate (s)
      return
   end if
   ttt = sqrt(ttt)
   cube = step * (sg_b+ttt+sss)/(sg_b-sg_a+2.0d0*ttt)
   if (cube.lt.0.0d0 .or. cube.gt.step) then
      f = f_b
      status = 'IntplnErr'
      deallocate (x_0)
      deallocate (s)
      return
   end if
   do i = 1, nloc
      iglob = glob(i)
      do j = 1, 3
         x(3*(iglob-1)+j) = x(3*(iglob-1)+j) - cube*s(3*(iglob-1)+j)
      end do
   end do
!
!     send s,x_0,x,g and p vector among the neighbors, the dd is going to change
!
   call sendvecmin(s)
   call sendvecmin(x_0)
   call sendvecmin(x)
   call sendvecmin(g)
   call sendvecmin(p)
!
!     get new function and gradient, then test for termination
!
   ncalls = ncalls + 1
   f_c = fgvalue (x,g)
   sg_c = 0.0d0
   do i = 1, nloc
      iglob = glob(i)
      do j = 1, 3
         sg_c = sg_c + s(3*(iglob-1)+j)*g(3*(iglob-1)+j)
      end do
   end do
   call MPI_ALLREDUCE(MPI_IN_PLACE,sg_c,1,MPI_REAL8,&
   &MPI_SUM,COMM_TINKER,ierr)
   if (abs(sg_c/sg_0) .le. cappa) then
      f = f_c
      if (status .eq. blank)  status = ' Success '
      deallocate (x_0)
      deallocate (s)
      return
   end if
!
!     get the next pair of bracketing points by replacing one
!     of the current brackets with the interpolated point
!
   if (f_c.le.f_a .or. f_c.le.f_b) then
      cubstp = min(abs(cube),abs(step-cube))
      if (cubstp.ge.stpmin .and. intpln.lt.intmax) then
!
!     if the current brackets have slopes of opposite sign,
!     then substitute the interpolated point for the bracket
!     point with slope of same sign as the interpolated point
!
         if (sg_a*sg_b .lt. 0.0d0) then
            if (sg_a*sg_c .lt. 0.0d0) then
               f_b = f_c
               sg_b = sg_c
               step = step - cube
            else
               f_a = f_c
               sg_a = sg_c
               step = cube
               do i = 1, nloc
                  iglob = glob(i)
                  do j = 1, 3
                     x(3*(iglob-1)+j) = x(3*(iglob-1)+j) +&
                     &cube*s(3*(iglob-1)+j)
                  end do
               end do
            end if
!
!     if current brackets have slope of same sign, then replace
!     the far bracket if the interpolated point has a slope of
!     the opposite sign or a lower function value than the near
!     bracket, otherwise replace the far bracket point
!
         else
            if (sg_a*sg_c.lt.0.0d0 .or. f_a.le.f_c) then
               f_b = f_c
               sg_b = sg_c
               step = step - cube
            else
               f_a = f_c
               sg_a = sg_c
               step = cube
               do i = 1, nloc
                  iglob = glob(i)
                  do j = 1, 3
                     x(3*(iglob-1)+j) = x(3*(iglob-1)+j) +&
                     &cube*s(3*(iglob-1)+j)
                  end do
               end do
            end if
         end if
         goto 30
      end if
   end if
!
!     interpolation has failed, reset to best current point
!
   f_1 = min(f_a,f_b,f_c)
   if (f_1 .eq. f_a) then
      sg_1 = sg_a
      do i = 1, nloc
         iglob = glob(i)
         do j = 1, 3
            x(3*(iglob-1)+j) = x(3*(iglob-1)+j) +&
            &(cube-step)*s(3*(iglob-1)+j)
         end do
      end do
   else if (f_1 .eq. f_b) then
      sg_1 = sg_b
      do i = 1, nloc
         iglob = glob(i)
         do j = 1, 3
            x(3*(iglob-1)+j) = x(3*(iglob-1)+j) + cube*s(3*(iglob-1)+j)
         end do
      end do
   else if (f_1 .eq. f_c) then
      sg_1 = sg_c
   end if
!
!     try to restart from best point with smaller stepsize
!
   if (f_1 .gt. f_0) then
!
!     send s,x_0,x,g and p vector among the neighbors, the dd is going to change
!
      call sendvecmin(s)
      call sendvecmin(x_0)
      call sendvecmin(x)
      call sendvecmin(g)
      call sendvecmin(p)
      ncalls = ncalls + 1
      f = fgvalue (x,g)
      status = 'IntplnErr'
      deallocate (x_0)
      deallocate (s)
      return
   end if
   f_0 = f_1
   sg_0 = sg_1
   if (sg_1 .gt. 0.0d0) then
      do i = 1, nloc
         iglob = glob(i)
         do j = 1, 3
            s(3*(iglob-1)+j) = -s(3*(iglob-1)+j)
         end do
      end do
      sg_0 = -sg_1
   end if
   step = max(cube,step-cube) / 10.0d0
   if (step .lt. stpmin)  step = stpmin
!
!     if already restarted once, then return with best point
!
   if (status .eq. ' ReSearch') then
!
!     send s,x_0,x,g and p vector among the neighbors, the dd is going to change
!
      call sendvecmin(s)
      call sendvecmin(x_0)
      call sendvecmin(x)
      call sendvecmin(g)
      call sendvecmin(p)
      ncalls = ncalls + 1
      f = fgvalue (x,g)
      status = 'BadIntpln'
      deallocate (x_0)
      deallocate (s)
      return
   else
      status = ' ReSearch'
      goto 10
   end if
end
