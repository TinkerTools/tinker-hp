c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine lbfgs  --  limited memory BFGS optimization  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "lbfgs" is a limited memory BFGS quasi-newton nonlinear
c     optimization routine
c
c     literature references:
c
c     J. Nocedal, "Updating Quasi-Newton Matrices with Limited
c     Storage", Mathematics of Computation, 35, 773-782 (1980)
c
c     D. C. Lui and J. Nocedal, "On the Limited Memory BFGS Method
c     for Large Scale Optimization", Mathematical Programming,
c     45, 503-528 (1989)
c
c     J. Nocedal and S. J. Wright, "Numerical Optimization",
c     Springer-Verlag, New York, 1999, Section 9.1
c
c     variables and parameters:
c
c     nvar      number of parameters in the objective function
c     x0        contains starting point upon input, upon return
c                 contains the best point found
c     minimum   during optimization contains best current function
c                 value; returns final best function value
c     grdmin    normal exit if rms gradient gets below this value
c     ncalls    total number of function/gradient evaluations
c
c     required external routines:
c
c     fgvalue    function to evaluate function and gradient values
c     optsave    subroutine to write out info about current status
c
c
#include "tinker_precision.h"
      subroutine lbfgs (n,x0,minimum,grdmin,fgvalue,optsave)
      use domdec
      use inform
      use iounit
      use keys
      use linmin
      use math
      use minima
      use interfaces,only: search
      use output
      use scales
      use usage
      use tinheader ,only:ti_p,re_p
      use mpi
      implicit none
      integer maxsav
      parameter (maxsav=20)
      integer i,j,k,l,m,ierr,n
      integer next
      integer msav,muse
      integer niter,ncalls
      integer nerr,maxerr
      integer iglob
      real(r_p) f,f_old,fgvalue
      real(r_p) f_move,x_move
      real(r_p) g_norm,g_rms
      real(r_p) minimum
      real(r_p) temp0,temp1,temp2
      real(r_p) grdmin
      real(r_p) angle,beta
      real(r_p) rms
      real(r_p) ys,yy,gamma
      real(r_p) x0(*)
      real(t_p) rho(maxsav)
      real(t_p) alpha(maxsav)
      real(r_p), allocatable :: x_old(:)
      real(r_p), allocatable :: g(:)
      real(r_p), allocatable :: g_old(:)
      real(r_p), allocatable :: p(:)
      real(r_p), allocatable :: q(:)
      real(r_p), allocatable :: r(:)
      real(r_p), allocatable :: h0(:)
      real(r_p), allocatable :: s(:,:)
      real(r_p), allocatable :: y(:,:)
      logical done
      character*9 blank,status
      character*20 keyword
      character*120 record
      character*120 string
      external fgvalue,optsave
c
c
c     initialize some values to be used below
cc
c      if (3*n .gt. maxvar) then
c         if (rank.eq.0) then
c         write (iout,10)
c   10    format (/,' LBFGS  --  Too many Parameters,',
c     &              ' Increase the Value of MAXVAR')
c         return
c         end if
c      end if
      ncalls = 0
c      rms = sqrt(real(nvar,t_p))
      rms = sqrt(real(3*n,t_p))
      if (coordtype .eq. 'CARTESIAN') then
         rms = rms / sqrt(3.0_re_p)
      end if
      blank = '         '
      done = .false.
      nerr = 0
      maxerr = 2
      if (rank.eq.0.and.tinkerdebug) write(*,'(x,a)') "***lbfgs"
c
c     set default values for variable scale factors
c
      if (.not. set_scale) then
         do i = 1, 3*n
            if (scale(i) .eq. 0.0_re_p)  scale(i) = 1.0_re_p
         end do
      end if
c
c     set default parameters for the optimization
c
      msav = min(3*n,maxsav)
      if (fctmin .eq. 0.0_re_p)  fctmin = -1000000000.0_re_p
      if (maxiter .eq. 0)  maxiter = 1000000
      if (nextiter .eq. 0)  nextiter = 1
      if (iprint .lt. 0)  iprint = 1
      if (iwrite .lt. 0)  iwrite = 1
c
c     set default parameters for the line search
c
      if (stpmax .eq. 0.0_re_p)  stpmax = 5.0_re_p
      stpmin = 1.0d-16
      cappa = 0.9_re_p
      slpmax = 10000.0_re_p
      angmax = 180.0_re_p
      intmax = 5
c
c     search the keywords for optimization parameters
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:120)
         if (keyword(1:14) .eq. 'LBFGS-VECTORS ') then
            read (string,*,err=20,end=20)  msav
         else if (keyword(1:17) .eq. 'STEEPEST-DESCENT ') then
            msav = 0
         else if (keyword(1:7) .eq. 'FCTMIN ') then
            read (string,*,err=20,end=20)  fctmin
         else if (keyword(1:8) .eq. 'MAXITER ') then
            read (string,*,err=20,end=20)  maxiter
         else if (keyword(1:8) .eq. 'STEPMAX ') then
            read (string,*,err=20,end=20)  stpmax
         else if (keyword(1:8) .eq. 'STEPMIN ') then
            read (string,*,err=20,end=20)  stpmin
         else if (keyword(1:6) .eq. 'CAPPA ') then
            read (string,*,err=20,end=20)  cappa
         else if (keyword(1:9) .eq. 'SLOPEMAX ') then
            read (string,*,err=20,end=20)  slpmax
         else if (keyword(1:7) .eq. 'ANGMAX ') then
            read (string,*,err=20,end=20)  angmax
         else if (keyword(1:7) .eq. 'INTMAX ') then
            read (string,*,err=20,end=20)  intmax
         else if (keyword(1:7) .eq. 'IWRITE ') then
            read (string,*,err=20,end=20)  iwrite
         end if
   20    continue
      end do
c
c     check the number of saved correction vectors
c
      if (msav.lt.0 .or. msav.gt.min(3*n,maxsav)) then
         msav = min(3*n,maxsav)
         write (iout,30)  msav
   30    format (/,' LBFGS  --  Number of Saved Vectors Set to',
     &              ' the Maximum of',i5)
      end if
c
c     print header information about the optimization method
c
      if (iprint .gt. 0 .and. rank.eq.0) then
         if (msav .eq. 0) then
            write (iout,40)
   40       format (/,' Steepest Descent Gradient Optimization :')
            write (iout,50)
   50       format (/,' SD Iter     F Value      G RMS     F Move',
     &                 '    X Move   Angle  FG Call  Comment',/)
         else if (rank.eq.0) then
            write (iout,60)
   60       format (/,' Limited Memory BFGS Quasi-Newton',
     &                 ' Optimization :')
            write (iout,70)
   70       format (/,' QN Iter     F Value      G RMS     F Move',
     &                 '    X Move   Angle  FG Call  Comment',/)
         end if
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (x_old(3*n))
      allocate (g(3*n))
      allocate (g_old(3*n))
      allocate (p(3*n))
      allocate (q(3*n))
      allocate (r(3*n))
      allocate (h0(3*n))
      allocate (s(3*n,maxsav))
      allocate (y(3*n,maxsav))
!$acc enter data create(beta) async
!$acc data create(x_old,g,g_old,p,q,r,h0,s,y,alpha,rho)
!$acc&     present(x0,scale) async
c
c     evaluate the function and get the initial gradient
c
      niter   = nextiter - 1
      maxiter = niter + maxiter
      ncalls  = ncalls + 1
c
      f = fgvalue (x0,g)
c
c
      f_old  = f
      m      = 0
      gamma  = 1.0_re_p
      g_norm = 0.0_re_p
      g_rms  = 0.0_re_p
!$acc parallel loop async default(present)
      do i = 1, nloc
         do j = 1, 3
           iglob  = glob(i)
           g_norm = g_norm + g(3*(iglob-1)+j)*g(3*(iglob-1)+j)
           g_rms  = g_rms + (g(3*(iglob-1)+j)*scale(3*(iglob-1)+j))**2
         end do
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,g_norm,1,MPI_RPREC,MPI_SUM,
     $     COMM_TINKER,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,g_rms,1,MPI_RPREC,MPI_SUM,
     $     COMM_TINKER,ierr)
      g_norm  = sqrt(g_norm)
      g_rms   = sqrt(g_rms) / rms
      f_move  = 0.5_re_p * stpmax * g_norm
c
c     print initial information prior to first iteration
c
      if (iprint .gt. 0 .and. rank.eq.0) then
         if (f.lt.1.0d8 .and. f.gt.-1.0d7 .and. g_rms.lt.1.0d5) then
            write (iout,80)  niter,f,g_rms,ncalls
   80       format (i6,f14.4,f11.4,29x,i7)
         else
            write (iout,90)  niter,f,g_rms,ncalls
   90       format (i6,d14.4,d11.4,29x,i7)
         end if
      end if
c
c     write initial intermediate prior to first iteration
c
      if (iwrite .gt. 0 .and. rank.eq.0)  call optsave (niter,x0)
c
c     tests of the various termination criteria
c
      if (niter .ge. maxiter) then
         status = 'IterLimit'
         done = .true.
      end if
      if (f .le. fctmin) then
         status = 'SmallFct '
         done = .true.
      end if
      if (g_rms .le. grdmin) then
         status = 'SmallGrad'
         done = .true.
      end if
c
c     start of a new limited memory BFGS iteration
c
      do while (.not. done)
         niter = niter + 1
         muse  = min(niter-1,msav)
         m     = m + 1
         if (m .gt. msav)  m = 1
c
c     estimate Hessian diagonal and compute the Hg product
c
!$acc parallel loop collapse(2) async
         do i = 1, nloc
            do j = 1, 3
              iglob = glob(i)
              h0(3*(iglob-1)+j) = gamma
              q(3*(iglob-1)+j)  = g(3*(iglob-1)+j)
            end do
         end do
         k = m
         do j = 1, muse
            k = k - 1
            if (k .eq. 0)  k = msav
!$acc serial async
            alpha(k) = 0.0_re_p
!$acc end serial
!$acc parallel loop gang vector async
            do i = 1, nloc
              iglob = glob(i)
              do l = 1, 3
               temp0 = s(3*(iglob-1)+l,k)*q(3*(iglob-1)+l)
!$acc atomic
               alpha(k) = alpha(k) + temp0
              end do
            end do
            if (nproc.gt.1) then
!$acc wait
!$acc host_data use_device(alpha)
               call MPI_ALLREDUCE(MPI_IN_PLACE,alpha(k),1,MPI_RPREC,
     $              MPI_SUM,COMM_TINKER,ierr)
!$acc end host_data
            end if
!$acc serial async
            alpha(k) = alpha(k) * rho(k)
!$acc end serial
!$acc parallel loop collapse(2) async
            do i = 1, nloc
              do l = 1, 3
                iglob = glob(i)
                q(3*(iglob-1)+l) = q(3*(iglob-1)+l) - 
     $                  alpha(k)*y(3*(iglob-1)+l,k)
              end do
            end do
         end do
!$acc parallel loop collapse(2) async
         do i = 1, nloc
           do j = 1, 3
             iglob = glob(i)
             r(3*(iglob-1)+j) = h0(3*(iglob-1)+j) * q(3*(iglob-1)+j)
           end do
         end do
         do j = 1, muse
!$acc serial async present(beta)
            beta = 0.0_re_p
!$acc end serial
!$acc parallel loop collapse(2) async present(beta)
            do i = 1, nloc
              do l = 1, 3
                iglob = glob(i)
                beta = beta + y(3*(iglob-1)+l,k)*r(3*(iglob-1)+l)
              end do
            end do
            if (nproc.gt.1) then
!$acc wait
!$acc host_data use_device(beta)
            call MPI_ALLREDUCE(MPI_IN_PLACE,beta,1,MPI_RPREC,
     $           MPI_SUM,COMM_TINKER,ierr)
!$acc end host_data
            end if
!$acc serial async present(beta)
            alpha(k) = alpha(k) * rho(k)
            beta = beta * rho(k)
!$acc end serial
!$acc parallel loop collapse(2) async present(beta)
            do i = 1, nloc
              do l = 1, 3
                iglob = glob(i)
                r(3*(iglob-1)+l) = r(3*(iglob-1)+l) + 
     $             s(3*(iglob-1)+l,k)*(alpha(k)-beta)
              end do
            end do
            k = k + 1
            if (k .gt. msav)  k = 1
         end do
c
c     set search direction and store current point and gradient
c
!$acc parallel loop collapse(2) async
         do i = 1, nloc
           do j = 1, 3
             iglob = glob(i)
                 p(3*(iglob-1)+j) = -r(3*(iglob-1)+j)
             x_old(3*(iglob-1)+j) = x0(3*(iglob-1)+j)
             g_old(3*(iglob-1)+j) = g (3*(iglob-1)+j)
           end do
         end do
         call sendvecmin(x_old)
         call sendvecmin(g_old)
c
c     perform line search along the new conjugate direction
c
         status = blank
         call search (n,f,g,x0,p,f_move,angle,ncalls,fgvalue,status)
c
c     update variables based on results of this iteration
c
         ys = 0.0_re_p
         yy = 0.0_re_p
!$acc parallel loop collapse(2) async
         do i = 1, nloc
           do j = 1, 3
             iglob =  glob(i)
             s(3*(iglob-1)+j,m) = x0(3*(iglob-1)+j) - 
     $              x_old(3*(iglob-1)+j)
             y(3*(iglob-1)+j,m) = g(3*(iglob-1)+j) - 
     $              g_old(3*(iglob-1)+j)
             ys = ys + y(3*(iglob-1)+j,m)*s(3*(iglob-1)+j,m)
             yy = yy + y(3*(iglob-1)+j,m)*y(3*(iglob-1)+j,m)
           end do
         end do
         if (nproc.gt.1) then
            call MPI_ALLREDUCE(MPI_IN_PLACE,ys,1,MPI_RPREC,
     $           MPI_SUM,COMM_TINKER,ierr)
            call MPI_ALLREDUCE(MPI_IN_PLACE,yy,1,MPI_RPREC,
     $           MPI_SUM,COMM_TINKER,ierr)
         end if
         gamma = abs(ys/yy)
!$acc serial async
         rho(m) = 1.0_re_p / ys
!$acc end serial
c
c     get the sizes of the moves made during this iteration
c
         f_move = f_old - f
         f_old  = f
         x_move = 0.0_re_p
!$acc parallel loop collapse(2) async
         do i = 1, nloc
           do j = 1, 3
             iglob = glob(i)
             if (use(iglob)) then
              x_move = x_move +((x0(3*(iglob-1)+j)-x_old(3*(iglob-1)+j))
     $                 /scale(3*(iglob-1)+j))**2
             end if
           end do
         end do
         if (nproc.gt.1) then
            call MPI_ALLREDUCE(MPI_IN_PLACE,x_move,1,MPI_RPREC,
     $           MPI_SUM,COMM_TINKER,ierr)
         end if
         x_move = sqrt(x_move) / rms
         if (coordtype .eq. 'INTERNAL') then
            x_move = radian * x_move
         end if
c
c     compute the rms gradient per optimization parameter
c
         g_rms = 0.0_re_p
!$acc parallel loop collapse(2) async
         do i = 1, nloc
           do j = 1, 3
             iglob = glob(i)
             g_rms = g_rms + (g(3*(iglob-1)+j)*scale(3*(iglob-1)+j))**2
           end do
         end do

         if (nproc.gt.1) then
            call MPI_ALLREDUCE(MPI_IN_PLACE,g_rms,1,MPI_RPREC,
     $           MPI_SUM,COMM_TINKER,ierr)
         end if
         g_rms = sqrt(g_rms) / rms
c
c     test for error due to line search problems
c
         if (status.eq.'BadIntpln' .or. status.eq.'IntplnErr') then
            nerr = nerr + 1
            if (nerr .ge. maxerr)  done = .true.
         else
            nerr = 0
         end if
c
c     test for too many total iterations
c
         if (niter .ge. maxiter) then
            status = 'IterLimit'
            done = .true.
         end if
c
c     test the normal termination criteria
c
         if (f .le. fctmin) then
            status = 'SmallFct '
            done = .true.
         end if
         if (g_rms .le. grdmin) then
            status = 'SmallGrad'
            done = .true.
         end if
c
c     print intermediate results for the current iteration
c
         if (iprint .gt. 0 .and. rank.eq.0) then
            if (done .or. mod(niter,iprint).eq.0) then
               if (f.lt.1.0d8 .and. f.gt.-1.0d7 .and.
     &             g_rms.lt.1.0d5 .and. f_move.lt.1.0d5) then
                  write (iout,100)  niter,f,g_rms,f_move,x_move,
     &                              angle,ncalls,status
  100             format (i6,f14.4,f11.4,f11.4,f10.4,f8.2,i7,3x,a9)
               else
                  write (iout,110)  niter,f,g_rms,f_move,x_move,
     &                              angle,ncalls,status
  110             format (i6,d14.4,d11.4,d11.4,f10.4,f8.2,i7,3x,a9)
               end if
            end if
         end if
c
c     write intermediate results for the current iteration
c
         if ((iwrite .gt. 0).and.(rank.eq.0)) then
            if (done .or. mod(niter,iwrite).eq.0) then
               call optsave (niter,x0)
            end if
         end if
      end do
c
c     perform deallocation of some local arrays
c
!$acc end data
!$acc exit data delete(beta) async
      deallocate (x_old)
      deallocate (g)
      deallocate (g_old)
      deallocate (p)
      deallocate (q)
      deallocate (r)
      deallocate (h0)
      deallocate (s)
      deallocate (y)
c
c     set final value of the objective function
c
      minimum = f
      if (iprint .gt. 0 .and. rank.eq.0) then
         if (status.eq.'SmallGrad' .or. status.eq.'SmallFct ') then
            write (iout,120)  status
  120       format (/,' LBFGS  --  Normal Termination due to ',a9)
         else
            write (iout,130)  status
  130       format (/,' LBFGS  --  Incomplete Convergence due to ',a9)
         end if
      end if
      return
      end
