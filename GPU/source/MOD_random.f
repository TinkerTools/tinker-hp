c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #############################################################
c     ##                                                         ##
c     ##  function random  --  portable random number generator  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "random" generates a random number on [0,1] via a long
c     period generator due to L'Ecuyer with Bays-Durham shuffle
c
c     literature references:
c
c     P. L'Ecuyer, Communications of the ACM, 31, 742-774 (1988)
c
c     W. H. Press, S. A. Teukolsky, W. T. Vetterling and B. P.
c     Flannery, Numerical Recipes (Fortran), 2nd Ed., Cambridge
c     University Press, 1992, Section 7.1
c
c
#include "tinker_macro.h"
      module random_mod
      use domdec
      use inform
      use iounit
      use keys
      use math
      use tinheader,only:ti_p
      use tinMemory,only:prmem_request
#ifdef _OPENACC
      !use iso_c_binding
      use curand
      use openacc
      !use curand_device
      !use cudafor
      !use openacc_curand
#endif
      implicit none
#include "tinker_macro.h"
      private
      integer im1,ia1,iq1,ir1
      integer im2,ia2,iq2,ir2
      integer big,nshuffle
      integer imm1,ndiv
      integer iyr
      integer:: seed,seed2,seed_save=0
      integer(8) npick,npickD,npickDr
      real(t_p) factor_u,store
      real(t_p),public,allocatable::samplevec(:)
      real(t_p),public :: pick
      real(r_p),public ::pick1
      logical first
      logical compute
#ifdef _OPENACC
c      type(curandStateXORWOW),pointer :: curng(:)
c      integer Ncurng
c!$acc declare device_resident(curng)
c      parameter (Ncurng=5120)
      type(curandGenerator) :: curng
      logical, public :: host_rand_platform=.false.
#else
      logical,parameter,public::host_rand_platform=.true.
#endif
      parameter (im1=2147483563)
      parameter (ia1=40014)
      parameter (iq1=53668)
      parameter (ir1=12211)
      parameter (im2=2147483399)
      parameter (ia2=40692)
      parameter (iq2=52774)
      parameter (ir2=3791)
      parameter (big=141803398)
      parameter (nshuffle=32)
      parameter (imm1=im1-1)
      parameter (ndiv=1+imm1/nshuffle)
      parameter (factor_u=1.0_ti_p/im1)
      integer ishuffle(nshuffle)
      data first  / .true. /
      data compute  / .true. /

      public :: init_rand_engine,random,randomvec,normal,normalvec
     &         ,ranvec,sphere,get_randseed,get_pickcount,normalarray
#ifdef _OPENACC
     &         ,init_curand_engine,randomgpu,normalgpu,rand_unitgpu
     &         ,normalgpuR4,reset_curand_seed,disp_ranvec

      interface normalgpu
         module procedure normalgpuR4
         module procedure normalgpuR8
      end interface
#endif

      interface normalarray
         module procedure normalarrayR4
         module procedure normalarrayR8
      end interface
      contains

      integer function get_randseed()
      get_randseed = seed_save
      end function

      subroutine get_pickcount(npic,npicD,npicDr)
      implicit none
      integer(8) npic, npicD, npicDr
      npic  = npick
      npicD = npickD
      npicDr= npickDr
      end subroutine

      subroutine init_rand_engine
      implicit none
      integer i,k,next
      integer year,month,day
      integer hour,minute,second
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     random number seed is first set to a big number,
c     then incremented by the seconds elapsed this decade
c
      first = .false.
      if (seed_save.eq.0) then
         seed = big
         call calendar (year,month,day,hour,minute,second)
         year = mod(year,10)
         seed = seed + 32140800*year + 2678400*(month-1)
         seed = seed + 86400*(day-1) + 3600*hour
         seed = seed + 60*minute + second + 10*rank
c
c        search the keywords for a random number seed
c
         do i = 1, nkey
            next = 1
            record = keyline(i)
            call gettext (record,keyword,next)
            call upcase (keyword)
            if (keyword(1:11) .eq. 'RANDOMSEED ') then
               string = record(next:240)
               read (string,*,err=10)  seed
               seed = max(1,seed)
            end if
   10       continue
         end do
         seed = seed + ranktot
      else
         seed = seed_save
      end if
c
c     print the value used for the random number seed
c
      if ((verbose).and.(rank.eq.0)) then
         write (iout,20)  seed
   20    format (/,' Random Number Generator Initialized',
     &              ' with SEED :',3x,i12)
      end if
c
c     Initiate counter
c
      npick=0
c
c     warm up and then load the shuffling table
c
      seed_save = seed
      seed2 = seed
      do i = nshuffle+8, 1, -1
         k = seed / iq1
         seed = ia1 * (seed-k*iq1) - k*ir1
         if (seed .lt. 0)  seed = seed + im1
         if (i .le. nshuffle)  ishuffle(i) = seed
      end do
      iyr = ishuffle(1)
      end subroutine


      function random ()
      implicit none
      real(t_p) random
      integer k,i

      if (first) call init_rand_engine
c
c     get a new random number value each call
c
      k = seed / iq1
      seed = ia1*(seed-k*iq1) - k*ir1
      if (seed .lt. 0)  seed = seed + im1
      k = seed2 / iq2
      seed2 = ia2*(seed2-k*iq2) - k*ir2
      if (seed2 .lt. 0)  seed2 = seed2 + im2
      i = 1 + iyr/ndiv
      iyr = ishuffle(i) - seed2
      ishuffle(i) = seed
      if (iyr .lt. 1)  iyr = iyr + imm1
      random = factor_u * iyr
c
c     print the value of the current random number
c
      !Incement counter
      npick = npick+1
c     if (debug) then
c        write (iout,30)  random
c  30    format (' RANDOM  --  The Random Number Value is',f12.8)
c     end if
      return
      end

      subroutine normalarrayR4(array,n,mean_,stddev_,stream_)
      implicit none
      real(4), intent(inout) ::  array(*)
      integer,   intent(in)::n
      integer,   intent(in),optional::stream_
      real(4), intent(in),optional::mean_,stddev_
      integer:: i

#ifdef _OPENACC
        if(.not. host_rand_platform) then
          call normalgpuR4(array,n,mean_,stddev_,stream_)
        endif
#endif
        if (host_rand_platform) then
          do i = 1, n
              array(i) = normal()
          end do
          if (present(stddev_)) array(1:n) = array(1:n) * stddev_
          if (present(mean_))   array(1:n) = array(1:n) + mean_
!$acc update device(array(1:n)) async
        end if

      end subroutine normalarrayR4

      subroutine normalarrayR8(array,n,mean_,stddev_,stream_)
      implicit none
      real(8), intent(inout) ::  array(*)
      integer,   intent(in)::n
      integer,   intent(in),optional::stream_
      real(8), intent(in),optional::mean_,stddev_
      integer:: i

#ifdef _OPENACC
        if (.not. host_rand_platform) then
           call normalgpuR8(array,n,mean_,stddev_,stream_)
        endif
#endif
        if (host_rand_platform) then
          do i = 1, n
              array(i) = normal()
          end do
          if (present(stddev_)) array(1:n) = array(1:n) * stddev_
          if (present(mean_))   array(1:n) = array(1:n) + mean_
!$acc update device(array(1:n)) async
        end if

      end subroutine normalarrayR8

#ifdef _OPENACC
c
c      All functions within this scope are specific to the GPU
c
c     Inititiate curand generator
c
c      subroutine init_curand_engine
c      implicit none
c      type(c_ptr) :: c_curng
c      type(c_devptr) :: cdev_curng
c      integer i
c      write(*,*),'***** Init curand generator engine with seed ',
c     &          ,seed, ' ***** '
c
c      if (first) call init_rand_engine
c
c      cdev_curng = acc_malloc(Ncurng*sizeof(curandStateXORWOW))
c      c_curng    = acc_hostptr(cdev_curng) 
cc
cc     Translate c to fortran pointer
cc
c      call c_f_pointer(c_curng,curng,(/Ncurng/))
c
c!$acc parallel loop deviceptr(curng)
c      do i = 1,Ncurng
c         call curand_init(seed_save,i,0,curng(i))
c      end do
c      end

      subroutine init_curand_engine(n,stream)
      implicit none
      integer,intent(in)::n
      integer,  optional::stream
      integer :: istat=0,siz

      if (first) call init_rand_engine
c
c     Allocate samplevec
c
      first = .false.
      siz   = 3*n
      if (btest(siz,0)) siz = siz + 1  !check if siz is odd
      call prmem_request(samplevec,siz)
c
c     Select the cuRand pseudo Generator type
c
c     istat = istat + curandCreateGenerator(curng,
c    &                   CURAND_RNG_PSEUDO_XORWOW)
      istat = istat + curandCreateGenerator(curng,
     &                   CURAND_RNG_PSEUDO_DEFAULT)
c     istat = istat + curandCreateGenerator(curng,
c    &                   CURAND_RNG_PSEUDO_MRGP32)

      istat = istat + curandSetPseudoRandomGeneratorSeed(curng,
     &                   seed_save)
c
c     Initiate counter
c
      npickD = 0
      npickDr= 0
c
c     Set cuda stream if necessary
c
      if (present(stream))
     &   istat = istat + curandSetStream(curng,
     &                       acc_get_cuda_stream(stream))

      if (istat.ne.0)
     &   print*,'curand initialisation failed ',istat

      end subroutine

      subroutine reset_curand_seed
      implicit none
      integer istat
      istat = curandSetPseudoRandomGeneratorSeed(curng,
     &                   seed_save)
      if (istat.ne.0)
     &   print*,'curand set seed failed ',istat
      end subroutine
c
c     get a sample following uniform distribution on GPU
c
      subroutine randomgpu(array,n,stream_)
      implicit none
      real(t_p), array(*)
      integer, intent(in):: n
      integer, intent(in),optional::stream_
      integer:: istat
      istat = 0

      if (present(stream_))
     &   istat = istat + curandSetStream(curng,
     &                       acc_get_cuda_stream(stream_))

!$acc host_data use_device(array)
      istat = istat + curandGenerate(curng,array,n)
!$acc end host_data
      if (istat.ne.0) print*,'ERROR in curandGenerate',istat

      !Increment counter
      npickDr= npickDr+ n
      end
c
c     get a sample following normal distribution on GPU
c
      subroutine normalgpuR8(array,n,mean_,stddev_,stream_)
      implicit none
      real(8), array(*)
      integer,   intent(in)::n
      integer,   intent(in),optional::stream_
      real(8), intent(in),optional::mean_,stddev_
      real(8) mean,stddev
      integer istat,num
      istat=0

      mean   = 0.0
      stddev = 1.0
      if (present(mean_)) mean = mean_
      if (present(stddev_)) stddev = stddev_
      if (present(stream_))
     &   istat = istat + curandSetStream(curng,
     &                       acc_get_cuda_stream(stream_))

      ! normal generator needs 'num' to be even to work with
      ! pseudo rand generator
      if (btest(n,0)) then
         num = n + 1  !first multiple of two following n
      else
         num = n
      end if

!$acc host_data use_device(array)
      istat = istat + curandGenerate(curng,array,num,mean,stddev)
!$acc end host_data

      if (istat.ne.0) print*,'ERROR in curandGenerateNormalR8',istat
      !Increment counter
      npickD = npickD + num
      end subroutine
c
c     get a sample following normal distribution on GPU with forced simple precision
c
      subroutine normalgpuR4(array,n,mean_,stddev_,stream_)
      implicit none
      real(4), array(*)
      integer,   intent(in)::n
      integer,   intent(in),optional::stream_
      real(4), intent(in),optional::mean_,stddev_
      real(4):: mean,stddev
      integer :: istat,num
      istat=0

      mean=0.0; stddev=1.0
      if (present(mean_)) mean = mean_
      if (present(stddev_)) stddev = stddev_
      if (present(stream_))
     &   istat = istat + curandSetStream(curng,
     &                       acc_get_cuda_stream(stream_))

      ! normal generator needs 'num' to be even to work with
      ! pseudo rand generator
      if (btest(n,0)) then
         num = n + 1  !first multiple of two following n
      else
         num = n
      end if

!$acc host_data use_device(array)
      istat = istat + curandGenerate(curng,array,num,mean,stddev)
!$acc end host_data

      if (istat.ne.0) print*,'ERROR in curandGenerateNormalR4',istat
      !Increment counter
      npickD = npickD + num
      end subroutine
c
c     Following the principle of ranvec but using the spherical approch
c     to minimize calls to cuRand library and also operates on an array
c
      subroutine rand_unitgpu (vector,n)
      implicit none
      integer i,j
      integer,  intent(in) :: n
      real(t_p),intent(out):: vector(3,*)
      real(t_p), theta,phi,pi
      real(t_p) x,y,s
      parameter(pi=3.141592653589793d0)
c
c     get a pair of random spherical coordinates
c
      call randomgpu(samplevec,2*n+2)

!$acc parallel loop collapse(2) default(present) async
      do i = 0, n-1; do j = 1, 2
         if (btest(j,0)) then
            samplevec(1+2*i) = -pi + 2*pi*samplevec(1+2*i) !theta angle
         else
            samplevec(2+2*i) = -0.5*pi + pi*samplevec(2+2*i) !phi angle
         end if
      end do; end do
c
c     construct the 3-dimensional random unit vector
c
!$acc parallel loop collapse(2) default(present) async
      do i = 0, n-1; do j= 1, 3
         if (j.eq.1) then
            vector(1,i+1) = sin(samplevec(1+2*i))*cos(samplevec(2+2*i))
         else if (j.eq.2) then
            vector(2,i+1) = sin(samplevec(1+2*i))*sin(samplevec(2+2*i))
         else
            vector(3,i+1) = cos(samplevec(1+2*i))
         end if
      end do; end do
c
c     print the components of the random unit vector
c
c     if (debug) then
c        write (iout,10)  vector(1),vector(2),vector(3)
c  10    format (' RANVEC  --  The Random Vector is',3f10.4)
c     end if
      end

c
c     Destroy rng
c
      subroutine destroy_curand_engine
      implicit none
      integer:: istat=0

      istat = istat + curandDestroyGenerator(curng)
      if (istat.ne.0)
     &   print*,'curand engine destruction failed'

      end
#endif
      subroutine disp_ranvec(n)
      implicit none
      integer i,n
 16   format(I4,10F9.6)
      !call randomgpu(samplevec(1),20)
      !call normalgpu(samplevec(1),10)
      !call randomgpu(samplevec(11),10)
      !call normalgpu(samplevec(21),20)
!$acc wait
!$acc update host(samplevec(1:500))
      do i = 0,4
         print 16, i*10, samplevec(i*10+1:i*10+10)
      end do
      end subroutine
c
c
c     ############################################################
c     ##                                                        ##
c     ##  function normal  --  random number from normal curve  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "normal" generates a random number from a normal Gaussian
c     distribution with a mean of zero and a variance of one
c
c
      function normal ()
      implicit none
      real(t_p) v1,v2,rsq
      real(t_p) factorn,normal
c
c
c     get a pair of random values from the distribution
c
      if (compute) then
   10    continue
         v1 = 2.0_ti_p * random () - 1.0_ti_p
         v2 = 2.0_ti_p * random () - 1.0_ti_p
         rsq = v1**2 + v2**2
         if (rsq .ge. 1.0_ti_p)  goto 10
         factorn = sqrt(-2.0_ti_p*log(rsq)/rsq)
         store = v1 * factorn
         normal = v2 * factorn
         compute = .false.
c
c     use the second random value computed at the last call
c
      else
         normal = store
         compute = .true.
      end if
c
c     print the value of the current random number
c
c     if (debug) then
c        write (iout,20)  normal
c  20    format (' NORMAL  --  The Random Number Value is',f12.8)
c     end if
      return
      end
c
c     get a sample following uniform distribution
c
      subroutine randomvec(array,n)
      implicit none
      real(t_p) array(*)
      integer n,i
      do i = 1,n; array(i) = random(); end do
      end subroutine
c
c     get a sample following normal distribution
c
      subroutine normalvec(array,n)
      implicit none
      real(t_p) array(n)
      integer i,n
      do i = 1,n; array(i) = normal(); end do;
      end subroutine
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine ranvec  --  unit vector in random direction  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "ranvec" generates a unit vector in 3-dimensional
c     space with uniformly distributed random orientation
c
c     literature references:
c
c     G. Marsaglia, Ann. Math. Stat., 43, 645 (1972)
c
c     R. C. Rapaport, The Art of Molecular Dynamics Simulation,
c     2nd Edition, Cambridge University Press, 2004, Section 18.4
c
c
      subroutine ranvec (vector)
      implicit none
      real(t_p) x,y,s
      real(t_p) vector(3)
c
c
c     get a pair of appropriate components in the plane
c
      s = 2.0_ti_p
      do while (s .ge. 1.0_ti_p)
         x = 2.0_ti_p * random () - 1.0_ti_p
         y = 2.0_ti_p * random () - 1.0_ti_p
         s = x**2 + y**2
      end do
c
c     construct the 3-dimensional random unit vector
c
      vector(3) = 1.0_ti_p - 2.0_ti_p*s
      s = 2.0_ti_p * sqrt(1.0_ti_p - s)
      vector(2) = s * y
      vector(1) = s * x
c
c     print the components of the random unit vector
c
c     if (debug) then
c        write (iout,10)  vector(1),vector(2),vector(3)
c  10    format (' RANVEC  --  The Random Vector is',3f10.4)
c     end if
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine sphere  --  uniform set of points on sphere  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "sphere" finds a specified number of uniformly distributed
c     points on a sphere of unit radius centered at the origin
c
c     literature reference:
c
c     E. B. Saff and A. B. J. Kuijlaars, "Distributing Many
c     Points on a Sphere", The Mathematical Intelligencer,
c     19, 5-11 (1997)
c
c
      subroutine sphere (ndot,dot)
      implicit none
      integer i,ndot
      real(t_p) theta,phi
      real(t_p) h,phiold
      real(t_p) tot,tot1
      real(t_p) dot(3,*)
c
c
c     find spherical coordinates then convert to Cartesian
c
      tot = real(ndot,t_p)
      tot1 = real(ndot-1,t_p)
      do i = 1, ndot
         h = -1.0_ti_p + 2.0_ti_p*real(i-1,t_p)/tot1
         theta = acos(h)
         if (i.eq.1 .or. i.eq.ndot) then
            phi = 0.0_ti_p
         else
            phi = mod(phiold+3.6_ti_p/sqrt(tot*(1.0_ti_p-h*h)),
     &                2.0_ti_p*pi)
         end if
         dot(1,i) = sin(theta) * cos(phi)
         dot(2,i) = sin(theta) * sin(phi)
         dot(3,i) = cos(theta)
         phiold = phi
      end do
      return
      end
      end module
