c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #######################################################################
c     ##                                                                   ##
c     ##  utils module  -- Tinker-HP utilities funcitons and subroutines   ##
c     ##                                                                   ##
c     #######################################################################
c
c
c     rpole_ind_extract  indice to extract rpole
c     rpole_der_ext      indice to extract rpole derivative term
c
#include "tinker_precision.h"

      module utils

      integer rpole_ind_extract(10)
      integer rpole_der_ext(9)
      real(t_p) rpole_scale(10)
      integer deriv1(10),deriv2(10),deriv3(10)

      interface front_convert_base
      module procedure front_convert_base3
      module procedure front_convert_base5
      end interface

      interface back_convert_base
      module procedure back_convert_base5
      end interface

      interface is_find
      module procedure is_find4
      module procedure is_find8
      end interface
!      This interface is broken with gcc
!      interface set_to_zero
!      module procedure set_to_zero1
!      module procedure set_to_zero1m
!      module procedure set_to_zero1_int
!      module procedure set_to_zero1_int1
!      module procedure set_to_zero1_int8
!      module procedure set_to_zero2
!      module procedure set_to_zero2_int
!      module procedure set_to_zero3
!      module procedure set_to_zero5
!      end interface

      interface set_value
      module procedure setr
      module procedure seti
      module procedure set_dip
      end interface
      
      interface atomic_Add
      module procedure atomic_Adds
      module procedure atomic_Adds2
      module procedure atomic_Addv
      end interface

      parameter(
     &  rpole_ind_extract =(/1,2,3,4,5,9,13,6,7,10/),
c          each index correspond to
c          qixx, qixy, qixz, qixy, qiyy, qiyz,qixz, qiyz, qizz
     &  rpole_der_ext = (/ 5,6,7,6,9,10,7,10,13 /),
     &  rpole_scale =(/1.0,1.0,1.0,1.0,1.0,1.0,1.0,2.0,2.0,2.0/),
c
c       indices into the electrostatic field array
c
     &  deriv1  =(/ 2, 5,  8,  9, 11, 16, 18, 14, 15, 20 /),
     &  deriv2  =(/ 3, 8,  6, 10, 14, 12, 19, 16, 20, 17 /),
     &  deriv3  =(/ 4, 9, 10,  7, 15, 17, 13, 20, 18, 19 /)
     &)

!$acc declare copyin(rpole_ind_extract,rpole_scale,
!$acc&               deriv1,deriv2,deriv3)

      contains
c
c     ####################################################################
c     ##                                                                ##
c     ##  subroutine imagevec2  --  compute the minimum image distance  ##
c     ##                                                                ##
c     ####################################################################
c
c
c     "imagevec" takes the components of pairwise distances between
c     two points in a periodic box and converts to the components
c     of the minimum image distances
c
#include "tinker_precision.h"
      subroutine imagevec2 (pos2,n)
      use sizes
      use boxes
      use cell

      implicit none
      integer n,i
      real(t_p) pos2(n,3)
      !DIR$ ASSUME_ALIGNED pos2:64
       do while (any(abs(pos2(1:n,1)).gt.xcell2))
          where (    abs(pos2(1:n,1)).gt.xcell2)
             pos2(:,1) = pos2(:,1) -sign(xcell,pos2(:,1))
          end where
       enddo
       do while (any(abs(pos2(1:n,2)).gt.ycell2))
          where (    abs(pos2(1:n,2)).gt.ycell2)
             pos2(:,2) = pos2(:,2) -sign(ycell,pos2(:,2))
          end where
       enddo
       do while (any(abs(pos2(1:n,3)).gt.zcell2))
          where (    abs(pos2(1:n,3)).gt.zcell2)
             pos2(:,3) = pos2(:,3) -sign(zcell,pos2(:,3))
          end where
       enddo
      return
      end
c
c     Convert forward a series of number to base n
c
      subroutine front_convert_base3 (n2,n1,n0,number)
!$acc routine
      implicit none
      integer  ,intent(in) :: n2,n1,n0
      integer*8,intent(out):: number
      integer*8,parameter  :: base=1000

!      if (n2.gt.base.or.n1.gt.base.or.n0.gt.base) then
!         print*, 'ERROR in front_convert_base3 :: ',
!     &           'call arguments are \n',n2,n1,n0 
!      end if
      number = int(n0,8) + int(n1,8)*base + int(n2,8)*base**2
      end subroutine
c
      subroutine front_convert_base5 (n4,n3,n2,n1,n0,number)
!$acc routine
      implicit none
      integer  ,intent(in) :: n4,n3,n2,n1,n0
      integer*8,intent(out):: number
      !real*8   ,parameter  :: base=1000
      integer*8,parameter  :: base=1000
      integer*8 :: n4_8,n3_8,n2_8,n1_8,n0_8
!      if (n4.gt.base.or.n3.gt.base.or.n2.gt.base.or.
!     &    n1.gt.base.or.n0.gt.base) then
!         print*, 'ERROR in front_convert_base5 :: ',
!     &           'call arguments are \n',n4,n3,n2,n1,n0 
!      end if
      n0_8   = int(n0,8)
      n1_8   = int(n1,8)
      n2_8   = int(n2,8)
      n3_8   = int(n3,8)
      n4_8   = int(n4,8)
      number = n0_8 + n1_8*base    + n2_8*base**2 
     &              + n3_8*base**3 + n4_8*base**4

      end subroutine
c
c     Decomposition of base-n number
c
c      subroutine back_convert_base5 (n4,n3,n2,n1,n0,number)
c!$acc routine
c      implicit none
c      integer*8,intent(in) :: number
c      integer  ,intent(out):: n4,n3,n2,n1,n0
c      integer*8,parameter  :: base=1000
c      real*8 :: cnum, quot, rest, quot_r
c
c      cnum=dble(number)
c
c      if (cnum.lt.0) then
c         print*,"ERROR ! Can't convert back a negative number",
c     &          number
c      end if
c
c      quot_r = cnum/base
c      quot   = aint(quot_r)
c      n0     = nint((quot_r - aint(quot_r))*base)
c      quot_r = quot/base
c      quot   = aint(quot_r)
c      n1     = nint((quot_r - aint(quot_r))*base)
c      quot_r = quot/base
c      quot   = aint(quot_r)
c      n2     = nint((quot_r - aint(quot_r))*base)
c      quot_r = quot/base
c      quot   = aint(quot_r)
c      n3     = nint((quot_r - aint(quot_r))*base)
c      quot_r = quot/base
c      quot   = aint(quot_r)
c      n4     = nint((quot_r - aint(quot_r))*base)
c      quot_r = quot/base
c      quot   = aint(quot_r)
c
c      if (quot.ne.0)
c     &   print*,"ERROR ! convert unfinished",cnum
c      end subroutine

      subroutine back_convert_base5 (n4,n3,n2,n1,n0,number)
!$acc routine
      implicit none
      integer*8,intent(in) :: number
      integer  ,intent(out):: n4,n3,n2,n1,n0
      integer*8,parameter  :: base=1000
      integer*8:: cnum

      cnum=number
!      if (cnum.lt.0) then
!         print*,"ERROR ! Can't convert back a negative number",
!     &          number
!      end if

      n0   = int(mod(cnum,base),kind(n0))
      cnum = cnum/base
      n1   = int(mod(cnum,base),kind(n1))
      cnum = cnum/base
      n2   = int(mod(cnum,base),kind(n2))
      cnum = cnum/base
      n3   = int(mod(cnum,base),kind(n3))
      cnum = cnum/base
      n4   = int(mod(cnum,base),kind(n4))
      cnum = cnum/base

!      if (cnum.ne.0)
!     &   print*,"ERROR ! convert unfinished",cnum
      end subroutine
c
c     Search for a number inside an array
c
      function is_find8(array,n,number) result(bool)
      implicit none
      integer*8,intent(in ):: array(*)
      integer*8,intent(in ):: number
      integer  ,intent(in ):: n
      logical:: bool
      integer i 

      bool = .false.
      do i = 1, n
         if (array(i).eq.number) then
            bool = .true.
            return
         end if
      end do
      end
      function is_find4(array,n,number) result(bool)
      implicit none
      integer  ,intent(in ):: array(*)
      integer  ,intent(in ):: number
      integer  ,intent(in ):: n
      logical:: bool
      integer i 

      bool = .false.
      do i = 1, n
         if (array(i).eq.number) then
            bool = .true.
            return
         end if
      end do
      end

      subroutine Tinker_Wait(sec,rank)
      implicit none
      integer,intent(in)::sec,rank
      integer i

      if (rank.eq.0) 
     &   write(0,'(A)',advance='no') 'Tinker wait countdown'
 13   format(I3)
      do i = sec,1,-1
         if (btest(i,0).and.rank.eq.0) then
            write(0,13,advance='no'),i
         end if
         if (i.eq.1.and.rank.eq.0) write(0,'(A)') '  Resume'
         call sleep(1)
      end do
      end subroutine
c
c     set_to_zero routines
c
      subroutine set_to_zero1(array1,n,queue)
      implicit none
      integer,intent(in):: n
      integer,intent(in)::queue
      real(t_p),intent(out):: array1(*)
      integer:: i
      real(t_p),parameter:: zero=0

!$acc parallel loop 
!$acc&         present(array1(1:n)) async(queue) 
      do i= 1, n
         array1(i) = zero
      end do
      end subroutine
      subroutine set_to_zero1m(array1,n,queue)
      implicit none
      integer  ,intent(in ):: n
      integer  ,intent(in ):: queue
      real(r_p),intent(out):: array1(*)
      integer:: i
      real(r_p),parameter:: zero=0

!$acc parallel loop 
!$acc&         present(array1(1:n)) async(queue) 
      do i= 1, n
         array1(i) = zero
      end do
      end subroutine
      subroutine set_to_zero1_int(array1,n,queue)
      implicit none
      integer,intent(in ):: n
      integer,intent(in ):: queue
      integer,intent(out):: array1(*)
      integer:: i
      integer,parameter:: zero=0

!$acc parallel loop 
!$acc&         present(array1(1:n)) async(queue) 
      do i= 1, n
         array1(i) = zero
      end do
      end subroutine
      subroutine set_to_zero1_int8(array1,n,queue)
      implicit none
      integer(8),intent(in ):: n
      integer   ,intent(in ):: queue
      integer   ,intent(out):: array1(*)
      integer:: i
      integer,parameter:: zero=0

!$acc parallel loop 
!$acc&         present(array1(1:n)) async(queue) 
      do i= 1, n
         array1(i) = zero
      end do
      end subroutine
      subroutine set_to_zero1_int1(array1,n,queue)
      implicit none
      integer   ,intent(in ):: n
      integer   ,intent(in ):: queue
      integer(1),intent(out):: array1(*)
      integer:: i
      integer,parameter:: zero=0

!$acc parallel loop 
!$acc&         present(array1(1:n)) async(queue) 
      do i= 1, n
         array1(i) = zero
      end do
      end subroutine
      subroutine set_to_zero2(array1,array2,n,queue)
      implicit none
      integer  ,intent(in ):: n
      integer  ,intent(in ):: queue
      real(t_p),intent(out):: array1(*),array2(*)
      integer:: i
      real(t_p),parameter:: zero=0.0

!$acc parallel loop 
!$acc&         present(array1(1:n),array2(1:n)) async(queue) 
      do i= 1, n
         array1(i) = zero
         array2(i) = zero
      end do
      end subroutine
      subroutine set_to_zero2m(array1,array2,n,queue)
      implicit none
      integer  ,intent(in ):: n
      integer  ,intent(in ):: queue
      real(r_p),intent(out):: array1(*),array2(*)
      integer:: i
      real(r_p),parameter:: zero=0.0

!$acc parallel loop 
!$acc&         present(array1(1:n),array2(1:n)) async(queue) 
      do i= 1, n
         array1(i) = zero
         array2(i) = zero
      end do
      end subroutine
      subroutine set_to_zero2_int(array1,array2,n,queue)
      implicit none
      integer,intent(in ):: n
      integer,intent(in ):: queue
      integer,intent(out):: array1(*),array2(*)
      integer:: i
      integer,parameter:: zero=0

!$acc parallel loop 
!$acc&         present(array1(1:n),array2(1:n)) async(queue) 
      do i= 1, n
         array1(i) = zero
         array2(i) = zero
      end do
      end subroutine
      subroutine set_to_zero3(array1,array2,array3,n,queue)
      implicit none
      integer  ,intent(in ):: n
      integer  ,intent(in ):: queue
      real(t_p),intent(out):: array1(*),array2(*),array3(*)
      integer:: i
      real(t_p),parameter:: zero=0

!$acc parallel loop 
!$acc&         present(array1(1:n),array2(1:n),array3(1:n))
!$acc&         async(queue) 
      do i= 1, n
         array1(i) = zero
         array2(i) = zero
         array3(i) = zero
      end do
      end subroutine
      subroutine set_to_zero5(array1,array2,array3,array4,array5,
     &                        n,queue)
      implicit none
      integer  ,intent(in ):: n
      integer  ,intent(in ):: queue
      real(t_p),intent(out):: array1(*),array2(*),array3(*),
     &                        array4(*),array5(*)
      integer:: i
      real(t_p),parameter:: zero=0

!$acc parallel loop 
!$acc&         present(array1(1:n),array2(1:n),array3(1:n),
!$acc&  array4(1:n),array5(1:n))
!$acc&         async(queue) 
      do i= 1, n
         array1(i) = zero
         array2(i) = zero
         array3(i) = zero
         array4(i) = zero
         array5(i) = zero
      end do
      end subroutine
c
      function comput_norm(array,siz,expo)
      use tinheader ,only: ti_p
      implicit none
      integer  ,intent(in):: siz,expo
      real(t_p),intent(in):: array(*)
      real(8) comput_norm
      integer i

      comput_norm = 0.0_ti_p
      if (expo.eq.0) then
         do i=1,siz
            if (comput_norm<abs(array(i)))
     &         comput_norm = abs(array(i))
         end do
      else if (btest(expo,0)) then  
         do i=1,siz
            comput_norm = comput_norm + abs(array(i))**expo
         end do
      else
         do i=1,siz
            comput_norm = comput_norm + array(i)**expo
         end do
      end if

      end function
      function comput_normr(array,siz,expo)
      use tinheader ,only: re_p
      implicit none
      integer  ,intent(in):: siz,expo
      real(r_p),intent(in):: array(*)
      real(8) comput_normr
      integer i

      comput_normr = 0.0_re_p
      if (expo.eq.0) then
         do i=1,siz
            if (comput_normr<abs(array(i)))
     &         comput_normr = abs(array(i))
         end do
      else if (btest(expo,0)) then  
         do i=1,siz
            comput_normr = comput_normr + abs(array(i))**expo
         end do
      else
         do i=1,siz
            comput_normr = comput_normr + array(i)**expo
         end do
      end if
      end function

c
      subroutine utils_amove(n,a,b,stream_)
      implicit none
      integer::j,stream
      integer,intent(in)::n
      integer,intent(in),optional::stream_
      real(t_p),intent(in ):: a(*)
      real(t_p),intent(out):: b(*)

      stream = -1
      if (present(stream_)) stream = stream_
!$acc parallel loop present(b(1:n),a(1:n)) async(stream)
      do j = 1, n
         b(j) = a(j)
      end do
      end
c
      subroutine amove2gpu(src1,src2,dst1,dst2,n,queue_)
#ifdef _OPENACC
      use openacc ,only: acc_async_sync
#endif
      implicit none
      real(t_p),intent(in)  :: src1(*),src2(*)
      real(t_p),intent(out) :: dst1(*),dst2(*)
      integer,intent(in)::n
      integer,intent(in),optional::queue_
      integer j
#ifdef _OPENACC
      integer queue

      if (present(queue_)) then
         queue=queue_
      else
         queue=acc_async_sync
      end if
#endif

!$acc parallel loop present(src1(1:n),src2(1:n),
!$acc&                      dst1(1:n),dst2(1:n))
!$acc&         async(queue)
      do j = 1, n
        dst1(j) = src1(j)
        dst2(j) = src2(j)
      end do
      end
c
c
      subroutine convolution_product (scalevec,n,grid,queue_)
      implicit none
      integer::i,j,queue=-1
      integer,intent(in)::n
      integer,intent(in),optional::queue_
      real(t_p),intent(in )::scalevec(*)
      real(t_p),intent(out)::grid(*)

      if (present(queue_)) queue = queue_

!$acc parallel loop present(scalevec,grid(1:n)) async(queue)
      do i = 1, n
         j = ishft(i,-1)+iand(i,1)
         grid(i) = scalevec(j)*grid(i)
      end do

      end subroutine
c
c     setter
c
      subroutine setr(index,value,array)
!$acc routine
      implicit none
      integer  ,intent(in )::index
      real(t_p),intent(in )::value
      real(t_p),intent(out)::array(*)

      array(index) = value

      end subroutine
c
      subroutine set_dip(index,value,array)
!$acc routine
      implicit none
      integer  ,intent(in )::index
      real(t_p),intent(in )::value(9)
      real(t_p),intent(out)::array(10,*)
      integer i

!$acc loop seq
      do i = 1,9
         array(i+1,index) = value(i)
      end do

      end subroutine
c
      subroutine seti(index,value,array)
!$acc routine
      implicit none
      integer,intent(in )::index
      integer,intent(in )::value
      integer,intent(out)::array(*)

      array(index) = value

      end subroutine
c
      subroutine atomic_Adds(loc,value,array)
!$acc routine
      implicit none
      integer  ,intent(in )::loc
      real(t_p),intent(in )::value
      real(t_p),intent(inout)::array(*)

!$acc atomic update
      array(loc) = array(loc) + value

      end subroutine
c
      subroutine atomic_Adds2(loc,value1,value2,array)
!$acc routine
      implicit none
      integer  ,intent(in )::loc
      real(t_p),intent(in )::value1,value2
      real(t_p),intent(inout)::array(*)

!$acc atomic update
      array(loc)   = array(loc) + value1
!$acc atomic update
      array(loc+1) = array(loc+1) + value2

      end subroutine
c
      subroutine atomic_Addv(loc,value,n,array)
!$acc routine
      implicit none
      integer  ,intent(in )::loc,n
      real(t_p),intent(in )::value(n)
      real(t_p),intent(inout)::array(*)
      integer i

      do i= 0, n-1
!$acc atomic update
         array(loc+i) = array(loc+i) + value(i+1)
      end do

      end subroutine
c
      subroutine atomic_Add_field(loc,value,array)
!$acc routine
      implicit none
      integer  ,intent(in )::loc
      real(t_p),intent(in )::value(6)
      real(t_p),intent(inout),contiguous::array(:,:,:)
      integer i,j

      do j = 1, 2
         do i = 1, 3
!$acc atomic update
            array(i,j,loc) = array(i,j,loc) + value(2*j+i)
         end do
      end do

      end subroutine

      end module


c
c     ###############################################################
c     #                                                             #
c     #  -- integrator workSpace module --                          #
c     #  holds data structure to be used inside integrator routines #
c     #                                                             #
c     ###############################################################
c

c     derivs  stores forces computes by gradient routine
c     etot    holds the system total energy at current timestep
c     epot    holds the system potential energy at current timestep
c     eksum   holds the system kinetic energy at current timestep
c     temp    holds the system temperature at current timestep
c     pres    holds the system pressure at current timestep

      module integrate_ws
      real(r_p),allocatable::derivs(:,:)
      real(r_p) etot,epot,eksum
      real(r_p) temp,pres
      real(r_p) ekin(3,3)
      real(r_p) stress(3,3)
      end module
