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
