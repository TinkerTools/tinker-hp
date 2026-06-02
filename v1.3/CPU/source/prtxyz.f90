!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ##############################################################
!     ##                                                          ##
!     ##  subroutine prtxyz  --  output of Cartesian coordinates  ##
!     ##                                                          ##
!     ##############################################################
!
!
!     "prtxyz" writes out a set of Cartesian coordinates
!     to an external disk file
!
!
!> @brief 
!> writes out a set of Cartesian coordinates
!> to an external disk file
!> @param[in] ixyz: unit of file
subroutine prtxyz (ixyz)
   use atmtyp
   use atoms
   use bound
   use boxes
   use couple
   use files
   use inform
   use titles
   implicit none
   integer i,k,ixyz
   integer size,crdsiz
   real*8 crdmin,crdmax
   logical opened
   character*2 atmc
   character*2 crdc
   character*2 digc
   character*25 fstr
   character*240 xyzfile
!
!
!     open the output unit if not already done
!
   inquire (unit=ixyz,opened=opened)
   if (.not. opened) then
      xyzfile = filename(1:leng)//'.xyz'
      call version (xyzfile,'new')
      open (unit=ixyz,file=xyzfile,status='new')
   end if
!
!     check for large systems needing extended formatting
!
   atmc = 'i6'
   if (n .ge. 100000)  atmc = 'i7'
   if (n .ge. 1000000)  atmc = 'i8'
   crdmin = 0.0d0
   crdmax = 0.0d0
   do i = 1, n
      crdmin = min(crdmin,xwrite(i),ywrite(i),zwrite(i))
      crdmax = max(crdmax,xwrite(i),ywrite(i),zwrite(i))
   end do
   crdsiz = 6
   if (crdmin .le. -1000.0d0)  crdsiz = 7
   if (crdmax .ge. 10000.0d0)  crdsiz = 7
   if (crdmin .le. -10000.0d0)  crdsiz = 8
   if (crdmax .ge. 100000.0d0)  crdsiz = 8
   crdsiz = crdsiz + max(6,digits)
   size = 0
   call numeral (crdsiz,crdc,size)
   if (digits .le. 6) then
      digc = '6 '
   else if (digits .le. 8) then
      digc = '8'
   else
      digc = '10'
   end if
!
!     write out the number of atoms and the title
!
   if (ltitle .eq. 0) then
      fstr = '('//atmc//')'
      write (ixyz,fstr(1:4))  n
   else
      fstr = '('//atmc//',2x,a)'
      write (ixyz,fstr(1:9))  n,title(1:ltitle)
   end if
!
!     write out the periodic cell lengths and angles
!
   if (use_bounds) then
      fstr = '(1x,6f'//crdc//'.'//digc//')'
      write (ixyz,fstr)  xbox,ybox,zbox,alpha,beta,gamma
   end if
!
!     write out the coordinate line for each atom
!
   fstr = '('//atmc//',2x,a3,3f'//crdc//&
   &'.'//digc//',i6,8'//atmc//')'
   do i = 1, n
      write (ixyz,fstr)  i,name(i),xwrite(i),ywrite(i),zwrite(i),&
      &type(i),(i12(k,i),k=1,n12(i))
   end do
!
!     close the output unit if opened by this routine
!
   if (.not. opened)  close (unit=ixyz)
   return
end
