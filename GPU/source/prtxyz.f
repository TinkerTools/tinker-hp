c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine prtxyz  --  output of Cartesian coordinates  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "prtxyz" writes out a set of Cartesian coordinates
c     to an external disk file
c
c
#include "tinker_precision.h"
      subroutine prtxyz (ixyz)
      use atmtyp
      use atoms    ,only: type
      use atomsMirror
      use bound
      use boxes
      use couple
      use files
      use inform
      use timestat ,only:timer_io,timer_enter,timer_exit,quiet_timers
      use tinheader
      use titles
      implicit none
      integer i,k,ixyz
      integer size,crdsiz
      real(r_p) crdmin,crdmax
      logical opened
      character*2 atmc
      character*2 crdc
      character*2 digc
      character*25 fstr
      character*240 xyzfile
c
      call timer_enter(timer_io)
c
c     open the output unit if not already done
c
      inquire (unit=ixyz,opened=opened)
      if (.not. opened) then
         xyzfile = filename(1:leng)//'.xyz'
         call version (xyzfile,'new')
         open (unit=ixyz,file=xyzfile,status='new')
      end if
c
c     check for large systems needing extended formatting
c
      atmc = 'i6'
      if (n .ge. 100000)  atmc = 'i7'
      if (n .ge. 1000000)  atmc = 'i8'
      crdmin = 0.0_ti_p
      crdmax = 0.0_ti_p
!$acc update host(x,y,z) async
!$acc parallel loop present(x,y,z) async
      do i = 1, n
         crdmin = min(crdmin,x(i),y(i),z(i))
         crdmax = max(crdmax,x(i),y(i),z(i))
      end do
!$acc wait
      crdsiz = 6
      if (crdmin .le. -1000.0_re_p)  crdsiz = 7
      if (crdmax .ge. 10000.0_re_p)  crdsiz = 7
      if (crdmin .le.-10000.0_re_p)  crdsiz = 8
      if (crdmax .ge.100000.0_re_p)  crdsiz = 8
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
c
c     write out the number of atoms and the title
c
      if (ltitle .eq. 0) then
         fstr = '('//atmc//')'
         write (ixyz,fstr(1:4))  n
      else
         fstr = '('//atmc//',2x,a)'
         write (ixyz,fstr(1:9))  n,title(1:ltitle)
      end if
c
c     write out the periodic cell lengths and angles
c
      if (use_bounds) then
         fstr = '(1x,6f'//crdc//'.'//digc//')'
         write (ixyz,fstr)  xbox,ybox,zbox,alpha,beta,gamma
      end if
c
c     write out the coordinate line for each atom
c
      fstr = '('//atmc//',2x,a3,3f'//crdc//
     &          '.'//digc//',i6,8'//atmc//')'
      do i = 1, n
         write (ixyz,fstr)  i,name(i),x(i),y(i),z(i),type(i),
     &                      (i12(k,i),k=1,n12(i))
      end do
c
c     close the output unit if opened by this routine
c
      if (.not. opened)  close (unit=ixyz)
      call timer_exit(timer_io)
      end
