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
#include "tinker_macro.h"
subroutine unwrap_pos(xr,yr,zr,wrapIdx,k)
   use cell
   implicit none
   integer,intent(in)::k
   integer(1) wrapIdx(*)
   real(r_p),intent(inout):: xr,yr,zr
   integer kk

   kk = 4*(k-1)
   xr = xr + int(wrapIdx(kk+1))*xcell
   yr = yr + int(wrapIdx(kk+2))*ycell
   zr = zr + int(wrapIdx(kk+3))*zcell
end subroutine

subroutine prtxyz (ixyz,x_frame,y_frame,z_frame,atom_count_frame)
   use atmtyp
   use atoms    ,only: type,pbcunwrap,pbcWrapIdx
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
   real(r_p), optional, intent(in) :: x_frame(:),y_frame(:),z_frame(:)
   integer, optional, intent(in) :: atom_count_frame
   integer atom_count
   integer size,crdsiz
   real(r_p) crdmin,crdmax,xi,yi,zi
   logical opened
   character*2 atmc
   character*2 crdc
   character*2 digc
   character*25 fstr
   character*240 xyzfile
!
   atom_count = n
   if (present(atom_count_frame)) atom_count = atom_count_frame
   call timer_enter(timer_io)
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
   if (atom_count .ge. 100000)  atmc = 'i7'
   if (atom_count .ge. 1000000)  atmc = 'i8'
   crdmin = 0.0_ti_p
   crdmax = 0.0_ti_p
   if (present(x_frame)) then
      do i = 1, atom_count
         crdmin = min(crdmin,x_frame(i),y_frame(i),z_frame(i))
         crdmax = max(crdmax,x_frame(i),y_frame(i),z_frame(i))
      end do
   else
!$acc update host(x,y,z) async
!$acc parallel loop present(x,y,z) async
   do i = 1, atom_count
      crdmin = min(crdmin,x(i),y(i),z(i))
      crdmax = max(crdmax,x(i),y(i),z(i))
   end do
!$acc wait
   end if
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
!
!     write out the number of atoms and the title
!
   if (ltitle .eq. 0) then
      fstr = '('//atmc//')'
      write (ixyz,fstr(1:4))  atom_count
   else
      fstr = '('//atmc//',2x,a)'
      write (ixyz,fstr(1:9))  atom_count,title(1:ltitle)
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
   if (pbcunwrap) then
      do i = 1, atom_count
         if (present(x_frame)) then
            xi = x_frame(i)
            yi = y_frame(i)
            zi = z_frame(i)
         else
            xi = x(i)
            yi = y(i)
            zi = z(i)
         end if
         call unwrap_pos(xi,yi,zi,pbcWrapIdx,i)
         write (ixyz,fstr)  i,name(i),xi,yi,zi,type(i),&
            &(i12(k,i),k=1,n12(i))
      end do
   else
      do i = 1, atom_count
         if (present(x_frame)) then
            write (ixyz,fstr)  i,name(i),x_frame(i),y_frame(i),z_frame(i),type(i),&
               &(i12(k,i),k=1,n12(i))
         else
            write (ixyz,fstr)  i,name(i),x(i),y(i),z(i),type(i),&
               &(i12(k,i),k=1,n12(i))
         end if
      end do
   end if
!
!     close the output unit if opened by this routine
!
   if (.not. opened)  close (unit=ixyz)
   call timer_exit(timer_io)
end
