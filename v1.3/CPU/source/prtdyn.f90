!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###############################################################
!     ##                                                           ##
!     ##  subroutine prtdyn  --  output of MD restart information  ##
!     ##                                                           ##
!     ###############################################################
!
!
!     "prtdyn" writes out the information needed to restart a
!     molecular dynamics trajectory to an external disk file
!
!
!> @brief 
!> writes out the information needed to restart a
!> molecular dynamics trajectory to an external disk file
!> @param[in] suffix: suffix of the written file
subroutine prtdyn(suffix)
   use atoms
   use boxes
   use files
   use group
   use mdstuf
   use moldyn
   use titles
   implicit none
   character(*), intent(in) :: suffix
   integer i,idyn
   integer freeunit
   logical exist
   character*2 atmc
   character*40 fstr
   character*240 dynfile

!
!     update an existing restart file or open a new one
!
   idyn = freeunit ()
   dynfile = filename(1:leng)//trim(suffix)//'.dyn'
   inquire (file=dynfile,exist=exist)
   if (exist) then
      open (unit=idyn,file=dynfile,status='old')
      rewind (unit=idyn)
   else
      open (unit=idyn,file=dynfile,status='new')
   end if
!
!     save the number of atoms and the title string
!
   fstr = '('' Number of Atoms and Title :'')'
   write (idyn,fstr(1:32))
   atmc = 'i6'
   if (n .ge. 100000)  atmc = 'i7'
   if (n .ge. 1000000)  atmc = 'i8'
   if (ltitle .eq. 0) then
      fstr = '('//atmc//')'
      write (idyn,fstr(1:4))  n
   else
      fstr = '('//atmc//',2x,a)'
      write (idyn,fstr(1:9))  n,title(1:ltitle)
   end if
!
!     save the periodic box edge lengths and angles
!
   fstr = '('' Periodic Box Dimensions :'')'
   write (idyn,fstr(1:30))
   fstr = '(3d26.16)'
   write (idyn,fstr(1:9))  xbox,ybox,zbox
   write (idyn,fstr(1:9))  alpha,beta,gamma
!
!     save the atomic positions, velocities and accelerations
!
   fstr = '('' Current Atomic Positions :'')'
   write (idyn,fstr(1:31))
   fstr = '(3d26.16)'
   do i = 1, n
      write (idyn,fstr(1:9))  x(i),y(i),z(i)
   end do
   fstr = '('' Current Atomic Velocities :'')'
   write (idyn,fstr(1:32))
   fstr = '(3d26.16)'
   do i = 1, n
      write (idyn,fstr(1:9))  v(1,i),v(2,i),v(3,i)
   end do
   fstr =  '('' Current Atomic Accelerations :'')'
   write (idyn,fstr(1:36))
   fstr = '(3d26.16)'
   do i = 1, n
      write (idyn,fstr(1:9))  a(1,i),a(2,i),a(3,i)
   end do
   fstr =  '('' Alternate Atomic Accelerations :'')'
   write (idyn,fstr(1:38))
   fstr = '(3d26.16)'
   do i = 1, n
      write (idyn,fstr(1:9))  aalt(1,i),aalt(2,i),aalt(3,i)
   end do
!      fstr =  '('' Alternate 2 Atomic Accelerations :'')'
!      write (idyn,fstr(1:40))
!      fstr = '(3d26.16)'
!      do i = 1, n
!         write (idyn,fstr(1:9))  aalt2(1,i),aalt2(2,i),aalt2(3,i)
!      end do
   fstr =  '('' pbc wrap index :'')'
   write (idyn,fstr(1:21))
   fstr = '(3i3)'
   do i = 1, n
      write (idyn,fstr(1:9))  pbcwrapindex(1,i),pbcwrapindex(2,i)&
      &,pbcwrapindex(3,i)
   end do
!
!     close the dynamics trajectory restart file
!
   close (unit=idyn)
   return
end
