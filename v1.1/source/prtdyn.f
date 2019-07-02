c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine prtdyn  --  output of MD restart information  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "prtdyn" writes out the information needed to restart a
c     molecular dynamics trajectory to an external disk file
c
c
      subroutine prtdyn
      use atoms
      use boxes
      use files
      use group
      use mdstuf
      use moldyn
      use titles
      implicit none
      integer i,idyn
      integer freeunit
      logical exist
      character*2 atmc
      character*39 fstr
      character*120 dynfile
c
c
c     update an existing restart file or open a new one
c
      idyn = freeunit ()
      dynfile = filename(1:leng)//'.dyn'
      inquire (file=dynfile,exist=exist)
      if (exist) then
         open (unit=idyn,file=dynfile,status='old')
         rewind (unit=idyn)
      else
         open (unit=idyn,file=dynfile,status='new')
      end if
c
c     save the number of atoms and the title string
c
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
c
c     save the periodic box edge lengths and angles
c
      fstr = '('' Periodic Box Dimensions :'')'
      write (idyn,fstr(1:30))
      fstr = '(3d26.16)'
      write (idyn,fstr(1:9))  xbox,ybox,zbox
      write (idyn,fstr(1:9))  alpha,beta,gamma
c
c     save the atomic positions, velocities and accelerations
c
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
c
c     close the dynamics trajectory restart file
c
      close (unit=idyn)
      return
      end
