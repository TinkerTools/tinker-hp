c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine readdyn  --  input of MD restart information  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "readdyn" get the positions, velocities and accelerations
c     for a molecular dynamics restart from an external disk file
c
c
      subroutine readdyn (idyn)
      use atoms
      use boxes
      use files
      use group
      use iounit
      use mdstuf
      use moldyn
      implicit none
      integer i,idyn,ndyn,ios
      logical exist,opened,quit
      character*240 dynfile
      character*240 record
c
c     open the input file if it has not already been done
c
      inquire (unit=idyn,opened=opened)
      if (.not. opened) then
         dynfile = filename(1:leng)//'.dyn'
         call version (dynfile,'old')
         inquire (file=dynfile,exist=exist)
         if (exist) then
            open (unit=idyn,file=dynfile,status='old')
            rewind (unit=idyn)
         else
            write (iout,10)
   10       format (/,' READDYN  --  Unable to Find the Dynamics',
     &                 ' Restart File')
            call fatal
         end if
      end if
c
c     initialize error handling during reading of the file
c
      i = 0
      quit = .true.
c
c     get the number of atoms and check for consistency
c
      read (idyn,20)
   20 format ()
      read (idyn,30)  record
   30 format (a240)
      read (record,*,err=250,end=250)  ndyn
      if (ndyn .ne. n) then
         write (iout,40)
   40    format (/,' READDYN  --  Restart File has Incorrect',
     &              ' Number of Atoms')
         call fatal
      end if
c
c     get the periodic box edge lengths and angles
c
      read (idyn,50)
   50 format ()
      read (idyn,60)  record
   60 format (a240)
      read (record,*,err=250,end=250)  xbox,ybox,zbox
      read (idyn,70)  record
   70 format (a240)
      read (record,*,err=250,end=250)  alpha,beta,gamma
      read (idyn,80)
   80 format ()
c
c     set the box volume and additional periodic box values
c
      call lattice
c
      quit = .true.
c
c     get the atomic positions, velocities and accelerations
c
      do i = 1, n
         read (idyn,160)  record
  160    format (a240)
         read (record,*,err=250,end=250)  x(i),y(i),z(i)
      end do
      read (idyn,170)
  170 format ()
      do i = 1, n
         read (idyn,180)  record
  180    format (a240)
         read (record,*,err=250,end=250)  v(1,i),v(2,i),v(3,i)
      end do
      read (idyn,190)
  190 format ()
      do i = 1, n
         read (idyn,200)  record
  200    format (a240)
         read (record,*,err=250,end=250)  a(1,i),a(2,i),a(3,i)
      end do
      read (idyn,210)
  210 format ()
      do i = 1, n
         read (idyn,220)  record
  220    format (a240)
         read (record,*,err=250,end=250)  aalt(1,i),aalt(2,i),
     &                                    aalt(3,i)
      end do
c      read (idyn,230)
c  230 format ()
c      do i = 1, n
c         read (idyn,240)  record
c  240    format (a240)
c         read (record,*,err=250,end=250)  aalt2(1,i),aalt2(2,i),
c     &                                    aalt2(3,i)
c      end do
      read (idyn,230,iostat=ios)
  230 format ()
      if (ios == 0) then
        do i = 1, n
           read (idyn,240)  record
  240       format (a240)
           read (record,*,err=250,end=250)  pbcwrapindex(1,i)
     &                     ,pbcwrapindex(2,i),pbcwrapindex(3,i)
        end do
      end if
      
  250 continue
      quit = .false.
      if (.not. opened)  close (unit=idyn)
c
c     report any error in reading the dynamics restart file
c
      if (quit) then
         write (iout,260)  i
  260    format (/,' READDYN  --  Error in Dynamics Restart',
     &              ' File at Atom',i6)
         call fatal
      end if
      return
      end
