!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ##############################################################
!     ##                                                          ##
!     ##  subroutine readxyz  --  input of Cartesian coordinates  ##
!     ##                                                          ##
!     ##############################################################
!
!
!     "readxyz" gets a set of Cartesian coordinates from
!     an external disk file
!
!
!> @brief 
!> gets a set of Cartesian coordinates from
!> an external disk file
!> @param[in] ixyz: unit of file
subroutine readxyz (ixyz)
   use sizes
   use atmtyp
   use atoms
   use bound
   use boxes
   use couple
   use files
   use inform
   use iounit
   use titles
   implicit none
   integer i,j,k,m
   integer ixyz,nmax
   integer next,size
   integer first,last
   integer nexttext
   integer trimtext
   integer, allocatable :: list(:)
   real*8 :: xlen,ylen,zlen
   real*8 :: aang,bang,gang
   logical exist,opened
   logical quit,reorder
   logical clash
   character*240 xyzfile
   character*240 record
   character*240 string
!
!
!     initialize the total number of atoms in the system
!
   n = 0
!
!     open the input file if it has not already been done
!
   inquire (unit=ixyz,opened=opened)
   if (.not. opened) then
      xyzfile = filename(1:leng)//'.xyz'
      call version (xyzfile,'old')
      inquire (file=xyzfile,exist=exist)
      if (exist) then
         open (unit=ixyz,file=xyzfile,status='old')
         rewind (unit=ixyz)
      else
         write (iout,10)
10       format (/,' READXYZ  --  Unable to Find the Cartesian',&
         &' Coordinates File')
         call fatal
      end if
   end if
!
!     read first line and return if already at end of file
!
   quit = .false.
   abort = .true.
   size = 0
   do while (size .eq. 0)
      read (ixyz,20,err=80,end=80)  record
20    format (a240)
      size = trimtext (record)
   end do
   abort = .false.
   quit = .true.
!
!     parse the title line to get the number of atoms
!
   i = 0
   next = 1
   call gettext (record,string,next)
   read (string,*,err=80,end=80)  n
!
!     allocate global arrays
!
   if (allocated(x)) deallocate(x)
   allocate (x(n))
   if (allocated(y)) deallocate(y)
   allocate (y(n))
   if (allocated(z)) deallocate(z)
   allocate (z(n))
   if (allocated(xold)) deallocate(xold)
   allocate (xold(n))
   if (allocated(yold)) deallocate(yold)
   allocate (yold(n))
   if (allocated(zold)) deallocate(zold)
   allocate (zold(n))
   if (allocated(type)) deallocate(type)
   allocate (type(n))
   if (allocated(i12))  deallocate(i12)
   allocate (i12(maxvalue,n))
   if (allocated(n12)) deallocate(n12)
   allocate (n12(n))
   if (allocated(tag)) deallocate(tag)
   allocate (tag(n))
   if (allocated(name)) deallocate(name)
   allocate (name(n))
   x = 0d0
   y = 0d0
   z = 0d0
   xold = 0d0
   yold = 0d0
   zold = 0d0
   type = 0
   n12 = 0
   i12 = 0
   tag = 0
   maxbnd = 4*n
   maxang = 6*n
   maxtors = 18*n
   maxbitor = 54*n

!
!     extract the title and determine its length
!
   string = record(next:240)
   first = nexttext (string)
   last = trimtext (string)
   if (last .eq. 0) then
      title = ' '
      ltitle = 0
   else
      title = string(first:last)
      ltitle = trimtext (title)
   end if
!
!     check for too few or too many total atoms in the file
!
   if (n .le. 0) then
      write (iout,30)
30    format (/,' READXYZ  --  The Coordinate File Does Not',&
      &' Contain Any Atoms')
      call fatal
   end if
!
!     initialize coordinates and connectivities for each atom
!
   do i = 1, n
      tag(i) = 0
      name(i) = '   '
      x(i) = 0.0d0
      y(i) = 0.0d0
      z(i) = 0.0d0
      type(i) = 0
      do j = 1, maxvalue
         i12(j,i) = 0
      end do
   end do
!
!     read the coordinates and connectivities for each atom
!
   do i = 1, n
      next = 1
      size = 0
      do while (size .eq. 0)
         read (ixyz,50,err=80,end=80)  record
50       format (a240)
         size = trimtext (record)
         if (i .eq. 1) then
            next = 1
            call getword (record,name(i),next)
            if (name(i) .ne. '   ')  goto 60
            read (record,*,err=60,end=60)  xlen,ylen,zlen,&
            &aang,bang,gang
            size = 0
            xbox = xlen
            ybox = ylen
            zbox = zlen
            alpha = aang
            beta = bang
            gamma = gang
            use_bounds = .true.
!               call lattice
60          continue
         end if
      end do
      read (record,*,err=80,end=80)  tag(i)
      call getword (record,name(i),next)
      string = record(next:240)
      read (string,*,err=70,end=70)  x(i),y(i),z(i),type(i),&
      &(i12(j,i),j=1,maxvalue)
70    continue
   end do
   quit = .false.
80 continue
   if (.not. opened)  close (unit=ixyz)
!
!     an error occurred in reading the coordinate file
!
   if (quit) then
      write (iout,90)  i
90    format (/,' READXYZ  --  Error in Coordinate File at Atom',i6)
      call fatal
   end if
!
!     for each atom, count and sort its attached atoms
!
   do i = 1, n
      n12(i) = 0
      do j = maxvalue, 1, -1
         if (i12(j,i) .ne. 0) then
            n12(i) = j
            goto 100
         end if
      end do
100   continue
      call sort (n12(i),i12(1,i))
   end do
!
!     perform dynamic allocation of some local arrays
!
   nmax = 0
   do i = 1, n
      nmax = max(tag(i),nmax)
      do j = 1, n12(i)
         nmax = max(i12(j,i),nmax)
      end do
   end do
   allocate (list(nmax))
!
!     check for scrambled atom order and attempt to renumber
!
   reorder = .false.
   do i = 1, n
      list(tag(i)) = i
      if (tag(i) .ne. i)  reorder = .true.
   end do
   if (reorder) then
      write (iout,110)
110   format (/,' READXYZ  --  Atom Labels not Sequential,',&
      &' Attempting to Renumber')
      do i = 1, n
         tag(i) = i
         do j = 1, n12(i)
            i12(j,i) = list(i12(j,i))
         end do
         call sort (n12(i),i12(1,i))
      end do
   end if
!
!     perform deallocation of some local arrays
!
   deallocate (list)
!
!     check for atom pairs with identical coordinates
!
   clash = .false.
   if (n .le. 10000)  call chkxyz (clash)
!
!     make sure that all connectivities are bidirectional
!
   do i = 1, n
      do j = 1, n12(i)
         k = i12(j,i)
         do m = 1, n12(k)
            if (i12(m,k) .eq. i)  goto 130
         end do
         write (iout,120)  k,i
120      format (/,' READXYZ  --  Check Connection of Atom',&
         &i6,' to Atom',i6)
         call fatal
130      continue
      end do
   end do
   return
end
