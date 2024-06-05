c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine unitcell  --  get periodic boundary conditions  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "unitcell" gets the periodic boundary box size and related
c     values from an external keyword file
c
c
      subroutine unitcell
      use bound
      use boxes
      use cutoff
      use domdec
      use keys
      use inform
      use iounit
      implicit none
      integer i,next
      real*8 boxmax,rmax
      logical nosymm
      character*20 keyword
      character*240 record
      character*240 string
c
      if (deb_Path) write(iout,*), 'unitcell '
c
c
c
c     set the default values for periodic boundary conditions
c
      use_bounds = .false.
c      use_replica = .false.
cc
cc     set the default values for the unitcell variables
cc
c      xbox = 0.0d0
c      ybox = 0.0d0
c      zbox = 0.0d0
c      alpha = 0.0d0
c      beta = 0.0d0
c      gamma = 0.0d0
      orthogonal = .false.
      monoclinic = .false.
      triclinic = .false.
      octahedron = .false.
      spacegrp = '          '
      nosymm = .false.
c
c     get keywords containing crystal lattice dimensions
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:240)
         if (keyword(1:7) .eq. 'X-AXIS ') then
            read (string,*,err=10,end=10)  xbox
         else if (keyword(1:7) .eq. 'Y-AXIS ') then
            read (string,*,err=10,end=10)  ybox
         else if (keyword(1:7) .eq. 'Z-AXIS ') then
            read (string,*,err=10,end=10)  zbox
         else if (keyword(1:7) .eq. 'A-AXIS ') then
            read (string,*,err=10,end=10)  xbox
         else if (keyword(1:7) .eq. 'B-AXIS ') then
            read (string,*,err=10,end=10)  ybox
         else if (keyword(1:7) .eq. 'C-AXIS ') then
            read (string,*,err=10,end=10)  zbox
         else if (keyword(1:6) .eq. 'ALPHA ') then
            read (string,*,err=10,end=10)  alpha
         else if (keyword(1:5) .eq. 'BETA ') then
            read (string,*,err=10,end=10)  beta
         else if (keyword(1:6) .eq. 'GAMMA ') then
            read (string,*,err=10,end=10)  gamma
         else if (keyword(1:11) .eq. 'OCTAHEDRON ') then
            octahedron = .true.
         else if (keyword(1:11) .eq. 'SPACEGROUP ') then
            call getword (record,spacegrp,next)
         else if (keyword(1:12) .eq. 'NO-SYMMETRY ') then
            nosymm = .true.
         end if
   10    continue
      end do
c
c     use periodic boundary conditions if a cell was defined
c
      boxmax = max(xbox,ybox,zbox)
      if (boxmax .ne. 0.0d0)  use_bounds = .true.
c
c     return error for truncated octahedron
c
      if (octahedron) then
        if (rank.eq.0) then
          write(iout,*) 'Truncated octahedron not implemented yet'
        end if
        call fatal
      end if
c
c     set unspecified periodic boundary box lengths and angles
c
      if (use_bounds) then
         if (xbox .eq. 0.0d0)  xbox = boxmax
         if (ybox .eq. 0.0d0)  ybox = xbox
         if (zbox .eq. 0.0d0)  zbox = xbox
         if (alpha .eq. 0.0d0)  alpha = 90.0d0
         if (beta .eq. 0.0d0)  beta = 90.0d0
         if (gamma .eq. 0.0d0)  gamma = 90.0d0
c
c     determine the general periodic boundary lattice type
c
         if (nosymm) then
            triclinic = .true.
         else if (alpha.eq.90.0d0 .and. beta.eq.90.0d0
     &               .and. gamma.eq.90.0d0) then
            orthogonal = .true.
         else if (alpha.eq.90.0d0 .and. gamma.eq.90.0d0) then
            monoclinic = .true.
         else
            triclinic = .true.
         end if
      end if
c
c     check for proper use of truncated octahedron boundary
c
      if (octahedron) then
         if (xbox.eq.ybox .and. xbox.eq.zbox .and. orthogonal) then
            orthogonal = .false.
            monoclinic = .false.
            triclinic = .false.
         else
            write (iout,20)
   20       format (/,' UNITCELL  --  Truncated Octahedron',
     &                 ' Incompatible with Defined Cell')
            call fatal
         end if
      end if
c
c     set the system extent for nonperiodic Ewald summation
c
      if (.not. use_bounds) then
         call extent (rmax)
         xbox = 2.0d0 * (rmax+max(ewaldcut,dewaldcut))
         ybox = xbox
         zbox = xbox
         alpha = 90.0d0
         beta = 90.0d0
         gamma = 90.0d0
         orthogonal = .true.
c         call lattice
c         boundary = 'NONE'
c         dens = 0.75d0
      end if
c
c     return error for non orthorombic unit cells
c
      if ((alpha.ne.90d0).or.(beta.ne.90d0).or.(gamma.ne.90d0)) then
        if (rank.eq.0) then
          write(iout,*) 'Only orthorombic unit cells are implemented'
        end if
        call fatal
      end if
c
      return
      end
