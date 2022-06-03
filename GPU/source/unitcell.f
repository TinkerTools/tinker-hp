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
#include "tinker_precision.h"
      subroutine unitcell
      use bound
      use boxes
      use cutoff
      use domdec
      use keys
      use iounit
      use tinheader
      implicit none
      integer i,next
      real(t_p) boxmax,rmax
      logical nosymm
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     set the default values for periodic boundary conditions
c
      use_bounds = .false.
c      use_replica = .false.
cc
cc     set the default values for the unitcell variables
cc
c      xbox = 0.0_re_p
c      ybox = 0.0_re_p
c      zbox = 0.0_re_p
c      alpha = 0.0_ti_p
c      beta = 0.0_ti_p
c      gamma = 0.0_ti_p
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
      if (boxmax .ne. 0.0_ti_p)  use_bounds = .true.
c
c     set unspecified periodic boundary box lengths and angles
c
      if (use_bounds) then
         if (xbox .eq. 0.0_ti_p)  xbox = boxmax
         if (ybox .eq. 0.0_ti_p)  ybox = xbox
         if (zbox .eq. 0.0_ti_p)  zbox = xbox
         if (alpha .eq. 0.0_ti_p)  alpha = 90.0_ti_p
         if (beta .eq. 0.0_ti_p)  beta = 90.0_ti_p
         if (gamma .eq. 0.0_ti_p)  gamma = 90.0_ti_p
c
c     determine the general periodic boundary lattice type
c
         if (nosymm) then
            triclinic = .true.
         else if (alpha.eq.90.0_ti_p .and. beta.eq.90.0_ti_p
     &               .and. gamma.eq.90.0_ti_p) then
            orthogonal = .true.
         else if (alpha.eq.90.0_ti_p .and. gamma.eq.90.0_ti_p) then
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
            triclinic  = .false.
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
         xbox = 2.0_re_p * (rmax+ewaldcut)
         ybox = xbox
         zbox = xbox
         alpha = 90.0_ti_p
         beta = 90.0_ti_p
         gamma = 90.0_ti_p
         orthogonal = .true.
c         call lattice
c         boundary = 'NONE'
c         dens = 0.75_ti_p
      end if
c
c     return error for non orthorombic unit cells
c
      if ((alpha.ne.90.0_ti_p).or.(beta.ne.90.0_ti_p).or.
     &   (gamma.ne.90.0_ti_p)) then
        if (rank.eq.0) then
          write(iout,*) 'Only orthorombic unit cells are implemented'
        end if
        call fatal
      end if
c
!$acc enter data copyin(use_bounds)
c
!$acc update device(orthogonal,xbox,ybox,zbox,
!$acc&     alpha,beta,gamma,octahedron)

      end
