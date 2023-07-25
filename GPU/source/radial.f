c
c
c     ##############################################################
c     ##  COPYRIGHT (C) 1995 by Yong Kong and Jay William Ponder  ##
c     ##                   All Rights Reserved                    ##
c     ##############################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  program radial  --  compute radial distribution function  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "radial" finds the radial distribution function for a specified
c     pair of atom types via analysis of a set of coordinate frames
c
c
#include "tinker_precision.h"
      module radial_inl
        contains
#include "image.f.inc"
      end module

      program radial
      use mpi
#ifdef _OPENACC
      use utilgpu,only: bind_gpu
#endif
      implicit none
      integer ierr,nthreadsupport
#ifdef _OPENACC
      call bind_gpu
#endif
c      call MPI_INIT_THREAD(MPI_THREAD_MULTIPLE,nthreadsupport,ierr)
      call MPI_INIT(ierr)
      call radial_bis
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_FINALIZE(ierr)
      end

      subroutine radial_bis
      use argue
      use atmtyp
      use atoms
      use bound
      use boxes
      use files
      use inform
      use iounit
      use cutoff
      use math
      use molcul
      use potent
      use mpi
      use radial_inl
      use dcdmod
      use domdec
      implicit none
      integer i,j,k
      integer nframe,iframe
      integer iarc,next
      integer molj,molk
      integer numj,numk
      integer typej,typek
      integer start,stop
      integer step,skip
      integer nbin,bin
      integer, allocatable :: hist(:)
      real(r_p) xj,yj,zj
      real(t_p) dx,dy,dz
      real(r_p) rjk,rmax,width
      real(r_p) rlower,rupper
      real(r_p) factor,pairs
      real(r_p) volume,expect
      real(r_p), allocatable :: gr(:)
      real(r_p), allocatable :: gs(:)
      logical exist,query
      logical intramol
      character*1 answer
      character*3 namej,namek
      character*6 labelj,labelk
      character*240 record
      character*240 string
      character*240 xyzfile
      character*240 dcdfile
      integer freeunit
      type(dcdinfo_t) :: dcdinfo

c
c
c     perform the standard initialization functions
c
      call initial
      call initmpi
      if(nproctot/=1 .and. ranktot==0) then
        write(0,*) "Error: radial is a serial program"
        call fatal
      endif
c
c     open the trajectory archive and read the initial frame
c
      call getxyz
c
c     get the unitcell parameters and number of molecules
c
      call unitcell
c
c     set cutoffs small to enforce use of minimum images
c
      use_vdw = .true.
      use_charge = .false.
      use_mpole = .false.
      use_ewald = .false.
      call cutoffs
      call lattice
      call drivermpi
      call mechanic
c
c     get numbers of the coordinate frames to be processed
c
      start = 1
      stop = 100000
      step = 1
      query = .true.
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=10,end=10)  start
         query = .false.
      end if
      call nextarg (string,exist)
      if (exist)  read (string,*,err=10,end=10)  stop
      call nextarg (string,exist)
      if (exist)  read (string,*,err=10,end=10)  step
   10 continue
      if (query) then
         write (iout,20)
   20    format (/,' Numbers of First & Last Frame and Step',
     &              ' Increment :  ',$)
         read (input,30)  record
   30    format (a240)
         read (record,*,err=40,end=40)  start,stop,step
   40    continue
      end if
c
c     get the names of the atoms to be used in rdf computation
c
      call nextarg (labelj,exist)
      call nextarg (labelk,exist)
      if (.not. exist) then
         write (iout,50)
   50    format (/,' Enter 1st & 2nd Atom Names or Type Numbers :  ',$)
         read (input,60)  record
   60    format (a240)
         next = 1
         call gettext (record,labelj,next)
         call gettext (record,labelk,next)
      end if
c
c     convert the labels to either atom names or type numbers
c
      namej = '   '
      typej = -1
      read (labelj,*,err=70,end=70)  typej
   70 continue
      if (typej .le. 0) then
         next = 1
         call gettext (labelj,namej,next)
      end if
      namek = '   '
      typek = -1
      read (labelk,*,err=80,end=80)  typek
   80 continue
      if (typek .le. 0) then
         next = 1
         call gettext (labelk,namek,next)
      end if
c
c     get maximum distance from input or minimum image convention
c
      if (.not. use_bounds) then
         rmax = -1.0d0
         query = .true.
         call nextarg (string,exist)
         if (exist) then
            read (string,*,err=90,end=90)  rmax
            query = .false.
         end if
   90    continue
         if (query) then
            write (iout,100)
  100       format (/,' Enter Maximum Distance to Accumulate',
     &                 ' [10.0 Ang] :  ',$)
            read (input,110)  rmax
  110       format (f20.0)
         end if
         if (rmax .le. 0.0d0)  rmax = 10.0d0
      else if (octahedron) then
         rmax = (sqrt(3.0d0)/4.0d0) * xbox
         rmax = 0.95d0 * rmax
      else
         rmax = min(xbox2*beta_sin*gamma_sin,ybox2*gamma_sin,
     &                         zbox2*beta_sin)
         rmax = 0.95d0 * rmax
      end if
c
c     get the desired width of the radial distance bins
c
      width = -1.0d0
      query = .true.
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=120,end=120)  width
         query = .false.
      end if
  120 continue
      if (query) then
         write (iout,130)
  130    format (/,' Enter Width of Distance Bins [0.01 Ang] :  ',$)
         read (input,140)  width
  140    format (f20.0)
      end if
      if (width .le. 0.0d0)  width = 0.01d0
c
c     decide whether to restrict to intermolecular atom pairs
c
      intramol = .false.
      call nextarg (answer,exist)
      if (.not. exist) then
         write (iout,150)
  150    format (/,' Include Intramolecular Pairs in Distribution',
     &              ' [N] :  ',$)
         read (input,160)  record
  160    format (a240)
         next = 1
         call gettext (record,answer,next)
      end if
      call upcase (answer)
      if (answer .eq. 'Y')  intramol = .true.
c
c     count the number of coordinate frames in the archive file
c
      abort = .false.
      if (dcdio) then
        dcdfile = filename(1:leng)//'.dcd'
        call dcdfile_open(dcdinfo,dcdfile)
        call dcdfile_read_header(dcdinfo,.false.)
        call dcdfile_read_next(dcdinfo)
        call dcdfile_skip_next(dcdinfo,0)
      else
        iarc = freeunit ()
        xyzfile = filename
        call suffix (xyzfile,'xyz','old')
        open (unit=iarc,file=xyzfile,status ='old')
        rewind (unit=iarc)
        call readxyz (iarc)
      end if
c      nframe = 0
c      do while (.not. abort)
c         if (.not.dcdio) then
c           call readxyz (iarc)
c         else 
c           call dcdfile_read_next
c           call dcdfile_skip_next(0)
c         end if
c         nframe = nframe + 1
c      end do
c      nframe = nframe - 1      
c      if (.not.dcdio) then
c        rewind (unit=iarc)
c      else 
c        call dcdfile_close
c        call dcdfile_open(dcdfile)
c        call dcdfile_read_header(.false.)
c      end if
c      stop = min(nframe,stop)
c      nframe = (stop-start)/step + 1
c      write (iout,170)  nframe
c  170 format (/,' Number of Coordinate Frames :',i14)
c
c     set the number of distance bins to be accumulated
c
      nbin = int(rmax/width)
      write (iout,180)  nbin
  180 format (' Number of Distance Bins :',i18)
c
c     perform dynamic allocation of some local arrays
c
      allocate (hist(nbin))
      allocate (gr(nbin))
      allocate (gs(nbin))
c
c     zero out the distance bins and distribution functions
c
      do i = 1, nbin
         hist(i) = 0
         gr(i) = 0.0d0
         gs(i) = 0.0d0
      end do
!$acc wait
!$acc data copy(hist) copyin(name)
c
c     get the archived coordinates for each frame in turn
c
      write (iout,190)
  190 format (/,' Reading the Coordinates Archive File :',/)
      nframe = 0
      iframe = start
      skip = start
      do while (iframe.ge.start .and. iframe.le.stop)
         do j = 1, skip-1
            if (.not.dcdio) then
              call readxyz (iarc)
            else 
              call dcdfile_read_next(dcdinfo)
              call dcdfile_skip_next(dcdinfo,0)
            end if
         end do
         iframe = iframe + step
         skip = step
         if (.not.dcdio) then
           call readxyz (iarc)
         else 
           call dcdfile_read_next(dcdinfo)
           call dcdfile_skip_next(dcdinfo,0)
         end if
         if (.not. abort) then
            nframe = nframe + 1
            if (mod(nframe,10) .eq. 0) then
               write (iout,200)  nframe
  200          format (4x,'Processing Coordinate Frame',i13)
            end if
!$acc parallel loop collapse(2) default(present)
            do j = 1, n ; do k = 1, n
              if (name(j).eq.namej .or. type(j).eq.typej) then
                if (name(k).eq.namek .or. type(k).eq.typek) then
                  xj = x(j)
                  yj = y(j)
                  zj = z(j)
                  molj = molcule(j)                
                  if (j .ne. k) then
                    molk = molcule(k)
                    if (intramol .or. molj.ne.molk) then
                      dx = x(k) - xj
                      dy = y(k) - yj
                      dz = z(k) - zj
                      call image_inl (dx,dy,dz)
                      rjk = sqrt(dx*dx + dy*dy + dz*dz)
                      bin = int(rjk/width) + 1
                      if (bin .le. nbin) then
!$acc atomic update
                        hist(bin) = hist(bin) + 1
                      endif
                    end if
                  end if
                end if
              end if
            enddo ; enddo
         end if
      end do

!$acc end data

c
c     ensure a valid frame is loaded and report total frames
c
      if(abort) then
        if (.not.dcdio) then
          rewind (unit=iarc)
          call readxyz (iarc)
        else 
          call dcdfile_close(dcdinfo)
          call dcdfile_open(dcdinfo,dcdfile)
          call dcdfile_read_header(dcdinfo,.false.)
          call dcdfile_read_next(dcdinfo)
          call dcdfile_skip_next(dcdinfo,0)
        end if
      endif
      close (unit=iarc)
      if (mod(nframe,10) .ne. 0) then
         write (iout,210)  nframe
  210    format (4x,'Processing Coordinate Frame',i13)
      end if
c
c     count the number of occurrences of each atom type
c
      numj = 0
      numk = 0
      do i = 1, n
         if (name(i).eq.namej .or. type(i).eq.typej)  numj = numj + 1
         if (name(i).eq.namek .or. type(i).eq.typek)  numk = numk + 1
      end do
c
c     normalize the distance bins to give radial distribution
c
      if (numj.ne.0 .and. numk.ne.0) then
         factor = (4.0d0/3.0d0) * pi * dble(nframe)
         if (use_bounds) then
            pairs = dble(numj) * dble(numk)
            volume = (gamma_sin*gamma_term) * xbox * ybox * zbox
            if (octahedron)  volume = 0.5d0 * volume
            !if (dodecadron)  volume = volume / root2
            factor = factor * pairs / volume
         end if
         do i = 1, nbin
            rupper = dble(i) * width
            rlower = rupper - width
            expect = factor * (rupper**3 - rlower**3)
            gr(i) = dble(hist(i)) / expect
         end do
      end if
c
c     find the 5th degree polynomial smoothed distribution function
c
      if (nbin .ge. 5) then
         gs(1) = (69.0d0*gr(1) + 4.0d0*gr(2) - 6.0d0*gr(3)
     &             + 4.0d0*gr(4) - gr(5)) / 70.0d0
         gs(2) = (2.0d0*gr(1) + 27.0d0*gr(2) + 12.0d0*gr(3)
     &             - 8.0d0*gr(4) + 2.0d0*gr(5)) / 35.0d0
         do i = 3, nbin-2
            gs(i) = (-3.0d0*gr(i-2) + 12.0d0*gr(i-1) + 17.0d0*gr(i)
     &                + 12.0d0*gr(i+1) - 3.0d0*gr(i+2)) / 35.0d0
         end do
         gs(nbin-1) = (2.0d0*gr(nbin-4) - 8.0d0*gr(nbin-3)
     &                  + 12.0d0*gr(nbin-2) + 27.0d0*gr(nbin-1)
     &                       + 2.0d0*gr(nbin)) / 35.0d0
         gs(nbin) = (-gr(nbin-4) + 4.0d0*gr(nbin-3) - 6.0d0*gr(nbin-2)
     &                + 4.0d0*gr(nbin-1) + 69.0d0*gr(nbin)) / 70.0d0
         do i = 1, nbin
            gs(i) = max(0.0d0,gs(i))
         end do
      end if
c
c     output the final radial distribution function results
c
      write (iout,220)  labelj,labelk
  220 format (/,' Pairwise Radial Distribution Function :'
     &        //,7x,'First Name or Type :  ',a6,
     &           5x,'Second Name or Type :  ',a6)
      write (iout,230)
  230 format (/,5x,'Bin',9x,'Counts',7x,'Distance',7x,'Raw g(r)',
     &           4x,'Smooth g(r)',/)
      do i = 1, nbin
         write (iout,240)  i,hist(i),(dble(i)-0.5d0)*width,gr(i),gs(i)
  240    format (i8,i15,3x,f12.4,3x,f12.4,3x,f12.4)
      end do

c
c     perform deallocation of some local arrays
c
c     deallocate (hist)
c     deallocate (gr)
c     deallocate (gs)
c
c     perform any final tasks before program exit
c
      if (dcdio) then
        call dcdfile_close(dcdinfo)
      else
        close (unit=iarc)
      end if
      call final
      end
