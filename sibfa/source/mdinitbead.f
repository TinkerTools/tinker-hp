c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine mdinitbead  --  initialize a pimd trajectory  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "mdinitbead" initializes the velocities and positions
c     for a molecular pimd trajectory, including restarts
c
c
      subroutine mdinitbead(ibead,dt,restart)
      use atmtyp
      use atoms
      use bath
      use beads
      use bound
      use couple
      use domdec
      use files
      use keys
      use freeze
      use inform
      use iounit
      use langevin
      use math
      use mdstuf
      use molcul
      use moldyn
      use mpole
      use neigh
      use units 
      use uprior
      use adqtb
      use usage
      implicit none
      integer i,j,k,idyn,iglob,nh
      integer, intent(in) :: ibead
      integer next
      integer lext,freeunit
      integer ierr
      real*8 e
      real*8 maxwell,speed
      real*8 hmax,hmass
      real*8 sum,dmass
      real*8 normal
      real*8 vec(3)
      real*8 dt
      real*8, allocatable :: derivs(:,:)
      logical exist,heavy,restart
      character*7 ext
      character*20 keyword
      character*120 dynfile
      character*120 record
      character*120 string
      character*3 numberbeads
     
c
c     try to restart using prior velocities and accelerations
c
      write(numberbeads, '(i3.3)') globbead(ibead)
      dynfile = filename(1:leng)//'_beads'//numberbeads//'.dyn'
      call version (dynfile,'old')
      inquire (file=dynfile,exist=exist)
      restart = exist
      if (exist) then
         idyn = freeunit ()
         open (unit=idyn,file=dynfile,status='old')
         rewind (unit=idyn)
         call readdyn (idyn)
         close (unit=idyn)
c
c     Do the domain decomposition
c
         call drivermpi
         call reinitnl(0)
         call reassignpme(.false.)
         call mechanicstep(0)
c         call nblist(0)
c          call reassignpi(0)
      else
c
c     set velocities and accelerations for cartesian dynamics
c
c         allocate (derivs(3,nbloc))
c         derivs = 0d0
c        ! call allocstep
c         call gradient (epotpi_loc,derivs)
c         call commforces(derivs)
c
         do i = 1, nloc
            iglob = glob(i)
            if (use(iglob)) then
               speed = maxwell (mass(iglob),nbeads*kelvin)
               call ranvec (vec)
               do j = 1, 3
                  v(j,iglob) = speed * vec(j)
               end do
            else
               do j = 1, 3
                  v(j,iglob) = 0.0d0
               end do
            end if
         end do
c         deallocate (derivs)
c         if (nuse .eq. n)  call mdrest (0)
      end if
c
c     check for any prior dynamics coordinate sets
c
      i = 0
      exist = .true.
      do while (exist)
         i = i + 1
         lext = 3
         call numeral (i,ext,lext)
         dynfile = filename(1:leng)//'_beads'//numberbeads//
     $                        '.'//ext(1:lext)
         inquire (file=dynfile,exist=exist)
         if (.not.exist .and. i.lt.100) then
            lext = 2
            call numeral (i,ext,lext)
            dynfile = filename(1:leng)//'_beads'//numberbeads//
     $                        '.'//ext(1:lext)
            inquire (file=dynfile,exist=exist)
         end if
         if (.not.exist .and. i.lt.10) then
            lext = 1
            call numeral (i,ext,lext)
            dynfile = filename(1:leng)//'_beads'//numberbeads//
     $                        '.'//ext(1:lext)
            inquire (file=dynfile,exist=exist)
         end if
      end do
      nprior = i - 1
      return
      end
c
