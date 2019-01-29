c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine optsave  --  save optimization info and results  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "optsave" is used by the optimizers to write imtermediate
c     coordinates and other relevant information; also checks for
c     user requested termination of an optimization
c
c
      subroutine optsave (ncycle,f,xx)
      use atoms
      use files
      use iounit
      use math
      use output
      use scales
      use usage
      implicit none
      integer i,iopt,iend
      integer ncycle,nvar
      integer lext,freeunit
      real*8 f,xx(*)
      logical exist
      character*7 ext
      character*120 optfile
      character*120 endfile
c
c
c     nothing to do if coordinate type is undefined
c
      if (coordtype .eq. 'NONE')  return
c
c     check scaling factors for optimization parameters
c
      if (.not. set_scale) then
         set_scale = .true.
         if (coordtype .eq. 'CARTESIAN') then
            do i = 1, 3*n
               scale(i) = 1.0d0
            end do
         end if
c         else if (coordtype .eq. 'INTERNAL') then
c            do i = 1, nomega
c               scale(i) = 1.0d0
c            end do
c         end if
      end if
c
c     transform optimization parameters back to coordinates
c
      if (coordtype .eq. 'CARTESIAN') then
         nvar = 0
         do i = 1, n
            if (use(i)) then
               nvar = nvar + 1
               x(i) = xx(nvar) / scale(nvar)
               nvar = nvar + 1
               y(i) = xx(nvar) / scale(nvar)
               nvar = nvar + 1
               z(i) = xx(nvar) / scale(nvar)
            end if
         end do
      end if
c      else if (coordtype .eq. 'INTERNAL') then
c         do i = 1, nomega
c            dihed(i) = xx(i) / scale(i)
c            ztors(zline(i)) = dihed(i) * radian
c         end do
c      end if
c
c     get name of archive or intermediate coordinates file
c
      iopt = freeunit ()
      if (cyclesave) then
         if (archive) then
            optfile = filename(1:leng)
            call suffix (optfile,'arc','old')
            inquire (file=optfile,exist=exist)
            if (exist) then
               call openend (iopt,optfile)
            else
               open (unit=iopt,file=optfile,status='new')
            end if
         else
            lext = 3
            call numeral (ncycle,ext,lext)
            optfile = filename(1:leng)//'.'//ext(1:lext)
            call version (optfile,'new')
            open (unit=iopt,file=optfile,status='new')
         end if
      else
         optfile = outfile
         call version (optfile,'old')
         open (unit=iopt,file=optfile,status='old')
         rewind (unit=iopt)
      end if
c
c     update intermediate file with desired coordinate type
c
      if (coordtype .eq. 'CARTESIAN') then
         call prtxyz (iopt)
      end if
      close (unit=iopt)
c
c     test for requested termination of the optimization
c
      endfile = 'tinker.end'
      inquire (file=endfile,exist=exist)
      if (.not. exist) then
         endfile = filename(1:leng)//'.end'
         inquire (file=endfile,exist=exist)
         if (exist) then
            iend = freeunit ()
            open (unit=iend,file=endfile,status='old')
            close (unit=iend,status='delete')
         end if
      end if
      if (exist) then
         write (iout,10)
   10    format (/,' OPTSAVE  --  Optimization Calculation Ending',
     &              ' due to User Request')
         call fatal
      end if
      return
      end
