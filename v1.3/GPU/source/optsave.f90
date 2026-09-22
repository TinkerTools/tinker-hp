!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ##################################################################
!     ##                                                              ##
!     ##  subroutine optsave  --  save optimization info and results  ##
!     ##                                                              ##
!     ##################################################################
!
!
!     "optsave" is used by the optimizers to write imtermediate
!     coordinates and other relevant information; also checks for
!     user requested termination of an optimization
!
!
#include "tinker_precision.h"
subroutine optsave (ncycle,xx)
   use atomsMirror
   use domdec
   use files
   use iounit
   use mdtraj, only: frame_view_t, write_snapshot
   use math
   use output
   use scales
   use status, only: failed, status_t
   use tinheader
   use usage
   use mpi
   implicit none
   integer i,iend,iglob
   integer ncycle
   integer lext,freeunit
   real(r_p) xx(*)
   logical exist
   character*7 ext
   character*240 optfile
   character*240 endfile
   type(frame_view_t) frame
   type(status_t) snapshot_status
!
!
!     nothing to do if coordinate type is undefined
!
   if (coordtype .eq. 'NONE')  return
!
!     check scaling factors for optimization parameters
!
   if (.not. set_scale) then
      set_scale = .true.
      if (coordtype .eq. 'CARTESIAN') then
!$acc parallel loop default(present) async
         do i = 1, 3*n
            scale(i) = 1.0_re_p
         end do
      end if
!         else if (coordtype .eq. 'INTERNAL') then
!            do i = 1, nomega
!               scale(i) = 1.0_re_p
!            end do
!         end if
   end if
!
!     transform optimization parameters back to coordinates
!
   if (coordtype .eq. 'CARTESIAN') then
!         nvar = 0
!$acc parallel loop default(present) async
      do i = 1, nloc
         iglob = glob(i)
         if (use(iglob)) then
!               nvar = nvar + 1
!               x(i) = xx(nvar) / scale(nvar)
            x(iglob) = xx(3*(iglob-1)+1) / scale(3*(iglob-1)+1)
            y(iglob) = xx(3*(iglob-1)+2) / scale(3*(iglob-1)+2)
            z(iglob) = xx(3*(iglob-1)+3) / scale(3*(iglob-1)+3)
!               nvar = nvar + 1
!               y(i) = xx(nvar) / scale(nvar)
!               nvar = nvar + 1
!               z(i) = xx(nvar) / scale(nvar)
         end if
      end do
   end if
!      else if (coordtype .eq. 'INTERNAL') then
!         do i = 1, nomega
!            dihed(i) = xx(i) / scale(i)
!            ztors(zline(i)) = dihed(i) * radian
!         end do
!      end if
!
!     get name of archive or intermediate coordinates file
!
   if (cyclesave) then
      if (archive) then
         optfile = filename(1:leng)
         call suffix (optfile,'arc','old')
      else
         lext = 3
         call numeral (ncycle,ext,lext)
         optfile = filename(1:leng)//'.'//ext(1:lext)
         call version (optfile,'new')
      end if
   else
      optfile = outfile
      call version (optfile,'old')
   end if
!
!     update intermediate file with desired coordinate type
!
   if (coordtype .eq. 'CARTESIAN') then
!$acc update host(x,y,z) async
!$acc wait
      frame = frame_view_t(atoms=n,host_global=.true.,x=x,y=y,z=z)
      call write_snapshot(trim(optfile),frame,snapshot_status,&
         &cyclesave .and. archive)
      if (failed(snapshot_status)) then
         write (iout,'(/,a)') trim(snapshot_status%message)
         call fatal
      end if
   end if
!
!     test for requested termination of the optimization
!
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
10    format (/,' OPTSAVE  --  Optimization Calculation Ending',&
         &' due to User Request')
      call fatal
   end if
   return
end
