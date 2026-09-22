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
#include "tinker_macro.h"
subroutine prtdyn
   use atomsMirror
   use boxes
   use files
   use inform   ,only:n_fwriten
   use moldyn   ,only:v,a,aalt
   use replicas
   use output   ,only:new_restart
   use restart  ,only:write_restart
   use status   ,only:failed,status_t
   use titles
   use timestat ,only:timer_io,timer_enter,timer_exit,quiet_timers
   implicit none
   integer :: lext
   character*10 :: ext
   character*240 :: dynfile
   character*3 :: numberreps
   type(status_t) :: status_value

!     if multiple replicas, the restart are numbered
   if (use_reps) write(numberreps,'(i3.3)') rank_reploc

   call timer_enter(timer_io)

!     Resolve naming here; the restart codec owns all DYN file I/O.
   if (new_restart.and.n_fwriten.ne.0) then
      if (n_fwriten.eq.-1) then
         if (use_reps) then
            dynfile = filename(1:leng)//'_reps'//numberreps//'_err.dyn'
         else
            dynfile = filename(1:leng)//'_err.dyn'
            call version(dynfile,'new')
         end if
      else
         lext = 6
         call numeral(n_fwriten,ext,lext)
         if (use_reps) then
            dynfile = filename(1:leng)//'_reps'//numberreps//'_.dyn'
         else
            dynfile = filename(1:leng)//'_'//ext(1:lext)//'.dyn'
         end if
         call version(dynfile,'new')
      end if
   else if (use_reps) then
      dynfile = filename(1:leng)//'_reps'//numberreps//'.dyn'
   else
      dynfile = filename(1:leng)//'.dyn'
   end if

   call write_restart(trim(dynfile),title,xbox,ybox,zbox,alpha,beta,gamma,&
      &x,y,z,v,a,aalt,status_value)
   if (failed(status_value)) then
      write (6,'(/,a)') trim(status_value%message)
      call fatal
   end if
   call timer_exit(timer_io)
end
