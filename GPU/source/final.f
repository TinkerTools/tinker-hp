c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine final  --  final actions before program exit  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "final" performs any final program actions, prints a status
c     message, and then pauses if necessary to avoid closing the
c     execution window
c
c
#include "tinker_precision.h"
      subroutine final
      use deriv
      use domdec
      use inform
      use iounit
      use moldyn
      use neigh
      use mpi
      use pme
      use potent
      use uprior
      use sizes
      use timestat
      use utilgpu
#ifdef _OPENACC
      use cudafor
#endif
      implicit none
      integer i
c
      if (tinkertime) call display_timers

      if (.not.use_pmecore.or.use_pmecore.and.rank>=ndir) call fft_final
! Not supported for synchronous computation
c#ifdef _OPENACC
c      if (cudaStreamDestroy(dir_stream).ne.0)
c     &   print*, 'error destroying stream'
c      if (cudaStreamDestroy(rec_stream).ne.0)
c     &   print*, 'error destroying stream'
c#endif
c
c     perform deallocation of associated pointers arrays
c
!!$acc exit data delete(nvlst,nelst,elst,vlst)
      if (allocated(nvlst))      deallocate (nvlst)
      if (allocated(vlst))       deallocate (vlst)
      if (allocated(nelst))      deallocate (nelst)
      if (allocated(elst))       deallocate (elst)
      if (allocated(thetai1))    deallocate (thetai1)
      if (allocated(thetai2))    deallocate (thetai2)
      if (allocated(thetai3))    deallocate (thetai3)
      if (allocated(qgridin_2d))   deallocate (qgridin_2d)
      if (allocated(qgridout_2d))  deallocate (qgridout_2d)
      if (allocated(qgrid2in_2d))  deallocate (qgrid2in_2d)
      if (allocated(qgrid2out_2d)) deallocate (qgrid2out_2d)
      if (allocated(qfac_2d))    deallocate (qfac_2d)
      if (allocated(udalt))      deallocate (udalt)
      if (allocated(upalt))      deallocate (upalt)
      if (allocated(udshortalt))  deallocate (udshortalt)
      if (allocated(upshortalt))  deallocate (upshortalt)
!!$acc exit data delete(demrec,deprec)
      if (allocated(demrec))     deallocate (demrec)
      if (allocated(deprec))     deallocate (deprec)
      if (allocated(bufbegpole)) deallocate (bufbegpole)
      if (allocated(bufbegrec))  deallocate (bufbegrec)
      if (allocated(domlen))     deallocate (domlen)
      if (allocated(repart))     deallocate (repart)
      if (allocated(zbegproc))   deallocate (zbegproc)
      if (allocated(zendproc))   deallocate (zendproc)
      if (allocated(ybegproc))   deallocate (ybegproc)
      if (allocated(yendproc))   deallocate (yendproc)
      if (allocated(zbegproc))   deallocate (zbegproc)
      if (allocated(zendproc))   deallocate (zendproc)
      if (allocated(p_send1))    deallocate (p_send1)
      if (allocated(p_recep1))   deallocate (p_recep1)
      if (allocated(p_send2))    deallocate (p_send2)
      if (allocated(p_recep2))   deallocate (p_recep2)
      if (allocated(ptorque_send))  deallocate (ptorque_send)
      if (allocated(ptorque_recep)) deallocate (ptorque_recep)
      if (allocated(pneig_send))    deallocate (pneig_send)
      if (allocated(pneig_recep))   deallocate (pneig_recep)
      if (allocated(glob))     deallocate (glob)
      if (allocated(loc))      deallocate (loc)
      if (allocated(globrec))  deallocate (globrec)
      if (allocated(locrec))   deallocate (locrec)
      if (allocated(globrec1)) deallocate (globrec1)
      if (allocated(locrec1))  deallocate (locrec1)
      if (allocated(v))        deallocate (v)
      if (allocated(a))        deallocate (a)
      if (allocated(aalt))     deallocate (aalt)
c
!!$acc exit data delete(deb,dea,deba,deub,deaa,deopb,deopd,desum)
      if (allocated(desum)) deallocate (desum)
      if (allocated(deb))   deallocate (deb)
      if (allocated(dea))   deallocate (dea)
      if (allocated(deba))  deallocate (deba)
      if (allocated(deub))  deallocate (deub)
      if (allocated(deaa))  deallocate (deaa)
      if (allocated(deopb)) deallocate (deopb)
      if (allocated(deopd)) deallocate (deopd)
!!$acc exit data delete(deid,deit,det,dept,debt,dett,deg,dex)
      if (allocated(deid))  deallocate (deid)
      if (allocated(deit))  deallocate (deit)
      if (allocated(det))   deallocate (det)
      if (allocated(dept))  deallocate (dept)
      if (allocated(debt))  deallocate (debt)
      if (allocated(dett))  deallocate (dett)
      if (allocated(dev))   deallocate (dev)
      if (allocated(dem))   deallocate (dem)
      if (allocated(dep))   deallocate (dep)
      if (allocated(desmd)) deallocate (desmd)
      if (allocated(deamdD)) deallocate (deamdD)
      if (allocated(deamdP)) deallocate (deamdP)
      if (allocated(deW1aMD)) deallocate (deW1aMD)
      if (allocated(deW2aMD)) deallocate (deW2aMD)
c
c     print a final status message before exiting TINKER
c
      if (debug) then
         write (iout,10)
   10    format (/,' TINKER is Exiting following Normal Termination',
     &              ' of the Program',/)
      end if
c
c     may need a pause to avoid closing the execution window
c
      if (holdup) then
         read (input,20)
   20    format ()
      end if
      return
      end
