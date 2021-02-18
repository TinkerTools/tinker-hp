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
      subroutine final
      use deriv
      use domdec
      use inform
      use iounit
      use moldyn
      use neigh
      use pme
      use uprior
      implicit none
c
c      call fft_final
c
c     perform deallocation of allocated arrays
c
      if (allocated(nvlst))  deallocate (nvlst)
      if (allocated(vlst))  deallocate (vlst)
      if (allocated(nelst))  deallocate (nelst)
      if (allocated(elst))  deallocate (elst)
      if (allocated(thetai1))  deallocate (thetai1)
      if (allocated(thetai2))  deallocate (thetai2)
      if (allocated(thetai3))  deallocate (thetai3)
      if (allocated(qgridin_2d))  deallocate (qgridin_2d)
      if (allocated(qgridout_2d))  deallocate (qgridout_2d)
      if (allocated(qgrid2in_2d))  deallocate (qgrid2in_2d)
      if (allocated(qgrid2out_2d))  deallocate (qgrid2out_2d)
      if (allocated(qfac_2d))  deallocate (qfac_2d)
      if (allocated(udalt))  deallocate (udalt)
      if (allocated(upalt))  deallocate (upalt)
      if (allocated(udshortalt))  deallocate (udshortalt)
      if (allocated(upshortalt))  deallocate (upshortalt)
      if (allocated(demrec)) deallocate (demrec)
      if (allocated(deprec)) deallocate (deprec)
      if (allocated(bufbegpole)) deallocate (bufbegpole)
      if (allocated(bufbegrec)) deallocate (bufbegrec)
      if (allocated(domlen)) deallocate (domlen)
      if (allocated(repart)) deallocate (repart)
      if (allocated(zbegproc)) deallocate (zbegproc)
      if (allocated(zendproc)) deallocate (zendproc)
      if (allocated(ybegproc)) deallocate (ybegproc)
      if (allocated(yendproc)) deallocate (yendproc)
      if (allocated(zbegproc)) deallocate (zbegproc)
      if (allocated(zendproc)) deallocate (zendproc)
      if (allocated(p_send1)) deallocate (p_send1)
      if (allocated(p_recep1)) deallocate (p_recep1)
      if (allocated(p_send2)) deallocate (p_send2)
      if (allocated(p_recep2)) deallocate (p_recep2)
      if (allocated(ptorque_send)) deallocate (ptorque_send)
      if (allocated(ptorque_recep)) deallocate (ptorque_recep)
      if (allocated(pneig_send)) deallocate (pneig_send)
      if (allocated(pneig_recep)) deallocate (pneig_recep)
      if (allocated(glob)) deallocate (glob)
      if (allocated(loc)) deallocate (loc)
      if (allocated(nelst)) deallocate (nelst)
      if (allocated(elst)) deallocate (elst)
      if (allocated(nvlst)) deallocate (nvlst)
      if (allocated(vlst)) deallocate (vlst)
      if (allocated(globrec)) deallocate (globrec)
      if (allocated(locrec)) deallocate (locrec)
      if (allocated(globrec1)) deallocate (globrec1)
      if (allocated(locrec1)) deallocate (locrec1)
      if (allocated(v)) deallocate (v)
      if (allocated(a)) deallocate (a)
      if (allocated(aalt)) deallocate (aalt)
c
      if (allocated(desum)) deallocate (desum)
      if (allocated(deb)) deallocate (deb)
      if (allocated(dea)) deallocate (dea)
      if (allocated(deba)) deallocate (deba)
      if (allocated(deub)) deallocate (deub)
      if (allocated(deaa)) deallocate (deaa)
      if (allocated(deopb)) deallocate (deopb)
      if (allocated(deopd)) deallocate (deopd)
      if (allocated(deid)) deallocate (deid)
      if (allocated(deit)) deallocate (deit)
      if (allocated(det)) deallocate (det)
      if (allocated(dept)) deallocate (dept)
      if (allocated(deat)) deallocate (deat)
      if (allocated(debt)) deallocate (debt)
      if (allocated(dett)) deallocate (dett)
      if (allocated(dev)) deallocate (dev)
      if (allocated(dem)) deallocate (dem)
      if (allocated(dep)) deallocate (dep)
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
