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
      implicit none
      include 'sizes.i'
      include 'deriv.i'
      include 'inform.i'
      include 'iounit.i'
      include 'moldyn.i'
      include 'neigh.i'
      include 'pme.i'
      include 'uprior.i'
      include 'openmp.i'
c
c
c      call fft_final
c
c     perform deallocation of associated pointers arrays
c
      if (associated(nvlst))  deallocate (nvlst)
      if (associated(vlst))  deallocate (vlst)
      if (associated(nelst))  deallocate (nelst)
      if (associated(elst))  deallocate (elst)
      if (associated(thetai1))  deallocate (thetai1)
      if (associated(thetai2))  deallocate (thetai2)
      if (associated(thetai3))  deallocate (thetai3)
      if (associated(qgridin_2d))  deallocate (qgridin_2d)
      if (associated(qgridout_2d))  deallocate (qgridout_2d)
      if (associated(qgrid2in_2d))  deallocate (qgrid2in_2d)
      if (associated(qgrid2out_2d))  deallocate (qgrid2out_2d)
      if (associated(qfac_2d))  deallocate (qfac_2d)
      if (associated(udalt))  deallocate (udalt)
      if (associated(upalt))  deallocate (upalt)
c      if (associated(iuse))  deallocate (iuse)
c      if (associated(use))  deallocate (use)
      if (associated(demrec)) deallocate (demrec)
      if (associated(deprec)) deallocate (deprec)
      if (associated(torquerec)) deallocate (torquerec)
      if (associated(torquedir)) deallocate (torquedir)
      if (associated(torquerecp)) deallocate (torquerecp)
      if (associated(torquedirp)) deallocate (torquedirp)
      if (associated(bufbegpole)) deallocate (bufbegpole)
      if (associated(bufbegrec)) deallocate (bufbegrec)
      if (associated(domlen)) deallocate (domlen)
      if (associated(repart)) deallocate (repart)
      if (associated(zbegproc)) deallocate (zbegproc)
      if (associated(zendproc)) deallocate (zendproc)
      if (associated(ybegproc)) deallocate (ybegproc)
      if (associated(yendproc)) deallocate (yendproc)
      if (associated(zbegproc)) deallocate (zbegproc)
      if (associated(zendproc)) deallocate (zendproc)
      if (associated(p_send1)) deallocate (p_send1)
      if (associated(p_recep1)) deallocate (p_recep1)
      if (associated(p_send2)) deallocate (p_send2)
      if (associated(p_recep2)) deallocate (p_recep2)
      if (associated(ptorque_send)) deallocate (ptorque_send)
      if (associated(ptorque_recep)) deallocate (ptorque_recep)
      if (associated(pneig_send)) deallocate (pneig_send)
      if (associated(pneig_recep)) deallocate (pneig_recep)
      if (associated(glob)) deallocate (glob)
      if (associated(loc)) deallocate (loc)
      if (associated(nelst)) deallocate (nelst)
      if (associated(elst)) deallocate (elst)
      if (associated(nvlst)) deallocate (nvlst)
      if (associated(vlst)) deallocate (vlst)
      if (associated(globrec)) deallocate (globrec)
      if (associated(locrec)) deallocate (locrec)
      if (associated(globrec1)) deallocate (globrec1)
      if (associated(locrec1)) deallocate (locrec1)
      if (associated(v)) deallocate (v)
      if (associated(a)) deallocate (a)
      if (associated(aalt)) deallocate (aalt)
      if (associated(desum)) deallocate (desum)
      if (associated(deb)) deallocate (deb)
      if (associated(dea)) deallocate (dea)
      if (associated(deba)) deallocate (deba)
      if (associated(deub)) deallocate (deub)
      if (associated(deaa)) deallocate (deaa)
      if (associated(deopb)) deallocate (deopb)
      if (associated(deopd)) deallocate (deopd)
      if (associated(deid)) deallocate (deid)
      if (associated(deit)) deallocate (deit)
      if (associated(det)) deallocate (det)
      if (associated(dept)) deallocate (dept)
      if (associated(debt)) deallocate (debt)
      if (associated(dett)) deallocate (dett)
      if (associated(dev)) deallocate (dev)
      if (associated(dem)) deallocate (dem)
      if (associated(dep)) deallocate (dep)
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
