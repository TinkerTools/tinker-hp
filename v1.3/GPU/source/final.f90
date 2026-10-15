!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###############################################################
!     ##                                                           ##
!     ##  subroutine final  --  final actions before program exit  ##
!     ##                                                           ##
!     ###############################################################
!
!
!     "final" performs any final program actions, prints a status
!     message, and then pauses if necessary to avoid closing the
!     execution window
!
!
#include "tinker_precision.h"
subroutine final(final_status)
   use analyz
   use angle
   use atmlst
   use atmtyp
   use atoms
   use beads, only: close_polymer_trajectories, polymer_ctr, polymer_main
   use charge
   use couple
   use dcdmod
   use deriv
   use disp
   use divcon
   use domdec
   use energi, only: delete_energi_on_device
   use fft
   !use fft    ,only: free_FFTgrid_p,is_fftInit
   use freeze
   use group
   use inform
   use iounit
   use kopbnd
   use langevin
   use mdstuf1
   use moldyn
   use molcul
   use mdstate
   use mlpot, only: release_mlpot_buffers
   use mpi
   use mdevice, only: default_work, direct_work, finalize_device, phase, phase_active, &
      &reciprocal_work, sync_work
   use mdtraj , only: close_trajectory, dynamic_trajectory
   use mpole
   use neigh
   use pme
   use polar
   use polpot
   use potent
   use sizes
   use reduction, only: delete_reduction_buffers
   use status, only: enrich_status, failed, set_status, StatusDeviceError, status_t
   use cflux, only: release_cflux_workspaces
   use tmatxb, only: release_tmatxb_workspace
   use tinMemory, only: s_shmem
   use timestat
   use uprior
   use vdw
   use virial, only: delete_vir_device
#ifdef USE_NVSHMEM
   use nvshmem, only: nvshmem_finalize
#endif
#ifdef PLUMED
   use plumed, only: lplumed, ncount, plumed_finalize
#endif
   use parallel, only: active_polar_tickets
#ifdef _OPENACC
   use cudafor
   use thrust, only: thrust_clear_cache_memory
   use tinMemory, only: sd_prmem
#endif
   implicit none
#ifdef COLVARS
   external colvars_finalize
   character(len=8) :: repeat_colvars_finalize
#endif
   integer i
   integer pme_cleanup_code
   integer final_fail_code, final_fail_any
   integer final_red_stat
   logical cleanup_failed, final_failed
   type(status_t), optional, intent(inout) :: final_status
   type(status_t) :: trajectory_status
   type(status_t) :: first_cleanup_status
   type(status_t) :: dev_clean_stat

   final_failed = .false.

#ifdef PLUMED
!     finalize adapters while borrowed MPI resources still exist
   if (lplumed) then
      if (ncount <= 0) then
         write (iout,'(/,a)') 'PLUMED finalization found no eplumed call'
         final_failed = .true.
      else
         write (iout,'(/,a,i0)') 'PLUMED eplumed calls: ', ncount
      end if
   end if
   call plumed_finalize()
   if (lplumed) then
      write (iout,'(/,a)') 'PLUMED finalization left an active adapter'
      final_failed = .true.
   end if
#endif
#ifdef COLVARS
   call colvars_finalize()
   repeat_colvars_finalize = ''
   call get_environment_variable('TINKERHP_COLVARS_REPEAT_FINALIZE', &
      &repeat_colvars_finalize)
   if (trim(repeat_colvars_finalize) == '1') call colvars_finalize()
#endif
!
   cleanup_failed = .false.
   call close_trajectory(dynamic_trajectory,trajectory_status)
   if (failed(trajectory_status)) then
      first_cleanup_status = trajectory_status
      cleanup_failed = .true.
   end if
   call close_polymer_trajectories(polymer_main,trajectory_status)
   if (failed(trajectory_status) .and. .not.cleanup_failed) then
      first_cleanup_status = trajectory_status
      cleanup_failed = .true.
   end if
   call close_polymer_trajectories(polymer_ctr,trajectory_status)
   if (failed(trajectory_status) .and. .not.cleanup_failed) then
      first_cleanup_status = trajectory_status
      cleanup_failed = .true.
   end if
   if (cleanup_failed) then
      write (iout,'(/,a)') trim(first_cleanup_status%message)
      final_failed = .true.
   end if
   call timer_exit(timer_prog)
   if (tinkertime.eq.all_time) call display_timers
   if (fw_mds) call mds_prt
!
!     Free ressources
!
   if (active_polar_tickets() /= 0) then
      write (iout,'(/,a)') 'Finalization found an active opaque MPI ticket'
      final_failed = .true.
   end if

   if (is_fftInit) then
      call fft_final
      call free_FFTgrid_p
      is_fftInit = .false.
   end if
!$acc wait(direct_work%queue)
!$acc wait(reciprocal_work%queue)
!$acc wait
   call release_mlpot_buffers
!$acc wait
   call delete_reduction_buffers(direct_work%queue)
   call release_cflux_workspaces(direct_work%queue)
   call release_polar_workspace(direct_work%queue)
   call release_tmatxb_workspace(direct_work%queue)
   call scaling_postop
   call dealloc_shared_scaling
   call mem_free_deriv
   call thrust_clear_cache_memory(sd_prmem)
   pme_cleanup_code = 0
   call release_pme_device(reciprocal_work%queue, pme_cleanup_code)
   if (pme_cleanup_code /= MPI_SUCCESS) then
      write (iout,'(/,a,i0)') 'PME finalization request cleanup failed: ', pme_cleanup_code
      final_failed = .true.
   end if
   if (pme_device_ready) then
      write (iout,'(/,a)') 'PME finalization left active device state'
      final_failed = .true.
   end if
   call delete_energi_on_device()
   call delete_vir_device()
!$acc wait
!
!     Free MPI shared memory segments before platform teardown
!
   if (use_vdw   ) call dealloc_shared_vdw
   if (use_charge) call dealloc_shared_chg
   if (use_mpole ) call dealloc_shared_mpole
   if (use_polar ) call dealloc_shared_polar
   call dealloc_shared_angles
   call dealloc_shared_angang
   if (use_opbend) call delete_data_kopbend
   if (use_tors.or.use_tortor) call delete_data_ktors
   if (use_chgflx) call dealloc_shared_chgflx
   if (use_repuls) call dealloc_shared_rep
   if (use_chgtrn) call dealloc_shared_chgct
   if (use_disp)   call dealloc_shared_disp
   call dealloc_shared_mutate
   call delete_data_kgeom
   call dealloc_shared_active
   call delete_data_bonds
   call dealloc_shared_attach
   call dealloc_shared_bitors
   call dealloc_shared_katom
   call dealloc_shared_mol
   call dealloc_shared_torsions
   call delete_data_kopdist
   call dealloc_polar_predictor
#ifdef USE_NVSHMEM
   call nvshmem_finalize()
#endif
   if (phase == phase_active) then
      call finalize_device(dev_clean_stat)
      if (failed(dev_clean_stat)) then
         write (iout,'(/,a)') trim(dev_clean_stat%message)
         final_failed = .true.
      end if
   end if
   if (phase == phase_active .or. direct_work%queue >= 0 .or. &
      &reciprocal_work%queue >= 0 .or. default_work%queue >= 0 .or. &
      &sync_work%queue >= 0 .or. direct_work%stream /= 0 .or. &
      &reciprocal_work%stream /= 0 .or. default_work%stream /= 0 .or. &
      &sync_work%stream /= 0) then
      write (iout,'(/,a)') 'Finalization left an active GPU execution channel'
      final_failed = .true.
   end if
! Not supported for synchronous computation
!#ifdef _OPENACC
!      if (cudaStreamDestroy(dir_stream).ne.0)
!     &   print*, 'error destroying stream'
!      if (cudaStreamDestroy(rec_stream).ne.0)
!     &   print*, 'error destroying stream'
!#endif
!
!     module "analyz"
!
   if (allocated(aesum)) deallocate (aesum)
   if (allocated(aem))   deallocate (aem)
   if (allocated(aep))   deallocate (aep)
   if (allocated(aev))   deallocate (aev)
   if (allocated(aeb))   deallocate (aeb)
   if (allocated(aea))   deallocate (aea)
   if (allocated(aeba))  deallocate (aeba)
   if (allocated(aub))   deallocate (aub)
   if (allocated(aeaa))  deallocate (aeaa)
   if (allocated(aeopb)) deallocate (aeopb)
   if (allocated(aeopd)) deallocate (aeopd)
   if (allocated(aeid))  deallocate (aeid)
   if (allocated(aeit))  deallocate (aeit)
   if (allocated(aet))   deallocate (aet)
   if (allocated(aept))  deallocate (aept)
   if (allocated(aebt))  deallocate (aebt)
   if (allocated(aett))  deallocate (aett)
   if (allocated(aeg))   deallocate (aeg)
   if (allocated(aex))   deallocate (aex)
   if (allocated(aec))   deallocate (aec)
   if (allocated(aeat))  deallocate (aeat)
   if (allocated(aer)) deallocate (aer)
   if (allocated(aedsp)) deallocate (aedsp)
   if (allocated(aect)) deallocate (aect)
!
!     module "angle"
!
   if (allocated(angleloc)) deallocate (angleloc)
!
!     module "atmlst"
!
   if (allocated(bndglob))    deallocate (bndglob)
   if (allocated(angleglob))  deallocate (angleglob)
   if (allocated(torsglob))   deallocate (torsglob)
   if (allocated(bitorsglob)) deallocate (bitorsglob)
   if (allocated(strbndglob)) deallocate (strbndglob)
   if (allocated(ureyglob))   deallocate (ureyglob)
   if (allocated(angangglob)) deallocate (angangglob)
   if (allocated(opbendglob)) deallocate (opbendglob)
   if (allocated(opdistglob)) deallocate (opdistglob)
   if (allocated(impropglob)) deallocate (impropglob)
   if (allocated(imptorglob)) deallocate (imptorglob)
   if (allocated(pitorsglob)) deallocate (pitorsglob)
   if (allocated(strtorglob)) deallocate (strtorglob)
   if (allocated(angtorglob)) deallocate (angtorglob)
   if (allocated(tortorglob)) deallocate (tortorglob)
   if (allocated(vdwglob))    deallocate (vdwglob)
   if (allocated(poleglob))   deallocate (poleglob)
   if (allocated(polerecglob))deallocate (polerecglob)
   if (allocated(dispglob))   deallocate (dispglob)
   if (allocated(disprecglob))deallocate (disprecglob)
   if (allocated(chgglob))    deallocate (chgglob)
   if (allocated(chgrecglob)) deallocate (chgrecglob)
   if (allocated(molculeglob))deallocate (molculeglob)
   if (allocated(npfixglob))  deallocate (npfixglob)
   if (allocated(ndfixglob))  deallocate (ndfixglob)
   if (allocated(nafixglob))  deallocate (nafixglob)
   if (allocated(ntfixglob))  deallocate (ntfixglob)
   if (allocated(ngfixglob))  deallocate (ngfixglob)
   if (allocated(nchirglob))  deallocate (nchirglob)
   if (allocated(ratglob))    deallocate (ratglob)
   if (allocated(chgglobnl))  deallocate (chgglobnl)
   if (allocated(vdwglobnl))  deallocate (vdwglobnl)
   if (allocated(poleglobnl)) deallocate (poleglobnl)
   if (allocated(dispglobnl)) deallocate (dispglobnl)
!
!     module "atmtyp"
!
!
!     module "atoms"
!
   if (associated(type)) deallocate (type)
   if (allocated(x))     deallocate (x)
   if (allocated(y))     deallocate (y)
   if (allocated(z))     deallocate (z)
   if (allocated(xold))  deallocate (xold)
   if (allocated(yold))  deallocate (yold)
   if (allocated(zold))  deallocate (zold)
!
!     module "charge"
!
   if (allocated(chgloc))    deallocate (chgloc)
   if (allocated(chglocnl))  deallocate (chglocnl)
   if (allocated(chgrecloc)) deallocate (chgrecloc)
!
!     module "couple"
!
   if (allocated(n12)) deallocate (n12)
   if (allocated(i12)) deallocate (i12)
!c
!c     module "dcdmod"
!c
!      if (allocated(titlesdcd)) deallocate (titlesdcd)
!
!     module "deriv"
!
!
!     module "disp"
!
   if (allocated(displocnl)) deallocate(displocnl)
   if (allocated(disprecloc))deallocate(disprecloc)
!
!     module "divcon"
!
   if (allocated(grplst))   deallocate (grplst)
   if (allocated(atmofst))  deallocate (atmofst)
   if (allocated(kofst))    deallocate (kofst)
   if (allocated(npergrp))  deallocate (npergrp)
   if (allocated(knblist))  deallocate (knblist)
   if (allocated(point))    deallocate (point)
   if (allocated(klst))     deallocate (klst)
   if (allocated(zmat))     deallocate (zmat)
   if (allocated(means))    deallocate (means)
   if (allocated(oldmeans)) deallocate (oldmeans)
   if (allocated(ytab))     deallocate (ytab)
!
!     module "domdec"
!
   if (allocated(domlen))          deallocate (domlen)
   if (allocated(domlenrec))       deallocate (domlenrec)
   if (allocated(domlenpole))      deallocate (domlenpole)
   if (allocated(domlenpolerec))   deallocate (domlenpolerec)
   if (allocated(p_recep1))        deallocate (p_recep1)
   if (allocated(p_recep2))        deallocate (p_recep2)
   if (allocated(p_send1))         deallocate (p_send1)
   if (allocated(p_send2))         deallocate (p_send2)
   if (allocated(p_recepshort1))   deallocate (p_recepshort1)
   if (allocated(p_recepshort2))   deallocate (p_recepshort2)
   if (allocated(p_sendshort1))    deallocate (p_sendshort1)
   if (allocated(p_sendshort2))    deallocate (p_sendshort2)
   if (allocated(pneig_recep))     deallocate (pneig_recep)
   if (allocated(pneig_send))      deallocate (pneig_send)
   if (allocated(precdir_recep))   deallocate (precdir_recep)
   if (allocated(precdir_send))    deallocate (precdir_send)
   if (allocated(precdir_recep1))  deallocate (precdir_recep1)
   if (allocated(precdir_send1))   deallocate (precdir_send1)
   if (allocated(precdir_recep2))  deallocate (precdir_recep2)
   if (allocated(precdir_send2))   deallocate (precdir_send2)
   if (allocated(pbig_recep))      deallocate (pbig_recep)
   if (allocated(pbig_send))       deallocate (pbig_send)
   if (allocated(pbigshort_recep)) deallocate (pbigshort_recep)
   if (allocated(pbigshort_send))  deallocate (pbigshort_send)
   if (allocated(loc))             deallocate (loc)
   if (allocated(glob))            deallocate (glob)
   if (allocated(locrec))          deallocate (locrec)
   if (allocated(globrec))         deallocate (globrec)
   if (allocated(prec_recep))      deallocate (prec_recep)
   if (allocated(prec_send))       deallocate (prec_send)
   if (allocated(repartrec))       deallocate (repartrec)
   if (allocated(repart))          deallocate (repart)
   if (allocated(bufbeg))          deallocate (bufbeg)
   if (allocated(bufbegpole))      deallocate (bufbegpole)
   if (allocated(bufbegrec))       deallocate (bufbegrec)
   if (allocated(buflen1))         deallocate (buflen1)
   if (allocated(buflen2))         deallocate (buflen2)
   if (allocated(bufbeg1))         deallocate (bufbeg1)
   if (allocated(bufbeg2))         deallocate (bufbeg2)
   if (allocated(buf1))            deallocate (buf1)
   if (allocated(buf2))            deallocate (buf2)
   if (allocated(zbegproc))        deallocate (zbegproc)
   if (allocated(zendproc))        deallocate (zendproc)
   if (allocated(ybegproc))        deallocate (ybegproc)
   if (allocated(yendproc))        deallocate (yendproc)
   if (allocated(xbegproc))        deallocate (xbegproc)
   if (allocated(xendproc))        deallocate (xendproc)
!
!     module "fft"
!
   if (allocated(istart1)) deallocate (istart1)
   if (allocated(iend1))   deallocate (iend1)
   if (allocated(jstart1)) deallocate (jstart1)
   if (allocated(jend1))   deallocate (jend1)
   if (allocated(kstart1)) deallocate (kstart1)
   if (allocated(kend1))   deallocate (kend1)
   if (allocated(isize1))  deallocate (isize1)
   if (allocated(jsize1))  deallocate (jsize1)
   if (allocated(ksize1))  deallocate (ksize1)
   if (allocated(istart2)) deallocate (istart2)
   if (allocated(iend2))   deallocate (iend2)
   if (allocated(jstart2)) deallocate (jstart2)
   if (allocated(jend2))   deallocate (jend2)
   if (allocated(kstart2)) deallocate (kstart2)
   if (allocated(kend2))   deallocate (kend2)
   if (allocated(isize2))  deallocate (isize2)
   if (allocated(jsize2))  deallocate (jsize2)
   if (allocated(ksize2))  deallocate (ksize2)
!
!     module "freeze"
!
   if (allocated(buflenrat1)) deallocate (buflenrat1)
   if (allocated(buflenrat2)) deallocate (buflenrat2)
   if (allocated(bufbegrat1)) deallocate (bufbegrat1)
   if (allocated(bufbegrat2)) deallocate (bufbegrat2)
   if (allocated(bufrat1))    deallocate (bufrat1)
   if (allocated(bufrat2))    deallocate (bufrat2)
!
!     module "group"
!
   if (allocated(kgrp))    deallocate (kgrp)
   if (allocated(grplist)) deallocate (grplist)
   if (allocated(igrp))    deallocate (igrp)
   if (allocated(grpmass)) deallocate (grpmass)
   if (allocated(wgrp))    deallocate (wgrp)
!
!     module "kopbnd"
!
   if (allocated(jopb)) deallocate (jopb)
!
!     module "langevin"
!
   if (allocated(Rn)) deallocate (Rn)
!
!     module "moldyn"
!
   if (allocated(v))     deallocate (v)
   if (allocated(a))     deallocate (a)
   if (allocated(aalt))  deallocate (aalt)
   if (allocated(aalt2)) deallocate (aalt2)
!
!     module "neigh"
!
   if (allocated(nvlst))       deallocate (nvlst)
   if (allocated(vlst))        deallocate (vlst)
   if (allocated(nelst))       deallocate (nelst)
   if (allocated(elst))        deallocate (elst)
   if (allocated(nshortvlst))  deallocate (nshortvlst)
   if (allocated(shortvlst))   deallocate (shortvlst)
   if (allocated(nshortelst))  deallocate (nshortelst)
   if (allocated(shortelst))   deallocate (shortelst)
   if (allocated(ineignl))     deallocate (ineignl)
   if (allocated(neigcell))    deallocate (neigcell)
   if (allocated(numneigcell)) deallocate (numneigcell)
   if (allocated(repartcell))  deallocate (repartcell)
   if (allocated(cell_len))    deallocate (cell_len)
   if (allocated(indcell))     deallocate (indcell)
   if (allocated(bufbegcell))  deallocate (bufbegcell)
   if (allocated(xbegcell))    deallocate (xbegcell)
   if (allocated(ybegcell))    deallocate (ybegcell)
   if (allocated(zbegcell))    deallocate (zbegcell)
   if (allocated(xendcell))    deallocate (xendcell)
   if (allocated(yendcell))    deallocate (yendcell)
   if (allocated(zendcell))    deallocate (zendcell)
!
!     module "mpole"
!
   if (allocated(poleloc))    deallocate (poleloc)
   if (allocated(polelocnl))  deallocate (polelocnl)
   if (allocated(polerecloc)) deallocate (polerecloc)
!     if (allocated(zaxis))      deallocate (zaxis)
!     if (allocated(yaxis))      deallocate (yaxis)
!     if (allocated(xaxis))      deallocate (xaxis)
   if (allocated(rpole))      deallocate (rpole)
!
!     module mdstuf1
!
   call gpuFreeMdstuf1Data
!
!     module "polar"
!
   if (allocated(uind)) deallocate (uind)
   if (allocated(uinp)) deallocate (uinp)
!
!     module "polpot"
!
   if (allocated(residue)) deallocate (residue)
   if (allocated(munp))    deallocate (munp)
   if (allocated(efres))   deallocate (efres)
!
!     module "uprior"
!
   if (allocated(udalt))      deallocate (udalt)
   if (allocated(upalt))      deallocate (upalt)
   if (allocated(udshortalt)) deallocate (udshortalt)
   if (allocated(upshortalt)) deallocate (upshortalt)
!
!     module "vdw"
!
   if (allocated(vdwlocnl)) deallocate (vdwlocnl)
   if (s_shmem /= 0) then
      write (iout,'(/,a,i0)') 'Finalization left MPI shared memory bytes: ', s_shmem
      final_failed = .true.
   end if
   final_fail_code = 0
   if (final_failed) final_fail_code = 1
   call MPI_Allreduce(final_fail_code, final_fail_any, 1, MPI_INTEGER, &
      &MPI_MAX, MPI_COMM_WORLD, final_red_stat)
   if (final_red_stat /= MPI_SUCCESS .or. final_fail_any /= 0) &
      &final_failed = .true.
   if (final_failed) then
      write (iout,'(/,a)') 'finalization owner cleanup incomplete'
   else
      write (iout,'(/,a)') 'finalization owner cleanup complete'
   end if
   if (present(final_status) .and. final_failed) then
      if (failed(final_status)) then
         call enrich_status(final_status, 'Finalization owner cleanup incomplete')
      else
         call set_status(final_status, StatusDeviceError, &
            &'Finalization owner cleanup incomplete')
      end if
   end if
!
!     print a final status message before exiting TINKER-HP
!
   if (.not.final_failed .and. .not.final_message_suppressed .and. &
      &rank.eq.0) write (iout,10)
10 format (/,' TINKER-HP is Exiting following Normal Termination')
!
!     may need a pause to avoid closing the execution window
!
   if (.not.final_message_suppressed .and. holdup) then
      read (input,20)
20    format ()
   end if
   final_message_suppressed = .false.
   return
end
