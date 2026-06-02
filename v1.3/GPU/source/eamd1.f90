!
!           ##          ##    ##    ###
!          #  #         # #  # #    #  ##
!         #    #        #  ##  #    #    #
!        ########       #      #    #    #
!       #        #      #      #    #    #
!      #          #     #      #    #  ##
!     #            #    #      #    ###
!
!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###############################################################
!     ##                                                           ##
!     ##  subroutine aMD                                           ##
!     ##                                                           ##
!     ###############################################################
!
!
!     "aMD" calculates the aMD potentials and their respective
!     factors applied on the corresponded derivatives
!
#include "tinker_macro.h"
subroutine aMD (derivs,epot)
   use deriv
   use domdec
   use energi
   use mamd
   use mpi
   use potent
   use virial
   implicit none
   integer i, j
   real(r_p) derivs(3,nloc)
   real(r_p) epot
!
!     => GaMD case
!
   if (use_gamd.or.use_amd_ene.or.use_amd_dih.or.use_amd_wat1)&
      &then
      call resetForcesAMD
   end if

   if (use_gamd) call GaMD (derivs,epot)
!
!     => Classical aMD case
!
   if (.not. use_gamd) then                      ! LOOP0
      if (use_amd_dih .or. use_amd_ene) then        ! LOOP1
!$acc wait
!$acc update host(epot,eDaMD,eW1aMD)
         ePaMD = epot
         amd_factor_dih = 1.0
         amd_factor_tot = 1.0
         call eamd1 (epot)
!$acc update device(epot)
         if (use_amd_dih .and. .not. use_amd_ene) then
            if (amd_factor_dih .lt. 1.0) then
!$acc parallel loop collapse(2) async &
!$acc         default(present)
               do i = 1, nloc
                  do j = 1, 3
                     derivs(j,i) = derivs(j,i) + deamdD(j,i) *&
                        &(amd_factor_dih - 1.0)
                  end do
               end do
!           do i = 1, 3
!              do j = 1, 3
!                 vir(i,j) = vir(i,j) + viramdD(i,j)*
!    $            (amd_factor_dih - 1.0)
!              end do
!           end do
            end if
         end if
         if (use_amd_ene) then
            if (amd_factor_tot .lt. 1.0) then
!$acc parallel loop collapse(2) async &
!$acc         default(present)
               do i = 1, nloc
                  do j = 1, 3
                     derivs(j,i) = derivs(j,i) * amd_factor_tot
                  end do
               end do
!           do i = 1, 3
!              do j = 1, 3
!                 vir(i,j) = vir(i,j) + vir(i,j)*(amd_factor_tot-1.0)
!              end do
!           end do
            end if
         end if
         if (use_amd_dih .and. use_amd_ene) then
            if (amd_factor_dih .lt. 1.0) then
!$acc parallel loop collapse(2) async &
!$acc         default(present)
               do i = 1, nloc
                  do j = 1, 3
                     derivs(j,i) = derivs(j,i) + deamdD(j,i) *&
                        &(amd_factor_dih - amd_factor_tot)
                  end do
               end do
!           do i = 1, 3
!              do j = 1, 3
!                 vir(i,j) = vir(i,j) - viramdD(i,j)*(amd_factor_dih
!    $            - 1.0) + viramdD(i,j)*(amd_factor_tot-1.0)
!              end do
!           end do
            end if
         end if
      end if       ! END LOOP1
   end if       ! END LOOP0
   return
end

subroutine commforces_energy_aMD(energy,derivs)
   use domdec
   use mamd
   use mpi
   use potent
   use utilgpu,only: rec_queue
   implicit none

   real(r_p),intent(inout):: energy
   real(r_p),intent(inout):: derivs(3,nbloc)

   integer i,j,k
   integer commloc,ierr
   integer :: reqsend(nproc),reqrec(nproc)
   real(r_p),allocatable ::buffer(:,:,:)
   integer tagmpi
   integer status(MPI_STATUS_SIZE)
!
!     Communication of the total et
!
   if (use_pmecore) then
      commloc = comm_dir
      if (ndir.eq.1.or.rank.gt.ndir-1) return
   else
      commloc = COMM_TINKER
      if (nproc.eq.1) return
   end if
!
   call MPI_ALLREDUCE(MPI_IN_PLACE,energy,1,MPI_RPREC,MPI_SUM,&
      &commloc,ierr)
!
   allocate (buffer(3,max(1,nloc),nbig_send))
!$acc data create(buffer) async
!
!     MPI : begin reception in buffer
!
!$acc host_data use_device(buffer)
   do i = 1, nbig_send
      tagmpi = nproc*rank + pbig_send(i) + 1
      call MPI_IRECV(buffer(1,1,i),3*nloc,MPI_RPREC,pbig_send(i),&
         &tagmpi,commloc,reqrec(i),ierr)
   end do
!$acc end host_data
!
!$acc wait
!$acc host_data use_device(derivs)
   do i = 1, nbig_recep
      tagmpi = nproc*pbig_recep(i) + rank + 1
      call MPI_ISEND(derivs(1,bufbeg(pbig_recep(i)+1)),&
         &3*domlen(pbig_recep(i)+1),MPI_RPREC,pbig_recep(i),&
         &tagmpi,commloc,reqsend(i),ierr)
   end do
!$acc end host_data
!
   do i = 1, nbig_send
      call MPI_WAIT(reqrec(i),status,ierr)
   end do
   do i = 1, nbig_recep
      call MPI_WAIT(reqsend(i),status,ierr)
   end do
!
!     MPI : move in global arrays
!
   do i = 1, nbig_send
!$acc parallel loop collapse(2) async default(present)
      do j = 1, nloc
         do k = 1, 3
            derivs(k,j) = derivs(k,j) + buffer(k,j,i)
         end do
      end do
   end do
!
!     Deallocation of arrays
!
!$acc end data
   deallocate(buffer)
end subroutine
!
!     ###############################################################
!     ##                                                           ##
!     ##  subroutine GaMD                                          ##
!     ##                                                           ##
!     ###############################################################
!
!
!     "GaMD" calculates the GaMD potentials and their respective
!     factors applied on the corresponded derivatives
!
subroutine GaMD (derivs,epot)
   use atoms
   use deriv
   use domdec
   use energi
   use inform ,only: deb_Path
   use mamd
   use mpi
   use potent
   use sizes  ,only:tinkerdebug
   use virial
   implicit none
   integer i, j
   integer iglob
   real(r_p) derivs(3,nloc)
   real(r_p) epot
!
   if (use_amd_dih .or. use_amd_ene .or. use_amd_wat1) then        ! LOOP1
      if (deb_Path) write(*,'(A)') "GaMD"
!$acc wait
!$acc update host(epot,eDaMD,eW1aMD)
      ePaMD = epot
      gamd_factor_dih = 1.0
      gamd_factor_tot = 1.0
      gamd_factor_wat1 = 1.0
      call eamd1 (epot)
!$acc update device(epot) async
      if (use_amd_dih .and. .not. use_amd_ene) then
         if (gamd_factor_dih .lt. 1.0) then
!$acc parallel loop collapse(2) async &
!$acc         default(present)
            do i = 1, nloc
               do j = 1, 3
                  derivs(j,i) = derivs(j,i) + deamdD(j,i) *&
                     &(gamd_factor_dih - 1.0)
               end do
            end do
!              do i = 1, 3
!                 do j = 1, 3
!                    vir(i,j) = vir(i,j) + viramdD(i,j)*
!    $               (gamd_factor_dih - 1.0)
!                 end do
!              end do
         end if
      end if
      if (use_amd_ene) then
         if (gamd_factor_tot .lt. 1.0) then
!$acc parallel loop collapse(2) async &
!$acc         default(present)
            do i = 1, nloc
               do j = 1, 3
                  derivs(j,i) = derivs(j,i) * gamd_factor_tot
               end do
            end do
!             do i = 1, 3
!                do j = 1, 3
!                   vir(i,j) = vir(i,j) + vir(i,j)*(gamd_factor_tot-1.0)
!                end do
!             end do
         end if
      end if
      if (use_amd_dih .and. use_amd_ene) then
         if (gamd_factor_dih .lt. 1.0) then
!$acc parallel loop collapse(2) async &
!$acc         default(present)
            do i = 1, nloc
               do j = 1, 3
                  derivs(j,i) = derivs(j,i) + deamdD(j,i) *&
                     &(gamd_factor_dih - gamd_factor_tot)
               end do
            end do
!             do i = 1, 3
!                do j = 1, 3
!                   vir(i,j) = vir(i,j) - viramdD(i,j)*(gamd_factor_dih
!    $              - 1.0) + viramdD(i,j)*(gamd_factor_tot-1.0)
!                end do
!             end do
         end if
      end if
      if (use_amd_wat1) then
         if (gamd_factor_wat1 .lt. 1.0) then
!$acc parallel loop collapse(2) async &
!$acc         default(present)
            do i = 1, nloc
               do j = 1, 3
                  iglob = glob(i)
                  if (type(iglob)==aMDwattype(1) .or. type(iglob)==&
                     &aMDwattype(2)) then
                     derivs(j,i) = derivs(j,i) + deW1aMD(j,i) *&
                        &(gamd_factor_wat1 - 1.0)
                  end if
               end do
            end do
!              do i = 1, 3
!                 do j = 1, 3
!                    vir(i,j) = vir(i,j) + viramdW1(i,j)*
!    $                         (gamd_factor_wat1 - 1.0)
!                 end do
!              end do
         end if
      end if
   end if               ! END LOOP1
   return
end
!
!     ###############################################################
!     ##                                                           ##
!     ##  subroutine eamd1  --  make aMD calculations              ##
!     ##                                                           ##
!     ###############################################################
!
!
!     "eamd1" assigns the biased forced on each concerned atoms and
!     stores the amd and derivatives energy components
!
!
subroutine eamd1 (epot)
   use atmtyp
   use atoms
   use domdec
   use deriv
   use energi
   use inform ,only: deb_Path
   use mamd
   use mdstuf
   use mpi
   use potent
   use usage
   use sizes  ,only: tinkerdebug
   use tinheader
   implicit none
   integer ierr,commloc
   integer iglob
   integer i,j,c
   real(r_p) dedx_amd, dedy_amd, dedz_amd
   real(r_p) dedx_gamd, dedy_gamd, dedz_gamd
   real(r_p) energy_amd
   real(r_p) epot
!     real(r_p) derivs(3,nbloc)
!
!     Initialization of the values
!
   if (deb_Path) write(*,'(2x,a)') "eamd1"
   energy_amd = 0.0d0
   amdboost = 0.0d0
!
!     ################################################
!     ##### FIRST CASE: aMD applied on dihedrals #####
!     ################################################
!
   if (.not. use_gamd) then                               ! LOOP0
      if (use_amd_dih) then                                  ! LOOP1
         call commforces_energy_aMD(eDaMD,deamdD)
!
!     Distribution of the aMD force on each dihedral atoms
!
         if (eDaMD .gt. amd_dih_ene) energy_amd = 0.0d0
         if (eDaMD .lt. amd_dih_ene) then
            amd_factor_dih = amd_dih_alpha/(amd_dih_alpha +&
               &amd_dih_ene - eDaMD )
            amd_factor_dih = amd_factor_dih**2
            energy_amd = ((amd_dih_ene - eDaMD)**2)/(amd_dih_alpha +&
               &(amd_dih_ene - eDaMD))
            amdboost = energy_amd
            amdboostavg_dih = amdboostavg_dih + energy_amd
            epot = epot + energy_amd
         end if
         if (eDaMD .ge. amd_dih_ene) then
            amd_factor_dih = 1.00
            amdboost = 0.0d0
         end if
!
!     aMD output data
!
         tamd_dih = tamd_dih + aMDdt_dih
         amdtpass_dih = amdtpass_dih + 1
         c = modulo(amdtpass_dih,amdoutputFreq_dih)
         if (rank.eq.0) then
            if (c == 0 .and. .not. use_amd_ene) then
               open (iamdout,file='aMD_output.dat',position='append',&
                  &action="write",status='old')
               write (iamdout,45)  tamd_dih, amdboost,&
                  &amdboostavg_dih/real(amdoutputFreq_dih)
45             format ('aMD-D ',3f10.4)
               amdboostavg_dih = 0.0d0
               close(iamdout)
            end if
         end if
!
!          tamd_dih = tamd_dih + aMDdt_dih
!          amdtpass_dih = amdtpass_dih + 1
!
      end if            ! END LOOP1
!
!     ###################################################
!     ##### SECOND CASE: aMD applied in a dual mode #####
!     ###################################################
!
      if (use_amd_dih .and. use_amd_ene) then     ! LOOP2
!
         ePaMD = ePaMD - eDaMD
!
!     Check if aMD is applied on potential energy
!
         if (ePaMD .ge. amd_ep_ene) then
            amd_factor_tot = 1.00
            energy_amd = 0.0d0
         end if
         if (ePaMD .lt. amd_ep_ene) then       ! LOOP1
            amd_factor_tot = amd_ep_alpha/(amd_ep_alpha&
               &+ amd_ep_ene - ePaMD)
            amd_factor_tot = amd_factor_tot**2
            energy_amd = ((amd_ep_ene - ePaMD)**2)/(amd_ep_alpha +&
               &(amd_ep_ene - ePaMD))
            amdboost = amdboost + energy_amd
            amdboostavg = amdboostavg + energy_amd
            epot = epot + energy_amd
         end if
!
!     aMD output data
!
         tamd = tamd + aMDdt
         amdtpass = amdtpass + 1
         c = modulo(amdtpass,amdoutputFreq)
         if (rank.eq.0) then
            if (c == 0) then
               open (iamdout,file='aMD_output.dat',position='append',&
                  &action="write",status='old')
               write (iamdout,46)  tamd, amdboost,&
                  &amdboostavg/real(amdoutputFreq)
46             format ('aMD-P ',3f10.4)
               amdboostavg = 0.0d0
               close(iamdout)
            end if
         end if
!
!         tamd = tamd + aMDdt
!         amdtpass = amdtpass + 1
!
      end if            ! END LOOP2
!
!     ################################################################
!     ##### THIRD CASE: aMD applied only on the potential energy #####
!     ################################################################
!
      if (use_amd_ene .and. .not. use_amd_dih) then     ! LOOP3
!
!     Check if aMD is applied on potential energy
!
         if (ePaMD .ge. amd_ep_ene) then
            amd_factor_tot = 1.00
            energy_amd = 0.0d0
         end if
         if (ePaMD .lt. amd_ep_ene) then       ! LOOP1
            amd_factor_tot = amd_ep_alpha/(amd_ep_alpha&
               &+ amd_ep_ene - ePaMD)
            amd_factor_tot = amd_factor_tot**2
            energy_amd = ((amd_ep_ene - ePaMD)**2)/(amd_ep_alpha +&
               &(amd_ep_ene - ePaMD))
            amdboost = amdboost + energy_amd
            amdboostavg = amdboostavg + energy_amd
            epot = epot + energy_amd
         end if
!
!     aMD output data
!
         tamd = tamd + aMDdt
         amdtpass = amdtpass + 1
         c = modulo(amdtpass,amdoutputFreq)
         if (rank.eq.0) then
            if (c == 0) then
               open (iamdout,file='aMD_output.dat',position='append',&
                  &action="write",status='old')
               write (iamdout,46)  tamd, amdboost,&
                  &amdboostavg/real(amdoutputFreq)
               amdboostavg = 0.0d0
               close(iamdout)
            end if
         end if
!
!         tamd = tamd + aMDdt
!         amdtpass = amdtpass + 1
!
      end if            ! END LOOP3
!
   end if            ! END LOOP0
!
!     ###############################################
!     #### FOURTH CASE: Use of the GaMD procedure ###
!     ###############################################
!
   if (use_gamd) then        ! LOOP1
!
      if (tamd .ne. 0) then
         tamd = tamd + aMDdt
         amdtpass_dih = amdtpass_dih + 1
         amdtpass = amdtpass + 1
      end if
!
!##################################################################
!     STEP 1: Preparation for the cMD part
!##################################################################
!
      if (amdtpass .lt. gamd_cmdprepsteps) then
!     No data are collected
         gamdstep = 1
         call GaMD_update_stat(2)
      end if
!
!##################################################################
!     STEP 2: cMD production
!
!     => At this point, communications are mandaroty for each step
!     => We will perform them here !
!##################################################################
!
!     Communication of the dihedral energy eDaMD
!
      if (use_amd_dih) then             ! LOOP COMM DIH
         call commforces_energy_aMD(eDaMD,deamdD)
      end if                            ! END LOOP COMM DIH
!
!     Communication of the bond/angle energy eW1aMD
!
      if (use_amd_wat1) then             ! LOOP COMM WAT1
         call commforces_energy_aMD(eW1aMD,deW1aMD)
      end if                            ! END LOOP COMM WAT1
!
!     Update of the stats
!
      if (amdtpass .lt. gamd_cmdsteps .and. amdtpass .ge.&
         &gamd_cmdprepsteps) then
         gamdstep = 2
         call GaMD_update_stat(1)
      end if
!
      if (amdtpass .eq. gamd_cmdsteps) call GaMD_calc_E
!
!##################################################################
!     STEP 3: Preparation for the GaMD equilibration
!##################################################################
!
      if (amdtpass .lt. gamd_eqprepsteps .and. amdtpass .ge.&
         &gamd_cmdsteps) then
!
         gamdstep = 3
!
!        Apply the boost potential without updating parameters gamdE,
!        gamdk and gamdk0
!
         call GaMD_calc_force (epot)
!
!        NOT Updated of the stats
!
         call GaMD_update_stat(2)
!
      end if
!
!##################################################################
!     STEP 4: GaMD equilibration production
!##################################################################
!
      if (amdtpass .lt. gamd_eqsteps .and. amdtpass .ge.&
         &gamd_eqprepsteps) then
!
!        Apply the boost potential with updating parameters gamdE,
!        gamdk and gamdk0
!
         gamdstep = 4
!
         call GaMD_calc_force (epot)
!
!        Update of all GaMD parameters
!
         call GaMD_update_stat(1)
         call GaMD_calc_E
!
      end if
!
!##################################################################
!     STEP 5: GaMD production
!##################################################################
!
      if (amdtpass .ge. gamd_eqsteps) then
!
!        Apply the boost potential without updating parameters gamdE,
!        gamdk and gamdk0
!
         gamdstep = 5
!
         call GaMD_calc_force (epot)
!
!        Update of GaMD results
!
         call GaMD_update_stat(2)
!
      end if
!
      if (tamd .eq. 0) then
         tamd = tamd + aMDdt
         amdtpass_dih = amdtpass_dih + 1
         amdtpass = amdtpass + 1
      end if
!
   end if                    ! END LOOP 1
!
   flush (iamdout)
   return
end
!
!     ###############################################################
!     ##                                                           ##
!     ##  subroutine GaMD_update_stat                              ##
!     ##                                                           ##
!     ###############################################################
!
!
!     "GaMD_update_stat" updates data such as Vmin, Vmax, Vavg, ..
!     which are used after to calculate E and k0
!
!
subroutine GaMD_update_stat(i)
   use mamd
   use mpi
   use domdec
   use deriv
   use energi
   use atmtyp
   use atoms
   use potent
   use usage
   use tinheader
   implicit none
   integer i,c
!
!     DIHEDRAL CASE:
!
   if (use_amd_dih) then
      if (i .eq. 1) then
         if (amdtpass .eq. gamd_cmdprepsteps) then
            VminD = eDaMD
            VmaxD = eDaMD
         end if
         if (eDaMD .lt. VminD) VminD = eDaMD
         if (eDaMD .gt. VmaxD) VmaxD = eDaMD
         c = modulo(amdtpass,amdoutputFreq_dih)
         if (c == 0) then
            VdiffD = eDaMD - VavgD
            VavgD = VavgD + (VdiffD/cptgamdD)
            M2D = M2D + VdiffD*(eDaMD - VavgD)
            sigmaVD = sqrt(M2D/cptgamdD)
            cptgamdD = cptgamdD + 1.0
         end if
      end if
      c = modulo(amdtpass,amdoutputFreq_dih)
      if (rank.eq.0) then
         if (c==0) then
            open(igamdrestart,file='GaMD_restart.dat',action="write")
            write (igamdrestart,55) tamd, VminD, VmaxD, VavgD,&
               &sigmaVD, M2D, gamdED, gamdk0D, gamdkD, amdtpass, cptgamdD
            close(igamdrestart)
55          format ("#GaMD Restart file for Tinker-HP DIHEDRAL CASE: ",/&
               &,f12.3, 10x, "#tamd", /, f12.3, 10x, "Vmin", /, f12.3, 10x,&
               &"Vmax", /, f12.3, 10x, "Vavg", /, f12.3, 10x, "sigmaV", /,&
               &e14.5, 10x, "M2", /, f20.3, 10x, "E", /, f8.6, 10x, "k0", /,&
               &f8.6, 10x, "k", /, i14, 10x, "amdtpass", /, f14.2, 10x,&
               &"cptgamdD")
         end if
         c = modulo(amdtpass,amdoutputFreq_dih)
         if (c == 0) then
            open (iamdout,file='aMD_output.dat',position='append',&
               &action="write",status='old')
            write (iamdout,66) tamd, gamdstep, amdboostD,&
               &amdboostavgD/real(amdoutputFreq_dih), VmaxD, VminD, VavgD,&
               &sigmaVD, gamdED, gamdk0D, gamdkD
66          format("GaMD-D ",f10.4,5x,i1,7f12.3,5x,f8.6,5x,f8.6)
            amdboostavgD = 0.0d0
            close(iamdout)
         end if
      end if
   end if
!
!     TOTAL POTENTIAL ENERGY CASE:
!
   if (use_amd_ene) then
      if (i .eq. 1) then
         if (amdtpass .eq. gamd_cmdprepsteps) then
            VminP = ePaMD
            VmaxP = ePaMD
         end if
         if (ePaMD .lt. VminP) VminP = ePaMD
         if (ePaMD .gt. VmaxP) VmaxP = ePaMD
         c = modulo(amdtpass,amdoutputFreq)
         if (c == 0) then
            VdiffP = ePaMD - VavgP
            VavgP = VavgP + (VdiffP/cptgamdP)
            M2P = M2P + VdiffP*(ePaMD - VavgP)
            sigmaVP = sqrt(M2P/cptgamdP)
            cptgamdP = cptgamdP + 1.0
         end if
      end if
      c = modulo(amdtpass,amdoutputFreq)
      if (rank.eq.0) then
         if (c==0) then
            if (.not. use_amd_dih) then
               open(igamdrestart,file='GaMD_restart.dat',action="write")
               write (igamdrestart,56) tamd, VminP, VmaxP, VavgP,&
                  &sigmaVP, M2P, gamdEP, gamdk0P, gamdkP, amdtpass, cptgamdP
               close(igamdrestart)
56             format ("#GaMD Restart file for Tinker-HP POTENTIAL ENERGY&
                  &      CASE: ",/&
                  &,f12.3, 10x, "#tamd", /, f12.3, 10x, "Vmin", /, f12.3, 10x,&
                  &"Vmax", /, f12.3, 10x, "Vavg", /, f12.3, 10x, "sigmaV", /,&
                  &f20.5, 10x, "M2", /, f20.3, 10x, "E", /, f8.6, 10x, "k0", /,&
                  &f8.6, 10x, "k", /, i14, 10x, "amdtpass", /, f14.2, 10x,&
                  &"cptgamdP")
            else
               open(igamdrestart,file='GaMD_restart.dat',position='append',&
                  &action="write",status='old')
               write (igamdrestart,57) tamd, VminP, VmaxP, VavgP,&
                  &sigmaVP, M2P, gamdEP, gamdk0P, gamdkP, amdtpass, cptgamdP
               close(igamdrestart)
57             format ("#GaMD Restart file for Tinker-HP POTENTIAL ENERGY&
                  &      CASE: ",/&
                  &,f12.3, 10x, "#tamd", /, f12.3, 10x, "Vmin", /, f12.3, 10x,&
                  &"Vmax", /, f12.3, 10x, "Vavg", /, f12.3, 10x, "sigmaV", /,&
                  &f20.5, 10x, "M2", /, f20.3, 10x, "E", /, f8.6, 10x, "k0", /,&
                  &f8.6, 10x, "k", /, i14, 10x, "amdtpass", /, f14.2, 10x,&
                  &"cptgamdP")
            end if
         end if
         c = modulo(amdtpass,amdoutputFreq)
         if (c == 0) then
            open (iamdout,file='aMD_output.dat',position='append',&
               &action="write",status='old')
            write (iamdout,67) tamd, gamdstep, amdboostP,&
               &amdboostavgP/real(amdoutputFreq), VmaxP, VminP, VavgP,&
               &sigmaVP, gamdEP, gamdk0P, gamdkP
67          format("GaMD-P ",f10.4,5x,i1,7f12.3,5x,f8.6,5x,f8.6)
            amdboostavgP = 0.0d0
            close(iamdout)
         end if
      end if
   end if
!
!     BOND/ANGLE WATER POTENTIAL ENERGY CASE:
!
   if (use_amd_wat1) then
      if (i .eq. 1) then
         if (amdtpass .eq. gamd_cmdprepsteps) then
            VminW1 = eW1aMD
            VmaxW1 = eW1aMD
         end if
         if (eW1aMD .lt. VminW1) VminW1 = eW1aMD
         if (eW1aMD .gt. VmaxW1) VmaxW1 = eW1aMD
         c = modulo(amdtpass,amdoutputFreq_dih)
         if (c == 0) then
            VdiffW1 = eW1aMD - VavgW1
            VavgW1 = VavgW1 + (VdiffW1/cptgamdW1)
            M2W1 = M2W1 + VdiffW1*(eW1aMD - VavgW1)
            sigmaVW1 = sqrt(M2W1/cptgamdW1)
            cptgamdW1 = cptgamdW1 + 1.0
         end if
      end if
      c = modulo(amdtpass,amdoutputFreq_dih)
      if (rank.eq.0) then
         if (c==0) then
            if (use_amd_dih) then
               open(igamdrestart,file='GaMD_restart.dat',position='append',&
                  &action="write",status='old')
               write (igamdrestart,58) tamd, VminW1, VmaxW1, VavgW1,&
                  &sigmaVW1, M2W1, gamdEW1, gamdk0W1,gamdkW1,amdtpass,cptgamdW1
               close(igamdrestart)
            else
               open(igamdrestart,file='GaMD_restart.dat',action="write")
               write (igamdrestart,58) tamd, VminW1, VmaxW1, VavgW1,&
                  &sigmaVW1, M2W1, gamdEW1, gamdk0W1,gamdkW1,amdtpass,cptgamdW1
               close(igamdrestart)
            end if
58          format ("#GaMD Restart file for Tinker-HP BOND/ANGLE WATER&
               &      CASE: ",/&
               &,f12.3, 10x, "#tamd", /, f12.3, 10x, "#Vmin", /, f12.3, 10x,&
               &"#Vmax", /, f12.3, 10x,"#Vavg", /, f12.3, 10x, "#sigmaV", /,&
               &f20.5, 10x, "#M2",/, f20.3, 10x,"#E", /,f8.6, 10x, "#k0", /,&
               &f12.10, 10x, "#k", /, i14, 10x, "#amdtpass", /, f14.2, 10x,&
               &"cptgamdW1")
         end if
         c = modulo(amdtpass,amdoutputFreq_dih)
         if (c == 0) then
            open (iamdout,file='aMD_output.dat',position='append',&
               &action="write",status='old')
            write (iamdout,68) tamd, gamdstep, amdboostW1,&
               &amdboostavgW1/real(amdoutputFreq_dih),VmaxW1,VminW1, VavgW1,&
               &sigmaVW1, gamdEW1, gamdk0W1, gamdkW1
68          format("GaMD-W ",f10.4,5x,i1,7f12.3,5x,f8.6,5x,f12.10)
            amdboostavgW1 = 0.0d0
            close(iamdout)
         end if
      end if
   end if
   return
end
!
!     ###############################################################
!     ##                                                           ##
!     ##  subroutine GaMD_calc_E                                   ##
!     ##                                                           ##
!     ###############################################################
!
!
!     "GaMD_calc_E" calculates each E and k according to the chosen
!     conditions we need then to apply as a boost
!
!
subroutine GaMD_calc_E
   use mamd
   use mpi
   use domdec
   use deriv
   use energi
   use atmtyp
   use atoms
   use potent
   use usage
   implicit none
   integer c
!
!     DIHEDRAL CASE:
!
   if (use_amd_dih) then
      if (gamd_ie .eq. 1) then
         gamdED = VmaxD
         gamdk0bisD = (gamd_sigma0D/sigmaVD)*(VmaxD-VminD)/&
            &(VmaxD-VavgD)
         if (gamdk0bisD .gt. 1.0) gamdk0D = 1.0
         if (gamdk0bisD .lt. 1.0) gamdk0D = gamdk0bisD
      end if
      if (gamd_ie .eq. 2) then
         gamdk0terD = (1.0-gamd_sigma0D/sigmaVD)*(VmaxD-VminD)/&
            &(VavgD-VminD)
         if (gamdk0terD .gt. 1 .and. gamdk0terD .le. 0) then
            gamdk0D = gamdk0terD
            gamdED = VminD + (VmaxD-VminD)/gamdk0D
         else
            gamdED = VmaxD
            gamdk0bisD = (gamd_sigma0D/sigmaVD)*(VmaxD-VminD)/&
               &(VmaxD-VavgD)
            if (gamdk0bisD .gt. 1.0) gamdk0D = 1.0
            if (gamdk0bisD .lt. 1.0) gamdk0D = gamdk0bisD
         end if
      end if
      gamdkD = gamdk0D/(VmaxD - VminD)
   end if
!
!     TOTENTIAL POTENTIAL ENERGY CASE:
!
   if (use_amd_ene) then
      if (gamd_ie .eq. 1) then
         gamdEP = VmaxP
         gamdk0bisP = (gamd_sigma0P/sigmaVP)*(VmaxP-VminP)/&
            &(VmaxP-VavgP)
         if (gamdk0bisP .gt. 1.0) gamdk0P = 1.0
         if (gamdk0bisP .lt. 1.0) gamdk0P = gamdk0bisP
      end if
      if (gamd_ie .eq. 2) then
         gamdk0terP = (1.0-gamd_sigma0P/sigmaVP)*(VmaxP-VminP)/&
            &(VavgP-VminP)
         if (gamdk0terP .gt. 1 .and. gamdk0terP .le. 0) then
            gamdk0P = gamdk0terP
            gamdEP = VminP + (VmaxP-VminP)/gamdk0P
         else
            gamdEP = VmaxP
            gamdk0bisP = (gamd_sigma0P/sigmaVP)*(VmaxP-VminP)/&
               &(VmaxP-VavgP)
            if (gamdk0bisP .gt. 1.0) gamdk0P = 1.0
            if (gamdk0bisP .lt. 1.0) gamdk0P = gamdk0bisP
         end if
      end if
      gamdkP = gamdk0P/(VmaxP - VminP)
   end if
!
!     BOND/ANGLE WATER POTENTIAL ENERGY CASE:
!
   if (use_amd_wat1) then
      if (gamd_ie .eq. 1) then
         gamdEW1 = VmaxW1
         gamdk0bisW1 = (gamd_sigma0W1/sigmaVW1)*(VmaxW1-VminW1)/&
            &(VmaxW1-VavgW1)
         if (gamdk0bisW1 .gt. 1.0) gamdk0W1 = 1.0
         if (gamdk0bisW1 .lt. 1.0) gamdk0W1 = gamdk0bisW1
      end if
      if (gamd_ie .eq. 2) then
         gamdk0terW1 = (1.0-gamd_sigma0W1/sigmaVW1)*(VmaxW1-VminW1)/&
            &(VavgW1-VminW1)
         if (gamdk0terW1 .gt. 1 .and. gamdk0terW1 .le. 0) then
            gamdk0W1 = gamdk0terW1
            gamdEW1 = VminW1 + (VmaxW1-VminW1)/gamdk0W1
         else
            gamdEW1 = VmaxW1
            gamdk0bisW1 = (gamd_sigma0W1/sigmaVW1)*(VmaxW1-VminW1)/&
               &(VmaxW1-VavgW1)
            if (gamdk0bisW1 .gt. 1.0) gamdk0W1 = 1.0
            if (gamdk0bisW1 .lt. 1.0) gamdk0W1 = gamdk0bisW1
         end if
      end if
      gamdkW1 = gamdk0W1/(VmaxW1 - VminW1)
   end if
   return
end
!
!     ###############################################################
!     ##                                                           ##
!     ##  subroutine GaMD_calc_force                               ##
!     ##                                                           ##
!     ###############################################################
!
!
!     "GaMD_calc_force" calculates the GaMD potentials and their respective
!     factors applied on the corresponded derivtaives
!
subroutine GaMD_calc_force (epot)
   use mamd
   use mpi
   use domdec
   use deriv
   use energi
   use atmtyp
   use atoms
   use potent
   use usage
   implicit none
   real(r_p) epot
!
!     Apply the boost potential without updating parameters gamdE,
!     gamdk and gamdk0
!
!     DIHEDRAL CASE:
!
   if (use_amd_dih) then
      if (eDaMD .lt. gamdED) then
         gamd_factor_dih = 1.0 - gamdkD*(gamdED - eDaMD)
         gamd_deltaV = 0.5*gamdkD*(gamdED - eDaMD)**2
         epot = epot + gamd_deltaV
         amdboostD = gamd_deltaV
         amdboostavgD = amdboostavgD + amdboostD
      else
         gamd_factor_dih = 1.0
         amdboostD = 0.0d0
      end if
      eDaMD = eDaMD + amdboostD
   end if
!
!     TOTAL POTENTIAL ENERGY CASE:
!
   if (use_amd_ene) then
      if (ePaMD .lt. gamdEP) then
         gamd_factor_tot = 1.0 - gamdkP*(gamdEP - ePaMD)
         gamd_deltaV = 0.5*gamdkP*(gamdEP - ePaMD)**2
         epot = epot + gamd_deltaV
         amdboostP = gamd_deltaV
         amdboostavgP = amdboostavgP + amdboostP
         if (gamd_factor_tot .lt. 0.001) write (iamdout,80)
80       format("aMDinfo: WARNING ! GaMD factor is lower than 0.001",&
            &" Simulation may become unstable ...")
      else
         gamd_factor_tot = 1.0
         amdboostP = 0.0d0
      end if
      ePaMD = ePaMD + amdboostP
   end if
!
!     BOND/ANGLE WATER POTENTIAL ENERGY CASE:
!
   if (use_amd_wat1) then
      if (eW1aMD .lt. gamdEW1) then
         gamd_factor_wat1 = 1.0 - gamdkW1*(gamdEW1 - eW1aMD)
         gamd_deltaV = 0.5*gamdkW1*(gamdEW1 - eW1aMD)**2
         epot = epot + gamd_deltaV
         amdboostW1 = gamd_deltaV
         amdboostavgW1 = amdboostavgW1 + amdboostW1
      else
         gamd_factor_wat1 = 1.0
         amdboostW1 = 0.0d0
      end if
      eW1aMD = eW1aMD + amdboostW1
   end if
   return
end
!
!     ###############################################################
!     ##                                                           ##
!     ##  subroutine aMD_reduction                                 ##
!     ##                                                           ##
!     ###############################################################
!
!
!     "aMD_reduction" desactivates and reactivates uses of aMD components
!     in function of the multitimstep algorithm used
!
subroutine aMD_reduction (i)
   use potent
   use mamd
   implicit none
   integer i, j
   real(r_p) epot
!
!     Case 1: desactivates use_amd_ene (gradfast)
!
   if (i == 1) then
      if (use_amd_ene) then
         save_amd_ene = use_amd_ene
         use_amd_ene = .false.
      end if
   end if
!
!     Case 2: reactivates use_amd_ene (gradfast)
!
   if (i==2) then
      if (save_amd_ene) use_amd_ene = save_amd_ene
   end if
!
!     Case 3: desactivates use_amd_dih (gradslow)
!
   if (i == 1) then
      if (use_amd_dih) then
         save_amd_dih = use_amd_dih
         use_amd_dih = .false.
      end if
   end if
!
!     Case 4: reactivates use_amd_dih (gradslow)
!
   if (i==4) then
      if (save_amd_dih) use_amd_dih = save_amd_dih
   end if
!
   return
end
