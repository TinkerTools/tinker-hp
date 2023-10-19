c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ######################################################################
c     ##                                                                  ##
c     ##  subroutine baoab  --  baoab Langevin molecular dynamics step  ##
c     ##                                                                  ##
c     ######################################################################
c
c
c     "baoab" performs a single molecular dynamics time step
c     via the respa-baoab recursion formula
c
c     literature references:
c
c     Efficient molecular dynamics using geodesic integration 
c     and solvent-solute splitting B. Leimkuhler and C. Matthews,
c     Proceedings of the Royal Society A, 472: 20160138, 2016
c
c     D. D. Humphreys, R. A. Friesner and B. J. Berne, "A Multiple-
c     Time-Step Molecular Dynamics Algorithm for Macromolecules",
c     Journal of Physical Chemistry, 98, 6885-6892 (1994)
c
c     X. Qian and T. Schlick, "Efficient Multiple-Time-Step Integrators
c     with Distance-Based Force Splitting for Particle-Mesh-Ewald
c     Molecular Dynamics Simulations", Journal of Chemical Physics,
c     115, 4019-4029 (2001)
c
c

      subroutine baoabrespa (istep,dt)
      use atmtyp
      use atoms
      use bath
      use cutoff
      use domdec
      use deriv
      use energi
      use freeze
      use inform
      use mdstuf
      use moldyn
      use mpi
      use timestat
      use tortor
      use units
      use usage
      use virial
      use spectra
      use utilbaoab
      use mdstuf1
      implicit none
      real*8, intent(in) :: dt
      integer, intent(in) :: istep
      integer i,iglob,stepfast
      real*8 dta,dta_2,dt_2
      real*8 time0,time1
      real*8 dip(3,nalt),dipind(3,nalt)

      if (istep.eq.1) then
        pres=atmsph
      end if
      time0 = mpi_wtime()
c
c     set some time values for the dynamics integration
c
      dt_2  = 0.5d0 * dt
      dta   = dt / dshort
      dta_2 = 0.5d0 * dta
c
c     find quarter step velocities and half step positions via baoab recursion
c
      do i = 1, nloc
        iglob = glob(i)
        if (use(iglob)) then
           v(:,iglob) = v(:,iglob) + dt_2*a(:,iglob)
        end if
      end do
c
      if (use_rattle) call rattle2(dt_2)
c     initialize virial from fast-evolving potential energy terms
      viralt(:,:)=0.0d0

      if(use_piston) then
        call apply_a_piston(dt_2,-1,.FALSE.)
        call apply_o_piston(dt)
      endif
      time1 = mpi_wtime()
      timeinte = timeinte + time1-time0 
c
c     respa inner loop
c
      do stepfast = 1, nalt
c
c     find fast-evolving velocities and positions via BAOAB recursion
c
        do i = 1, nloc
          iglob = glob(i)
          if (use(iglob)) then
             v(:,iglob) = v(:,iglob) + dta_2*aalt(:,iglob)
          end if
        end do
c
        if (use_rattle)  call rattle2 (dta_2)
c
        call apply_a_block(dta_2)
        call apply_o_block(dta)
        call apply_a_block(dta_2)
c
c       Reassign the particules that have changed of domain
c
c       -> real space
        call reassignrespa(stepfast,nalt)
c
c       communicate positions
c
        call commposrespa(stepfast.ne.nalt)
        !call commposrespa(.TRUE.)
c
        if(allocated(derivs)) deallocate(derivs)
        allocate(derivs(3,nbloc))
        derivs(:,:)=0.d0
c
        call mechanic_up_para_respa(istep,.true.)
        call allocsteprespa(.true.)
c
c     get the fast-evolving potential energy and atomic forces
c
        call gradfast (ealt,derivs)
c
c       communicate forces
c
        call commforcesrespa(derivs,.true.)
c
c       MPI : get total energy
c
        call reduceen(ealt)
c
c     use Newton's second law to get fast-evolving accelerations;
c     update fast-evolving velocities using the Verlet recursion
c
        do i = 1, nloc
          iglob = glob(i)
          if (use(iglob)) then
            aalt(:,iglob) = -convert *
     $         derivs(:,i) / mass(iglob)
            v(:,iglob) = v(:,iglob) + aalt(:,iglob)*dta_2
          end if
        end do
        deallocate(derivs)
c
        if (use_rattle)  call rattle2 (dta_2)
c
c     increment average virial from fast-evolving potential terms
        viralt(:,:)=viralt(:,:)+vir(:,:)/dshort

        if(ir .and. stepfast<nalt) then
          !call rotpole
          call compute_dipole(dip(:,stepfast)
     &       ,dipind(:,stepfast),.FALSE.)
        endif
      end do

      if(use_piston) call apply_a_piston(dt_2,-1,.FALSE.)
c
c     Reassign the particules that have changed of domain
c
c     -> real space
c
      time0 = mpi_wtime()
      call reassignrespa(nalt,nalt)
c
c     -> reciprocal space
      call reassignpme(.false.)
c
c     communicate positions
c
      call commposrespa(.false.)
      call commposrec
      time1 = mpi_wtime()
      timecommpos = timecommpos + time1 - time0
c
c
      time0 = mpi_wtime()
      call reinitnl(istep)
      time1 = mpi_wtime()
      timeinte = timeinte + time1-time0
      time0 = mpi_wtime()
c
      call mechanic_up_para_respa(istep,.false.)

      call allocsteprespa(.false.)
      time1 = mpi_wtime()
      timeinte = timeinte + time1 - time0
c
c     rebuild the neighbor lists
c
      time0 = mpi_wtime()
      if (use_list) call nblist(istep)
      time1 = mpi_wtime()
      timenl = timenl + time1 - time0
c
      time0 = mpi_wtime()
      if(allocated(derivs)) deallocate(derivs)
      allocate(derivs(3,nbloc))
      derivs = 0d0
      time1 = mpi_wtime()
      timeinte = timeinte + time1 - time0
c
c     get the slow-evolving potential energy and atomic forces
c
      time0 = mpi_wtime()
      call gradslow (epot,derivs)
      time1 = mpi_wtime()
      timegrad = timegrad + time1 - time0
c
c     communicate forces
c
      time0 = mpi_wtime()
      call commforcesrespa(derivs,.false.)
      time1 = mpi_wtime()
      timecommforces = timecommforces + time1-time0
c
c     MPI : get total energy
c
      time0 = mpi_wtime()
      call reduceen(epot)
      time1 = mpi_wtime()
      timered = timered + time1 - time0
c
c     use Newton's second law to get the next accelerations;
c     find the full-step velocities using the BAOAB recursion
c
      time0 = mpi_wtime()
      do i = 1, nloc
        iglob = glob(i)
        if (use(iglob)) then
          a(:,iglob) = -convert * derivs(:,i)/mass(iglob)
          v(:,iglob) = v(:,iglob) + dt_2*a(:,iglob)
        end if
      end do
c
c     find the constraint-corrected full-step velocities
c
      if (use_rattle)  call rattle2 (dt)
c
c     total potential and virial from sum of fast and slow parts
c
      epot = epot + ealt
      vir = vir + viralt

      time1 = mpi_wtime()
      timeinte = timeinte + time1 - time0
c
      time0 = mpi_wtime()
      call temper (dt,eksum,ekin,temp)
      call pressure (dt,ekin,pres,stress,istep)
      call pressure2 (epot,temp)
      if(ir) then
        call compute_dipole(dip(:,nalt)
     &       ,dipind(:,nalt),full_dipole)
        call save_dipole_respa(dip,dipind)
      endif
      if(use_piston) then
        call apply_b_piston(dt,pres,stress)
      endif
      time1 = mpi_wtime()
      timetp = timetp + time1-time0
c
c     total energy is sum of kinetic and potential energies
c
      time0 = mpi_wtime()
      etot = eksum + epot
c
c     compute statistics and save trajectory for this step
c
      call mdsave (istep,dt,epot,derivs)
      call mdrest (istep)
      call mdstat (istep,dt,etot,epot,eksum,temp,pres)
      time1 = mpi_wtime()
      timeinte = timeinte + time1-time0
      return
      end
