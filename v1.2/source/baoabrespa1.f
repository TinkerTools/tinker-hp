c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###################################################################################
c     ##                                                                               ##
c     ##  subroutine baoabrespa1  --  baoab r-RESPA1 Langevin molecular dynamics step  ##
c     ##                                                                               ##
c     ###################################################################################
c
c
c     "baoabrespa1" performs a single multiple time step molecular dynamics
c     step using the reversible reference system propagation algorithm
c     (r-RESPA) via a BAOAB core with the potential split into fast-
c     intermediate and slow-evolving portions
c
c     literature references:
c
c     Pushing the Limits of Multiple-Time-Step Strategies 
c     for Polarizable Point Dipole Molecular Dynamics
c     L Lagardere, F Aviat, JP Piquemal
c     The journal of physical chemistry letters 10, 2593-2599
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
c     Ruhong Zhou, Edward Harder, Huafeng Xu and B. J. Berne, 
c     "Efficient multiple time step method for use with Ewald and 
c     particle mesh Ewald for large biomolecular systems", 
c     J. Chem. Phys. 115, 2348-2358 (2001)
c
c
      subroutine baoabrespa1(istep,dt)
      use atmtyp
      use atoms
      use cutoff
      use domdec
      use deriv
      use energi
      use freeze
      use inform
      use moldyn
      use mpi
      use timestat
      use units
      use usage
      use virial
      use utilbaoab
      use bath, only: use_piston
      use mdstuf1
      implicit none
      integer i,j,iglob
      integer istep
      integer,save::n_call=0
      real*8 dt,dt_2
      real*8 dta,dta_2,dta2
      real*8    time0,time1

      if(istep==1) then
        if(use_piston .and. rank==0) then
          write(0,*) "BAOABRESPA1 is not compatible with "
     &                 //"LANGEVIN PISTON barostat"
          call fatal
        endif
      endif

      time0 = mpi_wtime()
c
c     set some time values for the dynamics integration
c
      dt_2  = 0.5d0 * dt
      dta   = dt / dinter
      dta_2 = 0.5d0 * dta
c      
      dta2 = dta / dshort
c
c     store the current atom positions, then find half-step
c     velocities via BAOAB recursion
c
      do i = 1, nloc
        iglob = glob(i)
        if (use(iglob)) then
          v(:,iglob) = v(:,iglob) + a(:,iglob)*dt_2
        end if
      end do
c
      if (use_rattle) call rattle2(dt_2)
      time1 = mpi_wtime()
      timeinte = timeinte + time1-time0 
c
c     find intermediate-evolving velocities and positions via BAOAB recursion
c
      call baoabrespaint1(dta,dta2)
c
c     Reassign the particules that have changed of domain
c
c     -> reciprocal space
c
      time0 = mpi_wtime()
      call reassignpme(.false.)
c
c     communicate positions
c
      time0 = mpi_wtime()
      call commposrec
      time1 = mpi_wtime()
      timecommpos = timecommpos + time1 - time0
c
      time0 = mpi_wtime()
      call reinitnl(istep)
      time1 = mpi_wtime()
      timeinte = timeinte + time1-time0
c
      call mechanicsteprespa1(istep,2)

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
      if (allocated(derivs)) deallocate(derivs)
      allocate(derivs(3,nbloc))
      derivs = 0d0
      time1 = mpi_wtime()
      timeinte = timeinte + time1-time0
c
c     get the slow-evolving potential energy and atomic forces
c
      time0 = mpi_wtime()
      call gradslow1 (epot,derivs)
      time1 = mpi_wtime()
      timegrad = timegrad + time1 - time0
c
c     communicate some forces
c
      time0 = mpi_wtime()
      call commforcesrespa1(derivs,2)
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
c     make half-step temperature and pressure corrections
c
      time0 = mpi_wtime()
      call temper2 (temp)
      time1 = mpi_wtime()
      timetp = timetp + time1-time0
c
c     use Newton's second law to get the slow accelerations;
c     find full-step velocities using BAOAB recursion
c
      time0 = mpi_wtime()
      do i = 1, nloc
        iglob = glob(i)
        if (use(iglob)) then
          a(:,iglob) = -convert * derivs(:,i) / mass(iglob)
          v(:,iglob) = v(:,iglob) + a(:,iglob)*dt_2
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
      do i = 1, 3
         do j = 1, 3
            vir(j,i) = vir(j,i) + viralt(j,i)
         end do
      end do
      time1 = mpi_wtime()
      timeinte = timeinte + time1-time0
c
c     make full-step temperature and pressure corrections
c
      time0 = mpi_wtime()
      call temper (dt,eksum,ekin,temp)
      call pressure (dt,ekin,pres,stress,istep)
      call pressure2 (epot,temp)
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
      
      end
c
c     subroutine baoabrespaint1 : 
c     find intermediate-evolving velocities and positions via BAOAB recursion
c
      subroutine baoabrespaint1(dta,dta2)
      use atmtyp
      use cutoff
      use deriv
      use domdec
      use deriv
      use energi
      use freeze
      use inform
      use moldyn
      use timestat
      use units
      use usage
      use virial
      use mpi
      use mdstuf1
      implicit none
      integer i,j,k,iglob,stepint
      real*8 dta,dta_2,dta2
      real*8 time0,time1
      time0 = mpi_wtime()
      dta_2 = 0.5d0 * dta
c
c     initialize virial from fast-evolving potential energy terms
c
      viralt(:,:) = 0.d0
      time1 = mpi_wtime()
      timeinte = timeinte + time1-time0 

      do stepint = 1, nalt
         time0=mpi_wtime()
         do i = 1, nloc
            iglob = glob(i)
            if (use(iglob)) then
              v(:,iglob) = v(:,iglob) + aalt(:,iglob)*dta_2
            end if
         end do
c
         if (use_rattle)  call rattle2 (dta_2)
         time1 = mpi_wtime()
         timeinte = timeinte + time1-time0 
c
c     find fast-evolving velocities and positions via BAOAB recursion
c
         call baoabrespafast1(dta2)
c
         time0 = mpi_wtime()
         call mechanicsteprespa1(-1,1)
         time1 = mpi_wtime()
         timeparam = timeparam + time1 - time0

         time0 = mpi_wtime()
         call allocsteprespa(.false.)
c
         if(allocated(derivs)) deallocate(derivs)
         allocate(derivs(3,nbloc))
         derivs(:,:)=0.d0

         if (allocated(desave)) deallocate(desave)
         allocate(desave(3,nbloc))
         time1 = mpi_wtime()
         timeinte = timeinte + time1 - time0
c
c     get the fast-evolving potential energy and atomic forces
c
        time0 = mpi_wtime()
        call gradint1(ealt,derivs)
        time1 = mpi_wtime()
        timegrad = timegrad + time1 - time0
c
c     communicate forces
c
        time0 = mpi_wtime()
        call commforcesrespa1(derivs,1)
        time1 = mpi_wtime()
        timecommforces = timecommforces + time1-time0
c
c     MPI : get total energy
c
        time0 = mpi_wtime()
        call reduceen(ealt)
        time1 = mpi_wtime()
        timered = timered + time1 - time0
c
c     use Newton's second law to get fast-evolving accelerations;
c     update fast-evolving velocities using the BAOAB recursion
c
        time0 = mpi_wtime()
        do i = 1, nloc
           iglob = glob(i)
           if (use(iglob)) then
              do j = 1, 3
                 aalt(j,iglob) = -convert *
     $              derivs(j,i) / mass(iglob)
                 v(j,iglob) = v(j,iglob) + aalt(j,iglob)*dta_2
              end do
           end if
        end do
        deallocate (derivs)
        if (use_rattle)  call rattle2 (dta)
c
c     increment average virial from fast-evolving potential terms
c
        ealt = ealt + ealt2
        do i = 1, 3
           do j = 1, 3
              viralt(j,i) = viralt(j,i) + (viralt2(j,i) 
     $         + vir(j,i))/dinter
           end do
        end do
        time1 = mpi_wtime()
        timeinte = timeinte + time1-time0
      end do
      end
c
c     subroutine baoabrespafast1 : 
c     find fast-evolving velocities and positions via BAOAB recursion
c
      subroutine baoabrespafast1(dta)
      use atmtyp
      use atoms
      use bath
      use cutoff
      use domdec
      use deriv
      use energi
      use freeze
      use inform
      use moldyn
      use timestat
      use units
      use usage
      use virial
      use mpi
      use utilbaoab
      use mdstuf1
      implicit none
      integer i,j,k,iglob,stepfast
      real*8 dta,dta_2
      real*8 time0,time1
      real*8 a1,a2

      time0 = mpi_wtime()
c
c     set time values and coefficients for BAOAB integration
c
      dta_2 = 0.5d0 * dta
c
c     initialize virial from fast-evolving potential energy terms
      viralt2(:,:)=0.d0
      time1 = mpi_wtime()
      timeinte = timeinte + time1-time0 

      do stepfast = 1, nalt2
         time0=mpi_wtime()
        do i = 1, nloc
           iglob = glob(i)
           if (use(iglob)) then
             v(:,iglob) = v(:,iglob) + aalt2(:,iglob)*dta_2
           end if
        end do
c
        if (use_rattle) then
           call rattle2 (dta_2)
        end if

        call apply_a_block(dta_2)
        call apply_o_block(dta)
        call apply_a_block(dta_2)
        time1 = mpi_wtime()
        timeinte = timeinte + time1-time0 
c
c       Reassign the particules that have changed of domain
c       
c       -> real space
         time0 = mpi_wtime()
        call reassignrespa(stepfast,nalt2)
        time1 = mpi_wtime()
        timereneig = timereneig + time1 - time0
c
c       communicate positions
c
         time0 = mpi_wtime()
        call commposrespa(stepfast.ne.nalt2)
        time1 = mpi_wtime()
        timecommpos = timecommpos + time1 - time0
c
         time0 = mpi_wtime()
         if (allocated(derivs)) deallocate(derivs)
         allocate(derivs(3,nbloc))
         derivs(:,:)=0.d0
         time1 = mpi_wtime()
         timeinte = timeinte + time1-time0
c
         if (stepfast.eq.nalt2) then
            time0 = mpi_wtime()
            call mechanicsteprespa1(-1,0)
            timeparam = timeparam + time1 - time0
            time0 = mpi_wtime()
            call allocsteprespa(.true.)
            time1 = mpi_wtime()
            timeinte = timeinte + time1 - time0
         endif
c
c     get the fast-evolving potential energy and atomic forces
c
        time0 = mpi_wtime()
        call gradfast1(ealt2,derivs)
        time1 = mpi_wtime()
        timegrad = timegrad + time1 - time0
c
c       communicate forces
c
        time0 = mpi_wtime()
        call commforcesrespa1(derivs,0)
        time1 = mpi_wtime()
        timecommforces = timecommforces + time1-time0
c
c       MPI : get total energy
c
        time0 = mpi_wtime()
        call reduceen(ealt2)
        time1 = mpi_wtime()
        timered = timered + time1 - time0
c
c     use Newton's second law to get fast-evolving accelerations;
c     update fast-evolving velocities using the BAOAB recursion
c
        time0 = mpi_wtime()
        do i = 1, nloc
          iglob = glob(i)
          if (use(iglob)) then
             aalt2(:,iglob) = -convert *
     $          derivs(:,i) / mass(iglob)
             v(:,iglob) = v(:,iglob) + aalt2(:,iglob)*dta_2
          end if
        end do
        deallocate(derivs)
c
        if (use_rattle)  call rattle2 (dta_2)
c
c     increment average virial from fast-evolving potential terms
c
        do i = 1, 3
           do j = 1, 3
              viralt2(j,i) = viralt2(j,i) + vir(j,i)/dshort
           end do
        end do
        time1 = mpi_wtime()
        timeinte = timeinte + time1 - time0
      end do
      end
