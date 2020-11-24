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
c     wfor Polarizable Point Dipole Molecular Dynamics
c     L Lagardère, F Aviat, JP Piquemal
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
      use freeze
      use moldyn
      use timestat
      use units
      use usage
      use virial
      use mpi
      implicit none
      integer i,j,iglob
      integer istep
      real*8 dt,dt_2
      real*8 dta,dta_2,dta2
      real*8 epot,etot
      real*8 eksum
      real*8 temp,pres
      real*8 ealt
      real*8 ekin(3,3)
      real*8 stress(3,3)
      real*8 viralt(3,3)
      real*8 time0,time1
      real*8, allocatable :: derivs(:,:)
c
c
c     set some time values for the dynamics integration
c
      dt_2 = 0.5d0 * dt
      dta = dt / dinter
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
            do j = 1, 3
               v(j,iglob) = v(j,iglob) + a(j,iglob)*dt_2
            end do
         end if
      end do
c
      if (use_rattle) call rattle2(dt_2)
c
c
c     find intermediate-evolving velocities and positions via BAOAB recursion
c
      call baoabrespaint1(ealt,viralt,dta,dta2)
c
c     Reassign the particules that have changed of domain
c
c     -> real space
c
      time0 = mpi_wtime()
c
c      call reassignrespa(.false.,nalt,nalt)
c
c     -> reciprocal space
c
      call reassignpme(.false.)
      time1 = mpi_wtime()
      timereneig = timereneig + time1 - time0
c
c     communicate positions
c
      time0 = mpi_wtime()
c      call commposrespa(.false.)
      call commposrec
      time1 = mpi_wtime()
      timecommstep = timecommstep + time1 - time0
c
c
      call reinitnl(istep)
c
      time0 = mpi_wtime()
      call mechanicsteprespa1(istep,2)
      time1 = mpi_wtime()
      timeparam = timeparam + time1 - time0

      time0 = mpi_wtime()
      call allocsteprespa(.false.)
      time1 = mpi_wtime()
      timeclear = timeclear + time1 - time0
c
c     rebuild the neighbor lists
c
      if (use_list) call nblist(istep)
c
      allocate (derivs(3,nbloc))
      derivs = 0d0
c
c     get the slow-evolving potential energy and atomic forces
c
      call gradslowbaoab1 (epot,derivs)
c
c     communicate some forces
c
      call commforcesrespa1(derivs,2)
c
c     MPI : get total energy
c
      call reduceen(epot)
c
c     make half-step temperature and pressure corrections
c
      call temper2 (temp)
c      call pressure2 (epot,temp)
c
c     use Newton's second law to get the slow accelerations;
c     find full-step velocities using BAOAB recursion
c
      do i = 1, nloc
         iglob = glob(i)
         if (use(iglob)) then
            do j = 1, 3
               a(j,iglob) = -convert * derivs(j,i) / mass(iglob)
               v(j,iglob) = v(j,iglob) + a(j,iglob)*dt_2
            end do
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (derivs)
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
c
c     make full-step temperature and pressure corrections
c
      call temper (dt,eksum,ekin,temp)
      call pressure (dt,ekin,pres,stress,istep)
      call pressure2 (epot,temp)
c
c     total energy is sum of kinetic and potential energies
c
      etot = eksum + epot
c
c     compute statistics and save trajectory for this step
c
      call mdstat (istep,dt,etot,epot,eksum,temp,pres)
      call mdsave (istep,dt,epot)
      call mdrest (istep)
      return
      end
c
c
c
c     #########################################################################
c     ##                                                                     ##
c     ##  subroutine gradfastbaoab1  --   fast energy & gradient components  ##
c     ##                                                                     ##
c     #########################################################################
c
c
c     "gradfast" calculates the potential energy and first derivatives
c     for the fast-evolving local valence potential energy terms
c
c
      subroutine gradfastbaoab1(energy,derivs)
      use cutoff
      use potent
      implicit none
      real*8 energy
      real*8 derivs(3,*)
      logical save_vdw,save_charge
      logical save_mpole,save_polar
      logical save_list
c
c
c     save the original state of slow-evolving potentials
c
      save_vdw = use_vdw
      save_charge = use_charge
      save_mpole = use_mpole
      save_polar = use_polar
      save_list = use_list
c
c     turn off slow-evolving nonbonded potential energy terms
c
      use_vdw = .false.
      use_charge = .false.
      use_mpole = .false.
      use_polar = .false.
      use_list = .false.
c
c     get energy and gradient for fast-evolving potential terms
c
      call gradient (energy,derivs)
c
c     restore the original state of slow-evolving potentials
c
      use_vdw = save_vdw
      use_charge = save_charge
      use_mpole = save_mpole
      use_polar = save_polar
      use_list = save_list
      return
      end
c
c     ###############################################################################
c     ##                                                                           ##
c     ##  subroutine gradintbaoab1  --  intermediate energy & gradient components  ##
c     ##                                                                           ##
c     ###############################################################################
c
c
c     "gradintbaoab1" calculates the potential energy and first derivatives
c     for the fast-evolving local valence potential energy terms
c
c
      subroutine gradintbaoab1 (energy,derivs)
      use cutoff
      use deriv
      use energi
      use polpot
      use potent
      use virial
      implicit none
      real*8 energy
      real*8 derivs(3,*)
      real*8 save_tcgomega
      logical save_bond,save_angle
      logical save_strbnd,save_urey
      logical save_angang,save_opbend
      logical save_opdist,save_improp
      logical save_imptor,save_tors
      logical save_pitors,save_strtor
      logical save_tortor,save_geom
      logical save_extra
      logical save_mrec
      logical save_prec,save_crec
      integer save_polalg,save_tcgorder
      logical save_tcgprec,save_tcgguess,save_tcgpeek
c
c
c     save the original state of fast-evolving potentials
c
      save_bond = use_bond
      save_angle = use_angle
      save_strbnd = use_strbnd
      save_urey = use_urey
      save_angang = use_angang
      save_opbend = use_opbend
      save_opdist = use_opdist
      save_improp = use_improp
      save_imptor = use_imptor
      save_tors = use_tors
      save_pitors = use_pitors
      save_strtor = use_strtor
      save_tortor = use_tortor
      save_geom = use_geom
      save_extra = use_extra
      save_crec = use_crec
      save_mrec = use_mrec
      save_prec = use_prec
      save_polalg = polalg
      save_tcgorder = tcgorder
      save_tcgprec = tcgprec
      save_tcgguess = tcgguess
      save_tcgpeek = tcgpeek
      save_tcgomega = tcgomega
c
c     turn off fast-evolving valence potential energy terms
c
      use_bond = .false.
      use_angle = .false.
      use_strbnd = .false.
      use_urey = .false.
      use_angang = .false.
      use_opbend = .false.
      use_opdist = .false.
      use_improp = .false.
      use_imptor = .false.
      use_tors = .false.
      use_pitors = .false.
      use_strtor = .false.
      use_tortor = .false.
      use_geom = .false.
      use_extra = .false.
      use_crec = .false.
      use_mrec = .false.
      use_prec = .false.
      use_cself = .false.
      use_mself = .false.
      use_pself = .false.
      use_cshortreal = .true.
      use_mpoleshortreal = .true.
      use_vdwshort = .true.
      use_polarshortreal = .true.
      polalg = polalgshort
      tcgorder = tcgordershort
      tcgprec = tcgprecshort
      tcgguess = tcgguessshort
      tcgpeek = tcgpeekshort
      tcgomega = tcgomegashort
c
c     get energy and gradient for slow-evolving potential terms
c
      call gradient (energy,derivs)
c
c     restore the original state of fast-evolving potentials
c
      use_bond = save_bond
      use_angle = save_angle
      use_strbnd = save_strbnd
      use_urey = save_urey
      use_angang = save_angang
      use_opbend = save_opbend
      use_opdist = save_opdist
      use_improp = save_improp
      use_imptor = save_imptor
      use_tors = save_tors
      use_pitors = save_pitors
      use_strtor = save_strtor
      use_tortor = save_tortor
      use_geom = save_geom
      use_extra = save_extra
      use_crec = save_crec
      use_mrec = save_mrec
      use_prec = save_prec
      use_cself = .true.
      use_mself = .true.
      use_pself = .true.
      use_cshortreal = .false.
      use_mpoleshortreal = .false.
      use_vdwshort = .false.
      use_polarshortreal = .false.
      polalg = save_polalg
      tcgorder = save_tcgorder
      tcgprec = save_tcgprec
      tcgguess = save_tcgguess
      tcgpeek = save_tcgpeek
c
c     also save the values of the short range real space forces
c
      desave = dep
      esave = ep
      return
      end
c
c
c     ########################################################################
c     ##                                                                    ##
c     ##  subroutine gradslowbaoab1  --  slow energy & gradient components  ##
c     ##                                                                    ##
c     ########################################################################
c
c
c     "gradslow" calculates the potential energy and first derivatives
c     for the slow-evolving nonbonded potential energy terms
c
c
      subroutine gradslowbaoab1 (energy,derivs)
      use cutoff
      use deriv
      use domdec
      use energi
      use moldyn
      use potent
      use virial
      implicit none
      real*8 energy
      real*8 derivs(3,*)
      logical save_bond,save_angle
      logical save_strbnd,save_urey
      logical save_angang,save_opbend
      logical save_opdist,save_improp
      logical save_imptor,save_tors
      logical save_pitors,save_strtor
      logical save_tortor,save_geom
      logical save_extra
c
c     save the original state of fast-evolving potentials
c
      save_bond = use_bond
      save_angle = use_angle
      save_strbnd = use_strbnd
      save_urey = use_urey
      save_angang = use_angang
      save_opbend = use_opbend
      save_opdist = use_opdist
      save_improp = use_improp
      save_imptor = use_imptor
      save_tors = use_tors
      save_pitors = use_pitors
      save_strtor = use_strtor
      save_tortor = use_tortor
      save_geom = use_geom
      save_extra = use_extra
c
c     turn off fast-evolving valence potential energy terms
c
      use_bond = .false.
      use_angle = .false.
      use_strbnd = .false.
      use_urey = .false.
      use_angang = .false.
      use_opbend = .false.
      use_opdist = .false.
      use_improp = .false.
      use_imptor = .false.
      use_tors = .false.
      use_pitors = .false.
      use_strtor = .false.
      use_tortor = .false.
      use_geom = .false.
      use_extra = .false.
      use_vdwlong = .true.
      use_clong = .true.
      use_mpolelong = .true.

c
c     get energy and gradient for slow-evolving potential terms
c
      call gradient (energy,derivs)
c
c     restore the original state of fast-evolving potentials
c
      use_bond = save_bond
      use_angle = save_angle
      use_strbnd = save_strbnd
      use_urey = save_urey
      use_angang = save_angang
      use_opbend = save_opbend
      use_opdist = save_opdist
      use_improp = save_improp
      use_imptor = save_imptor
      use_tors = save_tors
      use_pitors = save_pitors
      use_strtor = save_strtor
      use_tortor = save_tortor
      use_geom = save_geom
      use_extra = save_extra
      use_vdwlong = .false.
      use_clong = .false.
      use_mpolelong = .false.
c
c     substract the previously stored short range energy and forces
c
      esum = esum - esave
      energy = energy - esave
      desum = desum - desave
      derivs(:,1:nloc) = derivs(:,1:nloc) - desave(:,1:nloc)
      vir = vir - virsave/dinter
      virsave = 0d0
      return
      end
c
c     subroutine baoabrespaint1 : 
c     find intermediate-evolving velocities and positions via BAOAB recursion
c
      subroutine baoabrespaint1(ealt,viralt,dta,dta2)
      use atmtyp
      use atoms
      use cutoff
      use deriv
      use domdec
      use freeze
      use moldyn
      use timestat
      use units
      use usage
      use virial
      use mpi
      implicit none
      integer i,j,k,iglob
      real*8 dta,dta_2,dta2
      real*8 ealt,ealt2
      real*8 time0,time1
      real*8, allocatable :: derivs(:,:)
      real*8 viralt(3,3),viralt2(3,3)

      dta_2 = 0.5d0 * dta
c
c     initialize virial from intermediate potential energy terms
c
      do i = 1, 3
         do j = 1, 3
            viralt(j,i) = 0.0d0
         end do
      end do

      do k = 1, nalt
         do i = 1, nloc
            iglob = glob(i)
            if (use(iglob)) then
               do j = 1, 3
                  v(j,iglob) = v(j,iglob) + aalt(j,iglob)*dta_2
               end do
            end if
         end do
c
       if (use_rattle)  call rattle2 (dta_2)
c
c
c     find fast-evolving velocities and positions via BAOAB recursion
c
       call baoabrespafast1(ealt2,viralt2,dta2)
c
c
c      Reassign the particules that have changed of domain
c
c      -> real space
c
       time0 = mpi_wtime()
c
c       call reassignrespa(.false.,nalt,nalt)
c
       time1 = mpi_wtime()
       timereneig = timereneig + time1 - time0
cc
cc      communicate positions
cc
c       time0 = mpi_wtime()
cc       call commposshort(.false.)
c       time1 = mpi_wtime()
c       timecommstep = timecommstep + time1 - time0
c
       time0 = mpi_wtime()
       call mechanicsteprespa1(-1,1)
c       call mechanicsteprespa1(istep,1)
       time1 = mpi_wtime()
       timeparam = timeparam + time1 - time0
       time0 = mpi_wtime()
       call allocsteprespa(.false.)
       time1 = mpi_wtime()
       timeclear = timeclear + time1 - time0
c
       allocate (derivs(3,nbloc))
       derivs = 0d0

       if (allocated(desave)) deallocate (desave)
       allocate (desave(3,nbloc))
c
c     get the fast-evolving potential energy and atomic forces
c
       call gradintbaoab1(ealt,derivs)
c
c      communicate forces
c
       call commforcesrespa1(derivs,1)
c
c      MPI : get total energy
c
       call reduceen(ealt)
c
c     use Newton's second law to get fast-evolving accelerations;
c     update fast-evolving velocities using the BAOAB recursion
c
          do i = 1, nloc
             iglob = glob(i)
             if (use(iglob)) then
                do j = 1, 3
                   aalt(j,iglob) = -convert *
     $                derivs(j,i) / mass(iglob)
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
     $           + vir(j,i))/dinter
           end do
        end do
      end do
      return
      end
c
c     subroutine baoabrespafast1 : 
c     find fast-evolving velocities and positions via BAOAB recursion
c
      subroutine baoabrespafast1(ealt,viralt,dta)
      use atmtyp
      use atoms
      use bath
      use cutoff
      use domdec
      use freeze
      use langevin
      use moldyn
      use timestat
      use units
      use usage
      use virial
      use mpi
      implicit none
      integer i,j,k,iglob
      real*8 dta,dta_2
      real*8 ealt
      real*8 time0,time1
      real*8, allocatable :: derivs(:,:)
      real*8 viralt(3,3)
      real*8 a1,a2,normal
c
c     set time values and coefficients for BAOAB integration
c
      a1 = exp(-gamma*dta)
      a2 = sqrt((1-a1**2)*boltzmann*kelvin)

      dta_2 = 0.5d0 * dta
c
c     initialize virial from fast-evolving potential energy terms
c
      do i = 1, 3
         do j = 1, 3
            viralt(j,i) = 0.0d0
         end do
      end do

      do k = 1, nalt2
        do i = 1, nloc
           iglob = glob(i)
           if (use(iglob)) then
              do j = 1, 3
                 v(j,iglob) = v(j,iglob) + aalt2(j,iglob)*dta_2
              end do
           end if
        end do
c
        if (use_rattle)  call rattle2 (dta_2)
c
        do i = 1, nloc
          iglob = glob(i)
          if (use(iglob)) then
            xold(iglob) = x(iglob)
            yold(iglob) = y(iglob)
            zold(iglob) = z(iglob)
            x(iglob) = x(iglob) + v(1,iglob)*dta_2
            y(iglob) = y(iglob) + v(2,iglob)*dta_2
            z(iglob) = z(iglob) + v(3,iglob)*dta_2
          end if
        end do
c
        if (use_rattle) call rattle(dta_2)
        if (use_rattle) call rattle2(dta_2)
c
c     compute random part
c
        deallocate (Rn)
        allocate (Rn(3,nloc))
        do i = 1, nloc
          do j = 1, 3
            Rn(j,i) = normal()
          end do
        end do
        do i = 1, nloc
           iglob = glob(i)
           if (use(iglob)) then
              do j = 1, 3
                 v(j,iglob) = a1*v(j,iglob) + 
     $              a2*Rn(j,i)/sqrt(mass(iglob))
              end do
           end if
        end do
c
        if (use_rattle) call rattle2(dta)
c
        do i = 1, nloc
          iglob = glob(i)
          if (use(iglob)) then
            xold(iglob) = x(iglob)
            yold(iglob) = y(iglob)
            zold(iglob) = z(iglob)
            x(iglob) = x(iglob) + v(1,iglob)*dta_2
            y(iglob) = y(iglob) + v(2,iglob)*dta_2
            z(iglob) = z(iglob) + v(3,iglob)*dta_2
          end if
        end do
c
c       Reassign the particules that have changed of domain
c
c       -> real space
c
        time0 = mpi_wtime()
c
c       call reassignrespa(.true.,k,nalt2)
        call reassignrespa(k,nalt2)
c
        time1 = mpi_wtime()
        timereneig = timereneig + time1 - time0
c
c       communicate positions
c
        time0 = mpi_wtime()
        call commposrespa(k.ne.nalt2)
        time1 = mpi_wtime()
        timecommstep = timecommstep + time1 - time0
c
        allocate (derivs(3,nbloc))
        derivs = 0d0
c
        time0 = mpi_wtime()
        if (k.eq.nalt2) call mechanicsteprespa1(-1,0)
        time1 = mpi_wtime()
        timeparam = timeparam + time1 - time0
        time0 = mpi_wtime()
        if (k.eq.nalt2) call allocsteprespa(.true.)
        time1 = mpi_wtime()
        timeclear = timeclear + time1 - time0
c
c     get the fast-evolving potential energy and atomic forces
c
        call gradfastbaoab1(ealt,derivs)
c
c       communicate forces
c
        call commforcesrespa1(derivs,0)
c
c       MPI : get total energy
c
        call reduceen(ealt)
c
c     use Newton's second law to get fast-evolving accelerations;
c     update fast-evolving velocities using the BAOAB recursion
c
        do i = 1, nloc
           iglob = glob(i)
           if (use(iglob)) then
              do j = 1, 3
                 aalt2(j,iglob) = -convert *
     $              derivs(j,i) / mass(iglob)
                 v(j,iglob) = v(j,iglob) + aalt2(j,iglob)*dta_2
              end do
           end if
        end do
        deallocate (derivs)
c
        if (use_rattle)  call rattle2 (dta_2)
c
c     increment average virial from fast-evolving potential terms
c
        do i = 1, 3
           do j = 1, 3
              viralt(j,i) = viralt(j,i) + vir(j,i)/dshort
           end do
        end do
      end do
      return
      end
