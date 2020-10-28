c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine respa1  --  r-RESPA1 molecular dynamics step  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "respa1" performs a single multiple time step molecular dynamics
c     step using the reversible reference system propagation algorithm
c     (r-RESPA) via a Verlet core with the potential split into fast-
c     intermediate and slow-evolving portions
c
c     literature references:
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
#include "tinker_precision.h"
      module respa1_mod
         real(r_p),allocatable::derivs(:,:)
         real(r_p) ealt2
         real(r_p) viralt2(3,3)
      end module

      subroutine respa1(istep,dt)
      use atmtyp
      use atomsMirror
      use cutoff
      use domdec
      use deriv
      use energi  ,only: info_energy
      use freeze
      use inform
      use moldyn
      use timestat
      use respa1_mod
      use utils   ,only: set_to_zero1m
      use utilgpu ,only: prmem_requestm,rec_queue
      use units
      use usage
      use virial
      use mpi
      implicit none
      integer i,j,iglob
      integer istep
      real(r_p) dt,dt_2
      real(r_p) dta,dta_2,dta2
      real(r_p) epot,etot
      real(r_p) eksum
      real(r_p) temp,pres
      real(r_p) ealt
      real(r_p) ekin(3,3)
      real(r_p) stress(3,3)
      real(r_p) viralt(3,3)
      real(8) time0,time1
c
      if (istep.eq.1) then
!$acc enter data create(ekin,viralt,viralt2,stress)
!$acc&           create(ealt,ealt2,epot,etot,eksum,temp,pres)
      end if
c
c     set some time values for the dynamics integration
c
      dt_2 = 0.5_re_p * dt
      dta = dt / dinter
      dta_2 = 0.5_re_p * dta
c
      dta2 = dta / dshort

!$acc data present(ekin,viralt,stress,ealt,epot,etot,eksum,temp,pres)
!$acc&     present(use,glob,v,a,mass,vir,x,y,z)
c
c     store the current atom positions, then find half-step
c     velocities via velocity Verlet recursion
c
!$acc parallel loop collapse(2) async
      do i = 1, nloc
         do j = 1, 3
            iglob = glob(i)
            if (use(iglob)) then
               v(j,iglob) = v(j,iglob) + a(j,iglob)*dt_2
            end if
         end do
      end do
c
c     find intermediate-evolving velocities and positions via velocity Verlet recursion
c
      call respaint1(ealt,viralt,dta,dta2)
c
c     Reassign the particules that have changed of domain
c
c     -> real space
c      call reassignrespa(.false.,nalt,nalt)
c
c     -> reciprocal space
      call reassignpme(.false.)
c
c     communicate positions
c
c     call commposrespa(.false.)
      call commposrec
      call reCast_position
c
c
      call reinitnl(istep)
c
      call mechanicsteprespa1(istep,2)

      call allocsteprespa(.false.)
c
c     rebuild the neighbor lists
c
      if (use_list) call nblist(istep)
c
      call prmem_requestm(derivs,3,nbloc,async=.true.)
      call set_to_zero1m(derivs,3*nbloc,rec_queue)
c
c     get the slow-evolving potential energy and atomic forces
c
      call gradslow1 (epot,derivs)
c
c     communicate some forces
c
      call commforcesrespa1(derivs,2)
c
c     MPI : get total energy
c
      call reduceen(epot)
c
c     Debug print information
c
      if(deb_Energy) call info_energy(rank)
      if(deb_Force)  call info_forces(cNBond)
      if(deb_Atom)   call info_minmax_pva
c
c     make half-step temperature and pressure corrections
c
      call temper2 (temp)
c     call pressure2 (epot,temp)
c
c     use Newton's second law to get the slow accelerations;
c     find full-step velocities using velocity Verlet recursion
c
!$acc parallel loop collapse(2) async
!$acc&         present(derivs)
      do i = 1, nloc
         do j = 1, 3
            iglob = glob(i)
            if (use(iglob)) then
               a(j,iglob) = -convert * derivs(j,i) / mass(iglob)
               v(j,iglob) = v(j,iglob) + a(j,iglob)*dt_2
            end if
         end do
      end do
c
c     perform deallocation of some local arrays
c
c
c     find the constraint-corrected full-step velocities
c
      if (use_rattle)  call rattle2 (dt)
c
c     total potential and virial from sum of fast and slow parts
c
!$acc serial async
      epot = epot + ealt
      do i = 1, 3
         do j = 1, 3
            vir(j,i) = vir(j,i) + viralt(j,i)
         end do
      end do
!$acc end serial
c
c     make full-step temperature and pressure corrections
c
      call temper   (dt,eksum,ekin,temp)
      call pressure (dt,ekin,pres,stress,istep)
c
c     total energy is sum of kinetic and potential energies
c
!$acc serial async
      etot = eksum + epot
!$acc end serial
c
c     compute statistics and save trajectory for this step
c
      call mdsave (istep,dt,epot)
      call mdrestgpu (istep)
      call mdstat (istep,dt,etot,epot,eksum,temp,pres)
!$acc end data
      end
c
c
c
c     ###################################################################
c     ##                                                               ##
c     ##  subroutine gradfast1  --  fast energy & gradient components  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     "gradfast" calculates the potential energy and first derivatives
c     for the fast-evolving local valence potential energy terms
c
c
      subroutine gradfast1(energy,derivs)
      use cutoff
      use domdec,only:nbloc
      use potent
      implicit none
      real(r_p) energy
      real(r_p) derivs(3,nbloc)
      logical save_vdw,save_charge

      logical save_mpole,save_polar
      logical save_list
      logical save_smdvel, save_smdfor
c
c
c     save the original state of slow-evolving potentials
c
      save_vdw    = use_vdw
      save_charge = use_charge
      save_mpole  = use_mpole
      save_polar  = use_polar
      save_list   = use_list
      save_smdvel = use_smd_velconst
      save_smdfor = use_smd_forconst
c
c     turn off slow-evolving nonbonded potential energy terms
c
      use_vdw     = .false.
      use_charge  = .false.
      use_mpole   = .false.
      use_polar   = .false.
      use_list    = .false.
      use_smd_velconst = .false.
      use_smd_forconst = .false.
c
c     get energy and gradient for fast-evolving potential terms
c
      call gradient (energy,derivs)
c
c     restore the original state of slow-evolving potentials
c
      use_vdw    = save_vdw
      use_charge = save_charge
      use_mpole  = save_mpole
      use_polar  = save_polar
      use_list   = save_list
      use_smd_velconst = save_smdvel
      use_smd_forconst = save_smdfor
      return
      end
c
c     ##########################################################################
c     ##                                                                      ##
c     ##  subroutine gradint1  --  intermediate energy & gradient components  ##
c     ##                                                                      ##
c     ##########################################################################
c
c
c     "gradfast" calculates the potential energy and first derivatives
c     for the fast-evolving local valence potential energy terms
c
c
      subroutine gradint1 (energy,derivs)
      use cutoff
      use deriv
      use domdec,only:nbloc
      use energi
      use polpot
      use potent
      implicit none
      real(r_p) energy
      real(r_p) derivs(3,nbloc)
      real(t_p) save_tcgomega
      integer i,j
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
      logical save_smdvel, save_smdfor
c
c
c     save the original state of fast-evolving potentials
c
      save_bond   = use_bond
      save_angle  = use_angle
      save_strbnd = use_strbnd
      save_urey   = use_urey
      save_angang = use_angang
      save_opbend = use_opbend
      save_opdist = use_opdist
      save_improp = use_improp
      save_imptor = use_imptor
      save_tors   = use_tors
      save_pitors = use_pitors
      save_strtor = use_strtor
      save_tortor = use_tortor
      save_geom   = use_geom
      save_extra  = use_extra
      save_crec   = use_crec
      save_mrec   = use_mrec
      save_prec   = use_prec
      save_polalg = polalg
      save_tcgorder = tcgorder
      save_tcgprec  = tcgprec
      save_tcgguess = tcgguess
      save_tcgpeek  = tcgpeek
      save_tcgomega = tcgomega
      save_smdvel   = use_smd_velconst
      save_smdfor   = use_smd_forconst
c
c     turn off fast-evolving valence potential energy terms
c
      use_bond    = .false.
      use_angle   = .false.
      use_strbnd  = .false.
      use_urey    = .false.
      use_angang  = .false.
      use_opbend  = .false.
      use_opdist  = .false.
      use_improp  = .false.
      use_imptor  = .false.
      use_tors    = .false.
      use_pitors  = .false.
      use_strtor  = .false.
      use_tortor  = .false.
      use_geom    = .false.
      use_extra   = .false.
      use_mrec    = .false.
      use_prec    = .false.
      use_cself   = .false.
      use_mself   = .false.
      use_pself   = .false.
      use_cshortreal     = .true.
      use_mpoleshortreal = .true.
      use_vdwshort       = .true.
      use_polarshortreal = .true.
      polalg      = polalgshort
      tcgorder    = tcgordershort
      tcgprec     = tcgprecshort
      tcgguess    = tcgguessshort
      tcgpeek     = tcgpeekshort
      tcgomega    = tcgomegashort
      use_smd_velconst   = .false.
      use_smd_forconst   = .false.
c
c     get energy and gradient for slow-evolving potential terms
c
      call gradient (energy,derivs)
c
c     restore the original state of fast-evolving potentials
c
      use_bond    = save_bond
      use_angle   = save_angle
      use_strbnd  = save_strbnd
      use_urey    = save_urey
      use_angang  = save_angang
      use_opbend  = save_opbend
      use_opdist  = save_opdist
      use_improp  = save_improp
      use_imptor  = save_imptor
      use_tors    = save_tors
      use_pitors  = save_pitors
      use_strtor  = save_strtor
      use_tortor  = save_tortor
      use_geom    = save_geom
      use_extra   = save_extra
      use_mrec    = save_mrec
      use_prec    = save_prec
      use_cself   = .true.
      use_mself   = .true.
      use_pself   = .true.
      use_cshortreal     = .false.
      use_mpoleshortreal = .false.
      use_vdwshort       = .false.
      use_polarshortreal = .false.
      polalg      = save_polalg
      tcgorder    = save_tcgorder
      tcgprec     = save_tcgprec
      tcgguess    = save_tcgguess
      tcgpeek     = save_tcgpeek
      use_smd_velconst   = save_smdvel
      use_smd_forconst   = save_smdfor
c
c     also save the values of the short range real space forces
c
!$acc parallel loop collapse(2) present(desave,dep) async
      do i = 1,nbloc
         do j = 1,3
            desave(j,i) = dep(j,i)
         end do
      end do

!$acc serial present(esave,ep) async
      esave = ep
!$acc end serial
      end
c
c
c     ###################################################################
c     ##                                                               ##
c     ##  subroutine gradslow1  --  slow energy & gradient components  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     "gradslow" calculates the potential energy and first derivatives
c     for the slow-evolving nonbonded potential energy terms
c
c
      subroutine gradslow1 (energy,derivs)
      use cutoff
      use deriv
      use domdec
      use energi
      use moldyn
      use potent
      use tinheader ,only: re_p
      use virial
      implicit none
      integer i,j
      real(r_p) energy
      real(r_p) derivs(3,*)
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
      save_bond   = use_bond
      save_angle  = use_angle
      save_strbnd = use_strbnd
      save_urey   = use_urey
      save_angang = use_angang
      save_opbend = use_opbend
      save_opdist = use_opdist
      save_improp = use_improp
      save_imptor = use_imptor
      save_tors   = use_tors
      save_pitors = use_pitors
      save_strtor = use_strtor
      save_tortor = use_tortor
      save_geom   = use_geom
      save_extra  = use_extra
c
c     turn off fast-evolving valence potential energy terms
c
      use_bond    = .false.
      use_angle   = .false.
      use_strbnd  = .false.
      use_urey    = .false.
      use_angang  = .false.
      use_opbend  = .false.
      use_opdist  = .false.
      use_improp  = .false.
      use_imptor  = .false.
      use_tors    = .false.
      use_pitors  = .false.
      use_strtor  = .false.
      use_tortor  = .false.
      use_geom    = .false.
      use_extra   = .false.
      use_vdwlong = .true.
      use_clong   = .true.
      use_mpolelong = .true.

c
c     get energy and gradient for slow-evolving potential terms
c
      call gradient (energy,derivs)
c
c     restore the original state of fast-evolving potentials
c
      use_bond   = save_bond
      use_angle  = save_angle
      use_strbnd = save_strbnd
      use_urey   = save_urey
      use_angang = save_angang
      use_opbend = save_opbend
      use_opdist = save_opdist
      use_improp = save_improp
      use_imptor = save_imptor
      use_tors   = save_tors
      use_pitors = save_pitors
      use_strtor = save_strtor
      use_tortor = save_tortor
      use_geom   = save_geom
      use_extra  = save_extra
      use_vdwlong   = .false.
      use_clong     = .false.
      use_mpolelong = .false.
c
c     substract the previously stored short range energy and forces
c
!$acc serial present(esum,energy,esave,vir,virsave) async
      esum   = esum - esave
      energy = energy - esave
      do i = 1,3
         do j = 1,3
            vir    (j,i) = vir(j,i) - virsave(j,i)/dinter
            virsave(j,i) = 0.0_re_p
         end do
      end do
!$acc end serial

!$acc parallel loop collapse(2) default(present) async
      do i = 1,nbloc
         do j = 1,3
            if (i.le.nloc) derivs(j,i) = derivs(j,i) - desave(j,i)
            desum(j,i)= desum(j,i) - desave(j,i)
         end do
      end do

      end
c
c     subroutine respaint1 :
c     find intermediate-evolving velocities and positions via velocity Verlet recursion
c
      subroutine respaint1(ealt,viralt,dta,dta2)
      use atmtyp
      use cutoff
      use deriv
      use domdec
      use deriv
      use energi  ,only: info_energy
      use freeze
      use inform
      use moldyn
      use timestat
      use respa1_mod
      use utils   ,only: set_to_zero1m
      use utilgpu ,only:prmem_requestm,rec_queue
      use units
      use usage
      use virial
      use mpi
      implicit none
      integer i,j,k,iglob
      real(r_p) dta,dta_2,dta2
      real(r_p) ealt
      real(r_p) viralt(3,3)

      dta_2 = 0.5_re_p * dta
c
c     initialize virial from fast-evolving potential energy terms
c
!$acc data present(ealt,ealt2,vir,viralt2,viralt)
!$acc&     present(v,aalt,glob,use,mass)

!$acc parallel loop collapse(2) async
      do i = 1, 3
         do j = 1, 3
            viralt(j,i) = 0.0_re_p
         end do
      end do

      do k = 1, nalt
!$acc parallel loop collapse(2) async
         do i = 1, nloc
            do j = 1, 3
               iglob = glob(i)
               if (use(iglob)) then
                  v(j,iglob) = v(j,iglob) + aalt(j,iglob)*dta_2
               end if
            end do
         end do
c
c     find fast-evolving velocities and positions via velocity Verlet recursion
c
         call respafast1(ealt2,viralt2,dta2)
c
c
c     Reassign the particules that have changed of domain
c
c       -> real space
c       call reassignrespa(.false.,nalt,nalt)
c
c     Communicate positions
c
c       call commposrespa(.false.)
c
         call mechanicsteprespa1(-1,1)
c
         call prmem_requestm(derivs,3,nbloc,async=.true.)
         call set_to_zero1m (derivs,3*nbloc,rec_queue)

         call prmem_requestm(desave,3,nbloc,async=.true.)
c
c     get the fast-evolving potential energy and atomic forces
c
         call gradint1(ealt,derivs)
c
c      communicate forces
c
         call commforcesrespa1(derivs,1)
c
c      MPI : get total energy
c
         call reduceen(ealt)
c
c     Debug print information
c
         if(deb_Energy) call info_energy(rank)
         if(deb_Force ) call info_forces(cSNBond)
         if(deb_Atom  ) call info_minmax_pva
c
c     use Newton's second law to get fast-evolving accelerations;
c     update fast-evolving velocities using the Verlet recursion
c
!$acc parallel loop collapse(2) present(derivs) async
         do i = 1, nloc
            do j = 1, 3
               iglob = glob(i)
               if (use(iglob)) then
                  aalt(j,iglob) = -convert *
     $               derivs(j,i) / mass(iglob)
                  v(j,iglob) = v(j,iglob) + aalt(j,iglob)*dta_2
               end if
            end do
         end do
         if (use_rattle)  call rattle2 (dta)
c
c     increment average virial from fast-evolving potential terms
c
!$acc serial async
         ealt = ealt + ealt2
         do i = 1, 3
            do j = 1, 3
               viralt(j,i) = viralt(j,i) + (viralt2(j,i) +
     $          vir(j,i))/dinter
            end do
         end do
!$acc end serial
      end do
!$acc end data
      end
c
c     subroutine respafast1 :
c     find fast-evolving velocities and positions via velocity Verlet recursion
c
      subroutine respafast1(ealt,viralt,dta)
      use atmtyp
      use atomsMirror
      use cutoff
      use domdec
      use deriv
      use energi  ,only: info_energy
      use freeze
      use inform
      use moldyn
      use mpi
      use respa1_mod ,only: derivs
      use timestat
      use utils   ,only: set_to_zero1m
      use utilgpu ,only: prmem_requestm,rec_queue
      use units
      use usage
      use virial
      implicit none
      integer i,j,k,iglob

      real(r_p) dta,dta_2
      real(r_p) ealt
      real(r_p) viralt(3,3)

      dta_2 = 0.5_re_p * dta
c
c     initialize virial from fast-evolving potential energy terms
c
!$acc data present(ealt,viralt,vir)
!$acc&     present(use,v,x,y,z,aalt2,mass,glob,xold,yold,zold)

!$acc parallel loop collapse(2) async
      do i = 1, 3
         do j = 1, 3
            viralt(j,i) = 0.0_re_p
         end do
      end do

      do k = 1, nalt2
!$acc parallel loop async
         do i = 1, nloc
            iglob = glob(i)
            if (use(iglob)) then
!$acc loop seq
               do j = 1, 3
                  v(j,iglob) = v(j,iglob) + aalt2(j,iglob)*dta_2
               end do
               xold(iglob) = x(iglob)
               yold(iglob) = y(iglob)
               zold(iglob) = z(iglob)
               x(iglob) = x(iglob) + v(1,iglob)*dta
               y(iglob) = y(iglob) + v(2,iglob)*dta
               z(iglob) = z(iglob) + v(3,iglob)*dta
            end if
         end do
         if (use_rattle)  call rattle (dta)
c
c       Reassign the particules that have changed of domain
c
c       -> real space
        call reassignrespa(k,nalt2)
c
c       communicate positions
c
        call commposrespa(k.ne.nalt2)
        call reCast_position
c
        call prmem_requestm(derivs,3,nbloc,async=.true.)
        call set_to_zero1m(derivs,3*nbloc,rec_queue)
c
        if (k.eq.nalt2) call mechanicsteprespa1(-1,0)
        if (k.eq.nalt2) call allocsteprespa(.true.)
c
c       get the fast-evolving potential energy and atomic forces
c
        call gradfast1(ealt,derivs)
c
c       communicate forces
c
        call commforcesrespa1(derivs,0)
c
c       aMD/GaMD contributions
c
        call aMD (derivs,ealt)
c
c       MPI : get total energy
c
        call reduceen(ealt)
c
c     Debug print information
c
        if(deb_Energy) call info_energy(rank)
        if(deb_Force)  call info_forces(cBond)
        if(deb_Atom)   call info_minmax_pva
c
c     use Newton's second law to get fast-evolving accelerations;
c     update fast-evolving velocities using the Verlet recursion
c
!$acc parallel loop collapse(2) present(derivs) async
        do i = 1, nloc
           do j = 1, 3
              iglob = glob(i)
              if (use(iglob)) then
                 aalt2(j,iglob) = -convert * derivs(j,i) / mass(iglob)
                 v(j,iglob)  = v(j,iglob) + aalt2(j,iglob)*dta_2
              end if
           end do
        end do
        if (use_rattle)  call rattle2 (dta)
c
c       increment average virial from fast-evolving potential terms
c
!$acc parallel loop collapse(2) async
        do i = 1, 3
           do j = 1, 3
              viralt(j,i) = viralt(j,i) + vir(j,i)/dshort
           end do
        end do
      end do
!$acc end data
      end
