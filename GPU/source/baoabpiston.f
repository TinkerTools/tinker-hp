c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################################################
c     ##                                                                                           ##
c     ##  subroutine baoabpiston  --  BAOAB/Langein Piston NPT Langevin molecular dynamics step    ##
c     ##                                                                                           ##
c     ###############################################################################################
c
c
c     "baoabpiston" performs a single molecular dynamics time step
c     via the BAOAB recursion formula and using a Langevin Piston barostat
c
c     literature reference:
c
c     Efficient molecular dynamics using geodesic integration 
c     and solvent-solute splitting B. Leimkuhler and C. Matthews,
c     Proceedings of the Royal Society A, 472: 20160138, 2016
c
c     
c      Constant pressure molecular dynamics simulation: The Langevin piston method
c      J. Chem. Phys. 103, 4613 (1995)
c      Scott E. Feller, Yuhong Zhang,Richard W. Pastor, Bernie Brooks
c
#include "tinker_macro.h"
      subroutine baoabpiston (istep,dt)
      use atmtyp
      use atomsMirror
      use bath
      use cutoff
      use domdec
      use deriv    ,only:info_forces,cDef,ftot_l,comm_forces
      use energi   ,only: info_energy,calc_e,chk_energy_fluct
      use freeze
      use inform
      use langevin
      use mdstuf
      use mdstuf1
      use moldyn
      use mpi
      use random_mod
      use timestat
      use tinMemory,only: prmem_request,prmem_requestm
      use units
      use usage
      use utils    ,only: set_to_zero1m
      use utilgpu  ,only: rec_queue
      use virial
      implicit none
      integer i,j,istep,iglob,ierr
      real(r_p) dt,dt_2
     &         ,a1,a2,a1piston,a2piston,R
c
c     set time values and coefficients for BAOAB integration
c
      a1       = exp(-gamma*dt)
      a2       = sqrt((1-a1**2)*boltzmann*kelvin)
      a1piston = exp(-gammapiston*dt)
      a2piston = sqrt((1-a1piston**2)*boltzmann*kelvin)
      dt_2     = 0.5*dt
c
c     find quarter step velocities and half step positions via BAOAB recursion
c
      vextvol = vextvol + dt_2*aextvol
      call integrate_vel( a,dt_2 )
c
      if (use_rattle) call rattle2(dt)
c
      extvolold = extvol
      extvol    = extvol + vextvol*dt_2
      call rescale(istep)
c
      if (use_rattle) call save_atoms_pos
      call integrate_pos( dt_2 )
c
      if (use_rattle) then
         call rattle(dt_2)
         call rattle2(dt_2)
         call rattle2(dt)
      end if
c
c     compute random part
c
      call prmem_request(Rn,3,nloc+1,async=.true.)
#ifdef _OPENACC
      call normalgpu(Rn(1,1),3*nloc)
#endif
      if (host_rand_platform) then
         call normalvec(Rn,3*nloc)
!$acc update device(Rn) async
      end if
!$acc parallel loop collapse(2) present(glob,use,v,mass) async
      do i = 1, nloc; do j = 1, 3
         iglob = glob(i)
         if (use(iglob)) then
            v(j,iglob) = a1*v(j,iglob) + 
     $                   a2*Rn(j,i)/sqrt(mass(iglob))
         end if
      end do; end do
      if (rank.eq.0) then
        R       = real(normal(),r_p)
        vextvol = a1piston*vextvol + a2piston*R/sqrt(masspiston)
      end if
      call MPI_BCAST(vextvol,1,MPI_RPREC,0,COMM_TINKER,ierr)
c
      if (use_rattle) call rattle2(dt)
c
c     find full step positions via BAOAB recursion
c
      if (use_rattle) call save_atoms_pos
      call integrate_pos( dt_2 )
c
      if (use_rattle) call rattle(dt_2)
      if (use_rattle) call rattle2(dt_2)

      extvolold = extvol
      extvol    = extvol + vextvol*dt_2
      call rescale(istep)
c
c     Reassign the particules that have changed of domain
c
c     -> real space
      call reassign
c     -> reciprocal space
      call reassignpme(.false.)
c
c     communicate positions
c
      call commpos
      call commposrec
c
      if (.not.ftot_l) then
         call prmem_requestm(derivs,3,nbloc,async=.true.)
         call set_to_zero1m(derivs,3*nbloc,rec_queue)
      end if
c
      call reinitnl(istep)
c
      call mechanicstep(istep)
c
      call allocstep
c
c     rebuild the neighbor lists
c
      if (use_list) call nblist(istep)
c
c     get the potential energy and atomic forces
c
      call gradient (epot,derivs)
c
c     MPI : get total energy
c
      call reduceen(epot)
c
c     communicate forces
c
      call comm_forces( derivs )
c
c     Debug print information
c
      if (deb_Energy) call info_energy(rank)
      if (deb_Force)  call info_forces(cDef)
      if (deb_Atom)   call info_minmax_pva
      if (abort)      call emergency_save
      if (abort)      __TINKER_FATAL__
c
c     use Newton's second law to get the next accelerations;
c     find the full-step velocities using the BAOAB recursion
c
      call integrate_vel(derivs,a,dt_2)
c
      call temper   (dt,eksum,ekin,temp)
      call pressure (dt,ekin,pres,stress,istep)
c
c     find the constraint-corrected full-step velocities
c
      if (use_rattle)  call rattle2 (dt)
c
c     total energy is sum of kinetic and potential energies
c
      if (calc_e) then
!$acc serial present(epot,eksum,etot) async
         etot = eksum + epot
!$acc end serial
      end if
c
c     compute statistics and save trajectory for this step
c
      call mdstat (istep,dt,etot,epot,eksum,temp,pres)
      call mdsave (istep,dt,epot,derivs)
      call mdrestgpu (istep)
c
!$acc update host(pres) async
!$acc wait
      aextvol = convert*(pres-atmsph)/(prescon*masspiston)
      vextvol = vextvol + dt_2*aextvol

      end
