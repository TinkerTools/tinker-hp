c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine beeman  --  Beeman molecular dynamics step  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "beeman" performs a single molecular dynamics time step
c     via the Beeman multistep recursion formula; uses original
c     coefficients or Bernie Brooks' "Better Beeman" values
c
c     literature references:
c
c     D. Beeman, "Some Multistep Methods for Use in Molecular
c     Dynamics Calculations", Journal of Computational Physics,
c     20, 130-139 (1976)
c
c     B. R. Brooks, "Algorithms for Molecular Dynamics at Constant
c     Temperature and Pressure", DCRT Report, NIH, April 1988
c
c
#include "tinker_macro.h"
      module beeman_inl
      contains
#include "convert.f.inc"
      end module

      subroutine beeman (istep,dt)
      use atmtyp
      use atomsMirror
      use beeman_inl
      use cutoff
      use domdec
      use deriv     ,only:info_forces,cBond,cNBond,comm_forces,tdes_l
     &              ,dr_stride,de_tot,derivx
      use energi
      use freeze
      use inform
      use mdstuf
      use mdstuf1
      use moldyn
      use mpi
      use random_mod
      use timestat
      use tinheader ,only: re_p
      use units
      use usage
      use utilgpu
      implicit none
      integer i,j,istep,iglob
      real(r_p) dt,dt_x,factor
      real(r_p) part1,part2
      mdyn_rtyp dx
c
c     set time values and coefficients for Beeman integration
c
      factor = real(bmnmix,r_p)
      dt_x   = dt / factor
      part1  = 0.5_re_p*factor + 1.0_re_p
      part2  = part1 - 2.0_re_p
c
c     store the current atom positions, then find half-step
c     velocities and full-step positions via Beeman recursion
c
      if (use_rattle) call save_atoms_pos
!$acc parallel loop collapse(2) default(present) async
      do i = 1,nloc; do j = 1,3
         iglob = glob(i)
         if (use(iglob)) then
            v(j,iglob) = v(j,iglob) + (part1*a(j,iglob)-
     $           aalt(j,iglob))*dt_x
         end if
      end do; end do
      call integrate_pos(dt)
c
c     Reassign the particules that have changed of domain
c
c     -> real space
      call reassign
c
c     -> reciprocal space
      call reassignpme(.false.)
c
c     get constraint-corrected positions and half-step velocities
c
      if (use_rattle)  call rattle (dt)
c
c     communicate positions
c
      call commpos
      call commposrec
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
c     make half-step temperature and pressure corrections
c
      call temper2 (temp)
      call pressure2 (epot,temp)
c
c     communicate forces
c
      call commforces(derivs)
c
c     use Newton's second law to get the next accelerations;
c     find the full-step velocities using the Beeman recursion
c
      if (tdes_l) then
!$acc parallel loop collapse(2) default(present) async
      do i = 1,nbloc; do j = 1,3
         iglob = glob(i)
         if (i.le.nloc) then
         if (useAll.or.use(iglob)) then
            aalt(j,iglob) = a(j,iglob)
            dx            = derivx(i + (j-1)*dr_stride )
            a(j,iglob)    = -convert * mdr2md(dx)/mass(iglob)
            v(j,iglob)    = v(j,iglob) +
     &                (part2*a(j,iglob)+aalt(j,iglob))*dt_x
            derivx(i + (j-1)*dr_stride ) = 0
         end if
         else
            derivx(i + (j-1)*dr_stride ) = 0
         end if
      end do; end do
      else
!$acc parallel loop collapse(2) default(present) async
      do i = 1,nbloc; do j = 1,3
         iglob = glob(i)
         if (i.le.nloc) then
         if (useAll.or.use(iglob)) then
            aalt(j,iglob) = a(j,iglob)
            a(j,iglob)    = -convert * mdr2md(de_tot(j,i))/mass(iglob)
            v(j,iglob)    = v(j,iglob) +
     &                (part2*a(j,iglob)+aalt(j,iglob))*dt_x
            de_tot(j,i)   = 0
         end if
         else
            de_tot(j,i)   = 0
         end if
      end do; end do
      end if
c
c     Debug print information
c
      if (deb_Energy) call info_energy(rank)
      if (deb_Force)  call info_forces(cNBond)
      if (deb_Atom)   call info_minmax_pva
      if (abort)      __TINKER_FATAL__
c
c     find the constraint-corrected full-step velocities
c
      if (use_rattle)  call rattle2 (dt)
c
c     make full-step temperature and pressure corrections
c
      call temper (dt,eksum,ekin,temp)
      call pressure (dt,ekin,pres,stress,istep)
c
c     total energy is sum of kinetic and potential energies
c
      if (calc_e) then
!$acc serial async present(etot,eksum,esum)
      etot = eksum + esum
!$acc end serial
      end if
c
c     compute statistics and save trajectory for this step
c
      call mdstat (istep,dt,etot,epot,eksum,temp,pres)
      call mdsave (istep,dt,epot)
      call mdrestgpu (istep)
      end
