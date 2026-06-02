!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #############################################################
!     ##                                                         ##
!     ##  subroutine beeman  --  Beeman molecular dynamics step  ##
!     ##                                                         ##
!     #############################################################
!
!
!     "beeman" performs a single molecular dynamics time step
!     via the Beeman multistep recursion formula; uses original
!     coefficients or Bernie Brooks' "Better Beeman" values
!
!     literature references:
!
!     D. Beeman, "Some Multistep Methods for Use in Molecular
!     Dynamics Calculations", Journal of Computational Physics,
!     20, 130-139 (1976)
!
!     B. R. Brooks, "Algorithms for Molecular Dynamics at Constant
!     Temperature and Pressure", DCRT Report, NIH, April 1988
!
!
#include "tinker_macro.h"
module beeman_inl
contains
#include "convert.inc.f90"
end module

subroutine beeman (istep,dt)
   use atmtyp
   use atomsMirror
   use beeman_inl
   use cutoff
   use domdec
   use deriv     ,only:info_forces,cBond,cNBond,comm_forces,tdes_l&
      &,dr_stride,de_tot,derivx
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
!
!     set time values and coefficients for Beeman integration
!
   factor = real(bmnmix,r_p)
   dt_x   = dt / factor
   part1  = 0.5_re_p*factor + 1.0_re_p
   part2  = part1 - 2.0_re_p
!
!     store the current atom positions, then find half-step
!     velocities and full-step positions via Beeman recursion
!
   if (use_rattle) call save_atoms_pos
!$acc parallel loop collapse(2) default(present) async
   do i = 1,nloc; do j = 1,3
         iglob = glob(i)
         if (use(iglob)) then
            v(j,iglob) = v(j,iglob) + (part1*a(j,iglob)-&
               &aalt(j,iglob))*dt_x
         end if
      end do; end do
   call integrate_pos(dt)
!
!     Reassign the particules that have changed of domain
!
!     -> real space
   call reassign
!
!     -> reciprocal space
   call reassignpme(.false.)
!
!     get constraint-corrected positions and half-step velocities
!
   if (use_rattle)  call rattle (dt)
!
!     communicate positions
!
   call commpos
   call commposrec
!
   call reinitnl(istep)
!
   call mechanicstep(istep)
!
   call allocstep
!
!     rebuild the neighbor lists
!
   if (use_list) call nblist(istep)
!
!     get the potential energy and atomic forces
!
   call gradient (epot,derivs)
!
!     MPI : get total energy
!
   call reduceen(epot)
!
!     make half-step temperature and pressure corrections
!
   call temper2 (temp)
   call pressure2 (epot,temp)
!
!     communicate forces
!
   call comm_forces(derivs)
!
!     use Newton's second law to get the next accelerations;
!     find the full-step velocities using the Beeman recursion
!
   if (tdes_l) then
!$acc parallel loop collapse(2) default(present) async
      do i = 1,nbloc; do j = 1,3
            iglob = glob(i)
            if (i.le.nloc) then
               if (useAll.or.use(iglob)) then
                  aalt(j,iglob) = a(j,iglob)
                  dx            = derivx(i + (j-1)*dr_stride )
                  a(j,iglob)    = -convert * mdr2md(dx)/mass(iglob)
                  v(j,iglob)    = v(j,iglob) +&
                     &(part2*a(j,iglob)+aalt(j,iglob))*dt_x
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
                  v(j,iglob)    = v(j,iglob) +&
                     &(part2*a(j,iglob)+aalt(j,iglob))*dt_x
                  de_tot(j,i)   = 0
               end if
            else
               de_tot(j,i)   = 0
            end if
         end do; end do
   end if
!
!     Debug print information
!
   if (deb_Energy) call info_energy(rank)
   if (deb_Force)  call info_forces(cNBond)
   if (deb_Atom)   call info_minmax_pva
   if (abort)      __TINKER_FATAL__
!
!     find the constraint-corrected full-step velocities
!
   if (use_rattle)  call rattle2 (dt)
!
!     make full-step temperature and pressure corrections
!
   call temper (dt,eksum,ekin,temp)
   call pressure (dt,ekin,pres,stress,istep)
!
!     total energy is sum of kinetic and potential energies
!
   if (calc_e) then
!$acc serial async present(etot,eksum,esum)
      etot = eksum + esum
!$acc end serial
   end if
!
!     compute statistics and save trajectory for this step
!
   call mdstat (istep,dt,etot,epot,eksum,temp,pres)
   call mdsave (istep,dt,epot)
   call mdrestgpu (istep)
end
