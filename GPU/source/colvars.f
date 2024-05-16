c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #############################################################
c     ##                                                         ##
c     ##  colvars Routines  --  Colvars/Tinker-HP interface      ##
c     ##                                                         ##
c     #############################################################
#include "tinker_macro.h"
      module colvars_inl
      real(8),target:: kel_vin
      contains
#include "convert.f.inc"
      end module

      subroutine set_cvatoms_ids(ncvatoms_in,cvatoms_ids_in)
      use iso_c_binding
      use atoms, only: n
      use colvars
      use mpi
      use sizes ,only: tinkerdebug
      use domdec
      use molcul
      implicit none
      integer :: ncvatoms_in
      integer ierr,i,j,k,l,mol,o,p,init,stop
      integer, dimension(ncvatoms_in) :: cvatoms_ids_in
      ncvatoms = ncvatoms_in
      if (tinkerdebug.gt.0) print*, 'set_cvatoms_ids', ncvatoms
      if (allocated(cvatoms_ids)) deallocate(cvatoms_ids)
      allocate (cvatoms_ids(ncvatoms))
      do i = 1, ncvatoms
        cvatoms_ids(i) = cvatoms_ids_in(i)+1
      end do
c      cvatoms_ids(1:ncvatoms) = cvatoms_ids_in(1:ncvatoms)
      allocate (cv_pos(3,ncvatoms))
      allocate (decv(3,ncvatoms),decv_tot(3,ncvatoms))
      cv_pos   = 0d0
      decv     = 0d0
      decv_tot = 0d0
c
c     also build the extended list of atoms involved in the CV+ones
c     of the associated molecules
      ncvatomsmol = ncvatoms 
      if (allocated(cvatomsmol)) deallocate(cvatomsmol)
      allocate (cvatomsmol(n))
      cvatomsmol(1:ncvatoms) = cvatoms_ids(1:ncvatoms)
      do i = 1, ncvatoms
        j = cvatoms_ids(i)
        mol = molcule(j)
        init = imol(1,mol)
        stop = imol(2,mol)
        do k = init, stop
          l = kmol(k)
c
c       check wether atom is already in the list
c
          do o = 1, ncvatomsmol
            p = cvatomsmol(o)
            if (p.eq.l) goto 10
          end do
          ncvatomsmol = ncvatomsmol + 1
          cvatomsmol(ncvatomsmol) = l
 10       continue
        end do
      end do
      end subroutine

c     get a copy of values rather than the same address (they don't
change during MD)
c
      subroutine get_mpi(commcv,rankcv,nproccv)
      use domdec
      use iso_c_binding
      implicit none
      integer(c_int) :: commcv,rankcv,nproccv
      commcv = COMM_TINKER
      rankcv = rank
      nproccv = nproc
      return
      end subroutine

      subroutine get_sim_temp(res)
      use bath
      use iso_c_binding
      implicit none
      real(c_double) :: res
      res = kelvin
      return
      end subroutine

      subroutine get_sim_boltzmann(res)
      use units
      use iso_c_binding
      implicit none
      real(c_double) :: res
      res = gasconst
      return
      end subroutine

      subroutine get_sim_dt(res)
      use colvars
      use iso_c_binding
      implicit none
      real(c_double) :: res
c
c     In colvars, the time step is in fs
c
      res = dt_sim
      return
      end subroutine
c
      subroutine rand_colvars(res)
      use random_mod
      use colvars
      use iso_c_binding
      implicit none
      real(c_double) :: res
      res = normal()
      end subroutine
c
      subroutine get_system_size(system_size)
      use atoms ,only: n
      use iso_c_binding
      implicit none
      integer(c_int) :: system_size
      system_size = n
      end subroutine

      subroutine get_mass_atom(atom_number,atom_mass)
      use atmtyp
      use iso_c_binding
      implicit none
      integer(c_int) :: atom_number
      real(c_double) :: atom_mass
      atom_mass = mass(atom_number)
      end subroutine

      subroutine get_pos_atom(atom_number,xr,yr,zr)
      !use atoms
      use colvars ,only: cv_pos
      use iso_c_binding
      implicit none
      integer(c_int) :: atom_number
      real(c_double) :: xr,yr,zr
      xr = cv_pos(1,atom_number)
      yr = cv_pos(2,atom_number)
      zr = cv_pos(3,atom_number)
c12   format(A,I5,3F14.6)
c     print 12,'get_pos', atom_number,xr,yr,zr
      end subroutine

      subroutine get_forces_tot(atom_number,fx,fy,fz)
      !use atoms
      use colvars
      !use deriv
      use iso_c_binding
      implicit none
      integer(c_int) :: atom_number
      real(c_double) :: fx,fy,fz
      fx = -decv_tot(1,atom_number)
      fy = -decv_tot(2,atom_number)
      fz = -decv_tot(3,atom_number)
      end subroutine

      subroutine get_input_filename(input_name,len_name)
      use colvars
      use files
      use iso_c_binding
      implicit none
      type(c_ptr) :: input_name
      integer(c_int) :: len_name
      logical :: exist
      colvarsinput = filename(1:leng)//'.colvars'
      inquire (file=colvarsinput,exist=exist)
      if (exist) then
        input_name = c_loc(colvarsinput(1:1))
        len_name = leng+8
      else
        input_name = c_loc(colvarsinput(1:1))
        len_name = 0
      end if
      end subroutine

      subroutine get_output_filename(output_name,len_name)
      use colvars
      use files
      use replicas
      use iso_c_binding
      implicit none
      type(c_ptr) :: output_name
      integer(c_int) :: len_name
      integer :: lenadd
      character*3 numberreps
      logical :: exist
      if (use_reps) then
        write(numberreps, '(i3.3)') rank_reploc
        colvarsoutput =
     $ filename(1:leng)//'_reps'//numberreps
        lenadd = 8
      else
        colvarsoutput = filename(1:leng)
        lenadd = 0
      end if
      output_name = c_loc(colvarsoutput(1:1))
      len_name = leng+lenadd
      end subroutine

      subroutine get_restart_filename(input_name,len_name)
      use colvars
      use files
      use replicas
      use iso_c_binding
      implicit none
      character*3 numberreps
      logical :: exist
      integer :: lenadd
      type(c_ptr) :: input_name
      integer(c_int) :: len_name
      if (use_reps) then
        write(numberreps, '(i3.3)') rank_reploc
        colvarsrestart =
     $ filename(1:leng)//'_reps'//numberreps//'.colvars.state'
        lenadd = 8
      else
        colvarsrestart = filename(1:leng)//'.colvars.state'
        lenadd = 0
      end if
      inquire (file=colvarsrestart,exist=exist)
      if (exist) then
        input_name = c_loc(colvarsrestart(1:1))
        len_name = leng+lenadd
      else
        input_name = c_loc(colvarsrestart(1:1))
        len_name = 0
      end if
      end subroutine

      subroutine add_energy_tinker(energy)
      use energi
      use iso_c_binding
      implicit none
      real(c_double) :: energy
!$acc serial present(ex)
      ex = ex + energy
      !if(energy.ne.0) print*, 'energy',energy
!$acc end serial
      end subroutine

      subroutine add_forces_tinker(iatom,forces)
      use colvars
      use domdec
      use iso_c_binding
      implicit none
      real(c_double), dimension(3) :: forces
      integer(c_int) :: iatom
      decv(1,iatom) = -forces(1)
      decv(2,iatom) = -forces(2)
      decv(3,iatom) = -forces(3)
      end subroutine

      subroutine pbc_image(r1,r2,dr)
      use iso_c_binding
      implicit none
      real(c_double), dimension(3) :: r2,r1,dr
#if TINKER_SINGLE_PREC + TINKER_MIXED_PREC
      real(t_p) buff(3)
      buff(1) = (r2(1)-r1(1))
      buff(2) = (r2(2)-r1(2))
      buff(3) = (r2(3)-r1(3))
      call image(buff(1),buff(2),buff(3))
      dr(:)   = buff(:)
#else
      dr(1) = (r2(1)-r1(1))
      dr(2) = (r2(2)-r1(2))
      dr(3) = (r2(3)-r1(3))
      call image(dr(1),dr(2),dr(3))
#endif
      end subroutine

      subroutine get_pbc(a,b,c)
      use boxes
      use iso_c_binding
      implicit none
      real(c_double) :: a,b,c
      a = xbox
      b = ybox
      c = zbox
      end subroutine

      subroutine fatal_error()
      use iso_c_binding
      implicit none
      call fatal
      end subroutine

      subroutine set_use_colvars(use_colvars_in)
      use colvars
      use iso_c_binding
      implicit none
      logical(c_bool) :: use_colvars_in
      use_colvars = use_colvars_in
      if (use_colvars) then
!$acc wait
      end if
      end subroutine

      subroutine get_alch_lambda(res)
      use mutant
      use iso_c_binding
      implicit none
      real(c_double) :: res
      res = lambda
      end subroutine

      subroutine set_alch_lambda(lambda_in)
      use mutant
      use iso_c_binding
      implicit none
      real(c_double) :: lambda_in
      lambda = lambda_in
      end subroutine

      subroutine get_delambda(res)
      use deriv
      use domdec  ,only: rank
      use inform  ,only: iprint,verbose
      use iounit
      use iso_c_binding
      use moldyn  ,only: step_c
      use mutant
      implicit none
      real(c_double) :: res
      logical disp
      delambda = delambdae * dlambdaelambda + delambdav * dlambdavlambda
      res      = delambda

#if (TINKER_SINGLE_PREC + TINKER_MIXED_PREC)
      disp = mod(step_c,iprint/100).eq.0
#else
      disp = .true.
#endif
      if (rank.eq.0.and.verbose.and.disp) then
        write(iout,10) delambdae,dlambdaelambda,delambdav,dlambdavlambda
 10     format('values of del/dlel | dvl/dlvl: ',2F11.3,3x,2F11.3)
        write(iout,5) lambda,delambda
 5      format('get_delambda -- values of lambda | delambda ',2F9.3)
      end if
c     write(iout,20) dlambdaelambda
c20   format('value of dlambdaelambda: ',F15.3)
c     write(iout,30) delambdav
c30   format('value of delambdav: ',F15.3)
c     write(iout,40) dlambdavlambda
c40   format('value of dlambdavlambda: ',F15.3)
c     write(iout,50) delambda
c50   format('value of delambda: ',F15.3)
      end subroutine

      subroutine get_d2edlambda2(res)
      use iounit
      use deriv
      use iso_c_binding
      use mutant
      implicit none
      real(c_double) :: res
      write(iout,5) lambda
 5    format('get_d2edlambda2 -- value of lambda: ',F15.3)
      write(iout,10) d2edlambdae2
 10   format('value of d2edlambdae2: ',F15.3)
      write(iout,20) dlambdaelambda
 20   format('value of dlambdaelambda: ',F15.3)
      write(iout,30) d2edlambdav2
 30   format('value of d2edlambdav2: ',F15.3)
      write(iout,40) dlambdavlambda
 40   format('value of dlambdavlambda: ',F15.3)
      d2edlambda2 = d2edlambdae2 * dlambdaelambda**2
     $     + d2edlambdav2 * dlambdavlambda**2
      res = d2edlambda2
      end subroutine

      subroutine apply_force_delambda(force)
      use mutant
      use iso_c_binding
      implicit none
      real(c_double) :: force
c     flambdabias is the energy derivative with respect to dE/dlambda
c     which is equal to minus the force on dE/dlambda
      flambdabias = -force
      end subroutine

c
c     the master gets all the cv positions and total forces
c
      subroutine prepare_colvars(derivs)
      use atoms    ,only: pbcWrapIdx
      use atmtyp   ,only: mass
      use atomsMirror
      use boxes
      use colvars
      use colvars_inl
      use deriv
      use domdec
      use inform
      use mpi
      use potent
      use moldyn
      use units
      use mdstuf
      implicit none
      integer i,j,k,idx,iglob,iloc,ilocrec,ierr,offr
      real(r_p) derivs(3,nbloc),temp
c
      if (deb_Path) write(*,*) 'prepare_colvars'
      decv     = 0d0
      decv_tot = 0d0
      cv_pos   = 0d0
      offr     = dr_obnbr

!$acc parallel loop gang vector async copyin(cvatoms_ids)
!$acc&         copy(decv_tot,cv_pos)
!$acc& present(x,y,z,de_tot,derivs,de_buffr,repart,repartrec,loc,locrec
!$acc&        ,pbcWrapIdx,aalt,aalt2,mass)
      do i = 1, ncvatoms
        iglob = cvatoms_ids(i)
        iloc = loc(iglob)
        if (repart(iglob).eq.rank) then
          cv_pos(1,i) = x(iglob) + pbcWrapIdx(4*(iglob-1)+1)*xbox
          cv_pos(2,i) = y(iglob) + pbcWrapIdx(4*(iglob-1)+2)*ybox
          cv_pos(3,i) = z(iglob) + pbcWrapIdx(4*(iglob-1)+3)*zbox
        end if
        if (ftot_l) then
           if ((iloc.gt.0).and.(iloc.le.nbloc)) then
             decv_tot(1,i) = mdr2md(de_tot(1,iloc))
             decv_tot(2,i) = mdr2md(de_tot(2,iloc))
             decv_tot(3,i) = mdr2md(de_tot(3,iloc))
           end if
        else
           if ((iloc.gt.0).and.(iloc.le.nbloc)) then
             decv_tot(1,i) = derivs(1,iloc)
             decv_tot(2,i) = derivs(2,iloc)
             decv_tot(3,i) = derivs(3,iloc)
           end if
        end if
        if (mts .and. repart(iglob).eq.rank) then
           temp = -mass(iglob)/convert
           decv_tot(1,i) = decv_tot(1,i) + 
     &                      temp*(aalt(1,iglob)+aalt2(1,iglob))
           decv_tot(2,i) = decv_tot(2,i) +
     &                      temp*(aalt(2,iglob)+aalt2(2,iglob))
           decv_tot(3,i) = decv_tot(3,i) + 
     &                      temp*(aalt(3,iglob)+aalt2(3,iglob))
         end if
c
c       also add reciprocal part of the forces (in a separate array for now)
c
        if (repartrec(iglob).eq.rank) then
          ilocrec = locrec(iglob)
!$acc loop seq
          do k = 0,dr_nbnbr-1; do j = 1,3;
             idx = offr + j+(ilocrec-1)*3 + k*dr_stride
             decv_tot(j,i) = decv_tot(j,i)+ de_buffr(idx)
          end do; end do
        end if
      end do
!$acc wait
c
      if (rank.eq.0) then
        call MPI_REDUCE(MPI_IN_PLACE,cv_pos,3*ncvatoms,MPI_RPREC
     $                 ,MPI_SUM,0,COMM_TINKER,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,decv_tot,3*ncvatoms,MPI_RPREC
     $                 ,MPI_SUM,0,COMM_TINKER,ierr)
      else
        call MPI_REDUCE(cv_pos,cv_pos,3*ncvatoms,MPI_RPREC
     $                 ,MPI_SUM,0,COMM_TINKER,ierr)
        call MPI_REDUCE(decv_tot,decv_tot,3*ncvatoms,MPI_RPREC
     $                 ,MPI_SUM,0,COMM_TINKER,ierr)
      end if
      if (use_lambdadyn) then
        delambdae = delambdae + delambdaesave
        delambdaesave = 0d0
        delambdav = delambdav + delambdavsave
        delambdavsave = 0d0
        if (rank.eq.0) then
          call MPI_REDUCE(MPI_IN_PLACE,delambdav,1,MPI_RPREC
     $                   ,MPI_SUM,0,COMM_TINKER,ierr)
          call MPI_REDUCE(MPI_IN_PLACE,delambdae,1,MPI_RPREC
     $                   ,MPI_SUM,0,COMM_TINKER,ierr)
        else
          call MPI_REDUCE(delambdav,delambdav,1,MPI_RPREC
     $                   ,MPI_SUM,0,COMM_TINKER,ierr)
          call MPI_REDUCE(delambdae,delambdae,1,MPI_RPREC
     $                   ,MPI_SUM,0,COMM_TINKER,ierr)
        end if
      end if
c      if (use_osrw) then
c        if (rank.eq.0) then
c          call MPI_REDUCE(MPI_IN_PLACE,d2edlambdav2,1,MPI_RPREC
c     $                   ,MPI_SUM,0,COMM_TINKER,ierr)
c          call MPI_REDUCE(MPI_IN_PLACE,d2edlambdae2,1,MPI_RPREC
c     $                   ,MPI_SUM,0,COMM_TINKER,ierr)
c        else
c          call MPI_REDUCE(d2edlambdav2,d2edlambdav2,1,MPI_RPREC
c     $                   ,MPI_SUM,0,COMM_TINKER,ierr)
c          call MPI_REDUCE(d2edlambdae2,d2edlambdae2,1,MPI_RPREC
c     $                   ,MPI_SUM,0,COMM_TINKER,ierr)
c        end if
c      end if
      end subroutine

      subroutine prepare_colvars1(derivs)
      use atoms       ,only: pbcWrapIdx
      use atmtyp      ,only: mass
      use atomsMirror
      use colvars
      use colvars_inl
      use cell        ,only: xcell,ycell,zcell
      use deriv
      use domdec
      use inform
      use mpi
      use potent
      use mutant
      use moldyn
      use units
      use mdstuf
      implicit none
      integer i,j,k,idx,iglob,iloc,ilocrec,ierr,offr
      real(r_p) derivs(3,nbloc),temp
c
      if (deb_Path) write(*,*) 'prepare_colvars1'

      decv     = 0d0
      decv_tot = 0d0
      cv_pos   = 0d0
      offr     = dr_obnbr


!$acc parallel loop gang vector copyin(cvatoms_ids)
!$acc&         copy(decv_tot,cv_pos)
!$acc&         present(x,y,z,de_tot,derivs,de_buffr,pbcWrapIdx,
!$acc&         mass,aalt,aalt2) async
      do i = 1, ncvatoms
        iglob = cvatoms_ids(i)
        cv_pos(1,i) = x(iglob) + pbcWrapIdx(4*(iglob-1)+1)*xcell
        cv_pos(2,i) = y(iglob) + pbcWrapIdx(4*(iglob-1)+2)*ycell
        cv_pos(3,i) = z(iglob) + pbcWrapIdx(4*(iglob-1)+3)*zcell
        if (ftot_l) then
          decv_tot(1,i) = mdr2md(de_tot(1,iglob))
          decv_tot(2,i) = mdr2md(de_tot(2,iglob))
          decv_tot(3,i) = mdr2md(de_tot(3,iglob))
        else
          decv_tot(1,i) = derivs(1,iglob)
          decv_tot(2,i) = derivs(2,iglob)
          decv_tot(3,i) = derivs(3,iglob)
        end if

        if (mts) then
           temp = -mass(iglob)/convert
           decv_tot(1,i) = decv_tot(1,i) + 
     &                        temp*(aalt(1,iglob)+aalt2(1,iglob))
           decv_tot(2,i) = decv_tot(2,i) + 
     &                        temp*(aalt(2,iglob)+aalt2(2,iglob))
           decv_tot(3,i) = decv_tot(3,i) + 
     &                        temp*(aalt(3,iglob)+aalt2(3,iglob))
         end if
        !also add reciprocal part of the forces 
        !     (in a separate array for now)
!$acc loop seq
        do k = 0,dr_nbnbr-1
!$acc loop seq
           do j = 1,3
             idx = offr + j+(iglob-1)*3 + k*dr_stride
             decv_tot(j,i) = decv_tot(j,i)+ de_buffr(idx)
           end do
        end do
      end do
!$acc wait
c12   format(A,3f10.4)
c13   format(A,6f16.10)
c     print 13, 'pos ',cv_pos(1:6,1)
c     print 13, 'frc_',decv_tot(1:6,1)
      if (use_lambdadyn) then
        delambdae = delambdae + delambdaesave
        delambdaesave = 0d0
        delambdav = delambdav + delambdavsave
        delambdavsave = 0d0
      end if
      end subroutine
c
c     the master sends the new forces, he already has the bias
c
      subroutine distrib_colvars(derivs)
      use atoms
      use colvars
      use colvars_inl
      use deriv
      use domdec
      use inform
      use mpi
      use mutant
      use potent
      use virial
      implicit none
      real(r_p) derivs(3,nbloc)
      real(r_p) vxx,vxy,vxz,vyy,vyz,vzz
      integer   i,iglob,iloc,ierr
c
      call MPI_BCAST(decv,3*ncvatoms,MPI_RPREC,0,COMM_TINKER,ierr)
      if (deb_Path) print*, 'distrib_colvars'

      if (use_lambdadyn) then
         call MPI_BCAST(lambda,1,MPI_TPREC,0,COMM_TINKER,ierr)
         call def_lambdadyn
         if (use_osrw) then
            call MPI_BCAST(flambdabias,1,MPI_TPREC,0,COMM_TINKER,ierr)
            do i = 1, nloc
              dxdelambda(:,i)= ( dxdelambdae(:,i)+dxdelambdav(:,i) )
     $                        *dlambdavlambda
              derivs(:,i)    = derivs(:,i) + flambdabias*dxdelambda(:,i)
            end do
         end if
      end if
c
      if (ftot_l) then
!$acc parallel loop async copyin(cvatoms_ids,decv)
!$acc&         present(de_tot,repart,loc)
         do i = 1, ncvatoms
           iglob = cvatoms_ids(i)
           if (repart(iglob).eq.rank) then
             iloc = loc(iglob)
             de_tot(1,iloc) = de_tot(1,iloc) + rp2mdr(decv(1,i))
             de_tot(2,iloc) = de_tot(2,iloc) + rp2mdr(decv(2,i))
             de_tot(3,iloc) = de_tot(3,iloc) + rp2mdr(decv(3,i))
           end if
         end do
      else
!$acc parallel loop async copyin(cvatoms_ids,decv)
!$acc&         present(derivs,repart,loc,x,y,z,vir)
         do i = 1, ncvatoms
           iglob = cvatoms_ids(i)
           if (repart(iglob).eq.rank) then
             iloc = loc(iglob)
             derivs(1,iloc) = derivs(1,iloc) + decv(1,i)
             derivs(2,iloc) = derivs(2,iloc) + decv(2,i)
             derivs(3,iloc) = derivs(3,iloc) + decv(3,i)
c
c         add virial contribution from colvars
c
             vxx =  x(iglob)*decv(1,i)
             vxy =  y(iglob)*decv(1,i)
             vxz =  z(iglob)*decv(1,i)
             vyy =  y(iglob)*decv(2,i)
             vyz =  z(iglob)*decv(2,i)
             vzz =  z(iglob)*decv(3,i)
             vir(1,1) = vir(1,1) + vxx
             vir(2,1) = vir(2,1) + vxy
             vir(3,1) = vir(3,1) + vxz
             vir(1,2) = vir(1,2) + vxy
             vir(2,2) = vir(2,2) + vyy
             vir(3,2) = vir(3,2) + vyz
             vir(1,3) = vir(1,3) + vxz
             vir(2,3) = vir(2,3) + vyz
             vir(3,3) = vir(3,3) + vzz
           end if
         end do
      end if
!$acc wait
      end subroutine

      subroutine distrib_colvars1(derivs)
      use atoms
      use colvars
      use colvars_inl
      use deriv
      use domdec
      use inform
      use mpi
      use mutant
      use potent
      use virial
      implicit none
      real(r_p) derivs(3,nbloc)
      real(r_p) vxx,vxy,vxz,vyy,vyz,vzz
      integer   i,iglob,iloc,ierr
c
      if (deb_Path) print*, 'distrib_colvars1'

      if (use_lambdadyn) then
         call def_lambdadyn
         if (use_osrw) then
            do i = 1, nloc
              dxdelambda(:,i)= ( dxdelambdae(:,i)+dxdelambdav(:,i) )
     $                        *dlambdavlambda
              ! FIXME
              derivs(:,i)    = derivs(:,i) + flambdabias*dxdelambda(:,i)
            end do
         end if
      end if
c
      if (ftot_l) then
!$acc parallel loop async copyin(cvatoms_ids,decv) present(de_tot)
         do i = 1, ncvatoms
             iglob = cvatoms_ids(i)
             de_tot(1,iglob) = de_tot(1,iglob) + rp2mdr(decv(1,i))
             de_tot(2,iglob) = de_tot(2,iglob) + rp2mdr(decv(2,i))
             de_tot(3,iglob) = de_tot(3,iglob) + rp2mdr(decv(3,i))
         end do
c12      format(A,6f16.10)
c        print 12, 'decv',decv(1:6,1)
c        print 12, 'de_t',de_tot(642*3+1:642*3+6,1)
      else
!$acc parallel loop async copyin(cvatoms_ids,decv) present(derivs)
         do i = 1, ncvatoms
           iglob = cvatoms_ids(i)
           derivs(1,iglob) = derivs(1,iglob) + decv(1,i)
           derivs(2,iglob) = derivs(2,iglob) + decv(2,i)
           derivs(3,iglob) = derivs(3,iglob) + decv(3,i)
c
c         add virial contribution from colvars
c
           vxx =  x(iglob)*decv(1,i)
           vxy =  y(iglob)*decv(1,i)
           vxz =  z(iglob)*decv(1,i)
           vyy =  y(iglob)*decv(2,i)
           vyz =  z(iglob)*decv(2,i)
           vzz =  z(iglob)*decv(3,i)
           vir(1,1) = vir(1,1) + vxx
           vir(2,1) = vir(2,1) + vxy
           vir(3,1) = vir(3,1) + vxz
           vir(1,2) = vir(1,2) + vxy
           vir(2,2) = vir(2,2) + vyy
           vir(3,2) = vir(3,2) + vyz
           vir(1,3) = vir(1,3) + vxz
           vir(2,3) = vir(2,3) + vyz
           vir(3,3) = vir(3,3) + vzz
         end do
      end if
!$acc wait
c!$acc data copy(decv,decv_tot)
c      call minmaxone(decv    ,size(decv),'dcv  ')
c      call minmaxone(decv_tot,size(decv),'dcv_t')
c!$acc end data
      end subroutine

      subroutine colvars_init(dt)
      use colvars
      use domdec
      use mpi
      use tinheader ,only: re_p
      implicit none
      real(r_p),intent(in):: dt
      integer ierr

      use_colvars = .false.
#ifdef COLVARS
      dt_sim = dt*1000_re_p
c
c     only the master does colvars computations, but other ranks need to allocate coord arrays
c
      if (rank.eq.0) then
         call allocate_colvars
      end if
      if (nproc.gt.1) 
     &   call MPI_BCAST(use_colvars,1,MPI_LOGICAL,0,COMM_TINKER,ierr)
      if (use_colvars) then
!$acc wait
        if (nproc.gt.1) 
     &     call MPI_BCAST(ncvatoms,1,MPI_INT,0,COMM_TINKER,ierr)
        if (rank.gt.0) then
          allocate (cvatoms_ids(ncvatoms))
        end if
        if (nproc.gt.1) call MPI_BCAST
     &     (cvatoms_ids,ncvatoms,MPI_INT,0,COMM_TINKER,ierr)

        if (rank.gt.0) then
           allocate (cv_pos(3,ncvatoms))
           allocate (decv(3,ncvatoms),decv_tot(3,ncvatoms))
           cv_pos   = 0d0
           decv     = 0d0
           decv_tot = 0d0
        end if
      end if
#endif
      end subroutine

      subroutine colvars_run(derivs)
      use colvars
      use domdec
      use deriv
      implicit none
      real(r_p) derivs(3,nbloc)
#ifdef COLVARS
      if (nproc.eq.1) then
         call prepare_colvars1(derivs)
      else
         call prepare_colvars(derivs)
      end if
 
      !only the master does colvars computations
      if (rank.eq.0) call compute_colvars_tinker()

      if (nproc.eq.1) then
         call distrib_colvars1(derivs)
      else
         call distrib_colvars(derivs)
      end if
#else
 15   format(/,
     &"                ***  FATAL ERROR *** ",/,
     &" ----------------------------------------------------",/,
     &" Colvars Feature is unavailable with this build !",/,
     &" Please Rebuild Tinker-HP with",/,
     &"       **** COLVARS_SUPPORT=1  ****",/,
     &" ----------------------------------------------------")
      if (use_colvars) then
         if (rank.eq.0) write(0,15)
         call fatal
      end if
#endif
      end subroutine
c
c     retreive the use of replicas
c
      subroutine get_use_reps(doreplica)
      use replicas
      use iso_c_binding
      implicit none
      logical(c_bool) :: doreplica
      doreplica = use_reps
      end subroutine
c
c     retreive index of local replica
c
      subroutine get_index_rep(index)
      use replicas
      use iso_c_binding
      implicit none
      integer(c_int) :: index
      index = rank_reploc
      end subroutine
c
c     retreive number of replicas
c
      subroutine get_num_rep(inter_num)
      use replicas
      use iso_c_binding
      implicit none
      integer(c_int) :: inter_num
      inter_num = nreps
      end subroutine
c
c     retreive root 2 root communicator
c
      subroutine get_root2root(intercomm)
      use domdec
      use replicas
      use iso_c_binding
      implicit none
c     targeted integers, to make pgfortran and gfortran happy with their c_bindings
      integer(c_int) :: intercomm
      intercomm = COMM_ROOT2ROOT
      return
      end subroutine
