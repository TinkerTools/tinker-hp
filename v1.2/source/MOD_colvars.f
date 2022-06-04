c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##################################################################
c     ##                                                              ##
c     ##  module colvars  --  Colvars/Tinker-HP interface             ##
c     ##                                                              ##
c     ##################################################################
c
      module colvars
      use iso_c_binding
      use bath
      implicit none
      interface
        subroutine allocate_colvars () bind
     $       ( C, name = "allocate_colvars" )
        end subroutine allocate_colvars
        subroutine compute_colvars_tinker () bind
     $  ( C, name = "compute_colvars" )
        end subroutine compute_colvars_tinker
        subroutine delete_colvars () bind
     $  ( C, name = "delete_colvars" )
        end subroutine delete_colvars
      end interface
c
      logical :: use_colvars
      integer :: ncvatoms
      integer, allocatable :: cvatoms_ids(:)
      real*8, allocatable :: cv_pos(:,:),decv(:,:),decv_tot(:,:)
      real*8, target :: dt_sim
      real*8, target :: temp_rand
      save
      end

      subroutine set_cvatoms_ids(ncvatoms_in,cvatoms_ids_in)
      use iso_c_binding
      use colvars
      use mpi
      use domdec
      implicit none
      integer :: ncvatoms_in
      integer ierr
      integer, dimension(ncvatoms_in) :: cvatoms_ids_in
      ncvatoms = ncvatoms_in
      allocate (cvatoms_ids(ncvatoms))
      cvatoms_ids(1:ncvatoms) = cvatoms_ids_in(1:ncvatoms)
      allocate (cv_pos(3,ncvatoms))
      allocate (decv(3,ncvatoms),decv_tot(3,ncvatoms))
      cv_pos = 0d0
      decv = 0d0
      decv_tot = 0d0
      return
      end subroutine

      subroutine get_mpi(commcv,rankcv,nproccv)
      use domdec
      use iso_c_binding
      implicit none
c     targeted integers, to make pgfortran and gfortran happy with their c_bindings
      type(c_ptr) :: commcv,rankcv,nproccv
      commcv = c_loc(COMM_TINKER)
      rankcv = c_loc(rank)
      nproccv = c_loc(nproc)
      return
      end subroutine



      subroutine get_sim_temp(res)
      use bath
      use iso_c_binding
      implicit none
      real*8, target :: kel_vin
      type(c_ptr) :: res
      kel_vin = kelvin
      res = c_loc(kel_vin)
      return
      end subroutine

      subroutine get_sim_boltzmann(res)
      use units
      use iso_c_binding
      implicit none
c     targeted real, to make pgfortran and gfortran happy with their c_bindings
c     must use a real variable to get gasconst
c     because c_bindings doesn't like const arg to c_loc
      real*8, target :: gas_const = gasconst
      type(c_ptr) :: res
c      res = c_loc(boltzmann)
      res = c_loc(gas_const)
      return
      end subroutine

      subroutine get_sim_dt(res)
      use colvars
      use iso_c_binding
      implicit none
      type(c_ptr) :: res
c
c     In colvars, the time step is in fs
c
      res = c_loc(dt_sim)
      return
      end subroutine
c
      subroutine rand_colvars(res)
      use colvars
      use iso_c_binding
      implicit none
      real(c_double) :: res
      real*8 normal
      res = normal()
      return
      end subroutine
c
      subroutine get_system_size(system_size)
      use atoms
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
      use colvars
      use iso_c_binding
      implicit none
      integer(c_int) :: atom_number
      real(c_double) :: xr,yr,zr
      xr = cv_pos(1,atom_number)
      yr = cv_pos(2,atom_number)
      zr = cv_pos(3,atom_number)
      end subroutine

      subroutine get_forces_tot(atom_number,fx,fy,fz)
      use atoms
      use colvars
      use deriv
      use iso_c_binding
      implicit none
      integer(c_int) :: atom_number
      real(c_double) :: fx,fy,fz
      fx = -decv_tot(1,atom_number)
      fy = -decv_tot(2,atom_number)
      fz = -decv_tot(3,atom_number)
      end subroutine

      subroutine get_input_filename(input_name,len_name)
      use files
      use iso_c_binding
      implicit none
      type(c_ptr) :: input_name
      integer(c_int) :: len_name
      character*240 filenamecolvars
      logical :: exist
      filenamecolvars = filename(1:leng)//'.colvars'
      inquire (file=filenamecolvars,exist=exist)
      if (exist) then
        input_name = c_loc(filename)
        len_name = leng
      else
        input_name = c_loc(filename)
        len_name = 0
      end if
      end subroutine

      subroutine get_restart_filename(input_name,len_name)
      use files
      use iso_c_binding
      implicit none
      character*240 filenamerestart
      logical :: exist
      type(c_ptr) :: input_name
      integer(c_int) :: len_name
      filenamerestart = filename(1:leng)//'.colvars.state'
      inquire (file=filenamerestart,exist=exist)
      if (exist) then
        input_name = c_loc(filenamerestart)
        len_name = leng+14
      else
        input_name = c_loc(filenamerestart)
        len_name = 0
      end if
      end subroutine

      subroutine add_energy_tinker(energy)
      use energi
      use iso_c_binding
      implicit none
      real(c_double) :: energy
      esum = esum + energy
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
      dr(1) = (r2(1)-r1(1))
      dr(2) = (r2(2)-r1(2))
      dr(3) = (r2(3)-r1(3))
      call image(dr(1),dr(2),dr(3))
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
      use iounit
      use deriv
      use iso_c_binding
      use mutant
      implicit none
      real(c_double) :: res
      write(iout,5) lambda
 5    format('get_delambda -- value of lambda: ',F15.3)
      write(iout,10) delambdae
 10   format('value of delambdae: ',F15.3)
      write(iout,20) dlambdaelambda
 20   format('value of dlambdaelambda: ',F15.3)
      write(iout,30) delambdav
 30   format('value of delambdav: ',F15.3)
      write(iout,40) dlambdavlambda
 40   format('value of dlambdavlambda: ',F15.3)
      delambda = delambdae * dlambdaelambda + delambdav * dlambdavlambda
      write(iout,50) delambda
 50   format('value of delambda: ',F15.3)
      res = delambda
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
c       flambdabias is the energy derivative with respect to dE/dlambda
c       which is equal to minus the force on dE/dlambda
        flambdabias = -force
        end subroutine

c
c     the master gets all the cv positions and total forces
c
      subroutine prepare_colvars
      use atoms
      use colvars
      use deriv
      use domdec
      use mpi
      use potent
      implicit none
      integer i,iglob,iloc,ilocrec,ierr
c
c
      decv =0d0
      decv_tot = 0d0
      cv_pos = 0d0

      do i = 1, ncvatoms
        iglob = cvatoms_ids(i)
        iloc = loc(iglob)
        if (repart(iglob).eq.rank) then
          cv_pos(1,i) = x(iglob)
          cv_pos(2,i) = y(iglob)
          cv_pos(3,i) = z(iglob)
        end if
        if ((iloc.gt.0).and.(iloc.le.nbloc)) then
          decv_tot(1,i) = desum(1,iloc)
          decv_tot(2,i) = desum(2,iloc)
          decv_tot(3,i) = desum(3,iloc)
        end if
c
c       also add reciprocal part of the forces (in a separate array for now)
c
        ilocrec = locrec(iglob)
        if ((ilocrec.gt.0).and.(ilocrec.le.nlocrec2)) then
            decv_tot(1,i) = decv_tot(1,i) + decrec(1,ilocrec) + 
     $       demrec(1,ilocrec) + deprec(1,ilocrec) + dedsprec(1,ilocrec)
            decv_tot(2,i) = decv_tot(2,i) + decrec(2,ilocrec) + 
     $       demrec(2,ilocrec) + deprec(2,ilocrec) + dedsprec(2,ilocrec)
            decv_tot(3,i) = decv_tot(3,i) + decrec(3,ilocrec) + 
     $       demrec(3,ilocrec) + deprec(3,ilocrec) + dedsprec(3,ilocrec)
        end if
      end do
c
      if (rank.eq.0) then
        call MPI_REDUCE(MPI_IN_PLACE,cv_pos,3*ncvatoms,MPI_REAL8,
     $     MPI_SUM,0,COMM_TINKER,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,decv_tot,3*ncvatoms,MPI_REAL8,
     $     MPI_SUM,0,COMM_TINKER,ierr)
      else
        call MPI_REDUCE(cv_pos,cv_pos,3*ncvatoms,MPI_REAL8,
     $     MPI_SUM,0,COMM_TINKER,ierr)
        call MPI_REDUCE(decv_tot,decv_tot,3*ncvatoms,MPI_REAL8,
     $     MPI_SUM,0,COMM_TINKER,ierr)
      end if
      if (use_lambdadyn) then
        if (rank.eq.0) then
          call MPI_REDUCE(MPI_IN_PLACE,delambdav,1,MPI_REAL8,MPI_SUM,
     $   0,COMM_TINKER,ierr)
          call MPI_REDUCE(MPI_IN_PLACE,delambdae,1,MPI_REAL8,MPI_SUM,
     $   0,COMM_TINKER,ierr)
        else
          call MPI_REDUCE(delambdav,delambdav,1,MPI_REAL8,MPI_SUM,
     $   0,COMM_TINKER,ierr)
          call MPI_REDUCE(delambdae,delambdae,1,MPI_REAL8,MPI_SUM,
     $   0,COMM_TINKER,ierr)
        end if
      end if
      if (use_osrw) then
        if (rank.eq.0) then
          call MPI_REDUCE(MPI_IN_PLACE,d2edlambdav2,1,MPI_REAL8,MPI_SUM,
     $   0,COMM_TINKER,ierr)
          call MPI_REDUCE(MPI_IN_PLACE,d2edlambdae2,1,MPI_REAL8,MPI_SUM,
     $   0,COMM_TINKER,ierr)
        else
          call MPI_REDUCE(d2edlambdav2,d2edlambdav2,1,MPI_REAL8,MPI_SUM,
     $   0,COMM_TINKER,ierr)
          call MPI_REDUCE(d2edlambdae2,d2edlambdae2,1,MPI_REAL8,MPI_SUM,
     $   0,COMM_TINKER,ierr)
        end if
      end if
      end subroutine
c
c     the master sends the new forces, he already has the bias
c
      subroutine distrib_colvars(derivs)
      use colvars
      use deriv
      use domdec
      use mpi
      use mutant
      use potent
      implicit none
      real*8 derivs(3,*)
      integer i,iglob,iloc,ierr
c
      call MPI_BCAST(decv,3*ncvatoms,MPI_REAL8,0,COMM_TINKER,ierr)

      if (use_lambdadyn) then
        call MPI_BCAST(lambda,1,MPI_REAL8,0,COMM_TINKER,ierr)
        call def_lambdadyn
        if (use_osrw) then
          call MPI_BCAST(flambdabias,1,MPI_REAL8,0,COMM_TINKER,ierr)
          do i  = 1, nloc
            dxdelambda(:,i) = dxdelambdae(:,i)*dlambdaelambda +
     $         dxdelambdav(:,i)*dlambdavlambda
            derivs(:,i) = derivs(:,i) + flambdabias*dxdelambda(:,i)
          end do
        end if
      end if
c
      do i = 1, ncvatoms
        iglob = cvatoms_ids(i)
        if (repart(iglob).eq.rank) then
          iloc = loc(iglob)
          derivs(1,iloc) = derivs(1,iloc) + decv(1,i)
          derivs(2,iloc) = derivs(2,iloc) + decv(2,i)
          derivs(3,iloc) = derivs(3,iloc) + decv(3,i)
        end if
      end do
      return
      end subroutine
