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
c     colvarsinput   base filename used by colvars input
c     colvarsoutput   base filename used by colvars output
c     colvarsrestart   base filename used by colvars restart
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
      integer :: ncvatoms,ncvatomsmol
      integer, allocatable :: cvatoms_ids(:),cvatomsmol(:)
      real*8, allocatable :: cv_pos(:,:),decv(:,:),decv_tot(:,:)
      real*8, target :: dt_sim
      real*8, target :: temp_rand
      character*240, target:: colvarsinput,colvarsoutput
      character*240, target :: colvarsrestart
      save
      end

      subroutine set_cvatoms_ids(ncvatoms_in,cvatoms_ids_in)
      use iso_c_binding
      use atoms
      use colvars
      use domdec
      use mpi
      use molcul
      implicit none
      integer :: ncvatoms_in
      integer ierr,i,j,k,l,mol,o,p,init,stop
      integer, dimension(ncvatoms_in) :: cvatoms_ids_in
      ncvatoms = ncvatoms_in
      if (allocated(cvatoms_ids)) deallocate(cvatoms_ids)
      allocate (cvatoms_ids(ncvatoms))
      do i = 1, ncvatoms
        cvatoms_ids(i) = cvatoms_ids_in(i)+1
      end do
c      cvatoms_ids(1:ncvatoms) = cvatoms_ids_in(1:ncvatoms)
      allocate (cv_pos(3,ncvatoms))
      allocate (decv(3,ncvatoms),decv_tot(3,ncvatoms))
      cv_pos = 0d0
      decv = 0d0
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
c      write(*,*) 'ncvatomsmol = ',ncvatomsmol
c      do j = 1, ncvatomsmol
c        write(*,*) 'j = ',j,'cvatomsmol = ',cvatomsmol(j)
c      end do
c
      return
      end subroutine

c
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

c
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
c
c     TODO : gestion multi repliques
c
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
c
c     TODO : gestion multi repliques
c
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
      use deriv
      use iounit
      use inform
      use iso_c_binding
      use mutant
      implicit none
      real(c_double) :: res
      if (verbose) then
        write(iout,5) lambda
 5      format('get_delambda -- value of lambda: ',F15.3)
        write(iout,10) delambdae
 10     format('value of delambdae: ',F15.3)
        write(iout,20) dlambdaelambda
 20     format('value of dlambdaelambda: ',F15.3)
        write(iout,30) delambdav
 30     format('value of delambdav: ',F15.3)
        write(iout,40) dlambdavlambda
 40     format('value of dlambdavlambda: ',F15.3)
      end if
      delambda = delambdae * dlambdaelambda + delambdav * dlambdavlambda
      if (verbose) then
        write(iout,50) delambda
 50     format('value of delambda: ',F15.3)
      end if
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
c     flambdabias is the energy derivative with respect to dE/dlambda
c     which is equal to minus the force on dE/dlambda
      flambdabias = -force
      end subroutine

c
c     the master gets all the cv positions and total forces
c
      subroutine prepare_colvars
      use atoms
      use boxes
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
c
c   get the unwrapped coordinates
c
          cv_pos(1,i) = x(iglob) + pbcwrapindex(1,iglob)*xbox
          cv_pos(2,i) = y(iglob) + pbcwrapindex(2,iglob)*ybox
          cv_pos(3,i) = z(iglob) + pbcwrapindex(3,iglob)*zbox
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
      use atoms
      use colvars
      use deriv
      use domdec
      use mpi
      use mutant
      use potent
      use virial
      implicit none
      real*8 derivs(3,*)
      real*8 vxx,vxy,vxz,vyy,vyz,vzz
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
      return
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
