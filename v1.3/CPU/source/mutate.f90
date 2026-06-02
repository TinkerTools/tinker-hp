!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###############################################################
!     ##                                                           ##
!     ##  subroutine mutate  --  set parameters for hybrid system  ##
!     ##                                                           ##
!     ###############################################################
!
!
!     "mutate" constructs the hybrid hamiltonian for a specified
!     initial state, final state and mutation parameter "lambda"
!
!
!> @brief 
!> constructs the hybrid hamiltonian for a specified
!> initial state, final state and mutation parameter "lambda"
!> @param no params
subroutine mutate
   use atmtyp
   use atoms
   use bond
   use domdec
   use keys
   use inform
   use iounit
   use katoms
   use mutant
   use mpi
   use potent
   implicit none
   integer i,j,k,ihyb
   integer it0,it1,next
   integer ntbnd
   integer list(40)
   integer ierr
   integer, allocatable :: itbnd(:,:)
   character*20 keyword
   character*240 record
   character*240 string
!
   if (deb_Path) write(iout,*), 'mutate '
!
!     deallocate global pointers if necessary
!
   call dealloc_shared_mutate
!
!     allocate global pointers
!
   call alloc_shared_mutate
!
   scexp = 2.0d0
   scalpha = 0.3d0
   if (hostrank.eq.0) then
!
!     set defaults for lambda and soft core vdw parameters
!
      lambda = 1.0d0
      tlambda = 1.0d0
      vlambda = 1.0d0
      elambda = 1.0d0
      sck = 6.0d0
      sct = 1.0d0
      scs = 2.0d0
      vcouple = 0
!
!     perform dynamic allocation of some local arrays
!
      allocate (itbnd(2,nbond))
!
!     zero number of hybrid atoms and hybrid atom list
!
      nmut = 0
      do i = 1, n
         mut(i) = .false.
      end do
      ntbnd = 0
      do i = 1, 40
         list(i) = 0
      end do
      do i = 1, nbond
         itbnd(1,i) = 0
         itbnd(2,i) = 0
      end do
!
!     search keywords for free energy perturbation options
!
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:7) .eq. 'LAMBDA ') then
            string = record(next:240)
            read (string,*,err=30)  lambda
         else if (keyword(1:12) .eq. 'TORS-LAMBDA ') then
            string = record(next:240)
            read (string,*,err=30)  tlambda
         else if (keyword(1:11) .eq. 'VDW-LAMBDA ') then
            string = record(next:240)
            read (string,*,err=30)  vlambda
         else if (keyword(1:11) .eq. 'ELE-LAMBDA ') then
            string = record(next:240)
            read (string,*,err=30)  elambda
         else if (keyword(1:17) .eq. 'BOUND-VDW-LAMBDA ') then
            string = record(next:240)
            read (string,*,err=30)  bvlambda
            if (rank.eq.0) write(iout,35) bvlambda
35          format('Intervall bound for lambda_vdw is  ', F15.3)
         else if (keyword(1:17) .eq. 'BOUND-ELE-LAMBDA ') then
            string = record(next:240)
            read (string,*,err=30)  belambda
            if (rank.eq.0) write(iout,45) belambda
45          format('Intervall bound for lambda_elec is  ', F15.3)
         else if (keyword(1:11) .eq. 'VDW-SC-EXP ') then
            string = record(next:240)
            read (string,*,err=30)  scexp
         else if (keyword(1:13) .eq. 'VDW-SC-ALPHA ') then
            string = record(next:240)
            read (string,*,err=30)  scalpha
         else if (keyword(1:9) .eq. 'VDW-SC-K ') then
            string = record(next:240)
            read (string,*,err=30)  sck
         else if (keyword(1:9) .eq. 'VDW-SC-T ') then
            string = record(next:240)
            read (string,*,err=30)  sct
         else if (keyword(1:9) .eq. 'VDW-SC-S ') then
            string = record(next:240)
            read (string,*,err=30)  scs
!           else if (keyword(1:10) .eq. 'LAMBDADYN ') then
!              use_lambdadyn = .true.
         else if (keyword(1:15) .eq. 'VDW-ANNIHILATE ') then
            vcouple = 1
         else if (keyword(1:7) .eq. 'MUTATE ') then
            string = record(next:240)
            read (string,*,err=30)  ihyb,it0,it1
            nmut = nmut + 1
            imut(nmut) = ihyb
            mut(ihyb) = .true.
            type0(nmut) = it0
            type1(nmut) = it1
            class0(nmut) = atmcls(it0)
            class1(nmut) = atmcls(it1)
         else if (keyword(1:7) .eq. 'LIGAND ') then
            string = record(next:240)
            read (string,*,err=10,end=10)  (list(k),k=1,20)
10          continue
            k = 1
            do while (list(k) .ne. 0)
               if (list(k) .gt. 0) then
                  j = list(k)
                  nmut = nmut + 1
                  imut(nmut) = j
                  mut(j) = .true.
                  type0(nmut) = 0
                  type1(nmut) = type(j)
                  class0(nmut) = 0
                  class1(nmut) = class(j)
                  k = k + 1
               else
                  do j = abs(list(k)), abs(list(k+1))
                     nmut = nmut + 1
                     imut(nmut) = j
                     mut(j) = .true.
                     type0(nmut) = 0
                     type1(nmut) = type(j)
                     class0(nmut) = 0
                     class1(nmut) = class(j)
                  end do
                  k = k + 2
               end if
            end do
         else if (keyword(1:15) .eq. 'ROTATABLE-BOND ') then
            do k = 1, 40
               list(k) = 0
            end do
            string = record(next:240)
            read (string,*,err=20,end=20)  (list(k),k=1,40)
20          continue
            k = 1
            do while (list(k) .ne. 0)
               ntbnd = ntbnd + 1
               itbnd(1,ntbnd) = list(k)
               itbnd(2,ntbnd) = list(k+1)
               k = k + 2
            end do
         end if
30       continue
      end do
   end if
   call MPI_BCAST(nmut,1,MPI_INT,0,hostcomm,ierr)
   call MPI_BCAST(lambda,1,MPI_REAL8,0,hostcomm,ierr)
   call MPI_BCAST(elambda,1,MPI_REAL8,0,hostcomm,ierr)
   call MPI_BCAST(vlambda,1,MPI_REAL8,0,hostcomm,ierr)
   call MPI_BCAST(tlambda,1,MPI_REAL8,0,hostcomm,ierr)
   call MPI_BCAST(vcouple,1,MPI_INT,0,hostcomm,ierr)
   call MPI_BCAST(scexp,1,MPI_REAL8,0,hostcomm,ierr)
   call MPI_BCAST(scalpha,1,MPI_REAL8,0,hostcomm,ierr)
   call MPI_BCAST(sck,1,MPI_REAL8,0,hostcomm,ierr)
   call MPI_BCAST(sct,1,MPI_REAL8,0,hostcomm,ierr)
   call MPI_BCAST(scs,1,MPI_REAL8,0,hostcomm,ierr)
!
!     scale electrostatic parameter values based on lambda
!
   if (hostrank.eq.0) then
      if (elambda.ge.0.0d0 .and. elambda.lt.1.0d0)  call altelec
      if (tlambda.ge.0.0d0 .and. tlambda.lt.1.0d0) then
         if (ntbnd.ne.0) call alttors(ntbnd,itbnd)
      end if
   end if
!
!     write the status of the current free energy perturbation step
!
   if (nmut.ne.0 .and. .not.silent .and. rank .eq. 0) then
      write (iout,40)  tlambda
40    format (/,' Free Energy Perturbation :',f15.3,&
      &' Lambda for Torsional Angles')
      write (iout,50)  vlambda
50    format (/,' Free Energy Perturbation :',f15.3,&
      &' Lambda for van der Waals')
      write (iout,60)  elambda
60    format (' Free Energy Perturbation :',f15.3,&
      &' Lambda for Electrostatics')
   end if
!
!     perform deallocation of some local arrays
!
   if (hostrank.eq.0) deallocate (itbnd)
   return
end
!
!     ############################################################
!     ##                                                        ##
!     ##  subroutine alttors  --  mutated torsional parameters  ##
!     ##                                                        ##
!     ############################################################
!
!
!     "alttors" constructs mutated torsional parameters based
!     on the lambda mutation parameter "tlambda"
!
!
!> @brief 
!> constructs mutated torsional parameters based
!> on the lambda mutation parameter "tlambda"
!> @param[in] ntbnd: number of rotatable bonds
!> @param[in] itbnd: rotatable bonds
subroutine alttors (ntbnd,itbnd)
   use inform
   use iounit
   use mutant
   use potent
   use tors
   implicit none
   integer i,j
   integer ia,ib,ic,id
   integer kb,kc
   integer ntbnd
   integer itbnd(2,*)
!
   if (deb_Path) write(iout,*), 'alttors '
!
!     set torsional parameters across freely rotatable bonds
!
   if (use_tors) then
      do i = 1, ntors
         ia = itors(1,i)
         ib = itors(2,i)
         ic = itors(3,i)
         id = itors(4,i)
         if (mut(ia) .and. mut(ib) .and. mut(ic) .and. mut(id)) then
            do j = 1, ntbnd
               kb = itbnd(1,j)
               kc = itbnd(2,j)
               if ((kb.eq.ib .and. kc.eq.ic) .or.&
               &(kb.eq.ic .and. kc.eq.ib)) then
                  tors1(1,i) = tors1(1,i) * tlambda
                  tors2(1,i) = tors2(1,i) * tlambda
                  tors3(1,i) = tors3(1,i) * tlambda
                  tors4(1,i) = tors4(1,i) * tlambda
                  tors5(1,i) = tors5(1,i) * tlambda
                  tors6(1,i) = tors6(1,i) * tlambda
               end if
            end do
         end if
      end do
   end if
   return
end
!
!
!     ################################################################
!     ##                                                            ##
!     ##  subroutine altelec  --  mutated electrostatic parameters  ##
!     ##                                                            ##
!     ################################################################
!
!
!     "altelec" constructs the mutated electrostatic parameters
!     based on the lambda mutation parameter "elmd"
!
!
!> @brief 
!> "altelec" constructs the mutated electrostatic parameters
!> based on the lambda mutation parameter "elmd"
!> @param no params
subroutine altelec
   use angle
   use bond
   use charge
   use chgpen
   use chgtrn
   use cflux
   use domdec
   use inform
   use iounit
   use mpole
   use mutant
   use polar
   use potent
   implicit none
   integer i,j,k
   integer ia,ib,ic
!
   if (deb_Path) write(iout,*), 'altelec '
!
!
!
!     set electrostatic parameters for partial charge models
!
   if (use_charge) then
      do i = 1, nion
         k = iion(i)
         if (mut(k)) then
!               pchg(i) = pchg(i) * elambda
            pchg(i) = pchg_orig(i) * elambda
         end if
         pchg0(i) = pchg(i)
      end do
   end if
!
!     set electrostatic parameters for polarizable multipole models
!
   if (use_mpole.or.use_polar) then
      do i = 1, npole
         k = ipole(i)
         if (mut(k)) then
            do j = 1, 13
!                  pole(j,i) = pole(j,i) * elambda
               pole(j,i) = pole_orig(j,i) * elambda
            end do
            mono0(i) = pole(1,i)
            if (use_chgpen) then
               pcore(i) = pcore(i) * elambda
               pval(i) = pval(i) * elambda
               pval0(i) = pval(i)
            end if
         end if
      end do
   end if

   if (use_polar) then
      do i = 1, npole
         k = ipole(i)
         if (mut(k)) then
!               polarity(i) = polarity(i) * elambda
            polarity(i) = polarity_orig(i) * elambda
         end if
      end do
   end if
!
!     set scaled parameters for charge transfer models
!
   if (use_chgtrn) then
      do i = 1, npole
         k = ipole(i)
         if (mut(k)) then
            chgct(i) = chgct(i) * elambda
         end if
      end do
   end if
!
!     set scaled parameters for bond stretch charge flux
!
   if (use_chgflx) then
      do i = 1, nbond
         ia = ibnd(1,i)
         ib = ibnd(2,i)
         if (mut(ia) .and. mut(ib)) then
            bflx(i) = bflx(i) * elambda
         end if
      end do
   end if
!
!     set scaled parameters for angle bend charge flux
!
   if (use_chgflx) then
      do i = 1, nangle
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         if (mut(ia) .and. mut(ib) .and. mut(ic)) then
            aflx(1,i) = aflx(1,i) * elambda
            aflx(2,i) = aflx(2,i) * elambda
            abflx(1,i) = abflx(1,i) * elambda
            abflx(2,i) = abflx(2,i) * elambda
         end if
      end do
   end if

   return
end
!
!     subroutine def_lambdadyn_init: lambda dynamics initialization
!
!> @brief 
!> lambda dynamics initialization
!> @param[in] rank_: MPI rank
subroutine def_lambdadyn_init(rank_)
   use inform
   use iounit
   use keys
   use mutant
   implicit none
   integer,intent(in) :: rank_
   integer i,next
   character*20 keyword
   character*240 record
   character*240 string
!
   if (deb_Path) write(iout,*), 'def_lambdadyn_init '
!
   bvlambda = 0.5d0
   belambda = 0.5d0
!
!     search keywords for lambda dynamics options
!
   do i = 1, nkey
      next = 1
      record = keyline(i)
      call gettext (record,keyword,next)
      call upcase (keyword)
      if (keyword(1:17) .eq. 'BOUND-VDW-LAMBDA ') then
         string = record(next:240)
         read (string,*,err=20)  bvlambda
         write(iout,30) bvlambda
30       format('Intervall bound for lambda_vdw is  ', F15.3)
      else if (keyword(1:17) .eq. 'BOUND-ELE-LAMBDA ') then
         string = record(next:240)
         read (string,*,err=20)  belambda
         write(iout,40) belambda
40       format('Intervall bound for lambda_elec is  ', F15.3)
      end if
20    continue
   end do
   call def_lambdadyn(rank_)
   return
end

!
!     #############################################################################
!     ##                                                                         ##
!     ##  subroutine def_lambdadyn -- state weighting values for van der Waals   ##
!     ##                              and electrostatic interactions             ##
!     ##                                                                         ##
!     #############################################################################
!
!     "def_lambdadyn" constructs the state weighting values vlambda and elambda for
!     the van der Waals and electrostatic interactions respectively, defining them
!     as functions of the generic state weighting value lambda
!
!
!> @brief 
!> constructs the state weighting values vlambda and elambda for
!> the van der Waals and electrostatic interactions respectively, defining them
!> as functions of the generic state weighting value lambda
!> @param[in] rank_: MPI rank
subroutine def_lambdadyn(rank_)
   use atmtyp
   use atoms
   use domdec
   use deriv
   use keys
   use inform
   use iounit
   use mutant
   use mpi
   use potent
   implicit none
   integer,intent(in) :: rank_
   integer ierr
!
   if (deb_Path) write(iout,*), 'def_lambdadyn '
!

!     checks if the intervall bounds for vlambda and elambda are consistent
   if (rank_.eq.0) then
      if ( bvlambda .le. 0.0d0 .OR. bvlambda .gt. 1.0d0 ) then
         write(iout,*) 'Intervall bound for vlambda must be between 0 ',&
         &'and 1'
         call fatal
      else if ( belambda .lt. 0.0d0 .OR. belambda .ge. 1.0d0 ) then
         write(iout,*) 'Intervall bound for elambda must be between 0 ',&
         &'and 1'
         call fatal
      end if

      if (bvlambda.eq.0.0d0) then
10       format('Forbidden value of bvlambda: ',F15.3)
         call fatal
      end if
      if (belambda.eq.1.0d0) then
15       format('Forbidden value of belambda: ',F15.3)
         call fatal
      end if

      if (bvlambda .lt. belambda ) then
         write(iout,*) 'Intervall bound for vlambda cannot be ',&
         &'less than the intervall bound for elambda'
         call fatal
      end if
   end if

   if (lambda.lt.0.0d0) then
      vlambda = 0.0d0
      dlambdavlambda = 0.0d0
      elambda = 0.0d0
      dlambdaelambda = 0.0d0
   else if (lambda.lt.belambda) then
      vlambda = lambda/bvlambda
      dlambdavlambda = 1.0d0 / bvlambda
      elambda = 0.0d0
      dlambdaelambda = 0.0d0
   else if (lambda.lt.bvlambda) then
      vlambda = lambda / bvlambda
      dlambdavlambda = 1.0d0/bvlambda
      elambda = (1.0d0/(1-belambda))*lambda - belambda/(1-belambda)
      dlambdaelambda = 1.0d0/(1-belambda)
   else if (lambda.lt.1.0d0) then
      vlambda = 1.0d0
      dlambdavlambda = 0.0d0
      elambda = (1.0d0/(1-belambda)) * lambda - belambda/(1-belambda)
      dlambdaelambda = 1.0d0/(1-belambda)
   else
      vlambda = 1.0d0
      dlambdavlambda = 0.0d0
      elambda = 1.0d0
      dlambdaelambda = 0.0d0
   end if

   if ((rank_.eq.0).and.(verbose)) then
      write(iout,20) vlambda, dlambdavlambda
20    format('Value of vlambda is ', F15.3,&
      &' Value of dlambdavlambda is ',F15.3)
      write(iout,25) elambda, dlambdaelambda
25    format('Value of elambda is ', F15.3,&
      &' Value of dlambdaelambda is ',F15.3)
   end if
   if (hostrank.eq.0) call altelec
   call MPI_BARRIER(hostcomm,ierr)
!     Initialising double derivatives in case of OSRW
   if (use_OSRW) then
      d2edlambda2 = 0.0d0
      d2edlambdae2 = 0.0d0
      d2edlambdav2 = 0.0d0
   end if
   return
end
!
!     subroutine dealloc_shared_mutate : deallocate shared memory pointers for mutate
!     parameter arrays
!
!> @brief 
!> deallocate shared memory pointers for mutate
!> parameter arrays
!> @param no params
subroutine dealloc_shared_mutate
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER
   use mutant
   use mpi
   implicit none
   INTEGER(KIND=MPI_ADDRESS_KIND) :: windowsize
   INTEGER :: disp_unit,ierr
   TYPE(C_PTR) :: baseptr
!
   if (associated(imut)) then
      CALL MPI_Win_shared_query(winimut, 0, windowsize, disp_unit,&
      &baseptr, ierr)
      CALL MPI_Win_free(winimut,ierr)
   end if
   if (associated(mut)) then
      CALL MPI_Win_shared_query(winmut, 0, windowsize, disp_unit,&
      &baseptr, ierr)
      CALL MPI_Win_free(winmut,ierr)
   end if
   if (associated(type0)) then
      CALL MPI_Win_shared_query(wintype0, 0, windowsize, disp_unit,&
      &baseptr, ierr)
      CALL MPI_Win_free(wintype0,ierr)
   end if
   if (associated(type1)) then
      CALL MPI_Win_shared_query(wintype1, 0, windowsize, disp_unit,&
      &baseptr, ierr)
      CALL MPI_Win_free(wintype1,ierr)
   end if
   if (associated(class0)) then
      CALL MPI_Win_shared_query(winclass0, 0, windowsize, disp_unit,&
      &baseptr, ierr)
      CALL MPI_Win_free(winclass0,ierr)
   end if
   if (associated(class1)) then
      CALL MPI_Win_shared_query(winclass1, 0, windowsize, disp_unit,&
      &baseptr, ierr)
      CALL MPI_Win_free(winclass1,ierr)
   end if
   return
end
!
!     subroutine alloc_shared_mutate : allocate shared memory pointers for mutate
!     parameter arrays
!
!> @brief 
!> allocate shared memory pointers for mutate
!> parameter arrays
!> @param no params
subroutine alloc_shared_mutate
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER
   use atoms
   use domdec
   use mutant
   use mpi
   implicit none
   INTEGER(KIND=MPI_ADDRESS_KIND) :: windowsize
   INTEGER :: disp_unit,ierr
   TYPE(C_PTR) :: baseptr
   integer :: arrayshape(1)
!
!     imut
!
   arrayshape=(/n/)
   if (hostrank == 0) then
      windowsize = int(n,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
   else
      windowsize = 0_MPI_ADDRESS_KIND
   end if
   disp_unit = 1
!
!    allocation
!
   CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,&
   &hostcomm, baseptr, winimut, ierr)
   if (hostrank /= 0) then
      CALL MPI_Win_shared_query(winimut, 0, windowsize, disp_unit,&
      &baseptr, ierr)
   end if
!
!    association with fortran pointer
!
   CALL C_F_POINTER(baseptr,imut,arrayshape)
!
!     mut
!
   arrayshape=(/n/)
   if (hostrank == 0) then
      windowsize = int(n,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
   else
      windowsize = 0_MPI_ADDRESS_KIND
   end if
   disp_unit = 1
!
!    allocation
!
   CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,&
   &hostcomm, baseptr, winmut, ierr)
   if (hostrank /= 0) then
      CALL MPI_Win_shared_query(winmut, 0, windowsize, disp_unit,&
      &baseptr, ierr)
   end if
!
!    association with fortran pointer
!
   CALL C_F_POINTER(baseptr,mut,arrayshape)
!
!     type0
!
   arrayshape=(/n/)
   if (hostrank == 0) then
      windowsize = int(n,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
   else
      windowsize = 0_MPI_ADDRESS_KIND
   end if
   disp_unit = 1
!
!    allocation
!
   CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,&
   &hostcomm, baseptr, wintype0, ierr)
   if (hostrank /= 0) then
      CALL MPI_Win_shared_query(wintype0, 0, windowsize, disp_unit,&
      &baseptr, ierr)
   end if
!
!    association with fortran pointer
!
   CALL C_F_POINTER(baseptr,type0,arrayshape)
!
!     type1
!
   arrayshape=(/n/)
   if (hostrank == 0) then
      windowsize = int(n,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
   else
      windowsize = 0_MPI_ADDRESS_KIND
   end if
   disp_unit = 1
!
!    allocation
!
   CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,&
   &hostcomm, baseptr, wintype1, ierr)
   if (hostrank /= 0) then
      CALL MPI_Win_shared_query(wintype1, 0, windowsize, disp_unit,&
      &baseptr, ierr)
   end if
!
!    association with fortran pointer
!
   CALL C_F_POINTER(baseptr,type1,arrayshape)
!
!     class0
!
   arrayshape=(/n/)
   if (hostrank == 0) then
      windowsize = int(n,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
   else
      windowsize = 0_MPI_ADDRESS_KIND
   end if
   disp_unit = 1
!
!    allocation
!
   CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,&
   &hostcomm, baseptr, winclass0, ierr)
   if (hostrank /= 0) then
      CALL MPI_Win_shared_query(winclass0, 0, windowsize, disp_unit,&
      &baseptr, ierr)
   end if
!
!    association with fortran pointer
!
   CALL C_F_POINTER(baseptr,class0,arrayshape)
!
!     class1
!
   arrayshape=(/n/)
   if (hostrank == 0) then
      windowsize = int(n,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
   else
      windowsize = 0_MPI_ADDRESS_KIND
   end if
   disp_unit = 1
!
!    allocation
!
   CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,&
   &hostcomm, baseptr, winclass1, ierr)
   if (hostrank /= 0) then
      CALL MPI_Win_shared_query(winclass1, 0, windowsize, disp_unit,&
      &baseptr, ierr)
   end if
!
!    association with fortran pointer
!
   CALL C_F_POINTER(baseptr,class1,arrayshape)
   return
end
