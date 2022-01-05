c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine mutate  --  set parameters for hybrid system  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "mutate" constructs the hybrid hamiltonian for a specified
c     initial state, final state and mutation parameter "lambda"
c
c
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
c
c     deallocate global pointers if necessary
c
      call dealloc_shared_mutate
c
c     allocate global pointers
c
      call alloc_shared_mutate
c
      scexp = 5.0d0
      scalpha = 0.7d0
      if (hostrank.eq.0) then
c
c     set defaults for lambda and soft core vdw parameters
c
      lambda = 1.0d0
      tlambda = 1.0d0
      vlambda = 1.0d0
      elambda = 1.0d0
      vcouple = 0
c
c     perform dynamic allocation of some local arrays
c
      allocate (itbnd(2,nbond))
c
c     zero number of hybrid atoms and hybrid atom list
c
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
c
c     search keywords for free energy perturbation options
c
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
   10         continue
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
                       imut(nmut) = i
                       mut(j) = .true.
                       type0(nmut) = 0
                       type1(nmut) = type(i)
                       class0(nmut) = 0
                       class1(nmut) = class(i)
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
   20         continue
              k = 1
              do while (list(k) .ne. 0)
                 ntbnd = ntbnd + 1
                 itbnd(1,ntbnd) = list(k)
                 itbnd(2,ntbnd) = list(k+1)
                 k = k + 2
              end do
             end if
   30      continue
        end do
      end if
      call MPI_BCAST(nmut,1,MPI_INT,0,hostcomm,ierr)
      call MPI_BCAST(lambda,1,MPI_REAL8,0,hostcomm,ierr)
      call MPI_BCAST(elambda,1,MPI_REAL8,0,hostcomm,ierr)
      call MPI_BCAST(vlambda,1,MPI_REAL8,0,hostcomm,ierr)
      call MPI_BCAST(tlambda,1,MPI_REAL8,0,hostcomm,ierr)
      call MPI_BCAST(vcouple,1,MPI_INT,0,hostcomm,ierr)
c
c     scale electrostatic parameter values based on lambda
c
      if (hostrank.eq.0) then
        if (elambda.ge.0.0d0 .and. elambda.lt.1.0d0)  call altelec
        if (tlambda.ge.0.0d0 .and. tlambda.lt.1.0d0) then
          if (ntbnd.ne.0) call alttors(ntbnd,itbnd)
        end if
      end if
c
c     write the status of the current free energy perturbation step
c
      if (nmut.ne.0 .and. .not.silent .and. rank .eq. 0) then
         write (iout,40)  tlambda
   40    format (/,' Free Energy Perturbation :',f15.3,
     &              ' Lambda for Torsional Angles')
         write (iout,50)  vlambda
   50    format (/,' Free Energy Perturbation :',f15.3,
     &              ' Lambda for van der Waals')
         write (iout,60)  elambda
   60    format (' Free Energy Perturbation :',f15.3,
     &              ' Lambda for Electrostatics')
      end if
c
c     perform deallocation of some local arrays
c
      if (hostrank.eq.0) deallocate (itbnd)
      return
      end
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine alttors  --  mutated torsional parameters  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "alttors" constructs mutated torsional parameters based
c     on the lambda mutation parameter "tlambda"
c
c
      subroutine alttors (ntbnd,itbnd)
      use mutant
      use potent
      use tors
      implicit none
      integer i,j
      integer ia,ib,ic,id
      integer kb,kc
      integer ntbnd
      integer itbnd(2,*)
c
c
c     set torsional parameters across freely rotatable bonds
c
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
                  if ((kb.eq.ib .and. kc.eq.ic) .or.
     &                (kb.eq.ic .and. kc.eq.ib)) then
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
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine altelec  --  mutated electrostatic parameters  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "altelec" constructs the mutated electrostatic parameters
c     based on the lambda mutation parameter "elmd"
c
c
      subroutine altelec
      use sizes
      use charge
      use domdec
      use mpole
      use mutant
      use polar
      use potent
      implicit none
      integer i,j,k
c
c
c     set electrostatic parameters for partial charge models
c
      if (use_charge) then
         do i = 1, nion
            k = iion(i)
            if (mut(k)) then
               pchg(i) = pchg(i) * elambda
            end if
         end do
      end if
c
c     set electrostatic parameters for polarizable multipole models
c
      if (use_mpole) then
         do i = 1, npole
            k = ipole(i)
            if (mut(k)) then
               do j = 1, 13
                  pole(j,i) = pole(j,i) * elambda
               end do
            end if
         end do
         do i = 1, npole
            k = ipole(i)
            if (mut(k)) then
               polarity(i) = polarity(i) * elambda
            end if
         end do
      end if
      return
      end
c
c     subroutine dealloc_shared_mutate : deallocate shared memory pointers for mutate
c     parameter arrays
c
      subroutine dealloc_shared_mutate
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER
      use mutant
      use mpi
      implicit none
      INTEGER(KIND=MPI_ADDRESS_KIND) :: windowsize
      INTEGER :: disp_unit,ierr
      TYPE(C_PTR) :: baseptr
c
      if (associated(imut)) then
        CALL MPI_Win_shared_query(winimut, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winimut,ierr)
      end if
      if (associated(mut)) then
        CALL MPI_Win_shared_query(winmut, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winmut,ierr)
      end if
      if (associated(type0)) then
        CALL MPI_Win_shared_query(wintype0, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(wintype0,ierr)
      end if
      if (associated(type1)) then
        CALL MPI_Win_shared_query(wintype1, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(wintype1,ierr)
      end if
      if (associated(class0)) then
        CALL MPI_Win_shared_query(winclass0, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winclass0,ierr)
      end if
      if (associated(class1)) then
        CALL MPI_Win_shared_query(winclass1, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winclass1,ierr)
      end if
      return
      end
c
c     subroutine alloc_shared_mutate : allocate shared memory pointers for mutate
c     parameter arrays
c
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
c
c     imut
c
      arrayshape=(/n/)
      if (hostrank == 0) then
        windowsize = int(n,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winimut, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winimut, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,imut,arrayshape)
c
c     mut
c
      arrayshape=(/n/)
      if (hostrank == 0) then
        windowsize = int(n,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winmut, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winmut, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,mut,arrayshape)
c
c     type0
c
      arrayshape=(/n/)
      if (hostrank == 0) then
        windowsize = int(n,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, wintype0, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(wintype0, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,type0,arrayshape)
c
c     type1
c
      arrayshape=(/n/)
      if (hostrank == 0) then
        windowsize = int(n,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, wintype1, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(wintype1, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,type1,arrayshape)
c
c     class0
c
      arrayshape=(/n/)
      if (hostrank == 0) then
        windowsize = int(n,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winclass0, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winclass0, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,class0,arrayshape)
c
c     class1
c
      arrayshape=(/n/)
      if (hostrank == 0) then
        windowsize = int(n,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winclass1, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winclass1, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,class1,arrayshape)
      return
      end
