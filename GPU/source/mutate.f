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
#include "tinker_precision.h"
      subroutine mutate
      use atmtyp
      use atoms
      use charge
      use domdec
      use keys
      use inform
      use iounit
      use katoms
      use mutant
      use mpole
      use potent
      use tinheader ,only:ti_p,re_p
      implicit none
      integer i,j,k,ihyb
      integer it0,it1,next
      integer list(20)
      integer ierr
      character*20 keyword
      character*120 record
      character*120 string
c
c     allocate arrays
c
      if (rank.eq.0.and.tinkerdebug) print*,'mutate'
      call alloc_shared_mutate
c
c     set defaults for lambda and soft core vdw parameters
c
      lambda = 1.0_ti_p
      vlambda = 1.0_ti_p
      elambda = 1.0_ti_p
      scexp = 5.0_ti_p
      scalpha = 0.7_ti_p
c
c     zero number of hybrid atoms and hybrid atom list
c
      nmut = 0
      do i = 1, n
         mut(i) = .false.
      end do
      do i = 1, 20
         list(i) = 0
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
            string = record(next:120)
            read (string,*,err=20)  lambda
         else if (keyword(1:11) .eq. 'VDW-LAMBDA ') then
            string = record(next:120)
            read (string,*,err=20)  vlambda
         else if (keyword(1:11) .eq. 'ELE-LAMBDA ') then
            string = record(next:120)
            read (string,*,err=20)  elambda
         else if (keyword(1:7) .eq. 'MUTATE ') then
            string = record(next:120)
            read (string,*,err=20)  ihyb,it0,it1
            nmut = nmut + 1
            imut(nmut) = ihyb
            mut(ihyb) = .true.
            type0(nmut) = it0
            type1(nmut) = it1
            class0(nmut) = atmcls(it0)
            class1(nmut) = atmcls(it1)
         else if (keyword(1:7) .eq. 'LIGAND ') then
            string = record(next:120)
            read (string,*,err=10,end=10)  (list(k),k=1,20)
   10       continue
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
         end if
   20    continue
      end do
      call MPI_BARRIER(hostcomm,ierr)
c
c     scale electrostatic parameter values based on lambda
c
      if (hostrank.eq.0) then
        if (elambda.ge.0.0_ti_p .and. elambda.lt.1.0_ti_p)  call altelec
      end if
      call MPI_BARRIER(hostcomm,ierr)
c update from altelec call
      call upload_device_mutate
c
c     write the status of the current free energy perturbation step
c
      if (nmut.ne.0 .and. .not.silent .and. rank .eq. 0) then
         write (iout,30)  vlambda
   30    format (/,' Free Energy Perturbation :',f15.3,
     &              ' Lambda for van der Waals')
         write (iout,40)  elambda
   40    format (' Free Energy Perturbation :',f15.3,
     &              ' Lambda for Electrostatics')
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
            if (mut(i)) then
               pchg(i) = pchg(i) * elambda
            end if
         end do
      end if
c
c     set electrostatic parameters for polarizable multipole models
c
      if (use_mpole.or.use_polar) then
         do i = 1, npole
            k = ipole(i)
            if (mut(k)) then
               do j = 1, 13
                  pole(j,i) = pole(j,i) * elambda
               end do
            end if
         end do
      end if

      if (use_polar) then
         do i = 1, npole
            k = ipole(i) 
            if (mut(k)) then
               polarity(i) = polarity(i) * elambda
            end if
         end do
      end if
      end

      subroutine upload_device_mutate
      use charge
      use domdec,only: rank
      use inform,only: deb_Path
      use mpole
      use mutant
      use polar
      use potent
      implicit none
#ifdef _OPENACC
 12   format(2x,'upload_device_mutate')
      if(deb_Path) print 12
#endif

      if (use_vdw) then
!$acc update device(mut) async
      end if
      if (use_charge) then
!$acc update device(pchg) async
      end if
      if (use_mpole.or.use_polar) then
!$acc update device(pole) async
      end if
      if (use_polar) then
!$acc update device(polarity) async
      end if
      end subroutine

      subroutine delete_data_mutate
      use domdec,only:rank
      use mutant
      use tinMemory
      use sizes ,only:tinkerdebug
      implicit none

 12   format(2x,'delete_data_mutate')
      if (rank.eq.0.and.tinkerdebug) print 12

      call shmem_request(imut,  winimut,  [0])
      call shmem_request(mut,   winmut,   [0], config=mhostacc)
      call shmem_request(type0, wintype0, [0])
      call shmem_request(type1, wintype1, [0])
      call shmem_request(class0,winclass0,[0])
      call shmem_request(class1,winclass1,[0])
      end subroutine
c
c     subroutine alloc_shared_mutate : allocate shared memory pointers for mutate
c     parameter arrays
c
      subroutine alloc_shared_mutate
      use atoms
      use domdec
      use mutant
      use mpi
      use tinMemory
      implicit none
c
      if (associated(imut).and.n.eq.size(imut)) return

      call shmem_request(imut  ,winimut,  [n])
      call shmem_request(mut   ,winmut,   [n], config=mhostacc)
      call shmem_request(type0 ,wintype0, [n])
      call shmem_request(type1 ,wintype1, [n])
      call shmem_request(class0,winclass0,[n])
      call shmem_request(class1,winclass1,[n])
      end
