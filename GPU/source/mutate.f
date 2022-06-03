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
#include "tinker_types.h"
      subroutine mutate
      use atmtyp
      use atoms
      use bond
      use charge
      use domdec
      use keys
      use inform
      use iounit
      use katoms
      use mutant
      use mpi
      use mpole
      use potent
      use tinheader ,only:ti_p,re_p
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
c     allocate arrays
c
      if (deb_Path) write(*,*) 'mutate'
      call alloc_shared_mutate
c
c     set defaults for lambda and soft core vdw parameters
c
      lambda  = 1.0_ti_p
      tlambda = 1.0_ti_p
      vlambda = 1.0_ti_p
      elambda = 1.0_ti_p
      scexp   = 5.0_ti_p
      scalpha = 0.7_ti_p
      sck     = 6.0_ti_p
      sct     = 1.0_ti_p
      scs     = 2.0_ti_p
      vcouple = 0
c
c     perform dynamic allocation of some local arrays
c
      allocate (itbnd(2,nbond))
c
c     zero number of hybrid atoms, hybrid atom list and mutated torsions
c
      nmut = 0
      do i = 1, n
         mut(i) = .false.
         mutInt(i) = 0
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
           else if (keyword(1:17) .eq. 'BOUND-VDW-LAMBDA ') then
              string = record(next:120)
              read (string,*,err=30)  bvlambda
              if (rank.eq.0) write(iout,35) bvlambda
 35           format('Intervall bound for lambda_vdw is  ', F15.3)
           else if (keyword(1:17) .eq. 'BOUND-ELE-LAMBDA ') then
              string = record(next:120)
              read (string,*,err=30)  belambda
              if (rank.eq.0) write(iout,45) belambda
 45           format('Intervall bound for lambda_elec is  ', F15.3)
           else if (keyword(1:17) .eq. 'BOUND-POL-LAMBDA ') then
              string = record(next:120)
              read (string,*,err=30)  bplambda
           else if (keyword(1:11) .eq. 'VDW-SC-EXP ') then
              string = record(next:120)
              read (string,*,err=30)  scexp
           else if (keyword(1:13) .eq. 'VDW-SC-ALPHA ') then
              string = record(next:120)
              read (string,*,err=30)  scalpha
           else if (keyword(1:9) .eq. 'VDW-SC-K ') then
              string = record(next:120)
              read (string,*,err=30)  sck
           else if (keyword(1:9) .eq. 'VDW-SC-T ') then
              string = record(next:120)
              read (string,*,err=30)  sct
           else if (keyword(1:9) .eq. 'VDW-SC-S ') then
              string = record(next:120)
              read (string,*,err=30)  scs
           else if (keyword(1:10) .eq. 'LAMBDADYN ') then
              use_lambdadyn = .true.
         else if (keyword(1:15) .eq. 'VDW-ANNIHILATE ') then
            vcouple = 1
         else if (keyword(1:7) .eq. 'MUTATE ') then
            string = record(next:240)
            read (string,*,err=30)  ihyb,it0,it1
            nmut = nmut + 1
            imut(nmut) = ihyb
            mut(ihyb) = .true.
            mutInt(ihyb) = 1
            type0(nmut) = it0
            type1(nmut) = it1
            class0(nmut) = atmcls(it0)
            class1(nmut) = atmcls(it1)
         else if (keyword(1:7) .eq. 'LIGAND ') then
            string = record(next:240)
            read (string,*,err=10,end=10)  (list(k),k=1,40)
   10       continue
            k = 1
            do while (list(k) .ne. 0)
               if (list(k) .gt. 0) then
                  j = list(k)
                  nmut = nmut + 1
                  imut(nmut) = j
                  mut(j) = .true.
                  mutInt(j) = 1
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
                     mutInt(j) = 1
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
   20       continue
            k = 1
            do while (list(k) .ne. 0)
               ntbnd = ntbnd + 1
               itbnd(1,ntbnd) = list(k)
               itbnd(2,ntbnd) = list(k+1)
               k = k + 2
            end do
         end if
   30    continue
      end do
      if (nproc.gt.1) then
      call MPI_BCAST(nmut,1,MPI_INT,0,hostcomm,ierr)
      call MPI_BCAST(lambda,1,MPI_TPREC,0,hostcomm,ierr)
      call MPI_BCAST(elambda,1,MPI_TPREC,0,hostcomm,ierr)
      call MPI_BCAST(vlambda,1,MPI_TPREC,0,hostcomm,ierr)
      call MPI_BCAST(tlambda,1,MPI_TPREC,0,hostcomm,ierr)
      call MPI_BCAST(vcouple,1,MPI_INT,0,hostcomm,ierr)
      call MPI_BCAST(scexp,1,MPI_TPREC,0,hostcomm,ierr)
      call MPI_BCAST(scalpha,1,MPI_TPREC,0,hostcomm,ierr)
      call MPI_BCAST(sck,1,MPI_TPREC,0,hostcomm,ierr)
      call MPI_BCAST(sct,1,MPI_TPREC,0,hostcomm,ierr)
      call MPI_BCAST(scs,1,MPI_TPREC,0,hostcomm,ierr)
      end if
c
c     scale electrostatic and torsion parameter values based on lambda
c
      if (elambda.ge.0.0_ti_p .and. elambda.lt.1.0_ti_p)  then
         if (use_lambdadyn) then
            call altelec(3)
         else
            call altelec_definitive
         end if
      end if
      if (hostrank.eq.0) then
        if (tlambda.ge.0.0_ti_p .and. tlambda.lt.1.0_ti_p) then
          if (ntbnd.ne.0) call alttors(ntbnd,itbnd)
        end if
      end if
      if (nproc.gt.1) call MPI_BARRIER(hostcomm,ierr)
c update from altelec call
      call upload_device_mutate
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
      deallocate (itbnd)
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
      subroutine altchg(option)
      use angle
      use bond
      use cflux
      use charge
      use chgpen
      use chgtrn
      use domdec
      use mpole
      use mpi
      use mutant
      use polar
      use potent
      use sizes
      implicit none
      integer ,intent(in):: option
      integer i,j,k,ierr
      integer ia,ib,ic
      enum, bind(c)
      enumerator hrank0,transfert
      end enum

      if (btest(option,hrank0).and.hostrank.ne.0) return
c
c     set scaled parameters for charge transfer models
c
      if (use_chgtrn) then
         if (btest(option,transfert)) then
         end if
!$acc parallel loop default(present) async
         do i = 1, npole
            k = ipole(i)
            if (mut(k)) then
               chgct(i) = chgct(i) * elambda
            end if
         end do
         if (btest(option,transfert)) then
!$acc update host(chgct) async
         end if
      end if
c
c     set scaled parameters for bond stretch charge flux
c
      if (use_chgflx) then
         if (btest(option,transfert)) then
         end if
!$acc parallel loop default(present) async
         do i = 1, nbond
            ia = ibnd(1,i)
            ib = ibnd(2,i)
            if (mut(ia) .and. mut(ib)) then
               bflx(i) = bflx(i) * elambda
            end if
         end do
         if (btest(option,transfert)) then
!$acc update host(bflx) async
         end if
      end if
c
c     set scaled parameters for angle bend charge flux
c
      if (use_chgflx) then
         if (btest(option,transfert)) then
         end if
!$acc parallel loop default(present) async
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
         if (btest(option,transfert)) then
!$acc update host(aflx,abflx) async
         end if
      end if
      end subroutine
c
c     altered original state of elec entities
c
      subroutine altelec_definitive
      use sizes
      use charge
      use chgpen
      use chgtrn
      use domdec
      use mpole
      use inform ,only: deb_Path
      use mpi
      use mutant
      use polar
      use potent
      implicit none
      integer i,j,k,ierr
c
      if (hostrank.eq.0) then
         if (deb_Path) write(*,'(3x,A)') "altelec_definitive"
c
c     set electrostatic parameters for partial charge models
c
      if (use_charge) then
         do i = 1, nion
            k = iion(i)
            if (mut(k)) then
               pchg(i) = pchg(i) * elambda
            end if
            pchg0(i) = pchg(i)
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
               polarity(i) = polarity(i) * elambda
            end if
         end do
      end if

      call altchg(3)

      end if
      call MPI_BARRIER(hostcomm,ierr)
      end
c
c     altered a copy of elec entities
c
      subroutine altelec(option)
      use sizes
      use charge
      use chgpen
      use domdec
      use inform ,only: deb_Path
      use mpole
      use mpi
      use mutant
      use polar
      use potent
      implicit none
      integer,intent(in):: option
      integer i,j,k,ierr
      enum, bind(C)
      enumerator hrank0,transfert
      end enum
c
      if (btest(option,hrank0).and.hostrank.ne.0) goto 20
      if (deb_Path) write(*,'(3x,A,I4)') "altelec",option
c
c     set electrostatic parameters for partial charge models
c
      if (use_charge) then
         if (btest(option,transfert)) then
!$acc update device(mut) async
         end if
!$acc parallel loop async default(present)
         do i = 1, nion
            k = iion(i)
            if (mut(k)) then
               pchg(i) = pchg_orig(i) * elambda
            end if
         end do
         if (btest(option,transfert)) then
!$acc update host(pchg) async
         end if
      end if
c
c     set electrostatic parameters for polarizable multipole models
c
      if (use_mpole.or.use_polar) then
         if (btest(option,transfert)) then
!$acc update device(mut) async
         end if
!$acc parallel loop gang vector async
!$acc&         present(pole,pole_orig,mono0,pcore,pval,pval0)
         do i = 1, npole
            k = ipole(i)
            if (mut(k)) then
!$acc loop seq
               do j = 1, 13
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
         if (btest(option,transfert)) then
!$acc update host(pole,mono0,pcore,pval,pval0) async
         end if
      end if

      if (use_polar) then
         if (btest(option,transfert)) then
!$acc update device(mut) async
         end if
!$acc parallel loop gang vector async default(present)
         do i = 1, npole
            k = ipole(i) 
            if (mut(k)) then
               polarity(i) = polarity_orig(i) * elambda
            end if
         end do
         if (btest(option,transfert)) then
!$acc update host(polarity) async
         end if
      end if

      call altchg(option)

 20   continue
      if (nproc.gt.1) call MPI_BARRIER(hostcomm,ierr)
      end
c
c     subroutine def_lambdadyn_init: lambda dynamics initialization
c
      subroutine def_lambdadyn_init
      use iounit
      use keys
      use mutant
      use tinheader
      implicit none
      integer i,next
      character*20 keyword
      character*240 record
      character*240 string
c
      bvlambda = 0.5_ti_p
      belambda = 0.5_ti_p
      bplambda = 0.7_ti_p
c
c     search keywords for lambda dynamics options
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:17) .eq. 'BOUND-VDW-LAMBDA ') then
            string = record(next:120)
            read (string,*,err=20)  bvlambda
            write(iout,30) bvlambda
 30         format('Intervall bound for lambda_vdw is  ', F15.3)
         else if (keyword(1:17) .eq. 'BOUND-ELE-LAMBDA ') then
            string = record(next:120)
            read (string,*,err=20)  belambda
            write(iout,40) belambda
 40         format('Intervall bound for lambda_elec is  ', F15.3)
         else if (keyword(1:17) .eq. 'BOUND-POL-LAMBDA ') then
            string = record(next:120)
            read (string,*,err=20)  bplambda
         end if
   20    continue
      end do
      call def_lambdadyn
      end
c
c     #############################################################################
c     ##                                                                         ##
c     ##  subroutine def_lambdadyn -- state weighting values for van der Waals   ##
c     ##                              and electrostatic interactions             ##
c     ##                                                                         ##
c     #############################################################################
c
c     "def_lambdadyn" constructs the state weighting values vlambda and elambda for
c     the van der Waals and electrostatic interactions respectively, defining them
c     as functions of the generic state weighting value lambda 
c
c
      subroutine def_lambdadyn
      use atmtyp
      use atoms
      use domdec
      use deriv
      use keys
      use inform
      use iounit
      use mutant
      use mpi
      use moldyn
      use potent
      use tinheader
      implicit none
      integer ierr
      logical disp

c     checks if the intervall bounds for vlambda and elambda are consistent
      if (rank.eq.0) then
       if ( bvlambda .le. 0.0 .OR. bvlambda .gt. 1.0 ) then
         write(iout,*) 'Intervall bound for vlambda must be between 0 ',
     $ 'and 1'
         call fatal
       else if ( belambda .lt. 0.0 .OR. belambda .ge. 1.0 ) then
         write(iout,*) 'Intervall bound for elambda must be between 0 ',
     $ 'and 1'
         call fatal
       end if

       if (bvlambda.eq.0.0) then
 10       format('Forbidden value of bvlambda: ',F15.3)
          call fatal
       end if
       if (belambda.eq.1.0) then
 15       format('Forbidden value of belambda: ',F15.3)
          call fatal
       end if

       if (bvlambda .lt. belambda ) then
          write(iout,*) 'Intervall bound for vlambda cannot be ',
     $ 'less than the intervall bound for elambda'
          call fatal
       end if
      end if

      if (lambda.lt.0.0) then
         vlambda = 0.0
         dlambdavlambda = 0.0
         elambda = 0.0
         dlambdaelambda = 0.0
      else if (lambda.lt.belambda) then
         vlambda = lambda/bvlambda
         dlambdavlambda = 1.0 / bvlambda
         elambda = 0.0
         dlambdaelambda = 0.0
      else if (lambda.lt.bvlambda) then
         vlambda = lambda / bvlambda
         dlambdavlambda = 1.0/bvlambda
         elambda = (1.0/(1-belambda))*lambda - belambda/(1-belambda)
         dlambdaelambda = 1.0/(1-belambda)
      else if (lambda.lt.1.0) then
         vlambda = 1.0
         dlambdavlambda = 0.0
         elambda = (1.0/(1-belambda)) * lambda - belambda/(1-belambda)
         dlambdaelambda = 1.0/(1-belambda)
      else
         vlambda = 1.0
         dlambdavlambda = 0.0
         elambda = 1.0
         dlambdaelambda = 0.0
      end if

#if (TINKER_SINGLE_PREC + TINKER_MIXED_PREC)
      disp = rank.eq.0.and.verbose.and.mod(step_c,iprint/100).eq.0
#else
      disp = rank.eq.0.and.verbose
#endif
      if (disp) then
         write(iout,20) vlambda, dlambdavlambda
 20      format('Values of vlambda/dlambdavlambda are ', 2F15.3)
         write(iout,25) elambda, dlambdaelambda
 25      format('Values of elambda/dlambdaelambda are ', 2F15.3)
      end if
#ifdef _OPENACC
      call altelec(0)
#else
      call altelec(1)
#endif
c     Initialising double derivatives in case of OSRW
      if (use_OSRW) then
         d2edlambda2  = 0.0
         d2edlambdae2 = 0.0
         d2edlambdav2 = 0.0
      end if
      end

      subroutine upload_device_mutate
      use charge
      use cflux
      use chgpen
      use chgtrn
      use domdec,only: rank
      use inform,only: deb_Path
      use mpole
      use mutant
      use polar
      use potent
      use tors
#ifdef _OPENACC
      use utilcu ,only: cu_update_vcouple
#endif
      implicit none
#ifdef _OPENACC
 12   format(2x,'upload_device_mutate')
      if(deb_Path) print 12
#endif

      if (use_vdw) then
!$acc update device(mut,mutInt) async
#ifdef _OPENACC
        call cu_update_vcouple(vcouple)
#endif
      end if
      if (use_charge) then
!$acc update device(pchg,pchg0) async
      end if
      if (use_chgtrn) then
!$acc update device(chgct,pcore,pval,pval0) async
      end if
      if (use_chgflx) then
!$acc update device(bflx,aflx,abflx) async
      end if
      if (use_mpole.or.use_polar) then
!$acc update device(pole,mono0) async
      end if
      if (use_polar) then
!$acc update device(polarity) async
      end if
      if (use_tors) then
!$acc update device(tors1,tors2,tors3,tors4,tors5,tors6) async
      end if
      end subroutine

      subroutine dealloc_shared_mutate
      use domdec,only:rank
      use mutant
      use tinMemory
      use sizes ,only:tinkerdebug
      implicit none

 12   format(2x,'dealloc_shared_mutate')
      if (rank.eq.0.and.tinkerdebug) print 12

      call shmem_request(imut,  winimut,  [0])
      !TODO Remove mut from device
      call shmem_request(mut,   winmut,   [0], config=mhostacc)
      call shmem_request(mutInt,winmutInt,[0], config=mhostacc)
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
      !TODO Remove mut from device
      call shmem_request(mut   ,winmut,   [n], config=mhostacc)
      call shmem_request(mutInt,winmutInt,[n], config=mhostacc)
      call shmem_request(type0 ,wintype0, [n])
      call shmem_request(type1 ,wintype1, [n])
      call shmem_request(class0,winclass0,[n])
      call shmem_request(class1,winclass1,[n])
      end
