
c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #############################################################################
c     ##                                                                         ##
c     ##  subroutine epolar3tcg1_2  --  Ewald polarization analysis with tcg1    ##
c     ##                                                                         ##
c     #############################################################################
c
c
c     "epolar3tcg1_2" calculates the induced dipole polarization energy with the tcg1 method,
c     and partitions the energy among atoms, with separate cores for PME
c  
      subroutine epolar3tcg1_2
      use action
      use analyz
      use atmlst
      use atoms
      use boxes
      use chgpot
      use domdec
      use energi
      use ewald
      use inform
      use iounit
      use math
      use mpole
      use pme
      use polar
      use polpot
      use potent
      use mpi
      implicit none
c
      logical :: peek, precond, isguess
      integer :: i, iipole, nrhs, j, ierr

      integer, allocatable :: reqrec(:),reqsend(:)
      integer, allocatable :: reqrecdirrec(:),reqrecdirsend(:)
      integer, allocatable :: req2rec(:),req2send(:)
      real*8 :: e, term, sprod, omega, f

      real*8, dimension(2) :: n0, t1, t4
      real*8, allocatable, dimension(:,:) :: cphi
      real*8, allocatable, dimension(:,:,:) :: ef , res, resbis,
     $         mu_tcg1, r0, r0bis, Tr0, Tr0bis, aefi, aTr0, Taefi, 
     $         mupeek, ar0, oldr0, oldTr0, mu_tcg1rec,
     $         mu0, mu0bis, res2,res2bis
      real*8, allocatable, dimension(:,:,:) :: buffermpimu
      real*8, allocatable, dimension(:,:) :: buffermpi

      parameter (nrhs=2)
c
      if (deb_Path) write(iout,*), 'epolar3tcg1_2 '
c

      if (rank.le.ndir-1) then
        allocate(res(3,2,max(1,npolebloc)),res2(3,nrhs,max(1,npoleloc)))
        allocate (cphi(10,max(npoleloc,1)))
        cphi = 0d0
        allocate (buffermpi(10,max(npoleloc,1)))
        allocate (buffermpimu(3,nrhs,max(npoleloc,1)))
        allocate(ef(3,nrhs,max(1,npolebloc)), mu_tcg1(3,nrhs,npolebloc),
     $         r0(3,nrhs,npolebloc),Tr0(3,nrhs,npolebloc),
     $         atr0(3,nrhs,npolebloc),aefi(3,nrhs,npolebloc),
     $         taefi(3,nrhs,npolebloc),mupeek(3,nrhs,npolebloc),
     $         ar0(3,nrhs,npolebloc),oldr0(3,nrhs,npolebloc), 
     $         oldTr0(3,nrhs,npolebloc),mu0(3,nrhs,npolebloc))
        mu0 = 0d0
        mupeek = 0d0
        ef = 0d0
        aefi = 0d0
        taefi = 0d0
        r0 = 0d0
        ar0 = 0d0
        Tr0 = 0d0
        atr0 = 0d0
        mu_tcg1 = 0d0
        res = 0d0
        res2 = 0d0
      else
        allocate (resbis(3,2,max(1,npolerecloc)),
     $         res2bis(3,nrhs,max(1,npolerecloc)))
        allocate (mu_tcg1rec(3,nrhs,max(1,npolerecloc)))
        mu_tcg1rec = 0d0
        if (allocated(cphirec)) deallocate(cphirec)
        allocate(cphirec(10, npolerecloc))
        if (allocated(fphirec)) deallocate(fphirec)
        allocate(fphirec(20, npolerecloc))
        allocate (buffermpi(10,max(npolerecloc,1)))
        allocate (buffermpimu(3,nrhs,max(npolerecloc,1)))
        allocate(r0bis(3,nrhs,npolerecloc),Tr0bis(3,nrhs,npolerecloc),
     $         mu0bis(3,nrhs,npolerecloc))
        allocate (cphi(10,max(1,npoleloc)))
        cphi = 0.0d0
      end if
      allocate (reqrec(nproc))
      allocate (reqsend(nproc))
      allocate (req2rec(nproc))
      allocate (req2send(nproc))
      allocate (reqrecdirrec(nproc))
      allocate (reqrecdirsend(nproc))
c
c     set the energy unit conversion factor
c
      f = electric / dielec

      precond = tcgprec
      isguess = tcgguess
      peek    = tcgpeek

      omega = tcgomega

c
      if (rank.le.ndir-1) then
        call commdirdir(nrhs,0,ef,reqrec,reqsend)
        call commrecdirfields(0,cphirec,cphi,buffermpi,buffermpi,
     $   reqrecdirrec,reqrecdirsend)
        call efld0_direct(nrhs,ef)
        call commfield(nrhs,ef)
        call commrecdirfields(2,cphirec,cphi,buffermpi,buffermpi,
     $   reqrecdirrec,reqrecdirsend)

        term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
        do i = 1, npoleloc
          iipole = poleglob(i)
          do j = 1, 3
            ef(j,1,i)  = ef(j,1,i) 
     $        - cphi(j+1,i) 
     $        + term*rpole(j+1,iipole)
            ef(j,2,i)  = ef(j,2,i) 
     $        - cphi(j+1,i) 
     $        + term*rpole(j+1,iipole)
          end do
        end do
        call commdirdir(nrhs,1,ef,reqrec,reqsend)
        call commdirdir(nrhs,2,ef,reqrec,reqsend)
      else
        call efld0_recip(cphi)
        call commrecdirfields(1,cphirec,cphi,buffermpi,buffermpi,
     $   reqrecdirrec,reqrecdirsend)
        call commrecdirfields(2,cphirec,cphi,buffermpi,buffermpi,
     $   reqrecdirrec,reqrecdirsend)
      end if
c

      if (rank.le.ndir-1) then
        if (isguess) then
          call diagvec(nrhs, ef, mu0)
          call commrecdirdip(nrhs,1,mu0bis,mu0,buffermpimu,
     $     buffermpimu,req2rec,req2send)
          call commrecdirdip(nrhs,2,mu0bis,mu0,buffermpimu,
     $    buffermpimu,req2rec,req2send)
          call commrecdirsolv(nrhs,0,res2bis,res2,buffermpi,
     $     buffermpi,reqrecdirrec,reqrecdirsend)

          call tmatxb_pme(2, .true., mu0, res)
          call commfield(nrhs,res)

          call commrecdirsolv(nrhs,2,res2bis,res2,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)

          res(:,:,1:npoleloc) = res(:,:,1:npoleloc) 
     $    +  res2(:,:,1:npoleloc) - term*mu0(:,:,1:npoleloc)
        else
          res = 0d0
        end if

        r0(:,:,1:npoleloc) = ef(:,:,1:npoleloc) - res(:,:,1:npoleloc)
        oldr0 = r0
        if (precond) call diagvec(nrhs, r0, r0)
        call commdirdir(nrhs,0,r0,reqrec,reqsend)
        call commdirdir(nrhs,1,r0,reqrec,reqsend)
        call commdirdir(nrhs,2,r0,reqrec,reqsend)
         
        call commrecdirdip(nrhs,1,r0bis,r0,buffermpimu,buffermpimu,
     $   req2rec,req2send)
        call commrecdirdip(nrhs,2,r0bis,r0,buffermpimu,buffermpimu,
     $   req2rec,req2send)
        call commrecdirsolv(nrhs,0,res2bis,res2,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)

        call tmatxb_pme(2, .true., r0, tr0)
        call commfield(nrhs,tr0)

        call commrecdirsolv(nrhs,2,res2bis,res2,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)

        Tr0(:,:,1:npoleloc) = Tr0(:,:,1:npoleloc)+res2(:,:,1:npoleloc)
     $   - term*r0(:,:,1:npoleloc)


        oldTr0 = Tr0
        if (precond) call diagvec(nrhs, Tr0, Tr0)

        
        do i = 1,nrhs
          n0(i) = sprod(3*npoleloc, oldr0(:,i,:), r0(:,i,:))
          t1(i) = sprod(3*npoleloc, r0(:,i,:), oldTr0(:,i,:))
          call MPI_ALLREDUCE(MPI_IN_PLACE,n0(i),1,MPI_REAL8,MPI_SUM,
     $       comm_dir,ierr)
          call MPI_ALLREDUCE(MPI_IN_PLACE,t1(i),1,MPI_REAL8,MPI_SUM,
     $       comm_dir,ierr)
          t4(i) = n0(i)/t1(i)
          mu_tcg1(:,i,:) = t4(i) * r0(:,i,:)
        end do

        if (isguess) then
           mu_tcg1 = mu_tcg1 + mu0
        end if
c
c     if omegafitstep then store intermediate values
c
        if (omegafitstep) then
           residue = r0(:,1,:) - t4(1)*tr0(:,1,:)
           munp = mu_tcg1(:,1,:)
           efres = ef(:,2,:)
        end if

        if (peek .and. precond) then
           mupeek = omega*r0 - omega*t4(1)*tr0
           mu_tcg1 = mu_tcg1 + mupeek
        else if (peek .and. .not. precond) then
           call diagvec(nrhs, r0, ar0)
           call diagvec(nrhs, tr0, atr0)
           mupeek = omega*ar0 - omega*t4(1)*atr0
           mu_tcg1 = mu_tcg1 + mupeek
        end if
        e = sprod(3*npoleloc, ef(:,2,:), mu_tcg1(:,1,:))
        ep =  -.5d0*e*f
        nep = nloc
c
c     move dipoles in global variables
c
        do i = 1, npolebloc
          iipole = poleglob(i)
          do j = 1,3
            uind(j,iipole) = mu_tcg1(j,1,i)
            uinp(j,iipole) = mu_tcg1(j,2,i)
          end do
        end do

      else
        if (isguess) then
          call commrecdirdip(nrhs,0,mu0bis,mu0,buffermpimu,
     $     buffermpimu,req2rec,req2send)
          call commrecdirdip(nrhs,2,mu0bis,mu0,buffermpimu,
     $     buffermpimu,req2rec,req2send)
          call tmatxbrecip(mu0,mu0bis,2, res2, res2bis)

          call commrecdirsolv(nrhs,1,res2bis,res2,buffermpi,
     $     buffermpi,reqrecdirrec,reqrecdirsend)
          call commrecdirsolv(nrhs,2,res2bis,res2,buffermpi,
     $     buffermpi,reqrecdirrec,reqrecdirsend)
        end if

        call commrecdirdip(nrhs,0,r0bis,r0,buffermpimu,buffermpimu,
     $   req2rec,req2send)
        call commrecdirdip(nrhs,2,r0bis,r0,buffermpimu,buffermpimu,
     $   req2rec,req2send)
        call tmatxbrecip(r0,r0bis,2, res2, res2bis)

        call commrecdirsolv(nrhs,1,res2bis,res2,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,res2bis,res2,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)

      end if
c
      return
      end 
c
c
c     ###############################################################################
c     ##                                                                           ##
c     ##  subroutine epolar3tcg2_2  --  Ewald polarization analysis with tcg2      ##
c     ##                                                                           ##
c     ###############################################################################
c
c
c     "epolar3tcg2_2" calculates the induced dipole polarization energy with the tcg2 method,
c     and partitions the energy among atoms with separate cores for PME
c  
      subroutine epolar3tcg2_2
      use action
      use analyz
      use atmlst
      use atoms
      use boxes
      use chgpot
      use domdec
      use energi
      use ewald
      use inform
      use iounit
      use math
      use mpole
      use pme
      use polar
      use polpot
      use potent
      use mpi
      implicit none
      logical :: peek, precond, isguess
      integer :: i, iipole, nrhs, j, ierr

      integer, allocatable :: reqrec(:),reqsend(:)
      integer, allocatable :: reqrecdirrec(:),reqrecdirsend(:)
      integer, allocatable :: req2rec(:),req2send(:)
      real*8 :: e, term, sprod, sprodtemp, omega, f

      real*8, dimension(2) :: n0, t1, t4, t2, t3, g1, np1
      real*8, allocatable, dimension(:,:) :: cphi
      real*8, allocatable, dimension(:,:,:) :: ef, res, resbis,
     $         mu_tcg2, r0, r0bis, Tr0, Tr0bis, aefi, aTr0, Taefi, 
     $         mupeek, ar0, T2r0, T2r0bis, P2, oldr0, oldTr0,
     $         mu0, mu0bis, at2r0, oldT2r0,res2,res2bis
      real*8, allocatable, dimension(:,:,:) :: buffermpimu
      real*8, allocatable, dimension(:,:) :: buffermpi

      parameter (nrhs=2)
c
      if (deb_Path) write(iout,*), 'epolar3tcg2_2 '
c

      if (rank.le.ndir-1) then
        allocate(res(3,2,max(1,npolebloc)),cphi(10, max(1,npoleloc)), 
     $           res2(3,nrhs,npoleloc))
        allocate(ef(3,nrhs,npolebloc), mu_tcg2(3,nrhs,npolebloc),
     $           r0(3,nrhs,npolebloc),Tr0(3,nrhs,npolebloc),
     $           T2r0(3,nrhs,npolebloc),at2r0(3,nrhs,npolebloc),
     $           atr0(3,nrhs,npolebloc), aefi(3,nrhs,npolebloc),
     $           taefi(3,nrhs,npolebloc), mupeek(3,nrhs,npolebloc),
     $           ar0(3,nrhs,npolebloc), P2(3,nrhs,npolebloc),
     $           oldr0(3,nrhs,npolebloc), oldTr0(3,nrhs,npolebloc),
     $           mu0(3,nrhs,npolebloc))
        allocate (buffermpi(10,max(npoleloc,1)))
        allocate (buffermpimu(3,nrhs,max(npoleloc,1)))
        mu0 = 0d0
        mupeek = 0d0
        cphi = 0d0
        ef = 0d0
        aefi = 0d0
        taefi = 0d0
        r0 = 0d0
        ar0 = 0d0
        Tr0 = 0d0
        atr0 = 0d0
        mu_tcg2 = 0d0
        res = 0d0
        res2 = 0d0
      else
        allocate(resbis(3,2,npolerecloc),res2bis(3,nrhs,npolerecloc))
        if (allocated(fphirec)) deallocate(fphirec)
        allocate(fphirec(20, npolerecloc))
        if (allocated(cphirec)) deallocate(cphirec)
        allocate(cphirec(10, max(1,npolerecloc)))
        allocate(r0bis(3,nrhs,npolerecloc),Tr0bis(3,nrhs,npolerecloc),
     $           T2r0bis(3,nrhs,npolerecloc),mu0bis(3,nrhs,npolerecloc))
        allocate (buffermpi(10,max(npolerecloc,1)))
        allocate (buffermpimu(3,nrhs,max(npolerecloc,1)))
      end if

      allocate (reqrec(nproc))
      allocate (reqsend(nproc))
      allocate (req2rec(nproc))
      allocate (req2send(nproc))
      allocate (reqrecdirrec(nproc))
      allocate (reqrecdirsend(nproc))
c
c
c     set the energy unit conversion factor
c
      f = electric / dielec
      precond = tcgprec
      isguess = tcgguess
      peek    = tcgpeek

      omega = tcgomega


      if (rank.le.ndir-1) then 
        call commdirdir(nrhs,0,ef,reqrec,reqsend)
        call commrecdirfields(0,cphirec,cphi,buffermpi,buffermpi,
     $   reqrecdirrec,reqrecdirsend)
        call efld0_direct(nrhs,ef)
        call commfield(nrhs,ef)
        call commrecdirfields(2,cphirec,cphi,buffermpi,buffermpi,
     $   reqrecdirrec,reqrecdirsend)
        term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
        do i = 1, npoleloc
          iipole = poleglob(i)
          do j = 1, 3
            ef(j,1,i)  = ef(j,1,i) 
     $        - cphi(j+1,i) 
     $        + term*rpole(j+1,iipole)
            ef(j,2,i)  = ef(j,2,i) 
     $        - cphi(j+1,i) 
     $        + term*rpole(j+1,iipole)
          end do
        end do

        call commdirdir(nrhs,1,ef,reqrec,reqsend)
        call commdirdir(nrhs,2,ef,reqrec,reqsend)
      else
        call efld0_recip(cphi)
        call commrecdirfields(1,cphirec,cphi,buffermpi,buffermpi,
     $   reqrecdirrec,reqrecdirsend)
        call commrecdirfields(2,cphirec,cphi,buffermpi,buffermpi,
     $   reqrecdirrec,reqrecdirsend)
      end if

      if (rank.le.ndir-1) then
        if (isguess) then
          call diagvec(nrhs, ef, mu0)
          call commrecdirdip(nrhs,1,mu0bis,mu0,buffermpimu,
     $     buffermpimu,req2rec,req2send)
          call commrecdirdip(nrhs,2,mu0bis,mu0,buffermpimu,
     $     buffermpimu,req2rec,req2send)
          call commrecdirsolv(nrhs,0,res2bis,res2,buffermpi,
     $      buffermpi,reqrecdirrec,reqrecdirsend)

          call tmatxb_pme(2, .true., mu0, res)
          call commfield(nrhs,res)

          call commrecdirsolv(nrhs,2,res2bis,res2,buffermpi,
     $     buffermpi,reqrecdirrec,reqrecdirsend)

          res(:,:,1:npoleloc) = res(:,:,1:npoleloc) 
     $     +  res2(:,:,1:npoleloc) - term*mu0(:,:,1:npoleloc)
        else
          res = 0d0
        end if

        r0(:,:,1:npoleloc) = ef(:,:,1:npoleloc) - res(:,:,1:npoleloc)

        oldr0 = r0
        if (precond) call diagvec(nrhs, r0, r0)

        call commdirdir(nrhs,0,r0,reqrec,reqsend)
        call commdirdir(nrhs,1,r0,reqrec,reqsend)
        call commdirdir(nrhs,2,r0,reqrec,reqsend)

        call commrecdirdip(nrhs,1,r0bis,r0,buffermpimu,buffermpimu,
     $   req2rec,req2send)
        call commrecdirdip(nrhs,2,r0bis,r0,buffermpimu,buffermpimu,
     $   req2rec,req2send)

        call commrecdirsolv(nrhs,0,res2bis,res2,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)

        call tmatxb_pme(2, .true., r0, tr0)
        call commfield(nrhs,tr0)
        call commrecdirsolv(nrhs,2,res2bis,res2,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)

        Tr0(:,:,1:npoleloc) = tr0(:,:,1:npoleloc)+res2(:,:,1:npoleloc) 
     $   - term*r0(:,:,1:npoleloc)

        oldTr0 = Tr0
        if (precond) call diagvec(nrhs, Tr0, Tr0)

        do i = 1,nrhs
           n0(i) = sprod(3*npoleloc, oldr0(:,i,:), r0(:,i,:))
           t1(i) = sprod(3*npoleloc, r0(:,i,:), oldTr0(:,i,:))
           call MPI_ALLREDUCE(MPI_IN_PLACE,n0(i),1,MPI_REAL8,MPI_SUM,
     $        comm_dir,ierr)
           call MPI_ALLREDUCE(MPI_IN_PLACE,t1(i),1,MPI_REAL8,MPI_SUM,
     $        comm_dir,ierr)
           t4(i) = n0(i)/t1(i)
        end do
        if (peek .and. precond) then
           mupeek = omega*r0 - omega*t4(1)*tr0
        else if (peek .and. .not. precond) then
           call diagvec(nrhs, r0, ar0)
           call diagvec(nrhs, tr0, atr0)
           mupeek(:,:,1:npoleloc) = omega*ar0(:,:,1:npoleloc)
     $      - omega*t4(1)*atr0(:,:,1:npoleloc)
        end if
        call commdirdir(nrhs,0,Tr0,reqrec,reqsend)
        call commdirdir(nrhs,1,Tr0,reqrec,reqsend)
        call commdirdir(nrhs,2,Tr0,reqrec,reqsend)

        call commrecdirdip(nrhs,1,Tr0bis,Tr0,buffermpimu,buffermpimu,
     $   req2rec,req2send)
        call commrecdirdip(nrhs,2,Tr0bis,Tr0,buffermpimu,buffermpimu,
     $   req2rec,req2send)
        call commrecdirsolv(nrhs,0,res2bis,res2,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)

        call tmatxb_pme(2, .true., Tr0, T2r0)
        call commfield(nrhs,T2r0)

        call commrecdirsolv(nrhs,2,res2bis,res2,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)

        T2r0(:,:,1:npoleloc) = T2r0(:,:,1:npoleloc)+res2(:,:,1:npoleloc)
     $   - term*Tr0(:,:,1:npoleloc)

        oldT2r0 = T2r0
        if (precond) then
           call diagvec(nrhs, T2r0, T2r0)
        end if
        do i = 1,nrhs
           np1(i) = sprod(3*npoleloc,oldTr0(:,i,:), Tr0(:,i,:))
           call MPI_ALLREDUCE(MPI_IN_PLACE,np1,nrhs,MPI_REAL8,MPI_SUM,
     $        comm_dir,ierr)
           t2(i) = n0(i)*np1(i)/(t1(i)*t1(i))
           P2(:,i,:) = t2(i)*Tr0(:,i,:) - t4(i)*T2r0(:,i,:)
           sprodtemp = sprod(3*npoleloc, oldTr0(:,i,:),P2(:,i,:) )
           call MPI_ALLREDUCE(MPI_IN_PLACE,sprodtemp,1,MPI_REAL8,
     $       MPI_SUM,comm_dir,ierr)
           t3(i) = t1(i)*sprodtemp
           g1(i) = (t1(i)*t1(i) - n0(i)*np1(i))/t3(i)
           mu_tcg2(:,i,:) = (g1(i)*t2(i) + t4(i))*r0(:,i,:) 
     $                         - g1(i)*t4(i)*Tr0(:,i,:)
        end do
c
c     if omegafitstep then store intermediate values
c
        if (isguess) mu_tcg2(:,:,1:npoleloc) = mu_tcg2(:,:,1:npoleloc) 
     $    + mu0(:,:,1:npoleloc)
        if (omegafitstep) then
           residue = r0(:,1,:) + tr0(:,1,:)*(-g1(1)*t2(1) - t4(1))
     $               + T2r0(:,1,:)*(g1(1)*t4(1))
           munp = mu_tcg2(:,1,:)
           efres = ef(:,2,:)
        end if

        if (peek .and. .not. precond) then
           call diagvec(nrhs, T2r0, at2r0)
           mupeek = mupeek - omega*g1(1)*t2(1)*aTr0 
     $                 + omega*g1(1)*t4(1)*aT2r0
        else if (peek .and. precond) then
           mupeek = mupeek - omega*g1(1)*t2(1)*Tr0
     $                 + omega*g1(1)*t4(1)*T2r0
        end if
c
        if (peek) then
           mu_tcg2(:,:,1:npoleloc) = mu_tcg2(:,:,1:npoleloc) 
     $    + mupeek(:,:,1:npoleloc)
        end if

        e = sprod(3*npoleloc, ef(:,2,:), mu_tcg2(:,1,:))
        ep =  -.5d0*e*f
        nep = nloc
c
c     move dipoles in global variables
c
        do i = 1, npolebloc
          iipole = poleglob(i)
          do j = 1,3
            uind(j,iipole) = mu_tcg2(j,1,i)
            uinp(j,iipole) = mu_tcg2(j,2,i)
          end do
        end do

      else 
        if (isguess) then
          call commrecdirdip(nrhs,0,mu0bis,mu0,buffermpimu,
     $     buffermpimu,req2rec,req2send)
          call commrecdirdip(nrhs,2,mu0bis,mu0,buffermpimu,
     $     buffermpimu,req2rec,req2send)

          call tmatxbrecip(mu0,mu0bis,2, res2, res2bis)

          call commrecdirsolv(nrhs,1,res2bis,res2,buffermpi,
     $      buffermpi,reqrecdirrec,reqrecdirsend)
          call commrecdirsolv(nrhs,2,res2bis,res2,buffermpi,
     $     buffermpi,reqrecdirrec,reqrecdirsend)
        end if
        call commrecdirdip(nrhs,0,r0bis,r0,buffermpimu,buffermpimu,
     $   req2rec,req2send)
        call commrecdirdip(nrhs,2,r0bis,r0,buffermpimu,buffermpimu,
     $   req2rec,req2send)
        call tmatxbrecip(r0,r0bis,2, res2, res2bis)
        call commrecdirsolv(nrhs,1,res2bis,res2,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,res2bis,res2,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)

        call commrecdirdip(nrhs,0,Tr0bis,Tr0,buffermpimu,buffermpimu,
     $   req2rec,req2send)
        call commrecdirdip(nrhs,2,Tr0bis,Tr0,buffermpimu,buffermpimu,
     $   req2rec,req2send)

        call tmatxbrecip(Tr0,Tr0bis,2, res2, res2bis)

        call commrecdirsolv(nrhs,1,res2bis,res2,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,res2bis,res2,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
      end if
      return
      end 
