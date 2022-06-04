c      
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c      
c     evaluate induced dipole moments and the polarization energy
c     using either a (preconditioned) conjugate gradient algorithm or
c     Jacobi iterations coupled with DIIS extrapolation.
c
c     literature reference:
c     "Scalable Evaluation of Polarization Energy and Associated Forces
c     in Polarizable Molecular Dynamics: II. Toward Massively Parallel
c     Computations Using Smooth Particle Mesh Ewald",L. Lagardere et al.,
c     J. Chem. Theory Comput., 2015, 11 (6), pp 2589–2599
c
      subroutine newinduce_pme
      use atmlst
      use domdec
      use ewald
      use iounit
      use math
      use mpole
      use pme
      use polar
      use polpot
      use potent
      use units
      use uprior
      use mpi
      implicit none
c
c     with separate cores for reciprocal part
c
      integer i, j, k, nrhs
c
c     MPI
c
      integer iipole
      integer, allocatable :: reqrecdirsend(:),reqrecdirrec(:)
      integer, allocatable :: reqsend(:),reqrec(:)
      integer, allocatable :: req2send(:),req2rec(:)


c
      parameter (nrhs=2)
      real*8  wtime0, wtime1, wtime2, udsum, upsum
      real*8  term, xx(1)
      real*8, allocatable :: ef(:,:,:), mu(:,:,:), murec(:,:,:)
      real*8, allocatable :: cphi(:,:)
c
      real*8, allocatable :: buffermpi(:,:), buffermpimu(:,:,:)
c
      external tmatxb_pme
c
      if (.not.use_polar) return
c
 1000 format(' illegal polalg in newinduce.')
 1010 format(' time for the ',a,F14.5)
 1020 format(' total elapsed time in newinduce: ',F14.5)
c
c     allocate some memory and clear the arrays:
c
      if (rank.le.ndir-1) then
        allocate (mu(3,nrhs,max(1,npolebloc)))
        mu = 0.0d0
        allocate (murec(3,nrhs,max(npolerecloc,1)))
        murec = 0.0d0
        allocate (ef(3,nrhs,max(1,npolebloc)))
        ef = 0.0d0
        allocate (cphi(10,max(1,npoleloc)))
        cphi = 0.0d0
        if (allocated(cphirec)) deallocate (cphirec)
        allocate (cphirec(10,max(npolerecloc,1)))
        cphirec = 0.0d0
        allocate (buffermpi(10,max(npoleloc,1)))
        buffermpi = 0.0d0
        allocate (buffermpimu(3,nrhs,max(npoleloc,1)))
        buffermpimu = 0.0d0
      else 
        allocate (mu(3,nrhs,max(1,npolebloc)))
        mu = 0.0d0
        allocate (murec(3,nrhs,max(npolerecloc,1)))
        murec = 0.0d0
        if (allocated(fphirec)) deallocate (fphirec)
        allocate (fphirec(20,max(npolerecloc,1)))
        fphirec = 0.0d0
        allocate (cphi(10,max(1,npoleloc)))
        cphi = 0.0d0
        if (allocated(cphirec)) deallocate (cphirec)
        allocate (cphirec(10,max(npolerecloc,1)))
        cphirec = 0.0d0
        allocate (buffermpi(10,max(npolerecloc,1)))
        buffermpi = 0.0d0
        allocate (buffermpimu(3,nrhs,max(npolerecloc,1)))
        buffermpimu = 0.0d0
      end if
c
      allocate (reqrecdirsend(nproc))
      allocate (reqrecdirrec(nproc))
      allocate (reqsend(nproc))
      allocate (reqrec(nproc))
      allocate (req2send(nproc))
      allocate (req2rec(nproc))
c
c     compute the electric fields:
c
      wtime0 = mpi_wtime()
c
      if (rank.le.ndir-1) then
        call commdirdir(nrhs,0,mu,reqrec,reqsend)
c
        call commrecdirfields(0,cphirec,cphi,buffermpi,buffermpi,
     $   reqrecdirrec,reqrecdirsend)
      else 
c
c    compute the reciprocal space contribution (fields)
c
        call efld0_recip(cphi)
c
        call commrecdirfields(1,cphirec,cphi,buffermpi,buffermpi,
     $   reqrecdirrec,reqrecdirsend)
        call commrecdirfields(2,cphirec,cphi,buffermpi,buffermpi,
     $   reqrecdirrec,reqrecdirsend)
      end if
c
c    The real space processes compute the real fields and  add them to the recip ones
c
      if (rank.le.ndir-1) then
        term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
c
        call efld0_direct(nrhs,ef)
        call commfield(nrhs,ef)
c
        call commrecdirfields(2,cphirec,cphi,buffermpi,buffermpi,
     $   reqrecdirrec,reqrecdirsend)
c
c       Add direct and reciprocal fields
c
        do i = 1, npoleloc
          iipole = poleglob(i)
          do j = 1, 3
            ef(j,1,i)  = ef(j,1,i) - cphi(j+1,i) +
     $        term*rpole(j+1,iipole)
            ef(j,2,i)  = ef(j,2,i) - cphi(j+1,i) +
     $        term*rpole(j+1,iipole)
          end do
        end do
      end if
      wtime1 = mpi_wtime()
c
c     guess the dipoles:
c
c     predicted values for always stable predictor-corrector method
c
      if (rank.le.ndir-1) then
       if (use_pred .and. nualt.eq.maxualt) then
         call ulspred
         do i = 1, npoleloc
           iipole = poleglob(i)
           do j = 1, 3
             udsum = 0.0d0
             upsum = 0.0d0
             do k = 1, nualt - 1
               udsum = udsum + bpred(k)*udalt(k,j,iipole)
               upsum = upsum + bpred(k)*upalt(k,j,iipole)
             end do
             mu(j,1,i) = udsum
             mu(j,2,i) = upsum
           end do
         end do
       else if (polgsf.eq.0) then
         do i = 1, npoleloc
           iipole = poleglob(i)
           do k = 1, nrhs
             do j = 1, 3
               mu(j,k,i) = polarity(iipole)*ef(j,k,i)
             end do
           end do
         end do
       else
         do i = 1, npoleloc
           iipole = poleglob(i)
           do j = 1, 3
             mu(j,1,i) = uind(j,iipole)
             mu(j,2,i) = uinp(j,iipole)
           end do
         end do
       end if
c
       call commdirdir(nrhs,1,mu,reqrec,reqsend)
c
        call commrecdirdip(nrhs,1,murec,mu,buffermpimu,buffermpimu,
     $   req2rec,req2send)
c
       call commdirdir(nrhs,2,mu,reqrec,reqsend)
       call commrecdirdip(nrhs,2,murec,mu,buffermpimu,buffermpimu,
     $  req2rec,req2send)
      else 
c
        call commrecdirdip(nrhs,0,murec,mu,buffermpimu,buffermpimu,
     $   req2rec,req2send)
        call commrecdirdip(nrhs,2,murec,mu,buffermpimu,buffermpimu,
     $   req2rec,req2send)
      end if
c
c     now, call the proper solver.
c
      if (polalg.eq.1) then
        if (rank.le.ndir-1) then
          call inducepcg_pme(tmatxb_pme,nrhs,.true.,ef,mu,xx)
        else
          call inducepcg_pme(tmatxb_pme,nrhs,.true.,xx,xx,murec)
        end if
      else if (polalg.eq.2) then
        if (rank.le.ndir-1) then
          call inducejac_pme(tmatxb_pme,nrhs,.true.,ef,mu,xx)
       else
          call inducejac_pme(tmatxb_pme,nrhs,.true.,xx,xx,murec)
       end if
      else 
         if (rank.eq.0) write(iout,1000) 
         call fatal
      end if
      wtime2 = mpi_wtime()
      if (polprt.ge.1.and.rank.eq.0) then
        if (polprt.ge.2) then
          write(iout,1010) 'fields:  ', wtime1-wtime0
          write(iout,1010) 'dipoles: ', wtime2-wtime1
        end if
        write(iout,1020) wtime2 - wtime0
      end if
c
c     move the computed dipoles in the module.
c
      if (rank.le.ndir-1) then
        do i = 1, npolebloc
          iipole = poleglob(i)
          do j = 1, 3
            uind(j,iipole) = mu(j,1,i)
            uinp(j,iipole) = mu(j,2,i)
          end do
        end do
      else  
        do i = 1, npolerecloc
          iipole = polerecglob(i)
          do j = 1, 3
            uind(j,iipole) = murec(j,1,i)
            uinp(j,iipole) = murec(j,2,i)
          end do
        end do
      end if
      deallocate (buffermpi)
      deallocate (buffermpimu)
      if (rank.le.ndir-1) then
        deallocate (ef)
        deallocate (cphi)
        deallocate (mu)
      else
        deallocate (murec)
      end if
      deallocate (reqrec)
      deallocate (reqsend)
      deallocate (req2send)
      deallocate (req2rec)
      deallocate (reqrecdirrec)
      deallocate (reqrecdirsend)
c
c     update the lists of previous induced dipole values
c
      if ((use_pred).and..not.(use_lambdadyn)) then
        if (rank.le.ndir-1) then
           nualt = min(nualt+1,maxualt)
           do i = 1, npolebloc
              iipole = poleglob(i)
              do j = 1, 3
                 do k = nualt, 2, -1
                    udalt(k,j,iipole) = udalt(k-1,j,iipole)
                    upalt(k,j,iipole) = upalt(k-1,j,iipole)
                 end do
                 udalt(1,j,iipole) = uind(j,iipole)
                 upalt(1,j,iipole) = uinp(j,iipole)
               end do
           end do
         end if
      end if
      return
      end
c
      subroutine inducepcg_pme(matvec,nrhs,precnd,ef,mu,murec)
      use atmlst
      use domdec
      use ewald
      use iounit
      use math
      use mpole
      use polar
      use polpot
      use potent
      use timestat
      use units
      use mpi
      implicit none
c
c     solves the polarizable dipoles linear equations by preconditioned
c     conjugate gradient. A diagonal preconditioner is used when precnd
c     is true, otherwise the standard conjugate gradient algorithm is
c     recovered by setting the preconditioner to one.
c
      integer nrhs
      real*8  ef(3,nrhs,*), mu(3,nrhs,*), murec(3,nrhs,*)
      logical precnd
      real*8, allocatable :: res(:,:,:), h(:,:,:), pp(:,:,:), zr(:,:,:),
     $  diag(:)
      integer i, iglob, it, j, k, iipole
      real*8  ggold(2), ggnew(2), gnorm(2), gg(2), alphacg(2), ene(2)
      real*8  zero, pt5, one, resnrm
      save    zero, pt5, one
      external matvec
      real*8, allocatable :: dipfield(:,:,:),dipfieldbis(:,:,:)
      real*8 term
c
c     MPI
c
      real*8, allocatable :: buffermpi1(:,:,:),buffermpi2(:,:,:)
      real*8, allocatable :: buffermpimu(:,:,:)
      integer ierr
      integer req1, req2, req3, req4
      integer,allocatable :: reqrecdirrec(:),reqrecdirsend(:)
      integer,allocatable :: reqendrec(:),reqendsend(:)
      integer,allocatable :: req2endrec(:),req2endsend(:)
      integer,allocatable :: reqrec(:),reqsend(:)
      integer,allocatable :: req2rec(:),req2send(:)
      integer status(MPI_STATUS_SIZE)
c
 1000 format(' cgiter converged after ',I3,' iterations.',/,
     $       ' final energy        = ',2D14.7,/,
     $       ' final residual norm = ',2D14.7)
 1010 format(' energy and residual norm at iteration ',I3,':',4D12.2)
 1020 format(' Conjugate gradient solver: induced dipoles',/,
     $  ' ipole       mux         muy         muz')
 1021 format(' Conjugate gradient solver: induced p-dipoles',/,
     $  ' ipole       mux         muy         muz')
 1030 format(i6,2x,f10.7,2x,f10.7,2x,f10.7)
 1040 format(' Using a diagonal preconditioner.')

      zero = 0.0d0
      one  = 1.0d0
      pt5  = 0.5d0
c
c
c     allocate some memory and setup the preconditioner:
c
      if (rank.le.ndir-1) then
        allocate (buffermpimu(3,nrhs,max(npoleloc,1)))
        buffermpimu = 0.0d0
        allocate (dipfield(3,nrhs,max(1,npoleloc)))
        dipfield = 0d0
        allocate (dipfieldbis(3,nrhs,max(1,npolerecloc)))
        dipfieldbis = 0d0
        allocate (res(3,nrhs,max(1,npoleloc)))
        allocate (zr(3,nrhs,max(1,npoleloc)))
        allocate (diag(npoleloc))
        allocate (h(3,nrhs,max(1,npolebloc)))
      else 
        allocate (buffermpimu(3,nrhs,max(1,npolerecloc)))
        buffermpimu = 0.0d0
        allocate (dipfieldbis(3,nrhs,max(1,npolerecloc)))
        dipfieldbis = 0d0
        allocate (dipfield(3,nrhs,max(1,npoleloc)))
        dipfield = 0d0
        allocate (diag(npoleloc))
      end if
      allocate (pp(3,nrhs,max(1,npolebloc)))
      allocate (buffermpi1(3,nrhs,max(npoleloc,1)))
      buffermpi1 = 0.0d0
      allocate (buffermpi2(3,nrhs,max(1,npolerecloc)))
      buffermpi2 = 0.0d0
      allocate (reqrecdirrec(nproc))
      allocate (reqrecdirsend(nproc))
      allocate (reqrec(nproc))
      allocate (reqsend(nproc))
      allocate (req2rec(nproc))
      allocate (req2send(nproc))
      allocate (reqendrec(nproc))
      allocate (reqendsend(nproc))
      allocate (req2endsend(nproc))
      allocate (req2endrec(nproc))
      if (precnd) then
        do i = 1, npoleloc
          iipole = poleglob(i)
          if (polarity(iipole).eq.0) then
            diag(i) = tinypol
          else
            diag(i) = polarity(iipole)
          end if
        end do
        if (polprt.ge.2.and.rank.eq.0) write (iout,1040)
      else
        diag = 1d0
      end if
c
c     initialize
c
      if (rank.le.ndir-1) then
        pp = 0.0d0
        h  = 0.0d0
        res  = 0.0d0
      end if
c
c     now, compute the initial direction
c
      ggold = 0.0d0
c
c     MPI : begin reception
c
      if (rank.le.ndir-1) then
        call commdirdir(nrhs,0,pp,reqrec,reqsend)
c
        call commrecdirsolv(nrhs,0,dipfieldbis,dipfield,buffermpi1,
     $      buffermpi2,reqrecdirrec,reqrecdirsend)
c
c     The PME processes compute the reciprocal matrix vector product
c
      else 
        call commrecdirdip(nrhs,0,murec,pp,buffermpimu,
     $   buffermpimu,req2rec,req2send)
c
        call tmatxbrecip(mu,murec,nrhs,dipfield,dipfieldbis)
        call commrecdirsolv(nrhs,1,dipfieldbis,dipfield,buffermpi1,
     $      buffermpi2,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,dipfieldbis,dipfield,buffermpi1,
     $      buffermpi2,reqrecdirrec,reqrecdirsend)
      end if
      if (rank.le.ndir-1) then
c
        call matvec(nrhs,.true.,mu,h)
        call commfield(nrhs,h)
c
        call commrecdirsolv(nrhs,2,dipfieldbis,dipfield,buffermpi1,
     $    buffermpi2,reqrecdirrec,reqrecdirsend)
        term = (4.0d0/3.0d0) * aewald**3 / sqrtpi

        do i = 1, npoleloc
          do k = 1, nrhs
            do j = 1, 3
              h(j,k,i) = h(j,k,i) + dipfield(j,k,i)-term*mu(j,k,i)
              res(j,k,i) = ef(j,k,i)-h(j,k,i)
              zr(j,k,i) = diag(i)*res(j,k,i)
              pp(j,k,i) = zr(j,k,i)
            end do
          end do
        end do

        ggold = 0d0
        do i = 1, npoleloc
          do k = 1, nrhs
            do j = 1, 3
              ggold(k) = ggold(k) + res(j,k,i)*zr(j,k,i)
            end do
          end do
        end do


        call MPI_IALLREDUCE(MPI_IN_PLACE,ggold(1),nrhs,MPI_REAL8,
     $     MPI_SUM,comm_dir,req1,ierr)
        call commdirdir(nrhs,1,pp,reqrec,reqsend)
c
        call commrecdirdip(nrhs,1,murec,pp,buffermpimu,
     $     buffermpimu,req2rec,req2send)
c
        call commdirdir(nrhs,2,mu,reqrec,reqsend)
c
        call commrecdirdip(nrhs,2,murec,pp,buffermpimu,
     $     buffermpimu,req2rec,req2send)
        call MPI_WAIT(req1,status,ierr)
c
      else 
        call commrecdirdip(nrhs,2,murec,pp,buffermpimu,
     $     buffermpimu,req2rec,req2send)
      end if
c
c     now, start the main loop:
c
      do it = 1, politer
        do k = 1, nrhs
          gg(k) = zero
          ggnew(k) = zero
        end do
c
c     MPI : begin reception
c
        if (rank.le.ndir-1) then
          call commdirdir(nrhs,0,pp,reqrec,reqsend)
        else
          call commrecdirdip(nrhs,0,murec,pp,buffermpimu,
     $    buffermpimu,req2rec,req2send)
        end if
c
        if (rank.le.ndir-1) then
          call commrecdirsolv(nrhs,0,dipfieldbis,dipfield,buffermpi1,
     $      buffermpi2,reqrecdirrec,reqrecdirsend)
c
        else 
          call tmatxbrecip(mu,murec,nrhs,dipfield,dipfieldbis)
          call commrecdirsolv(nrhs,1,dipfieldbis,dipfield,buffermpi1,
     $      buffermpi2,reqrecdirrec,reqrecdirsend)
          call commrecdirsolv(nrhs,2,dipfieldbis,dipfield,buffermpi1,
     $      buffermpi2,reqrecdirrec,reqrecdirsend)
        end if
        if (rank.le.ndir-1) then
c
          call matvec(nrhs,.true.,pp,h)
          call commfield(nrhs,h)
c
          call commrecdirsolv(nrhs,2,dipfieldbis,dipfield,buffermpi1,
     $      buffermpi2,reqrecdirrec,reqrecdirsend)
c
          do i = 1, npoleloc
            do k = 1, nrhs
              do j = 1, 3
                h(j,k,i) = h(j,k,i) + dipfield(j,k,i)-term*pp(j,k,i)
              end do
            end do
          end do
          gg = 0d0
          do i = 1, npoleloc
            do k = 1, nrhs
              do j = 1, 3
                gg(k) = gg(k) + pp(j,k,i)*h(j,k,i)
              end do
            end do
          end do

c          do k = 1, nrhs
c            gg(k) = sum(pp(:,k,1:npoleloc)*h(:,k,1:npoleloc))
c          end do
          call MPI_IALLREDUCE(MPI_IN_PLACE,gg(1),nrhs,MPI_REAL8,MPI_SUM,
     $      comm_dir,req2,ierr)
          call MPI_WAIT(req2,status,ierr)
          do k = 1, nrhs
            if (gg(k).eq.zero) return
            alphacg(k)  = ggold(k)/gg(k)
            ggnew(k)  = zero
            ene(k)    = zero
          end do

          do i = 1, npoleloc
            do k = 1, nrhs
              do j = 1, 3
                mu(j,k,i) = mu(j,k,i) + alphacg(k)*pp(j,k,i)
                res(j,k,i) = res(j,k,i)-alphacg(k)*h(j,k,i)
                zr(j,k,i) = diag(i)*res(j,k,i)
              end do
            end do
          end do

          do i = 1, npoleloc
            do k = 1, nrhs
              do j = 1, 3
                ggnew(k) = ggnew(k) + res(j,k,i)*zr(j,k,i)
                ene(k) = ene(k) - pt5*mu(j,k,i)*(res(j,k,i)+ef(j,k,i))
              end do
            end do
          end do
        end if
        call MPI_IALLREDUCE(MPI_IN_PLACE,ggnew(1),nrhs,MPI_REAL8,
     $      MPI_SUM,COMM_TINKER,req3,ierr)
        call MPI_WAIT(req3,status,ierr)
        if (rank.le.ndir-1) then
          call MPI_IALLREDUCE(MPI_IN_PLACE,ene(1),nrhs,MPI_REAL8,
     $      MPI_SUM,comm_dir,req4,ierr)
          call MPI_WAIT(req4,status,ierr)

          do i = 1, npoleloc
            do k = 1, nrhs
              do j = 1, 3
                pp(j,k,i) = zr(j,k,i)+ ggnew(k)/ggold(k)*pp(j,k,i)
              end do
            end do
          end do
c
          call commdirdir(nrhs,1,pp,reqrec,reqsend)
c
          call commrecdirdip(nrhs,1,murec,pp,buffermpimu,
     $    buffermpimu,req2rec,req2send)
c
          call commdirdir(nrhs,2,pp,reqrec,reqsend)
          call commrecdirdip(nrhs,2,murec,pp,buffermpimu,
     $    buffermpimu,req2rec,req2send)
        else 
          call commrecdirdip(nrhs,2,murec,pp,buffermpimu,
     $    buffermpimu,req2rec,req2send)
        end if
c
        ggold = ggnew
        resnrm = zero
        do k = 1, nrhs
          gnorm(k) = sqrt(ggnew(k)/dble(3*npolar))
          resnrm   = max(resnrm,gnorm(k))
        end do
        if (resnrm.lt.poleps) then
          if (polprt.ge.1.and.rank.eq.0) write(iout,1000) it,
     $      (ene(k)*coulomb, k = 1, nrhs), (gnorm(k), k = 1, nrhs)
          goto 10
        end if
        if (polprt.ge.2.and.rank.eq.0) write(iout,1010)
     $    it, (ene(k)*coulomb, gnorm(k), k = 1, nrhs)
      end do
 10   continue
c
      if (rank.le.ndir-1) then
        call commdirdir(nrhs,0,mu,reqendrec,reqendsend)
        call commdirdir(nrhs,1,mu,reqendrec,reqendsend)
        call commrecdirdip(nrhs,1,murec,mu,buffermpimu,
     $    buffermpimu,req2endrec,req2endsend)
c
        call commdirdir(nrhs,2,mu,reqendrec,reqendsend)
        call commrecdirdip(nrhs,2,murec,mu,buffermpimu,
     $    buffermpimu,req2endrec,req2endsend)
      else 
        call commrecdirdip(nrhs,0,murec,mu,buffermpimu,
     $    buffermpimu,req2endrec,req2endsend)
        call commrecdirdip(nrhs,2,murec,mu,buffermpimu,
     $    buffermpimu,req2endrec,req2endsend)
      end if
c
c     finalize and debug printing:
c
      if (polprt.ge.3) then
        write(iout,1020)
        do i = 1, npoleloc
          iglob = glob(i)
          write(iout,1030) iglob, (mu(j,1,i), j = 1, 3)
        end do
      end if
      if (polprt.ge.4) then
        write(iout,1021)
        do i = 1, npoleloc
          iglob = glob(i)
          write(iout,1030) iglob, (mu(j,2,i), j = 1, 3)
        end do
      end if
      deallocate (buffermpimu)
      deallocate (reqrecdirsend)
      deallocate (reqrecdirrec)
      deallocate (reqsend)
      deallocate (reqrec)
      deallocate (req2send)
      deallocate (req2rec)
      deallocate (reqendsend)
      deallocate (reqendrec)
      deallocate (req2endsend)
      deallocate (req2endrec)
      deallocate (diag)
      deallocate (pp)
      if (rank.le.ndir-1) then
        deallocate (buffermpi1)
        deallocate (res)
        deallocate (zr)
        deallocate (dipfield)
        deallocate (h)
      else
        deallocate (buffermpi2)
        deallocate (dipfieldbis)
      end if
      return
      end
c
      subroutine inducejac_pme(matvec,nrhs,dodiis,ef,mu,murec)
      use atmlst
      use domdec
      use ewald
      use iounit
      use math
      use mpole
      use polar
      use polpot
      use potent
      use timestat
      use mpi
      implicit none
c
c     solves the polarizable dipoles linear equations by jacobi
c     iterations coupled to the DIIS extrapolation.
c
      integer nrhs, info
      real*8  ef(3,nrhs,*), mu(3,nrhs,*), murec(3,nrhs,*)
      real*8, allocatable :: munew(:,:,:), h(:,:,:)
      real*8, allocatable :: xdiis(:,:), ediis(:,:), bmat(:,:),
     $   bloc(:,:), cex(:)
      integer i, j, k, iipole, it, ind, ndismx, nmat, lenb
      logical dodiis
      parameter (ndismx=25)
      integer ipiv(ndismx+1)
      real*8  zero, one, rnorm(2), rr, xx(1)
      real*8 term
      save    zero, one, xx
      external matvec
      real*8, allocatable :: dipfield(:,:,:),dipfieldbis(:,:,:)
c
c     MPI
c
      real*8, allocatable :: buffermpi1(:,:,:),buffermpi2(:,:,:)
      real*8, allocatable :: buffermpimu(:,:,:)
      integer ierr
      integer reqnorm,reqdiis(2*ndismx+1)
      integer,allocatable :: reqrecdirrec(:),reqrecdirsend(:)
      integer,allocatable :: reqrec(:),reqsend(:)
      integer,allocatable :: req2rec(:),req2send(:)
      integer status(MPI_STATUS_SIZE)
c
 1000 format(' itsolv converged after ',I3,' iterations.',/,
     $       ' final residual norm = ',3D14.7)
 1010 format(' residual norm at iteration ',I3,':',3D12.2)
 1020 format(' Jacobi/DIIS solver: induced dipoles',/,
     $  ' ipole       mux         muy         muz')
 1021 format(' Jacobi/DIIS solver: induced p-dipoles',/,
     $  ' ipole       mux         muy         muz')
 1030 format(i6,2x,f10.7,2x,f10.7,2x,f10.7)
c
c
      zero  = 0.0d0
      one   = 1.0d0
      xx(1) = 0.0d0

      if (rank.le.ndir-1) then
        allocate (buffermpimu(3,nrhs,max(npoleloc,1)))
        buffermpimu = 0.0d0
        allocate (munew(3,nrhs,max(1,npolebloc)))
        munew = 0.0d0
        allocate (h(3,nrhs,max(1,npolebloc)))
        h = 0.0d0
        if (dodiis) then
          nmat = 1
          lenb = ndismx + 1
          allocate (xdiis(3*nrhs*max(1,npoleloc),ndismx))
          allocate (ediis(3*nrhs*max(1,npoleloc),ndismx))
          allocate (bmat(lenb,lenb))
          allocate (bloc(lenb,lenb))
          allocate (cex(lenb))
          bmat = 0.0d0
        end if
      else
        allocate (munew(3,nrhs,max(1,npolebloc)))
        munew = 0.0d0
        allocate (buffermpimu(3,nrhs,max(1,npolerecloc)))
        buffermpimu = 0.0d0
      end if
      allocate (dipfieldbis(3,nrhs,max(1,npolerecloc)))
      dipfieldbis = 0d0
      allocate (dipfield(3,nrhs,max(1,npoleloc)))
      dipfield = 0d0
      allocate (buffermpi1(3,nrhs,max(npoleloc,1)))
      buffermpi1 = 0d0
      allocate (buffermpi2(3,nrhs,max(npolerecloc,1)))
      buffermpi2 = 0d0
      allocate (reqrecdirsend(nproc))
      allocate (reqrecdirrec(nproc))
      allocate (reqrec(nproc))
      allocate (reqsend(nproc))
      allocate (req2rec(nproc))
      allocate (req2send(nproc))
c
c     main loop:
c
      do it = 1, politer
        rnorm = 0.0d0
c
c     MPI : begin reception
c
        if (rank.le.ndir-1) then
          call commdirdir(nrhs,0,mu,reqrec,reqsend)
        else 
          call commrecdirdip(nrhs,0,dipfieldbis,dipfield,
     $      buffermpimu,buffermpimu,req2rec,req2send)
        end if
c
        if (rank.le.ndir-1) then
           call commrecdirsolv(nrhs,0,dipfieldbis,dipfield,buffermpi1,
     $     buffermpi2,reqrecdirrec,reqrecdirsend)
c
c     The PME processes compute the reciprocal matrix vector product
c
c        else
        else 
          call tmatxbrecip(mu,murec,nrhs,dipfield,dipfieldbis)
          call commrecdirsolv(nrhs,1,dipfieldbis,dipfield,buffermpi1,
     $      buffermpi2,reqrecdirrec,reqrecdirsend)
          call commrecdirsolv(nrhs,2,dipfieldbis,dipfield,buffermpi1,
     $      buffermpi2,reqrecdirrec,reqrecdirsend)
        end if
c
        if (rank.le.ndir-1) then
c
c     jacobi step:
c
c    The real space processes extract the recip fields, compute the real fields
c    and add them
c
         term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
         call matvec(nrhs,.false.,mu,h)
         call commfield(nrhs,h)
c
          call commrecdirsolv(nrhs,2,dipfieldbis,dipfield,buffermpi1,
     $      buffermpi2,reqrecdirrec,reqrecdirsend)
c
         do i = 1, npoleloc
           do k = 1, nrhs
             do j = 1, 3
               h(j,k,i) = h(j,k,i) + dipfield(j,k,i)-term*mu(j,k,i)
             end do
           end do
         end do
c
          do i = 1, npoleloc
            iipole = poleglob(i)
            do k = 1, nrhs
              do j = 1, 3
                munew(j,k,i) = polarity(iipole)*(ef(j,k,i) - h(j,k,i))
              end do
            end do
          end do
c
          do i = 1, npoleloc
            do k = 1, nrhs
              do j = 1,3 
                rnorm(k) = rnorm(k) + (munew(j,k,i)-mu(j,k,i))**2
              end do
            end do
          end do

c
        end if
        call MPI_IALLREDUCE(MPI_IN_PLACE,rnorm(1),nrhs,MPI_REAL8,
     $    MPI_SUM,COMM_TINKER,reqnorm,ierr)
        if (rank.le.ndir-1) then
          if (dodiis) then
            ind = 0
            do i = 1, npoleloc
              do k = 1, nrhs
                do j = 1, 3
                  ind = ind + 1
                  xdiis(ind,nmat) = munew(j,k,i)
                  ediis(ind,nmat) = munew(j,k,i) - mu(j,k,i)
                end do
              end do
            end do
c
c           Compute Pulay's Matrix and extrapolate
c
            call diis(ndismx,3*nrhs*npoleloc,xdiis,ediis,
     $        bmat,nmat,reqdiis,comm_dir)
            do i = 1, 2*nmat-3
              call MPI_WAIT(reqdiis(i),status,ierr)
            end do
            bloc = bmat
            cex = 0.0d0
            cex(1) = one
            call dgesv(nmat,1,bloc,ndismx+1,ipiv,cex,nmat,info)
            munew = 0.0d0
            call extrap(3*nrhs*npoleloc,nmat-1,xdiis,cex,munew)
          end if
c
          call commdirdir(nrhs,1,munew,reqrec,reqsend)
          call commrecdirdip(nrhs,1,murec,munew,buffermpimu,
     $    buffermpimu,req2rec,req2send)
c
          call commdirdir(nrhs,2,mu,reqrec,reqsend)
          call commrecdirdip(nrhs,2,murec,munew,buffermpimu,
     $    buffermpimu,req2rec,req2send)
          mu(:,:,1:npoleloc) = munew(:,:,1:npoleloc)
        else 
          call commrecdirdip(nrhs,2,murec,munew,buffermpimu,
     $    buffermpimu,req2rec,req2send)
        end if
c
c     compute the norm of the increment.
c
        call MPI_WAIT(reqnorm,status,ierr)
        rr = zero
        do k = 1, nrhs
         rnorm(k) = sqrt(rnorm(k)/dble(3*npolar))
c          rnorm(k) = sqrt(rnorm(k))
          rr = max(rnorm(k),rr)
        end do
        if (polprt.ge.2.and.rank.eq.0)
     $         write(6,1010) it, (rnorm(k), k = 1, nrhs)
        if (rr.lt.poleps) then
          if (polprt.ge.1.and.rank.eq.0)
     $       write(6,1000) it, (rnorm(k), k = 1, nrhs)
          goto 10
        end if
      end do
 10   continue
c
      if (polprt.ge.3) then
        write(iout,1020)
        do i = 1, npoleloc
          iipole = poleglob(i)
          write(iout,1030) iipole, (mu(j,1,i), j = 1, 3)
        end do
      end if
      if (polprt.ge.4) then
        write(iout,1021)
        do i = 1, npoleloc
          iipole = poleglob(i)
          write(iout,1030) iipole, (mu(j,2,i), j = 1, 3)
        end do
      end if
c
c     free the memory.
c
      deallocate (reqsend)
      deallocate (reqrec)
      deallocate (req2send)
      deallocate (req2rec)
      deallocate (reqrecdirsend)
      deallocate (reqrecdirrec)
      deallocate (buffermpimu)
      if (rank.le.ndir-1) then
        deallocate (munew)
        deallocate (h)
        deallocate (dipfield)
        if (dodiis) then
          deallocate (xdiis)
          deallocate (ediis)
          deallocate (bmat)
          deallocate (bloc)
          deallocate (cex)
        end if
      else
        deallocate (dipfieldbis)
      end if
      return
      end
c
c
      subroutine efld0_direct(nrhs,ef)
c
c     Compute the direct space contribution to the permanent electric field.
c      Also compute the "p" field, which is used to
c     compute the energy according to the AMOEBA force field.
c
      use atmlst
      use atoms
      use bound
      use chgpen
      use couple
      use cutoff
      use domdec
      use ewald
      use iounit
      use math
      use mpole
      use neigh
      use polar
      use polgrp
      use polpot
      use potent
      use shunt
      use mpi
      implicit none
      integer i,iglob,kglob,nrhs,ipoleloc,nnelst
      real*8  ef(3,nrhs,npolebloc)
      integer ii, j, k, kkk, iipole, kkpole,kbis
      real*8 xr,yr,zr
      real*8 r,r2,rr1,rr2
      real*8 rr3,rr5,rr7
      real*8 rr3i,rr5i,rr7i
      real*8 rr3k,rr5k,rr7k
      real*8 ci,dix,diy,diz
      real*8 qixx,qiyy,qizz
      real*8 qixy,qixz,qiyz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkyy,qkzz
      real*8 qkxy,qkxz,qkyz
      real*8 dir,dkr
      real*8 qix,qiy,qiz,qir
      real*8 qkx,qky,qkz,qkr
      real*8 fid(3), fip(3)
      real*8 fkd(3), fkp(3)
      real*8 cutoff2
      real*8 corei,corek
      real*8 vali,valk
      real*8 alphai,alphak
      real*8 dmp3,dmp5,dmp7
      real*8 dmpi(7),dmpk(7)
      real*8 dmpik(7),dmpe(7)
      real*8 scalek
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: pscale(:)
      logical shortrange
      character*11 mode
      character*80 :: RoutineName
c
 1000 format(' Warning, system moved too much since last neighbor list
     $  update, try lowering nlupdate')
c
      shortrange = use_polarshortreal
      if (shortrange) then 
         RoutineName='efld0_shortreal'
         mode = 'SHORTEWALD'
      else
         RoutineName='efld0_direct'
         mode = 'EWALD'
      endif
c
      call switch (mode)
      cutoff2 = cut2
c
      allocate (dscale(n))
      allocate (pscale(n))
      dscale = 1.0d0
      pscale = 1.0d0

c
      do ii = 1, npolelocnl
        iipole = poleglobnl(ii)
        iglob = ipole(iipole)
        i = loc(iglob)
        ipoleloc = poleloc(iipole)
        if ((i.eq.0).or.(i.gt.nbloc)) then
          write(iout,1000)
          cycle
        end if
        ci  = rpole(1,iipole)
        dix = rpole(2,iipole)
        diy = rpole(3,iipole)
        diz = rpole(4,iipole)
        qixx = rpole(5,iipole)
        qixy = rpole(6,iipole)
        qixz = rpole(7,iipole)
        qiyy = rpole(9,iipole)
        qiyz = rpole(10,iipole)
        qizz = rpole(13,iipole)
        if (use_chgpen) then
           corei = pcore(iipole)
           vali = pval(iipole)
           alphai = palpha(iipole)
        end if
c
c     set exclusion coefficients for connected atoms
c
        if (dpequal) then
           do j = 1, n12(iglob)
              pscale(i12(j,iglob)) = p2scale
              do k = 1, np11(iglob)
                 if (i12(j,iglob) .eq. ip11(k,iglob))
     &              pscale(i12(j,iglob)) = p2iscale
              end do
              dscale(i12(j,iglob)) = pscale(i12(j,iglob))
           end do
           do j = 1, n13(iglob)
              pscale(i13(j,iglob)) = p3scale
              do k = 1, np11(iglob)
                 if (i13(j,iglob) .eq. ip11(k,iglob))
     &              pscale(i13(j,iglob)) = p3iscale
              end do
              dscale(i13(j,iglob)) = pscale(i13(j,iglob))
           end do
           do j = 1, n14(iglob)
              pscale(i14(j,iglob)) = p4scale
              do k = 1, np11(iglob)
                  if (i14(j,iglob) .eq. ip11(k,iglob))
     &              pscale(i14(j,iglob)) = p4iscale
              end do
              dscale(i14(j,iglob)) = pscale(i14(j,iglob))
           end do
           do j = 1, n15(iglob)
              pscale(i15(j,iglob)) = p5scale
              do k = 1, np11(iglob)
                 if (i15(j,iglob) .eq. ip11(k,iglob))
     &              pscale(i15(j,iglob)) = p5iscale
              end do
              dscale(i15(j,iglob)) = pscale(i15(j,iglob))
           end do
        else
           do j = 1, n12(iglob)
              pscale(i12(j,iglob)) = p2scale
              do k = 1, np11(iglob)
                 if (i12(j,iglob) .eq. ip11(k,iglob))
     &              pscale(i12(j,iglob)) = p2iscale
              end do
           end do
           do j = 1, n13(iglob)
              pscale(i13(j,iglob)) = p3scale
              do k = 1, np11(iglob)
                 if (i13(j,iglob) .eq. ip11(k,iglob))
     &              pscale(i13(j,iglob)) = p3iscale
              end do
           end do
           do j = 1, n14(iglob)
              pscale(i14(j,iglob)) = p4scale
              do k = 1, np11(iglob)
                  if (i14(j,iglob) .eq. ip11(k,iglob))
     &              pscale(i14(j,iglob)) = p4iscale
              end do
           end do
           do j = 1, n15(iglob)
              pscale(i15(j,iglob)) = p5scale
              do k = 1, np11(iglob)
                 if (i15(j,iglob) .eq. ip11(k,iglob))
     &              pscale(i15(j,iglob)) = p5iscale
              end do
           end do
           do j = 1, np11(iglob)
              dscale(ip11(j,iglob)) = d1scale
           end do
           do j = 1, np12(iglob)
              dscale(ip12(j,iglob)) = d2scale
           end do
           do j = 1, np13(iglob)
              dscale(ip13(j,iglob)) = d3scale
           end do
           do j = 1, np14(iglob)
              dscale(ip14(j,iglob)) = d4scale
           end do
        end if
c
         if (shortrange) then
           nnelst = nshortelst(ii)
         else
           nnelst = nelst(ii)
         end if
        do kkk = 1, nnelst
          if (shortrange) then
            kkpole = shortelst(kkk,ii)
          else
            kkpole = elst(kkk,ii)
          end if
          kbis = poleloc(kkpole)
          kglob = ipole(kkpole)
          if ((kbis.eq.0).or.(kbis.gt.npolebloc)) then
            write(iout,1000)
            cycle
          end if
          xr = x(kglob) - x(iglob)
          yr = y(kglob) - y(iglob)
          zr = z(kglob) - z(iglob)
          if (use_bounds)  call image (xr,yr,zr)
          r2 = xr*xr + yr* yr + zr*zr
          if (r2.le.cutoff2) then
            r = sqrt(r2)
            rr1 = 1.0d0 / r
            rr2 = rr1 * rr1
            rr3 = rr2 * rr1
            rr5 = 3.0d0 * rr2 * rr3
            rr7 = 5.0d0 * rr2 * rr5
            ck   = rpole(1,kkpole)
            dkx  = rpole(2,kkpole)
            dky  = rpole(3,kkpole)
            dkz  = rpole(4,kkpole)
            qkxx = rpole(5,kkpole)
            qkxy = rpole(6,kkpole)
            qkxz = rpole(7,kkpole)
            qkyy = rpole(9,kkpole)
            qkyz = rpole(10,kkpole)
            qkzz = rpole(13,kkpole)
c
c     intermediates involving moments and separation distance
c
            dir = dix*xr + diy*yr + diz*zr
            qix = qixx*xr + qixy*yr + qixz*zr
            qiy = qixy*xr + qiyy*yr + qiyz*zr
            qiz = qixz*xr + qiyz*yr + qizz*zr
            qir = qix*xr + qiy*yr + qiz*zr
            dkr = dkx*xr + dky*yr + dkz*zr
            qkx = qkxx*xr + qkxy*yr + qkxz*zr
            qky = qkxy*xr + qkyy*yr + qkyz*zr
            qkz = qkxz*xr + qkyz*yr + qkzz*zr
            qkr = qkx*xr + qky*yr + qkz*zr
c
c     calculate real space Ewald error function damping
c
            call dampewald (7,r,r2,1.0d0,dmpe)
c
c     find the field components for Thole polarization damping
c
            if (use_thole) then
               call dampthole (iipole,kkpole,7,r,dmpik)
               scalek = dscale(kglob)
               dmp3 = dmpe(3) - (1.0d0-scalek*dmpik(3))*rr3
               dmp5 = dmpe(5) - (1.0d0-scalek*dmpik(5))*rr5
               dmp7 = dmpe(7) - (1.0d0-scalek*dmpik(7))*rr7
               fid(1) = -xr*(dmp3*ck-dmp5*dkr+dmp7*qkr)
     &                     - dmp3*dkx + 2.0d0*dmp5*qkx
               fid(2) = -yr*(dmp3*ck-dmp5*dkr+dmp7*qkr)
     &                     - dmp3*dky + 2.0d0*dmp5*qky
               fid(3) = -zr*(dmp3*ck-dmp5*dkr+dmp7*qkr)
     &                     - dmp3*dkz + 2.0d0*dmp5*qkz
               fkd(1) = xr*(dmp3*ci+dmp5*dir+dmp7*qir)
     &                     - dmp3*dix - 2.0d0*dmp5*qix
               fkd(2) = yr*(dmp3*ci+dmp5*dir+dmp7*qir)
     &                     - dmp3*diy - 2.0d0*dmp5*qiy
               fkd(3) = zr*(dmp3*ci+dmp5*dir+dmp7*qir)
     &                     - dmp3*diz - 2.0d0*dmp5*qiz
               scalek = pscale(kglob)
               dmp3 = dmpe(3) - (1.0d0-scalek*dmpik(3))*rr3
               dmp5 = dmpe(5) - (1.0d0-scalek*dmpik(5))*rr5
               dmp7 = dmpe(7) - (1.0d0-scalek*dmpik(7))*rr7
               fip(1) = -xr*(dmp3*ck-dmp5*dkr+dmp7*qkr)
     &                     - dmp3*dkx + 2.0d0*dmp5*qkx
               fip(2) = -yr*(dmp3*ck-dmp5*dkr+dmp7*qkr)
     &                     - dmp3*dky + 2.0d0*dmp5*qky
               fip(3) = -zr*(dmp3*ck-dmp5*dkr+dmp7*qkr)
     &                     - dmp3*dkz + 2.0d0*dmp5*qkz
               fkp(1) = xr*(dmp3*ci+dmp5*dir+dmp7*qir)
     &                     - dmp3*dix - 2.0d0*dmp5*qix
               fkp(2) = yr*(dmp3*ci+dmp5*dir+dmp7*qir)
     &                     - dmp3*diy - 2.0d0*dmp5*qiy
               fkp(3) = zr*(dmp3*ci+dmp5*dir+dmp7*qir)
     &                     - dmp3*diz - 2.0d0*dmp5*qiz
c
c     find the field components for charge penetration damping
c
            else if (use_chgpen) then
               corek = pcore(kkpole)
               valk = pval(kkpole)
               alphak = palpha(kkpole)
               call dampdir (r,alphai,alphak,dmpi,dmpk)
               scalek = dscale(kglob)
               rr3i = dmpe(3) - (1.0d0-scalek*dmpi(3))*rr3
               rr5i = dmpe(5) - (1.0d0-scalek*dmpi(5))*rr5
               rr7i = dmpe(7) - (1.0d0-scalek*dmpi(7))*rr7
               rr3k = dmpe(3) - (1.0d0-scalek*dmpk(3))*rr3
               rr5k = dmpe(5) - (1.0d0-scalek*dmpk(5))*rr5
               rr7k = dmpe(7) - (1.0d0-scalek*dmpk(7))*rr7
               rr3 = dmpe(3) - (1.0d0-scalek)*rr3
               fid(1) = -xr*(rr3*corek + rr3k*valk
     &                     - rr5k*dkr + rr7k*qkr)
     &                     - rr3k*dkx + 2.0d0*rr5k*qkx
               fid(2) = -yr*(rr3*corek + rr3k*valk
     &                     - rr5k*dkr + rr7k*qkr)
     &                     - rr3k*dky + 2.0d0*rr5k*qky
               fid(3) = -zr*(rr3*corek + rr3k*valk
     &                     - rr5k*dkr + rr7k*qkr)
     &                     - rr3k*dkz + 2.0d0*rr5k*qkz
               fkd(1) = xr*(rr3*corei + rr3i*vali
     &                     + rr5i*dir + rr7i*qir)
     &                     - rr3i*dix - 2.0d0*rr5i*qix
               fkd(2) = yr*(rr3*corei + rr3i*vali
     &                     + rr5i*dir + rr7i*qir)
     &                     - rr3i*diy - 2.0d0*rr5i*qiy
               fkd(3) = zr*(rr3*corei + rr3i*vali
     &                     + rr5i*dir + rr7i*qir)
     &                     - rr3i*diz - 2.0d0*rr5i*qiz
               scalek = pscale(kglob)
               rr3 = rr2 * rr1
               rr3i = dmpe(3) - (1.0d0-scalek*dmpi(3))*rr3
               rr5i = dmpe(5) - (1.0d0-scalek*dmpi(5))*rr5
               rr7i = dmpe(7) - (1.0d0-scalek*dmpi(7))*rr7
               rr3k = dmpe(3) - (1.0d0-scalek*dmpk(3))*rr3
               rr5k = dmpe(5) - (1.0d0-scalek*dmpk(5))*rr5
               rr7k = dmpe(7) - (1.0d0-scalek*dmpk(7))*rr7
               rr3 = dmpe(3) - (1.0d0-scalek)*rr3
               fip(1) = -xr*(rr3*corek + rr3k*valk
     &                     - rr5k*dkr + rr7k*qkr)
     &                     - rr3k*dkx + 2.0d0*rr5k*qkx
               fip(2) = -yr*(rr3*corek + rr3k*valk
     &                     - rr5k*dkr + rr7k*qkr)
     &                     - rr3k*dky + 2.0d0*rr5k*qky
               fip(3) = -zr*(rr3*corek + rr3k*valk
     &                     - rr5k*dkr + rr7k*qkr)
     &                     - rr3k*dkz + 2.0d0*rr5k*qkz
               fkp(1) = xr*(rr3*corei + rr3i*vali
     &                     + rr5i*dir + rr7i*qir)
     &                     - rr3i*dix - 2.0d0*rr5i*qix
               fkp(2) = yr*(rr3*corei + rr3i*vali
     &                     + rr5i*dir + rr7i*qir)
     &                     - rr3i*diy - 2.0d0*rr5i*qiy
               fkp(3) = zr*(rr3*corei + rr3i*vali
     &                     + rr5i*dir + rr7i*qir)
     &                     - rr3i*diz - 2.0d0*rr5i*qiz
            end if
c
c     increment the field at each site due to this interaction
c
            do j = 1, 3
               ef(j,1,ipoleloc) = ef(j,1,ipoleloc) + fid(j)
               ef(j,1,kbis) = ef(j,1,kbis) + fkd(j)
               ef(j,2,ipoleloc) = ef(j,2,ipoleloc) + fip(j)
               ef(j,2,kbis) = ef(j,2,kbis) + fkp(j)
            end do
          end if
        end do
c
c     reset exclusion coefficients for connected atoms
c
         if (dpequal) then
            do j = 1, n12(iglob)
               pscale(i12(j,iglob)) = 1.0d0
               dscale(i12(j,iglob)) = 1.0d0
            end do
            do j = 1, n13(iglob)
               pscale(i13(j,iglob)) = 1.0d0
               dscale(i13(j,iglob)) = 1.0d0
            end do
            do j = 1, n14(iglob)
               pscale(i14(j,iglob)) = 1.0d0
               dscale(i14(j,iglob)) = 1.0d0
            end do
            do j = 1, n15(iglob)
               pscale(i15(j,iglob)) = 1.0d0
               dscale(i15(j,iglob)) = 1.0d0
            end do
         else
            do j = 1, n12(iglob)
               pscale(i12(j,iglob)) = 1.0d0
            end do
            do j = 1, n13(iglob)
               pscale(i13(j,iglob)) = 1.0d0
            end do
            do j = 1, n14(iglob)
               pscale(i14(j,iglob)) = 1.0d0
            end do
            do j = 1, n15(iglob)
               pscale(i15(j,iglob)) = 1.0d0
            end do
            do j = 1, np11(iglob)
               dscale(ip11(j,iglob)) = 1.0d0
            end do
            do j = 1, np12(iglob)
               dscale(ip12(j,iglob)) = 1.0d0
            end do
            do j = 1, np13(iglob)
               dscale(ip13(j,iglob)) = 1.0d0
            end do
            do j = 1, np14(iglob)
               dscale(ip14(j,iglob)) = 1.0d0
            end do
         end if
      end do
c
      deallocate (dscale)
      deallocate (pscale)
      return
      end
c
c
      subroutine tmatxb_pme(nrhs,dodiag,mu,efi)
c
c     Compute the direct space contribution to the electric field due to the current value
c     of the induced dipoles
c
      use atmlst
      use atoms
      use bound
      use chgpen
      use couple
      use domdec
      use ewald
      use math
      use mpole
      use neigh
      use polar
      use polgrp
      use polpot
      use potent
      use shunt
      use mpi
      implicit none
      integer i,nrhs,iglob,kglob,iipole,kkpole,kkpoleloc,nnelst
      integer ipoleloc
      real*8  mu(3,nrhs,*), efi(3,nrhs,npolebloc)
      real*8 xr,yr,zr
      real*8 r,r2,rr1,rr2
      real*8 rr3,rr5
      real*8 rr3ik,rr5ik
      real*8 scalek
      real*8 dmp3,dmp5
      real*8 fid(3),fkd(3)
      real*8 fip(3),fkp(3)
      real*8 dmpik(7),dmpe(7)
      real*8 dlocal(6)
      real*8 duix,duiy,duiz
      real*8 puix,puiy,puiz
      real*8 dukx,duky,dukz
      real*8 pukx,puky,pukz
      real*8 alphai,alphak
      real*8 pol
      real*8, allocatable :: uscale(:)
      real*8, allocatable :: wscale(:)
      logical dodiag
      logical shortrange
      integer j, ii, kkk, irhs
      real*8  cutoff2
      character*11 mode
      character*80 :: RoutineName

      shortrange = use_polarshortreal
      if (shortrange) then 
         RoutineName='tmatxb_shortreal'
         mode = 'SHORTEWALD'
      else
         RoutineName='tmatxb_pme'
         mode = 'EWALD'
      endif
c
c     initialize the result vector
c
      efi = 0d0
c
c     perform dynamic allocation of some local arrays
c
      allocate (uscale(n))
      allocate (wscale(n))
c
c     set arrays needed to scale connected atom interactions
c
      uscale = 1.0d0
      wscale = 1.0d0
c
c     gather some parameters, then set up the damping factors.
c
      call switch (mode)

      cutoff2 = cut2
c
      do ii = 1, npolelocnl
        iipole   = poleglobnl(ii)
        iglob    = ipole     (iipole)
        i        = loc       (iglob)
        ipoleloc = poleloc(iipole)
        if (i.eq.0) cycle
        duix = mu(1,1,ipoleloc)
        duiy = mu(2,1,ipoleloc)
        duiz = mu(3,1,ipoleloc)
        puix = mu(1,2,ipoleloc)
        puiy = mu(2,2,ipoleloc)
        puiz = mu(3,2,ipoleloc)
        if (use_chgpen) then
           alphai = palpha(iipole)
        end if
c
c
c     set exclusion coefficients for connected atoms
c
        do j = 1, n12(iglob)
           wscale(i12(j,iglob)) = w2scale
        end do
        do j = 1, n13(iglob)
           wscale(i13(j,iglob)) = w3scale
        end do
        do j = 1, n14(iglob)
           wscale(i14(j,iglob)) = w4scale
        end do
        do j = 1, n15(iglob)
           wscale(i15(j,iglob)) = w5scale
        end do
        do j = 1, np11(iglob)
           uscale(ip11(j,iglob)) = u1scale
        end do
        do j = 1, np12(iglob)
           uscale(ip12(j,iglob)) = u2scale
        end do
        do j = 1, np13(iglob)
           uscale(ip13(j,iglob)) = u3scale
        end do
        do j = 1, np14(iglob)
           uscale(ip14(j,iglob)) = u4scale
        end do
        if (shortrange) then
          nnelst = nshortelst(ii)
        else
          nnelst = nelst(ii)
        end if
        do kkk = 1, nnelst
          if (shortrange) then
            kkpole = shortelst(kkk,ii)
          else
            kkpole = elst(kkk,ii)
          end if
          kglob = ipole(kkpole)
          kkpoleloc = poleloc(kkpole)
          if (kkpoleloc.eq.0) cycle
          xr = x(kglob) - x(iglob)
          yr = y(kglob) - y(iglob)
          zr = z(kglob) - z(iglob)
          if (use_bounds)  call image (xr,yr,zr)
          r2 = xr*xr + yr* yr + zr*zr
          if (r2 .le. off2) then
c
c     compute the distances and the scaling factors according to
c     Thole's model.
c
            r = sqrt(r2)
            rr1 = 1.0d0 / r
            rr2 = rr1 * rr1
            rr3 = rr2 * rr1
            rr5 = 3.0d0 * rr2 * rr3
            dukx = mu(1,1,kkpoleloc)
            duky = mu(2,1,kkpoleloc)
            dukz = mu(3,1,kkpoleloc)
            pukx = mu(1,2,kkpoleloc)
            puky = mu(2,2,kkpoleloc)
            pukz = mu(3,2,kkpoleloc)
            if (use_chgpen) then
               alphak = palpha(kkpole)
            end if
c
c     calculate real space Ewald error function damping
c
            call dampewald (7,r,r2,1.0d0,dmpe)
c
c     find the field components for Thole polarization damping
c
            if (use_thole) then
               call dampthole2 (iipole,kkpole,5,r,dmpik)
               scalek = uscale(kglob)
               dmp3 = dmpe(3) - (1.0d0-scalek*dmpik(3))*rr3
               dmp5 = dmpe(5) - (1.0d0-scalek*dmpik(5))*rr5
               dlocal(1) = -dmp3 + dmp5*xr*xr
               dlocal(2) = dmp5*xr*yr
               dlocal(3) = dmp5*xr*zr
               dlocal(4) = -dmp3 + dmp5*yr*yr
               dlocal(5) = dmp5*yr*zr
               dlocal(6) = -dmp3 + dmp5*zr*zr
c
c     find the field components for charge penetration damping
c
            else if (use_chgpen) then
               call dampmut (r,alphai,alphak,dmpik)
               scalek = wscale(kglob)
               rr3ik = dmpe(3) - (1.0d0-scalek*dmpik(3))*rr3
               rr5ik = dmpe(5) - (1.0d0-scalek*dmpik(5))*rr5
               dlocal(1) = -rr3ik + rr5ik*xr*xr
               dlocal(2) = rr5ik*xr*yr
               dlocal(3) = rr5ik*xr*zr
               dlocal(4) = -rr3ik + rr5ik*yr*yr
               dlocal(5) = rr5ik*yr*zr
               dlocal(6) = -rr3ik + rr5ik*zr*zr
            end if
            fid(1) = dlocal(1)*dukx+dlocal(2)*duky+dlocal(3)*dukz
            fid(2) = dlocal(2)*dukx+dlocal(4)*duky+dlocal(5)*dukz
            fid(3) = dlocal(3)*dukx+dlocal(5)*duky+dlocal(6)*dukz
            fkd(1) = dlocal(1)*duix+dlocal(2)*duiy+dlocal(3)*duiz
            fkd(2) = dlocal(2)*duix+dlocal(4)*duiy+dlocal(5)*duiz
            fkd(3) = dlocal(3)*duix+dlocal(5)*duiy+dlocal(6)*duiz

            fip(1) = dlocal(1)*pukx+dlocal(2)*puky+dlocal(3)*pukz
            fip(2) = dlocal(2)*pukx+dlocal(4)*puky+dlocal(5)*pukz
            fip(3) = dlocal(3)*pukx+dlocal(5)*puky+dlocal(6)*pukz
            fkp(1) = dlocal(1)*puix+dlocal(2)*puiy+dlocal(3)*puiz
            fkp(2) = dlocal(2)*puix+dlocal(4)*puiy+dlocal(5)*puiz
            fkp(3) = dlocal(3)*puix+dlocal(5)*puiy+dlocal(6)*puiz
            do j = 1, 3
               efi(j,1,ipoleloc) = efi(j,1,ipoleloc) - fid(j)
               efi(j,1,kkpoleloc) = efi(j,1,kkpoleloc) - fkd(j)
               efi(j,2,ipoleloc) = efi(j,2,ipoleloc) - fip(j)
               efi(j,2,kkpoleloc) = efi(j,2,kkpoleloc) - fkp(j)
            end do
          end if
        end do
c
c     reset interaction scaling coefficients for connected atoms
c
        do j = 1, n12(iglob)
           wscale(i12(j,iglob)) = 1.0d0
        end do
        do j = 1, n13(iglob)
           wscale(i13(j,iglob)) = 1.0d0
        end do
        do j = 1, n14(iglob)
           wscale(i14(j,iglob)) = 1.0d0
        end do
        do j = 1, n15(iglob)
           wscale(i15(j,iglob)) = 1.0d0
        end do
        do j = 1, np11(iglob)
           uscale(ip11(j,iglob)) = 1.0d0
        end do
        do j = 1, np12(iglob)
           uscale(ip12(j,iglob)) = 1.0d0
        end do
        do j = 1, np13(iglob)
           uscale(ip13(j,iglob)) = 1.0d0
        end do
        do j = 1, np14(iglob)
           uscale(ip14(j,iglob)) = 1.0d0
        end do
      end do
      if(dodiag) then
c
c     if dodiag is true, also compute the "self-induced" field,
c     i.e., the diagonal portion of the matrix/vector product.
c
        do i = 1, npoleloc
          iipole = poleglob(i)
c
c     if no polarisability, take a negligeable value to allow convergence
c
          if (polarity(iipole).eq.0.0d0) then
             pol = tinypol ** -1
          else
             pol  = polarity(iipole) ** -1
          endif
          do irhs = 1, nrhs
            do j = 1, 3
              efi(j,irhs,i) = efi(j,irhs,i) +
     $           mu(j,irhs,i)*pol
            end do
          end do
        end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (uscale)
      deallocate (wscale)
      return
      end

      subroutine efld0_recip(cphi)
c
c     Compute the reciprocal space contribution to the electric field due to the permanent 
c     multipoles
c
      use atmlst
      use bound
      use boxes
      use domdec
      use ewald
      use fft
      use math
      use mpole
      use pme
      use potent
      use timestat
      use mpi
      implicit none
      integer ierr
      integer status(MPI_STATUS_SIZE),tag
      integer i,j,k,iglob,iipole,iloc
      integer k1,k2,k3
      integer m1,m2,m3
      integer nf1,nf2,nf3
      real*8 r1,r2,r3
      real*8 h1,h2,h3
      real*8 volterm,denom
      real*8 hsq,expterm
      real*8 term,pterm
      real*8 cmp(10),fmp(10),cphi(10,*)
      real*8, allocatable :: qgridmpi(:,:,:,:,:)
      integer, allocatable :: reqrec(:),reqsend(:)
      integer, allocatable :: reqbcastrec(:),reqbcastsend(:)
      integer nprocloc,commloc,rankloc
c
      if (use_pmecore) then
        nprocloc = nrec
        commloc  = comm_rec
        rankloc  = rank_bis
      else
        nprocloc = nproc
        commloc  = COMM_TINKER
        rankloc  = rank
      end if
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return

      allocate (reqrec(nprocloc))
      allocate (reqsend(nprocloc))
      allocate (reqbcastrec(nprocloc))
      allocate (reqbcastsend(nprocloc))
      allocate (qgridmpi(2,n1mpimax,n2mpimax,n3mpimax,nrec_recep))
c
c     MPI : Begin reception
c
      do i = 1, nrec_recep
        tag = nprocloc*rankloc + prec_recep(i) + 1
        call MPI_IRECV(qgridmpi(1,1,1,1,i),2*n1mpimax*n2mpimax*
     $   n3mpimax,MPI_REAL8,prec_recep(i),tag,
     $   commloc,reqrec(i),ierr)
      end do
c
c     zero out the PME charge grid
      qgridin_2d = 0d0
c
c     fill the pme grid, loop over the multipoles sites
c
      do i = 1, npolerecloc
         iipole = polerecglob(i)
         iglob = ipole(iipole)
         cmp(1) = rpole(1,iipole)
         cmp(2) = rpole(2,iipole)
         cmp(3) = rpole(3,iipole)
         cmp(4) = rpole(4,iipole)
         cmp(5) = rpole(5,iipole)
         cmp(6) = rpole(9,iipole)
         cmp(7) = rpole(13,iipole)
         cmp(8) = 2.0d0 * rpole(6,iipole)
         cmp(9) = 2.0d0 * rpole(7,iipole)
         cmp(10) = 2.0d0 * rpole(10,iipole)
c
c     compute B-spline coefficients
c
         call bspline_fill_site(iglob,i)
c
c     convert Cartesian multipoles to fractional coordinates
c
        call cmp_to_fmp_site(cmp,fmp)
c
c     assign PME grid
c
        call grid_mpole_site(iglob,i,fmp)
      end do
c
c     MPI : begin sending
c
      do i = 1, nrec_send
        tag = nprocloc*prec_send(i) + rankloc + 1
        call MPI_ISEND(qgridin_2d(1,1,1,1,i+1),
     $   2*n1mpimax*n2mpimax*n3mpimax,MPI_REAL8,
     $   prec_send(i),tag,commloc,reqsend(i),ierr)
      end do
c
      do i = 1, nrec_recep
        call MPI_WAIT(reqrec(i),status,ierr)
      end do
      do i = 1, nrec_send
        call MPI_WAIT(reqsend(i),status,ierr)
      end do
c
c     do the reduction 'by hand'
c
      do i = 1, nrec_recep
        call aadd(2*n1mpimax*n2mpimax*n3mpimax,
     $   qgridin_2d(1,1,1,1,1),
     $   qgridmpi(1,1,1,1,i),qgridin_2d(1,1,1,1,1))
      end do
c
c     Perform 3-D FFT forward transform
c
      call fft2d_frontmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
     $ n3mpimax)
c
c     make the scalar summation over reciprocal lattice
c
      if ((istart2(rankloc+1).eq.1).and.(jstart2(rankloc+1).eq.1).and.
     $   (kstart2(rankloc+1).eq.1)) then
           qfac_2d(1,1,1) = 0.0d0
      end if
      pterm = (pi/aewald)**2
      volterm = pi * volbox
      nf1 = (nfft1+1) / 2
      nf2 = (nfft2+1) / 2
      nf3 = (nfft3+1) / 2
      do k3 = kstart2(rankloc+1),kend2(rankloc+1)
        do k2 = jstart2(rankloc+1),jend2(rankloc+1)
          do k1 = istart2(rankloc+1),iend2(rankloc+1)
            m1 = k1 - 1
            m2 = k2 - 1
            m3 = k3 - 1
            if (k1 .gt. nf1)  m1 = m1 - nfft1
            if (k2 .gt. nf2)  m2 = m2 - nfft2
            if (k3 .gt. nf3)  m3 = m3 - nfft3
            if ((m1.eq.0).and.(m2.eq.0).and.(m3.eq.0)) goto 10
            r1 = dble(m1)
            r2 = dble(m2)
            r3 = dble(m3)
            h1 = recip(1,1)*r1 + recip(1,2)*r2 + recip(1,3)*r3
            h2 = recip(2,1)*r1 + recip(2,2)*r2 + recip(2,3)*r3
            h3 = recip(3,1)*r1 + recip(3,2)*r2 + recip(3,3)*r3
            hsq = h1*h1 + h2*h2 + h3*h3
            term = -pterm * hsq
            expterm = 0.0d0
            if ((term .gt. -50.0d0)) then
               denom = volterm*hsq*bsmod1(k1)*bsmod2(k2)*bsmod3(k3)
               expterm = exp(term) / denom
               if (.not. use_bounds) then
                  expterm = expterm * (1.0d0-cos(pi*xbox*sqrt(hsq)))
               else if (octahedron) then
                  if (mod(m1+m2+m3,2) .ne. 0)  expterm = 0.0d0
               end if
            qfac_2d(k1-istart2(rankloc+1)+1,k2-jstart2(rankloc+1)+1,
     $        k3-kstart2(rankloc+1)+1) = expterm
            end if
 10         continue
          end do
        end do
      end do
c
c     account for the zeroth grid point for a finite system
c
      if ((istart2(rankloc+1).eq.1).and.(jstart2(rankloc+1).eq.1).and.
     $   (kstart2(rankloc+1).eq.1)) then
        if (.not. use_bounds) then
           expterm = 0.5d0 * pi / xbox
           qfac_2d(1,1,1) = expterm
        end if
      end if
c
c     complete the transformation of the charge grid
c
      do k = 1, ksize2(rankloc+1)
         do j = 1, jsize2(rankloc+1)
           do i = 1, isize2(rankloc+1)
             term = qfac_2d(i,j,k)
             qgridout_2d(1,i,j,k) = term*qgridout_2d(1,i,j,k)
             qgridout_2d(2,i,j,k) = term*qgridout_2d(2,i,j,k)
           end do
         end do
      end do
c
c     perform 3-D FFT backward transform
c
      call fft2d_backmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
     $ n3mpimax)
c
c     MPI : Begin reception
c
      do i = 1, nrec_send
        tag = nprocloc*rankloc + prec_send(i) + 1
        call MPI_IRECV(qgridin_2d(1,1,1,1,i+1),
     $   2*n1mpimax*n2mpimax*n3mpimax,MPI_REAL8,
     $   prec_send(i),tag,commloc,reqbcastrec(i),ierr)
      end do
c
c     MPI : begin sending
c
      do i = 1, nrec_recep
        tag = nprocloc*prec_recep(i) + rankloc + 1
        call MPI_ISEND(qgridin_2d(1,1,1,1,1),
     $   2*n1mpimax*n2mpimax*n3mpimax,MPI_REAL8,
     $   prec_recep(i),tag,commloc,reqbcastsend(i),ierr)
      end do
c
      do i = 1, nrec_send
        call MPI_WAIT(reqbcastrec(i),status,ierr)
      end do
      do i = 1, nrec_recep
        call MPI_WAIT(reqbcastsend(i),status,ierr)
      end do
c
c     get field
c
      do i = 1, npolerecloc
        iipole = polerecglob(i)
        iglob = ipole(iipole)
        iloc = poleloc(iipole)
        call fphi_mpole_site(iglob,i)
        call fphi_to_cphi_site(fphirec(1,i),cphirec(1,i))
        if (.not.(use_pmecore).and.(repart(iglob).eq.rankloc)) then
          call amove(10,cphirec(1,i),cphi(1,iloc))
        end if
      end do
c
      deallocate (qgridmpi)
      deallocate (reqbcastrec)
      deallocate (reqbcastsend)
      deallocate (reqrec)
      deallocate (reqsend)
      return
      end
c
c
      subroutine tmatxbrecip(mu,murec,nrhs,dipfield,dipfieldbis)
c
c     Compute the reciprocal space contribution to the electric field due to the current  
c     value of the induced dipoles
c
      use atmlst
      use boxes
      use domdec
      use ewald
      use fft
      use math
      use mpole
      use pme
      use potent
      use mpi
      implicit none
      integer ierr,iglob
      integer status(MPI_STATUS_SIZE),tag
      integer nrhs,iipole
      integer i,j,k,iloc
      real*8 fuind(3),fuinp(3)
      real*8 term
      real*8 a(3,3)
      real*8 fdip_phi1(10), fdip_phi2(10), fdip_sum_phi(20)
      real*8 dipfield(3,nrhs,*),dipfieldbis(3,nrhs,*)
      real*8 mu(3,nrhs,*),murec(3,nrhs,*)
      real*8, allocatable :: qgridmpi(:,:,:,:,:)
      integer, allocatable :: reqbcastrec(:),reqbcastsend(:)
      integer, allocatable :: reqrec(:),reqsend(:)
      integer nprocloc,commloc,rankloc
c
      if (use_pmecore) then
        nprocloc = nrec
        commloc  = comm_rec
        rankloc  = rank_bis
      else
        nprocloc = nproc
        commloc  = COMM_TINKER
        rankloc  = rank
      end if

      allocate (reqbcastrec(nprocloc))
      allocate (reqbcastsend(nprocloc))
      allocate (reqrec(nprocloc))
      allocate (reqsend(nprocloc))
      allocate (qgridmpi(2,n1mpimax,n2mpimax,n3mpimax,nrec_recep))
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
c
c     MPI : Begin reception
c
      do i = 1, nrec_recep
        tag = nprocloc*rankloc + prec_recep(i) + 1
        call MPI_IRECV(qgridmpi(1,1,1,1,i),2*n1mpimax*n2mpimax*n3mpimax,
     $   MPI_REAL8,prec_recep(i),tag,
     $   commloc,reqrec(i),ierr)
      end do
c
c     zero out the PME charge grid
      qgrid2in_2d = 0d0
c
c     fill the pme grid, loop over the induced dipoles sites
c
      do j = 1, 3
        a(1,j) = dble(nfft1) * recip(j,1)
        a(2,j) = dble(nfft2) * recip(j,2)
        a(3,j) = dble(nfft3) * recip(j,3)
      end do
      do i = 1, npolerecloc
        iipole = polerecglob(i)
        iglob = ipole(iipole)
        iloc  = poleloc(iipole)
c
c       Convert cartesian dipoles to fractional coordinates
c
        if (repart(iglob).ne.rank) then
          do k = 1, 3
             fuind(k) = a(k,1)*murec(1,1,i) + a(k,2)*murec(2,1,i)
     &                       + a(k,3)*murec(3,1,i)
             fuinp(k) = a(k,1)*murec(1,2,i) + a(k,2)*murec(2,2,i)
     &                       + a(k,3)*murec(3,2,i)
          end do
        else
          do k = 1, 3
             fuind(k) = a(k,1)*mu(1,1,iloc) + a(k,2)*mu(2,1,iloc)
     &                       + a(k,3)*mu(3,1,iloc)
             fuinp(k) = a(k,1)*mu(1,2,iloc) + a(k,2)*mu(2,2,iloc)
     &                       + a(k,3)*mu(3,2,iloc)
          end do
        end if
c
c     assign PME grid
c
        call grid_uind_site(iglob,i,fuind,fuinp,qgrid2in_2d)
      end do
c
c     MPI : begin sending
c
      do i = 1, nrec_send
        tag = nprocloc*prec_send(i) + rankloc + 1
        call MPI_ISEND(qgrid2in_2d(1,1,1,1,i+1),
     $   2*n1mpimax*n2mpimax*n3mpimax,MPI_REAL8,
     $   prec_send(i),tag,commloc,reqsend(i),ierr)
      end do
c
      do i = 1, nrec_recep
        call MPI_WAIT(reqrec(i),status,ierr)
      end do
      do i = 1, nrec_send
        call MPI_WAIT(reqsend(i),status,ierr)
      end do
c
c     do the reduction 'by hand'
c
      do i = 1, nrec_recep
        call aadd(2*n1mpimax*n2mpimax*n3mpimax,
     $   qgrid2in_2d(1,1,1,1,1),qgridmpi(1,1,1,1,i),
     $   qgrid2in_2d(1,1,1,1,1))
      end do
c
c     Perform 3-D FFT forward transform
c
      call fft2d_frontmpi(qgrid2in_2d,qgrid2out_2d,n1mpimax,n2mpimax,
     $ n3mpimax)
c
c     complete the transformation of the charge grid
c
      do k = 1, ksize2(rankloc+1)
         do j = 1, jsize2(rankloc+1)
           do i = 1, isize2(rankloc+1)
             term = qfac_2d(i,j,k)
             qgrid2out_2d(1,i,j,k) = term*qgrid2out_2d(1,i,j,k)
             qgrid2out_2d(2,i,j,k) = term*qgrid2out_2d(2,i,j,k)
           end do
         end do
      end do
c
c     perform 3-D FFT backward transform
c
      call fft2d_backmpi(qgrid2in_2d,qgrid2out_2d,n1mpimax,n2mpimax,
     $ n3mpimax)
c
c     MPI : Begin reception
c
      do i = 1, nrec_send
        tag = nprocloc*rankloc + prec_send(i) + 1
        call MPI_IRECV(qgrid2in_2d(1,1,1,1,i+1),
     $   2*n1mpimax*n2mpimax*n3mpimax,MPI_REAL8,
     $   prec_send(i),tag,commloc,reqbcastrec(i),ierr)
      end do
c
c     MPI : begin sending
c
      do i = 1, nrec_recep
        tag = nprocloc*prec_recep(i) + rankloc + 1
        call MPI_ISEND(qgrid2in_2d,
     $   2*n1mpimax*n2mpimax*n3mpimax,MPI_REAL8,
     $   prec_recep(i),tag,commloc,reqbcastsend(i),ierr)
      end do
c
      do i = 1, nrec_send
        call MPI_WAIT(reqbcastrec(i),status,ierr)
      end do
      do i = 1, nrec_recep
        call MPI_WAIT(reqbcastsend(i),status,ierr)
      end do
c
c     get fields
c
      do i = 1, npolerecloc
        iipole = polerecglob(i)
        iglob = ipole(iipole)
        iloc  = poleloc(iipole)
        call fphi_uind_site(iglob,i,fdip_phi1,fdip_phi2,fdip_sum_phi)
        if (.not.(use_pmecore).and.(repart(iglob).eq.rankloc)) then
c
c     convert the dipole fields from fractional to Cartesian
c
          do k = 1, 3
             dipfield(k,1,iloc) = a(k,1)*fdip_phi1(2)
     &                           + a(k,2)*fdip_phi1(3)
     &                           + a(k,3)*fdip_phi1(4)
             dipfield(k,2,iloc) = a(k,1)*fdip_phi2(2)
     &                           + a(k,2)*fdip_phi2(3)
     $                           + a(k,3)*fdip_phi2(4)
          end do
        else
          do k = 1, 3
             dipfieldbis(k,1,i) = a(k,1)*fdip_phi1(2)
     &                           + a(k,2)*fdip_phi1(3)
     &                           + a(k,3)*fdip_phi1(4)
             dipfieldbis(k,2,i) = a(k,1)*fdip_phi2(2)
     &                           + a(k,2)*fdip_phi2(3)
     $                           + a(k,3)*fdip_phi2(4)
          end do
        end if
      end do
      deallocate (qgridmpi)
      deallocate (reqrec)
      deallocate (reqsend)
      deallocate (reqbcastsend)
      deallocate (reqbcastrec)
      return
      end
