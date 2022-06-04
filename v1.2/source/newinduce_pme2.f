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
      subroutine newinduce_pme2
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
c     without separate cores for reciprocal part
c
      integer i, j, k, nrhs
c
c     MPI
c
      integer iglob, iipole
      integer, allocatable :: reqrecdirsend(:),reqrecdirrec(:)
      integer, allocatable :: reqsend(:),reqrec(:)
      integer, allocatable :: req2send(:),req2rec(:)


c
      parameter (nrhs=2)
      real*8  wtime0, wtime1, wtime2, udsum, upsum
      real*8  term
      real*8, allocatable :: ef(:,:,:), mu(:,:,:), murec(:,:,:)
      real*8, allocatable :: cphi(:,:)
c
      real*8, allocatable :: buffermpi1(:,:),buffermpi2(:,:)
      real*8, allocatable :: buffermpimu1(:,:,:),buffermpimu2(:,:,:)
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
      allocate (mu(3,nrhs,max(1,npolebloc)))
      mu = 0d0
      allocate (murec(3,nrhs,max(1,npolerecloc)))
      murec = 0d0
c
      allocate (buffermpi1(10,max(npoleloc,1)))
      buffermpi1 = 0d0
      allocate (buffermpi2(10,max(npolerecloc,1)))
      buffermpi2 = 0d0
      allocate (buffermpimu1(3,nrhs,max(npoleloc,1)))
      buffermpimu1 = 0d0
      allocate (buffermpimu2(3,nrhs,max(npolerecloc,1)))
      buffermpimu2 = 0d0
c
      allocate (ef(3,nrhs,max(1,npolebloc)))
      ef = 0d0
      allocate (cphi(10,max(npoleloc,1)))
      cphi = 0d0
      if (allocated(cphirec)) deallocate (cphirec)
      allocate (cphirec(10,max(npolerecloc,1)))
      cphirec = 0d0
      if (allocated(fphirec)) deallocate(fphirec)
      allocate (fphirec(20,max(npolerecloc,1)))
      fphirec = 0d0
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
c
c    compute the reciprocal space contribution (fields)
c
      call efld0_recip(cphi)
c
      call commrecdirfields(0,cphirec,cphi,buffermpi1,buffermpi2,
     $ reqrecdirrec,reqrecdirsend)
      call commrecdirfields(1,cphirec,cphi,buffermpi1,buffermpi2,
     $ reqrecdirrec,reqrecdirsend)
      call commrecdirfields(2,cphirec,cphi,buffermpi1,buffermpi2,
     $ reqrecdirrec,reqrecdirsend)
c
      call efld0_direct(nrhs,ef)
c
      call commfield(nrhs,ef)
c
      call commdirdir(nrhs,0,mu,reqrec,reqsend)
c
c     Add direct and reciprocal fields
c
      do i = 1, npoleloc
        iipole = poleglob(i)
        do j = 1, 3
          ef(j,1,i)  = ef(j,1,i) - cphi(j+1,i) 
          ef(j,2,i)  = ef(j,2,i) - cphi(j+1,i) 
        end do
      end do
c
      term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
      do i = 1, npoleloc
        iipole = poleglob(i)
        do j = 1, 3
          ef(j,1,i)  = ef(j,1,i) + term*rpole(j+1,iipole)
          ef(j,2,i)  = ef(j,2,i) + term*rpole(j+1,iipole)
        end do
      end do
      wtime1 = mpi_wtime()
c
c     guess the dipoles:
c
c     predicted values for always stable predictor-corrector method
c
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
      call commdirdir(nrhs,2,mu,reqrec,reqsend)
      call commrecdirdip(nrhs,0,murec,mu,buffermpimu1,buffermpimu2,
     $ req2rec,req2send)
      call commrecdirdip(nrhs,1,murec,mu,buffermpimu1,buffermpimu2,
     $ req2rec,req2send)
      call commrecdirdip(nrhs,2,murec,mu,buffermpimu1,buffermpimu2,
     $ req2rec,req2send)
c
c     now, call the proper solver.
c
      if (polalg.eq.1) then
        call inducepcg_pme2(tmatxb_pme,nrhs,.true.,ef,mu,murec)
      else if (polalg.eq.2) then
        call inducejac_pme2(tmatxb_pme,nrhs,.true.,ef,mu,murec)
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
c     move the computed dipoles in the common block.
c
      do i = 1, npolebloc
        iipole = poleglob(i)
        do j = 1, 3
          uind(j,iipole) = mu(j,1,i)
          uinp(j,iipole) = mu(j,2,i)
        end do
      end do
      do i = 1, npolerecloc
        iipole = polerecglob(i)
        iglob = ipole(iipole)
        do j = 1, 3
          if (repart(iglob).ne.rank) then
            uind(j,iipole) = murec(j,1,i)
            uinp(j,iipole) = murec(j,2,i)
          else
            uind(j,iipole) = mu(j,1,poleloc(iipole))
            uinp(j,iipole) = mu(j,2,poleloc(iipole))
          end if
        end do
      end do
c
c     update the lists of previous induced dipole values
c
      if ((use_pred).and..not.(use_lambdadyn)) then
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
      deallocate (reqrec)
      deallocate (reqsend)
      deallocate (req2rec)
      deallocate (req2send)
      deallocate (reqrecdirrec)
      deallocate (reqrecdirsend)
      deallocate (buffermpi1)
      deallocate (buffermpi2)
      deallocate (buffermpimu1)
      deallocate (buffermpimu2)
      deallocate (ef)
c      deallocate (fphi)
      deallocate (mu)
      deallocate (murec)
      deallocate (cphi)
c      deallocate (cphirec)
      return
      end
c
      subroutine inducepcg_pme2(matvec,nrhs,precnd,ef,mu,murec)
      use atmlst
      use domdec
      use ewald
      use iounit
      use math
      use mpole
      use polar
      use polpot
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
      integer i, it, j, k
      real*8  ggold(2), ggnew(2), gnorm(2), gg(2), alphacg(2), ene(2)
      real*8  zero, pt5, one, resnrm, term
      save    zero, pt5, one
      data    zero/0.0d0/, pt5/0.50d0/, one/1.0d0/
      external matvec
      real*8, allocatable :: dipfield(:,:,:),dipfieldbis(:,:,:)
c
c     MPI
c
      real*8, allocatable :: buffermpi1(:,:,:),buffermpi2(:,:,:)
      real*8, allocatable :: buffermpimu1(:,:,:),buffermpimu2(:,:,:)
      integer iglob, iipole, ierr
      integer req1, req2, req3, req4
      integer status(MPI_STATUS_SIZE)
      integer,allocatable :: reqrecdirrec(:),reqrecdirsend(:)
      integer,allocatable :: reqrec(:),reqsend(:)
      integer,allocatable :: req2rec(:),req2send(:)
      integer,allocatable :: reqendrec(:),reqendsend(:)
      integer,allocatable :: req2endrec(:),req2endsend(:)
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
c
      allocate (buffermpi1(3,nrhs,max(npoleloc,1)))
      buffermpi1 = 0d0
      allocate (buffermpimu1(3,nrhs,max(npoleloc,1)))
      buffermpimu1 = 0d0
      allocate (buffermpi2(3,nrhs,max(npolerecloc,1)))
      buffermpi2 = 0d0
      allocate (buffermpimu2(3,nrhs,max(npolerecloc,1)))
      buffermpimu2 = 0d0
c
      allocate (dipfield(3,nrhs,max(1,npoleloc)))
      dipfield = 0d0
      allocate (dipfieldbis(3,nrhs,max(1,npolerecloc)))
      dipfieldbis = 0d0
      allocate (reqrecdirrec(nproc))
      allocate (reqrecdirsend(nproc))
      allocate (reqrec(nproc))
      allocate (reqsend(nproc))
      allocate (req2rec(nproc))
      allocate (req2send(nproc))
      allocate (reqendsend(nproc))
      allocate (reqendrec(nproc))
      allocate (req2endsend(nproc))
      allocate (req2endrec(nproc))
c
c     allocate some memory and setup the preconditioner:
c
      allocate (res(3,nrhs,max(1,npoleloc)))
      allocate (h(3,nrhs,max(1,npolebloc)))
      allocate (pp(3,nrhs,max(1,npolebloc)))
      allocate (zr(3,nrhs,max(1,npoleloc)))
      allocate (diag(npoleloc))
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
        do i = 1, npoleloc
          diag(i) = one
        end do
      end if
c
c     initialize
c
      res = 0d0
      pp = 0d0
      zr = 0d0
      h = 0d0

c
c     now, compute the initial direction
c
      ggold = 0d0
c
      call tmatxbrecip(mu,murec,nrhs,dipfield,dipfieldbis)
      call matvec(nrhs,.true.,mu,h)
      call commfield(nrhs,h)
c
      call commrecdirsolv(nrhs,0,dipfieldbis,dipfield,buffermpi1,
     $    buffermpi2,reqrecdirrec,reqrecdirsend)
      call commdirdir(nrhs,0,pp,reqrec,reqsend)
c
      call commrecdirsolv(nrhs,1,dipfieldbis,dipfield,buffermpi1,
     $    buffermpi2,reqrecdirrec,reqrecdirsend)
c
      call commrecdirdip(nrhs,0,murec,pp,buffermpimu1,
     $ buffermpimu2,req2rec,req2send)
c
      call commrecdirsolv(nrhs,2,dipfieldbis,dipfield,buffermpi1,
     $  buffermpi2,reqrecdirrec,reqrecdirsend)

      term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
      do i = 1, npoleloc
        do k = 1, nrhs
          do j = 1, 3
            h(j,k,i) = h(j,k,i) + dipfield(j,k,i)- term*mu(j,k,i)
            res(j,k,i) = ef(j,k,i)-h(j,k,i)
            zr(j,k,i) = diag(i)*res(j,k,i)
          end do
        end do
      end do
      pp(:,:,1:npoleloc) = zr(:,:,1:npoleloc)
      ggold = 0d0
      do i = 1, npoleloc
        do k = 1, nrhs
          do j = 1, 3
            ggold(k) = ggold(k) + res(j,k,i)*zr(j,k,i)
          end do
        end do
      end do

      call MPI_IALLREDUCE(MPI_IN_PLACE,ggold(1),nrhs,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,req1,ierr)
c
c     MPI : begin sending
c
      call commdirdir(nrhs,1,pp,reqrec,reqsend)
      call commrecdirdip(nrhs,1,murec,pp,buffermpimu1,
     $   buffermpimu2,req2rec,req2send)
      call commdirdir(nrhs,2,pp,reqrec,reqsend)
      call commrecdirdip(nrhs,2,murec,pp,buffermpimu1,
     $  buffermpimu2,req2rec,req2send)
      call MPI_WAIT(req1,status,ierr)
c
c     now, start the main loop:
c
      do it = 1, politer
        do k = 1, nrhs
          gg(k) = zero
        end do
c
        call tmatxbrecip(pp,murec,nrhs,dipfield,dipfieldbis)
        call matvec(nrhs,.true.,pp,h)
        call commfield(nrhs,h)
c
c     Begin the reception of the reciprocal fields
c
        call commrecdirsolv(nrhs,0,dipfieldbis,dipfield,buffermpi1,
     $    buffermpi2,reqrecdirrec,reqrecdirsend)
c
c     MPI : begin reception
c
        call commdirdir(nrhs,0,pp,reqrec,reqsend)
        call commrecdirsolv(nrhs,1,dipfieldbis,dipfield,buffermpi1,
     $    buffermpi2,reqrecdirrec,reqrecdirsend)
c
c     Begin reception of mu for PME
c
        call commrecdirdip(nrhs,0,murec,pp,buffermpimu1,
     $  buffermpimu2,req2rec,req2send)
c
c       Wait for the reciprocal fields
c
        call commrecdirsolv(nrhs,2,dipfieldbis,dipfield,buffermpi1,
     $    buffermpi2,reqrecdirrec,reqrecdirsend)

        term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
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
              gg(k) = gg(k)+pp(j,k,i)*h(j,k,i)
            end do
          end do
        end do
        call MPI_IALLREDUCE(MPI_IN_PLACE,gg(1),nrhs,MPI_REAL8,MPI_SUM,
     $    COMM_TINKER,req2,ierr)
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
              res(j,k,i) = res(j,k,i) - alphacg(k)*h(j,k,i)
            end do
          end do
        end do

        do i = 1, npoleloc
          do k = 1, nrhs
            do j = 1, 3
              zr(j,k,i) = diag(i)*res(j,k,i)
            end do
          end do
        end do
        ggnew = 0d0
        ene = 0d0
        do i = 1, npoleloc
          do k = 1, nrhs
            do j = 1, 3
              ggnew(k) = ggnew(k) + res(j,k,i)*zr(j,k,i)
              ene(k) = ene(k)-pt5*(mu(j,k,i)*(res(j,k,i)+ef(j,k,i)))
            end do
          end do
        end do
        call MPI_IALLREDUCE(MPI_IN_PLACE,ggnew(1),nrhs,MPI_REAL8,
     $    MPI_SUM,COMM_TINKER,req3,ierr)
        call MPI_IALLREDUCE(MPI_IN_PLACE,ene(1),nrhs,MPI_REAL8,MPI_SUM,
     $    COMM_TINKER,req4,ierr)
        call MPI_WAIT(req3,status,ierr)
        call MPI_WAIT(req4,status,ierr)
        resnrm = zero
        do k = 1, nrhs
          gnorm(k) = sqrt(ggnew(k)/dble(3*npolar))
          resnrm   = max(resnrm,gnorm(k))
        end do
        if (polprt.ge.2.and.rank.eq.0) write(iout,1010)
     $    it, (ene(k)*coulomb, gnorm(k), k = 1, nrhs)

        do i = 1, npoleloc
          do k = 1, nrhs
            do j = 1, 3
              pp(j,k,i) = zr(j,k,i)+ggnew(k)/ggold(k)*pp(j,k,i)
            end do
          end do
        end do
c
        call commdirdir(nrhs,1,pp,reqrec,reqsend)
c
        call commrecdirdip(nrhs,1,murec,pp,buffermpimu1,
     $  buffermpimu2,req2rec,req2send)

        call commdirdir(nrhs,2,pp,reqrec,reqsend)
        call commrecdirdip(nrhs,2,murec,pp,buffermpimu1,
     $  buffermpimu2,req2rec,req2send)
c
        ggold = ggnew
        if (resnrm.lt.poleps) then
          if (polprt.ge.1.and.rank.eq.0) write(iout,1000) it,
     $      (ene(k)*coulomb, k = 1, nrhs), (gnorm(k), k = 1, nrhs)
          goto 10
        end if
      end do
 10   continue
c
c     MPI : begin reception
c
      call commdirdir(nrhs,0,mu,reqendrec,reqendsend)
c
c     Begin reception of mu for PME
c
      call commrecdirdip(nrhs,0,murec,mu,buffermpimu1,
     $  buffermpimu2,req2endrec,req2endsend)
c
c     MPI : begin sending
c
      call commdirdir(nrhs,1,mu,reqendrec,reqendsend)
c
c
      call commrecdirdip(nrhs,1,murec,mu,buffermpimu1,
     $  buffermpimu2,req2endrec,req2endsend)
      call commdirdir(nrhs,2,mu,reqendrec,reqendsend)
      call commrecdirdip(nrhs,2,murec,mu,buffermpimu1,
     $  buffermpimu2,req2endrec,req2endsend)
c
c     finalize and debug printing:
c
      if (polprt.ge.3) then
        write(iout,1020)
        do i = 1, npoleloc
          iipole = poleglob(i)
          iglob = ipole(iipole)
          write(iout,1030) iglob, (mu(j,1,i), j = 1, 3)
        end do
      end if
      if (polprt.ge.4) then
        write(iout,1021)
        do i = 1, npoleloc
          iipole = poleglob(i)
          iglob = ipole(iipole)
          write(iout,1030) iglob, (mu(j,2,i), j = 1, 3)
        end do
      end if
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
      deallocate (buffermpi1)
      deallocate (buffermpimu1)
      deallocate (buffermpimu2)
      deallocate (buffermpi2)
      deallocate (dipfield)
      deallocate (dipfieldbis)
      deallocate (res)
      deallocate (h)
      deallocate (pp)
      deallocate (zr)
      deallocate (diag)
      return
      end
c
      subroutine inducejac_pme2(matvec,nrhs,dodiis,ef,mu,murec)
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
      integer i, j, k, it, ind, ndismx, nmat, lenb
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
      real*8, allocatable :: buffermpimu1(:,:,:),buffermpimu2(:,:,:)
      integer iglob, iipole, ierr
      integer reqnorm,reqdiis(2*ndismx+1)
      integer status(MPI_STATUS_SIZE)
      integer,allocatable :: reqrecdirrec(:),reqrecdirsend(:)
      integer,allocatable :: reqrec(:),reqsend(:)
      integer,allocatable :: req2rec(:),req2send(:)
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
     
      allocate (buffermpi1(3,nrhs,max(npoleloc,1)))
      buffermpi1 = 0d0
      allocate (buffermpimu1(3,nrhs,max(npoleloc,1)))
      buffermpimu1 = 0d0
      allocate (buffermpi2(3,nrhs,max(npolerecloc,1)))
      buffermpi2 = 0d0
      allocate (buffermpimu2(3,nrhs,max(npolerecloc,1)))
      buffermpimu2 = 0d0
c
      allocate (munew(3,nrhs,max(1,npolebloc)))
      munew = 0d0
      allocate (dipfield(3,nrhs,max(1,npoleloc)))
      dipfield = 0d0
      allocate (dipfieldbis(3,nrhs,max(1,npolerecloc)))
      dipfieldbis = 0d0
      allocate (reqrecdirsend(nproc))
      allocate (reqrecdirrec(nproc))
      allocate (reqrec(nproc))
      allocate (reqsend(nproc))
      allocate (req2rec(nproc))
      allocate (req2send(nproc))
c
      allocate (h(3,nrhs,max(1,npolebloc)))
      h = 0d0
      if (dodiis) then
        nmat = 1
        lenb = ndismx + 1
        allocate (xdiis(3*nrhs*max(npoleloc,1),ndismx))
        allocate (ediis(3*nrhs*max(npoleloc,1),ndismx))
        allocate (bmat(lenb,lenb))
        allocate (bloc(lenb,lenb))
        allocate (cex(lenb))
        bmat = 0d0
      end if
c
c     main loop:
c
      do it = 1, politer
        h = 0d0
        rnorm = 0d0
c
        call tmatxbrecip(mu,murec,nrhs,dipfield,dipfieldbis)
        call matvec(nrhs,.false.,mu,h)
        call commfield(nrhs,h)
c
        call commrecdirsolv(nrhs,0,dipfieldbis,dipfield,buffermpi1,
     $    buffermpi2,reqrecdirrec,reqrecdirsend)
c
        call commdirdir(nrhs,0,mu,reqrec,reqsend)
        call commrecdirsolv(nrhs,1,dipfieldbis,dipfield,buffermpi1,
     $    buffermpi2,reqrecdirrec,reqrecdirsend)
c
        call commrecdirdip(nrhs,0,murec,munew,buffermpimu1,
     $  buffermpimu2,req2rec,req2send)
c
c     jacobi step:
c
        call commrecdirsolv(nrhs,2,dipfieldbis,dipfield,buffermpi1,
     $    buffermpi2,reqrecdirrec,reqrecdirsend)
c        
        term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
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
        rnorm = 0d0
        do i = 1, npoleloc
          do k = 1, nrhs
            do j = 1, 3
              rnorm(k) = rnorm(k)+(munew(j,k,i)-mu(j,k,i))**2
            end do
          end do
        end do
c
        call MPI_IALLREDUCE(MPI_IN_PLACE,rnorm(1),nrhs,MPI_REAL8,
     $    MPI_SUM,COMM_TINKER,reqnorm,ierr)
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
c         Compute Pulay's Matrix and extrapolate
c
          call diis(ndismx,3*nrhs*npoleloc,xdiis,ediis,
     $      bmat,nmat,reqdiis,COMM_TINKER)
c
          do i = 1, 2*nmat-3
            call MPI_WAIT(reqdiis(i),status,ierr)
          end do
          bloc = bmat
          cex = 0d0
          cex(1) = one
          call dgesv(nmat,1,bloc,ndismx+1,ipiv,cex,nmat,info)
          munew = 0d0
          call extrap(3*nrhs*npoleloc,nmat-1,xdiis,cex,munew)
        end if
c
        call commdirdir(nrhs,1,munew,reqrec,reqsend)
        call commrecdirdip(nrhs,1,murec,munew,buffermpimu1,
     $  buffermpimu2,req2rec,req2send)
        call commdirdir(nrhs,2,mu,reqrec,reqsend)
        call commrecdirdip(nrhs,2,murec,munew,buffermpimu1,
     $  buffermpimu2,req2rec,req2send)

        do i = 1, npoleloc
          do k = 1, nrhs
            do j = 1, 3
              mu(j,k,i) = munew(j,k,i)
            end do
          end do
        end do
c
c     compute the norm of the increment.
c
        call MPI_WAIT(reqnorm,status,ierr)
        rr = zero
        do k = 1, nrhs
         rnorm(k) = sqrt(rnorm(k)/dble(3*npolar))
          rr = max(rnorm(k),rr)
        end do
        if (rr.lt.poleps) then
          if (polprt.ge.1.and.rank.eq.0)
     $       write(6,1000) it, (rnorm(k), k = 1, nrhs)
          goto 10
        end if
        if (polprt.ge.1.and.rank.eq.0)
     $       write(6,1010) it, (rnorm(k), k = 1, nrhs)
      end do
 10   continue
      if (polprt.ge.3) then
        write(iout,1020)
        do i = 1, npoleloc
          iipole = poleglob(i)
          iglob = ipole(iipole)
          write(iout,1030) iglob, (mu(j,1,i), j = 1, 3)
        end do
      end if
      if (polprt.ge.4) then
        write(iout,1021)
        do i = 1, npoleloc
          iipole = poleglob(i)
          iglob = poleglob(iipole)
          write(iout,1030) iglob, (mu(j,2,i), j = 1, 3)
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
      deallocate (buffermpi1)
      deallocate (buffermpimu2)
      deallocate (buffermpi2)
      deallocate (munew)
      deallocate (dipfield)
      deallocate (dipfieldbis)
      deallocate (h)
      if (dodiis) then
        deallocate (xdiis)
        deallocate (ediis)
        deallocate (bmat)
        deallocate (bloc)
        deallocate (cex)
      end if
      return
      end
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine ulspred  --  induced dipole prediction coeffs  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ulspred" uses standard extrapolation or a least squares fit
c     to set coefficients of an induced dipole predictor polynomial
c
c     literature references:
c
c     J. Kolafa, "Time-Reversible Always Stable Predictor-Corrector
c     Method for Molecular Dynamics of Polarizable Molecules", Journal
c     of Computational Chemistry, 25, 335-342 (2004)
c
c     W. Wang and R. D. Skeel, "Fast Evaluation of Polarizable Forces",
c     Journal of Chemical Physics, 123, 164107 (2005)
c
c
      subroutine ulspred
      use mpole
      use uprior
      implicit none
      integer i,j,k,m
      real*8 coeff,udk,upk
      real*8 amax,apmax
      real*8 b(maxualt)
      real*8 bp(maxualt)
      real*8 a(maxualt*(maxualt+1)/2)
      real*8 ap(maxualt*(maxualt+1)/2)
      real*8 c(maxualt,maxualt)
      real*8 cp(maxualt,maxualt)
c
c
c     set the Gear predictor binomial coefficients
c
      if (polpred .eq. 'GEAR') then
         do i = 1, nualt
            coeff = gear(i)
            bpred(i) = coeff
            bpredp(i) = coeff
            bpreds(i) = coeff
            bpredps(i) = coeff
         end do
c
c     set always stable predictor-corrector (ASPC) coefficients
c
      else if (polpred .eq. 'ASPC') then
         do i = 1, nualt
            coeff = aspc(i)
            bpred(i) = coeff
            bpredp(i) = coeff
            bpreds(i) = coeff
            bpredps(i) = coeff
         end do
c
c     derive normal equations corresponding to least squares fit
c
      else
         do k = 1, nualt
            b(k) = 0.0d0
            bp(k) = 0.0d0
            do m = k, nualt
               c(k,m) = 0.0d0
               cp(k,m) = 0.0d0
            end do
         end do
         do i = 1, npole
            do j = 1, 3
               do k = 1, nualt
                  udk = udalt(k,j,i)
                  upk = upalt(k,j,i)
                  do m = k, nualt
                     c(k,m) = c(k,m) + udk*udalt(m,j,i)
                     cp(k,m) = cp(k,m) + upk*upalt(m,j,i)
                  end do
               end do
            end do
         end do
         i = 0
         do k = 2, nualt
            b(k-1) = c(1,k)
            bp(k-1) = cp(1,k)
            do m = k, nualt
               i = i + 1
               a(i) = c(k,m)
               ap(i) = cp(k,m)
            end do
         end do
c
c     check for nonzero coefficients and solve normal equations
c
         k = nualt - 1
         amax = 0.0d0
         apmax = 0.0d0
         do i = 1, k*(k+1)/2
            amax = max(amax,a(i))
            apmax = max(apmax,ap(i))
         end do
         if (amax .ne. 0.0d0)  call cholesky (k,a,b)
         if (apmax .ne. 0.0d0)  call cholesky (k,ap,bp)
c
c     transfer the final solution to the coefficient vector
c
         do k = 1, nualt-1
            bpred(k) = b(k)
            bpredp(k) = bp(k)
            bpreds(k) = b(k)
            bpredps(k) = bp(k)
         end do
         bpred(nualt) = 0.0d0
         bpredp(nualt) = 0.0d0
         bpreds(nualt) = 0.0d0
         bpredps(nualt) = 0.0d0
      end if
      return
      end
c
!===================================================
!     sub diagvec
!===================================================
! Performs product of vector a with polarisabilities
!
      subroutine diagvec(nrhs, A, B)
      use atmlst
      use mpole
      use polar
      implicit none

      integer, intent(in) :: nrhs
      real*8, dimension(3,nrhs,npolebloc) :: A
      real*8, dimension(3,nrhs,npolebloc) :: B
      integer :: i,iipole, irhs, j

      do i = 1, npolebloc
         iipole = poleglob(i)
         do irhs = 1, nrhs
            do j = 1,3
               B(j,irhs,i) = A(j,irhs,i)*polarity(iipole)
            end do
         end do
      end do

      return
      end
