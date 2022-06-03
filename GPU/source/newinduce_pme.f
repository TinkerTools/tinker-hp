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
#include "tinker_precision.h"
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
      use tinheader ,only:ti_p,re_p
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
      integer ierr, iipole, proc,tag
      integer, allocatable :: reqrecdirsend(:),reqrecdirrec(:)
      integer, allocatable :: reqsend(:),reqrec(:)
      integer, allocatable :: req2send(:),req2rec(:)
      integer status(MPI_STATUS_SIZE)
c
      parameter (nrhs=2)
      real(t_p)  wtime0, wtime1, wtime2, udsum, upsum
      real(t_p)  term, xx(1)
      real(t_p), allocatable::ef(:,:,:),mu(:,:,:),murec(:,:,:)
      real(t_p), allocatable::cphi(:,:)
c
      real(t_p), allocatable::buffermpi(:,:),buffermpimu(:,:,:)
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
        mu = 0.0_ti_p
        allocate (murec(3,nrhs,max(npolerecloc,1)))
        murec = 0.0_ti_p
        allocate (ef(3,nrhs,max(1,npolebloc)))
        ef = 0.0_ti_p
        allocate (cphi(10,max(1,npoleloc)))
        cphi = 0.0_ti_p
        if (allocated(cphirec)) deallocate (cphirec)
        allocate (cphirec(10,max(npolerecloc,1)))
        cphirec = 0.0_ti_p
        allocate (buffermpi(10,max(npoleloc,1)))
        buffermpi = 0.0_ti_p
        allocate (buffermpimu(3,nrhs,max(npoleloc,1)))
        buffermpimu = 0.0_ti_p
      else
        allocate (mu(3,nrhs,max(1,npolebloc)))
        mu = 0.0_ti_p
        allocate (murec(3,nrhs,max(npolerecloc,1)))
        murec = 0.0_ti_p
        if (allocated(fphirec)) deallocate (fphirec)
        allocate (fphirec(20,max(npolerecloc,1)))
        fphirec = 0.0_ti_p
        allocate (cphi(10,max(1,npoleloc)))
        cphi = 0.0_ti_p
        if (allocated(cphirec)) deallocate (cphirec)
        allocate (cphirec(10,max(npolerecloc,1)))
        cphirec = 0.0_ti_p
        allocate (buffermpi(10,max(npolerecloc,1)))
        buffermpi = 0.0_ti_p
        allocate (buffermpimu(3,nrhs,max(npolerecloc,1)))
        buffermpimu = 0.0_ti_p
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
        term = (4.0_ti_p/3.0_ti_p) * aewald**3 / sqrtpi
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
             udsum = 0.0_ti_p
             upsum = 0.0_ti_p
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
c     move the computed dipoles in the common block.
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
      if (use_pred) then
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
      use tinheader ,only:ti_p,re_p
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
      real(t_p)  ef(3,nrhs,*), mu(3,nrhs,*), murec(3,nrhs,*)
      logical precnd
      real(t_p), allocatable :: res(:,:,:), h(:,:,:),
     &                pp(:,:,:), zr(:,:,:), diag(:)
      integer i, iglob, it, j, k, iipole
      real(t_p)  ggold(2), ggnew(2), gnorm(2), gg(2), 
     &         alphacg(2), ene(2)
      real(t_p)  zero, pt5, one, resnrm
      save    zero, pt5, one
      data    zero/0.0_ti_p/, pt5/0.50_ti_p/, one/1.0_ti_p/
      external matvec
      real(t_p), allocatable :: dipfield(:,:,:),
     &                          dipfieldbis(:,:,:)
      real(t_p) term
      real(t_p) time0,time1
c
c     MPI
c
      real(t_p), allocatable :: buffermpi1(:,:,:),
     &                          buffermpi2(:,:,:)
      real(t_p), allocatable :: buffermpimu(:,:,:)
      integer ierr, tag, proc
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
c
c
c     allocate some memory and setup the preconditioner:
c
      if (rank.le.ndir-1) then
        allocate (buffermpimu(3,nrhs,max(npoleloc,1)))
        buffermpimu = 0.0_ti_p
        allocate (dipfield(3,nrhs,max(1,npoleloc)))
        dipfield = 0_ti_p
        allocate (dipfieldbis(3,nrhs,max(1,npolerecloc)))
        dipfieldbis = 0_ti_p
        allocate (res(3,nrhs,max(1,npoleloc)))
        allocate (zr(3,nrhs,max(1,npoleloc)))
        allocate (diag(npoleloc))
        allocate (h(3,nrhs,max(1,npolebloc)))
      else 
        allocate (buffermpimu(3,nrhs,max(1,npolerecloc)))
        buffermpimu = 0.0_ti_p
        allocate (dipfieldbis(3,nrhs,max(1,npolerecloc)))
        dipfieldbis = 0_ti_p
        allocate (dipfield(3,nrhs,max(1,npoleloc)))
        dipfield = 0_ti_p
        allocate (diag(npoleloc))
      end if
      allocate (pp(3,nrhs,max(1,npolebloc)))
      allocate (buffermpi1(3,nrhs,max(npoleloc,1)))
      buffermpi1 = 0.0_ti_p
      allocate (buffermpi2(3,nrhs,max(1,npolerecloc)))
      buffermpi2 = 0.0_ti_p
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
          term   = polarity(iipole)
          diag(i) = merge(tinypol,term,term.eq.0.0)
        end do
        if (polprt.ge.2.and.rank.eq.0) write (iout,1040)
      else
        diag = 1_ti_p
      end if
c
c     initialize
c
      if (rank.le.ndir-1) then
        pp = 0.0_ti_p
        h  = 0.0_ti_p
        res  = 0.0_ti_p
      end if
c
c     now, compute the initial direction
c
      ggold = 0.0_ti_p
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
        time0 = mpi_wtime()
        call tmatxbrecip(mu,murec,nrhs,dipfield,dipfieldbis)
        time1 = mpi_wtime()
        timerecdip = timerecdip + time1-time0
        call commrecdirsolv(nrhs,1,dipfieldbis,dipfield,buffermpi1,
     $      buffermpi2,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,dipfieldbis,dipfield,buffermpi1,
     $      buffermpi2,reqrecdirrec,reqrecdirsend)
      end if
      if (rank.le.ndir-1) then
c
        time0 = mpi_wtime()
        call matvec(nrhs,.true.,mu,h)
        time1 = mpi_wtime()
        timerealdip = timerealdip + time1-time0
        call commfield(nrhs,h)
c
        call commrecdirsolv(nrhs,2,dipfieldbis,dipfield,buffermpi1,
     $    buffermpi2,reqrecdirrec,reqrecdirsend)
        term = (4.0_ti_p/3.0_ti_p) * aewald**3 / sqrtpi
        h(:,:,1:npoleloc) = h(:,:,1:npoleloc)+dipfield(:,:,1:npoleloc)
     $  - term*mu(:,:,1:npoleloc)
        res(:,:,1:npoleloc) = ef(:,:,1:npoleloc)-h(:,:,1:npoleloc)
        do k = 1, nrhs
          do j = 1, 3
            zr(j,k,1:npoleloc) = diag(1:npoleloc)*res(j,k,1:npoleloc)
          end do
        end do
        pp(:,:,1:npoleloc) = zr(:,:,1:npoleloc)
        do k = 1, nrhs
          ggold(k) = sum(res(:,k,:)*zr(:,k,:))
        end do
        call MPI_IALLREDUCE(MPI_IN_PLACE,ggold(1),nrhs,MPI_TPREC,
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
          h(:,:,1:npoleloc) = h(:,:,1:npoleloc)+dipfield(:,:,1:npoleloc)
     $     - term*pp(:,:,1:npoleloc)
          do k = 1, nrhs
            gg(k) = sum(pp(:,k,1:npoleloc)*h(:,k,1:npoleloc))
          end do
          call MPI_IALLREDUCE(MPI_IN_PLACE,gg(1),nrhs,MPI_TPREC,MPI_SUM,
     $      comm_dir,req2,ierr)
          call MPI_WAIT(req2,status,ierr)
          do k = 1, nrhs
            if (gg(k).eq.zero) return
            alphacg(k)  = ggold(k)/gg(k)
            ggnew(k)  = zero
            ene(k)    = zero
          end do
          do k = 1, nrhs
            mu(:,k,1:npoleloc) = mu(:,k,1:npoleloc) + alphacg(k)*
     $       pp(:,k,1:npoleloc)
            res(:,k,1:npoleloc) = res(:,k,1:npoleloc) - alphacg(k)*
     $       h(:,k,1:npoleloc)
          end do
          do k = 1, nrhs
            do j = 1, 3
              zr(j,k,1:npoleloc) = diag(1:npoleloc)*res(j,k,1:npoleloc)
            end do
          end do
          do k = 1, nrhs
            ggnew(k) = sum(res(:,k,:)*zr(:,k,:))
            ene(k) = -pt5*sum(mu(:,k,1:npoleloc)*(res(:,k,1:npoleloc)+
     $        ef(:,k,1:npoleloc)))
          end do
        end if
        call MPI_IALLREDUCE(MPI_IN_PLACE,ggnew(1),nrhs,MPI_TPREC,
     $      MPI_SUM,MPI_COMM_WORLD,req3,ierr)
        call MPI_WAIT(req3,status,ierr)
        if (rank.le.ndir-1) then
          call MPI_IALLREDUCE(MPI_IN_PLACE,ene(1),nrhs,MPI_TPREC,
     $      MPI_SUM,comm_dir,req4,ierr)
          call MPI_WAIT(req4,status,ierr)
          do k = 1, nrhs
            pp(:,k,1:npoleloc) = zr(:,k,1:npoleloc)+ggnew(k)/ggold(k)*
     $         pp(:,k,1:npoleloc)
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
          gnorm(k) = sqrt(ggnew(k)/real(3*npolar,t_p))
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
      use tinheader ,only:ti_p,re_p
      use mpi
      implicit none
c
c     solves the polarizable dipoles linear equations by jacobi
c     iterations coupled to the DIIS extrapolation.
c
      integer nrhs, info
      real(t_p)  ef(3,nrhs,*), mu(3,nrhs,*), murec(3,nrhs,*)
      real(t_p), allocatable :: munew(:,:,:), h(:,:,:)
      real(t_p), allocatable :: xdiis(:,:), ediis(:,:),
     $   bmat(:,:), bloc(:,:), cex(:)
      integer i, j, k, iipole, it, ind, ndismx, nmat, lenb
      logical dodiis
      parameter (ndismx=25)
      integer ipiv(ndismx+1)
      real(t_p)  zero, one, rnorm(2), rr, xx(1)
      real(t_p) term
      real(t_p) time0,time1,time2,time3
      save    zero, one, xx
      data    zero/0.0_ti_p/, one/1.0_ti_p/, xx/0.0_ti_p/
      external matvec
      real(t_p), allocatable :: dipfield(:,:,:),
     &                       dipfieldbis(:,:,:)
c
c     MPI
c
      real(t_p), allocatable :: buffermpi1(:,:,:),
     &                          buffermpi2(:,:,:)
      real(t_p), allocatable :: buffermpimu(:,:,:)
      integer ierr, tag, proc
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
      if (rank.le.ndir-1) then
        allocate (buffermpimu(3,nrhs,max(npoleloc,1)))
        buffermpimu = 0.0_ti_p
        allocate (munew(3,nrhs,max(1,npolebloc)))
        munew = 0.0_ti_p
        allocate (h(3,nrhs,max(1,npolebloc)))
        h = 0.0_ti_p
        if (dodiis) then
          nmat = 1
          lenb = ndismx + 1
          allocate (xdiis(3*nrhs*max(1,npoleloc),ndismx))
          allocate (ediis(3*nrhs*max(1,npoleloc),ndismx))
          allocate (bmat(lenb,lenb))
          allocate (bloc(lenb,lenb))
          allocate (cex(lenb))
          bmat = 0.0_ti_p
        end if
      else
        allocate (munew(3,nrhs,max(1,npolebloc)))
        munew = 0.0_ti_p
        allocate (buffermpimu(3,nrhs,max(1,npolerecloc)))
        buffermpimu = 0.0_ti_p
      end if
      allocate (dipfieldbis(3,nrhs,max(1,npolerecloc)))
      dipfieldbis = 0_ti_p
      allocate (dipfield(3,nrhs,max(1,npoleloc)))
      dipfield = 0_ti_p
      allocate (buffermpi1(3,nrhs,max(npoleloc,1)))
      buffermpi1 = 0_ti_p
      allocate (buffermpi2(3,nrhs,max(npolerecloc,1)))
      buffermpi2 = 0_ti_p
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
        rnorm = 0.0_ti_p
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
          time0 = mpi_wtime()
          call tmatxbrecip(mu,murec,nrhs,dipfield,dipfieldbis)
          time1 = mpi_wtime()
          if (it.eq.1) timerecdip = timerecdip + time1-time0
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
         term = (4.0_ti_p/3.0_ti_p) * aewald**3 / sqrtpi
         time0 = mpi_wtime()
         call matvec(nrhs,.false.,mu,h)
         time1 = mpi_wtime()
         if (it.eq.1) timerealdip = timerealdip + time1-time0
         call commfield(nrhs,h)
c
          call commrecdirsolv(nrhs,2,dipfieldbis,dipfield,buffermpi1,
     $      buffermpi2,reqrecdirrec,reqrecdirsend)
c
          h(1:3,1:nrhs,1:npoleloc) = h(:,:,1:npoleloc)+
     $      dipfield(:,:,1:npoleloc)
     $     - term*mu(:,:,1:npoleloc)
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
          do k = 1, nrhs
            rnorm(k)=sum((munew(1:3,k,1:npoleloc)-
     $        mu(1:3,k,1:npoleloc))**2)
          end do
c
        end if
        call MPI_IALLREDUCE(MPI_IN_PLACE,rnorm(1),nrhs,MPI_TPREC,
     $    MPI_SUM,MPI_COMM_WORLD,reqnorm,ierr)
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
            cex = 0.0_ti_p
            cex(1) = one
#ifdef _OPENACC
            call fatal_device("inducejac_pme")
#else
            call dgesv(nmat,1,bloc,ndismx+1,ipiv,cex,nmat,info)
#endif
            munew = 0.0_ti_p
            call extrap(3*nrhs*npoleloc,nmat-1,xdiis,cex,munew)
          end if
        time2 = mpi_wtime()
c
          call commdirdir(nrhs,1,munew,reqrec,reqsend)
          call commrecdirdip(nrhs,1,murec,munew,buffermpimu,
     $    buffermpimu,req2rec,req2send)
c
          time0 = mpi_wtime()
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
         rnorm(k) = sqrt(rnorm(k)/real(3*npolar,t_p))
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
      use tinheader ,only:ti_p,re_p
      use mpi
      implicit none
      integer ierr
      integer status(MPI_STATUS_SIZE),tag
      integer i,j,k,iglob,iipole,iloc
      integer k1,k2,k3
      integer m1,m2,m3
      integer nf1,nf2,nf3
      real(t_p) r1,r2,r3
      real(t_p) h1,h2,h3
      real(t_p) volterm,denom
      real(t_p) hsq,expterm
      real(t_p) term,pterm
      real(t_p) cmp(10),fmp(10),cphi(10,*)
      real(t_p), allocatable :: qgridmpi(:,:,:,:,:)
      integer, allocatable :: reqrec(:),reqsend(:)
      integer, allocatable :: reqbcastrec(:),reqbcastsend(:)
      integer nprocloc,commloc,rankloc,proc
      real(t_p) time0,time1
c
      if (use_pmecore) then
        nprocloc = nrec
        commloc  = comm_rec
        rankloc  = rank_bis
      else
        nprocloc = nproc
        commloc  = MPI_COMM_WORLD
        rankloc  = rank
      end if

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
     $   n3mpimax,MPI_TPREC,prec_recep(i),tag,
     $   commloc,reqrec(i),ierr)
      end do
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
c
c     zero out the PME charge grid
      qgridin_2d = 0_ti_p
c
c     fill the pme grid, loop over the multipoles sites
c
      time0 = mpi_wtime()
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
         cmp(8) = 2.0_ti_p * rpole(6,iipole)
         cmp(9) = 2.0_ti_p * rpole(7,iipole)
         cmp(10) = 2.0_ti_p * rpole(10,iipole)
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
      time1 = mpi_wtime()
      timegrid1 = timegrid1 + time1-time0
c
c     MPI : begin sending
c
      time0 = mpi_wtime()
      do i = 1, nrec_send
        tag = nprocloc*prec_send(i) + rankloc + 1
        call MPI_ISEND(qgridin_2d(1,1,1,1,i+1),
     $   2*n1mpimax*n2mpimax*n3mpimax,MPI_TPREC,
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
      time1 = mpi_wtime()
      timerecreccomm = timerecreccomm + time1 - time0
c
c     Perform 3-D FFT forward transform
c
      time0 = mpi_wtime()
      call fft2d_frontmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
     $ n3mpimax)
      time1 = mpi_wtime()
      timeffts = timeffts + time1-time0
c
c     make the scalar summation over reciprocal lattice
c
      time0 = mpi_wtime()
      if ((istart2(rankloc+1).eq.1).and.(jstart2(rankloc+1).eq.1).and.
     $   (kstart2(rankloc+1).eq.1)) then
           qfac_2d(1,1,1) = 0.0_ti_p
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
            r1 = real(m1,t_p)
            r2 = real(m2,t_p)
            r3 = real(m3,t_p)
            h1 = recip(1,1)*r1 + recip(1,2)*r2 + recip(1,3)*r3
            h2 = recip(2,1)*r1 + recip(2,2)*r2 + recip(2,3)*r3
            h3 = recip(3,1)*r1 + recip(3,2)*r2 + recip(3,3)*r3
            hsq = h1*h1 + h2*h2 + h3*h3
            term = -pterm * hsq
            expterm = 0.0_ti_p
            if ((term .gt. -50.0_ti_p)) then
               denom = volterm*hsq*bsmod1(k1)*bsmod2(k2)*bsmod3(k3)
               expterm = exp(term) / denom
               if (.not. use_bounds) then
                  expterm = expterm * (1.0_ti_p-cos(pi*xbox*sqrt(hsq)))
               else if (octahedron) then
                  if (mod(m1+m2+m3,2) .ne. 0)  expterm = 0.0_ti_p
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
           expterm = 0.5_ti_p * pi / xbox
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
      time1 = mpi_wtime()
      timescalar = timescalar + time1-time0
c
c     perform 3-D FFT backward transform
c
      time0 = mpi_wtime()
      call fft2d_backmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
     $ n3mpimax)
      time1 = mpi_wtime()
      timeffts = timeffts + time1-time0
c
c     MPI : Begin reception
c
      time0 = mpi_wtime()
      do i = 1, nrec_send
        tag = nprocloc*rankloc + prec_send(i) + 1
        call MPI_IRECV(qgridin_2d(1,1,1,1,i+1),
     $   2*n1mpimax*n2mpimax*n3mpimax,MPI_TPREC,
     $   prec_send(i),tag,commloc,reqbcastrec(i),ierr)
      end do
c
c     MPI : begin sending
c
      do i = 1, nrec_recep
        tag = nprocloc*prec_recep(i) + rankloc + 1
        call MPI_ISEND(qgridin_2d(1,1,1,1,1),
     $   2*n1mpimax*n2mpimax*n3mpimax,MPI_TPREC,
     $   prec_recep(i),tag,commloc,reqbcastsend(i),ierr)
      end do
c
      do i = 1, nrec_send
        call MPI_WAIT(reqbcastrec(i),status,ierr)
      end do
      do i = 1, nrec_recep
        call MPI_WAIT(reqbcastsend(i),status,ierr)
      end do
      time1 = mpi_wtime()
      timerecreccomm = timerecreccomm + time1-time0
c
c     get field
c
      time0 = mpi_wtime()
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
      time1 = mpi_wtime()
      timegrid2 = timegrid2 + time1-time0
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
      use tinheader ,only:ti_p,re_p
      implicit none
      integer ierr,iglob
      integer status(MPI_STATUS_SIZE),tag
      integer nrhs,iipole
      integer i,j,k,iloc
      real(t_p) fuind(3),fuinp(3)
      real(t_p) term
      real(t_p) a(3,3)
      real(t_p) fdip_phi1(10), fdip_phi2(10), fdip_sum_phi(20)
      real(t_p) dipfield(3,nrhs,*),dipfieldbis(3,nrhs,*)
      real(t_p) mu(3,nrhs,*),murec(3,nrhs,*)
      real(t_p), allocatable :: qgridmpi(:,:,:,:,:)
      integer, allocatable :: reqbcastrec(:),reqbcastsend(:)
      integer, allocatable :: reqrec(:),reqsend(:)
      integer nprocloc,commloc,rankloc,proc
c
      if (use_pmecore) then
        nprocloc = nrec
        commloc  = comm_rec
        rankloc  = rank_bis
      else
        nprocloc = nproc
        commloc  = MPI_COMM_WORLD
        rankloc  = rank
      end if

      allocate (reqbcastrec(nprocloc))
      allocate (reqbcastsend(nprocloc))
      allocate (reqrec(nprocloc))
      allocate (reqsend(nprocloc))
      allocate (qgridmpi(2,n1mpimax,n2mpimax,n3mpimax,nrec_recep))
c
c     MPI : Begin reception
c
      do i = 1, nrec_recep
        tag = nprocloc*rankloc + prec_recep(i) + 1
        call MPI_IRECV(qgridmpi(1,1,1,1,i),2*n1mpimax*n2mpimax*n3mpimax,
     $   MPI_TPREC,prec_recep(i),tag,
     $   commloc,reqrec(i),ierr)
      end do
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
c
c     zero out the PME charge grid
      qgrid2in_2d = 0_ti_p
c
c     fill the pme grid, loop over the induced dipoles sites
c
      do j = 1, 3
        a(1,j) = real(nfft1,t_p) * recip(j,1)
        a(2,j) = real(nfft2,t_p) * recip(j,2)
        a(3,j) = real(nfft3,t_p) * recip(j,3)
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
     $   2*n1mpimax*n2mpimax*n3mpimax,MPI_TPREC,
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
     $   2*n1mpimax*n2mpimax*n3mpimax,MPI_TPREC,
     $   prec_send(i),tag,commloc,reqbcastrec(i),ierr)
      end do
c
c     MPI : begin sending
c
      do i = 1, nrec_recep
        tag = nprocloc*prec_recep(i) + rankloc + 1
        call MPI_ISEND(qgrid2in_2d,
     $   2*n1mpimax*n2mpimax*n3mpimax,MPI_TPREC,
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
