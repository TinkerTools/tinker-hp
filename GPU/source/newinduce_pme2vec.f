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
      subroutine newinduce_pme2vec
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
      use timestat
      use tinheader ,only:ti_p,re_p
      implicit none
c
c     without separate cores for reciprocal part
c
      integer i, j, k, nrhs, proc
c
c     MPI
c
      integer iglob, ierr, iipole
      integer, allocatable :: reqrecdirsend(:),reqrecdirrec(:)
      integer, allocatable :: reqsend(:),reqrec(:)
      integer, allocatable :: req2send(:),req2rec(:)
      integer tag,status(MPI_STATUS_SIZE)
      real(t_p) time0,time1
c
      parameter (nrhs=2)
      real(t_p)  wtime0, wtime1, wtime2, udsum, upsum
      real(t_p)  term, xx(1)
      real(t_p), allocatable :: ef(:,:,:), mu(:,:,:), murec(:,:,:)
      real(t_p), allocatable :: cphi(:,:)
c
      real(t_p), allocatable :: buffermpi1(:,:),buffermpi2(:,:)
      real(t_p), allocatable :: buffermpimu1(:,:,:),buffermpimu2(:,:,:)
c
      external tmatxb_pmevec
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
      mu = 0_ti_p
      allocate (murec(3,nrhs,max(1,npolerecloc)))
      murec = 0_ti_p
c
      allocate (buffermpi1(10,max(npoleloc,1)))
      buffermpi1 = 0_ti_p
      allocate (buffermpi2(10,max(npolerecloc,1)))
      buffermpi2 = 0_ti_p
      allocate (buffermpimu1(3,nrhs,max(npoleloc,1)))
      buffermpimu1 = 0_ti_p
      allocate (buffermpimu2(3,nrhs,max(npolerecloc,1)))
      buffermpimu2 = 0_ti_p
c
      allocate (ef(3,nrhs,max(1,npolebloc)))
cnew  allocate (ef(max(1,npolebloc),3,nrhs))
      ef = 0_ti_p
      allocate (cphi(10,max(npoleloc,1)))
      cphi = 0_ti_p
      if (allocated(cphirec)) deallocate (cphirec)
      allocate (cphirec(10,max(npolerecloc,1)))
      cphirec = 0_ti_p
      if (allocated(fphirec)) deallocate(fphirec)
      allocate (fphirec(20,max(npolerecloc,1)))
      fphirec = 0_ti_p
c
      allocate (reqrecdirsend(nproc))
      allocate (reqrecdirrec(nproc))
      allocate (reqsend(nproc))
      allocate (reqrec(nproc))
      allocate (req2send(nproc))
      allocate (req2rec(nproc))
c
c    compute the electric fields:
c
c
c    compute the reciprocal space contribution (fields)
c
      time0 = mpi_wtime()
      call efld0_recip(cphi)
      time1 = mpi_wtime()
      timerecdip = timerecdip + time1-time0
c
      call commrecdirfields(0,cphirec,cphi,buffermpi1,buffermpi2,
     $ reqrecdirrec,reqrecdirsend)
      call commrecdirfields(1,cphirec,cphi,buffermpi1,buffermpi2,
     $ reqrecdirrec,reqrecdirsend)
      call commrecdirfields(2,cphirec,cphi,buffermpi1,buffermpi2,
     $ reqrecdirrec,reqrecdirsend)
c
      time0 = mpi_wtime()
      call efld0_directvec(nrhs,ef)
      time1 = mpi_wtime()
      timerealdip = timerealdip + time1-time0
c
      call commfieldvec(nrhs,ef)
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
cnew       ef(i,j,1)  = ef(i,j,1) - cphi(j+1,i) 
cnew       ef(i,j,2)  = ef(i,j,2) - cphi(j+1,i) 
        end do
      end do
c
      term = (4.0_ti_p/3.0_ti_p) * aewald**3 / sqrtpi
      do i = 1, npoleloc
        iipole = poleglob(i)
        do j = 1, 3
           ef(j,1,i)  = ef(j,1,i) + term*rpole(j+1,iipole)
           ef(j,2,i)  = ef(j,2,i) + term*rpole(j+1,iipole)
cnew      ef(i,j,1)  = ef(i,j,1) + term*rpole(j+1,iipole)
cnew      ef(i,j,2)  = ef(i,j,2) + term*rpole(j+1,iipole)
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
cnew          mu(j,k,i) = polarity(iipole)*ef(i,j,k)
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
        call inducepcg_pme2vec(tmatxb_pmevec,nrhs,.true.,ef,mu,murec)
      else if (polalg.eq.2) then
        call inducejac_pme2vec(tmatxb_pmevec,nrhs,.true.,ef,mu,murec)
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
      if (use_pred) then
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
      subroutine inducepcg_pme2vec(matvec,nrhs,precnd,ef,mu,murec)
      use atmlst
      use domdec
      use ewald
      use iounit
      use math
      use mpole
      use polar
      use polpot
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
cnew  real(t_p)  ef(npolebloc,3,nrhs), mu(3,nrhs,*), murec(3,nrhs,*)
      logical precnd
      real(t_p), allocatable :: res(:,:,:), h(:,:,:), pp(:,:,:),
     $  zr(:,:,:), diag(:)
      integer i, it, j, k
      real(t_p)  ggold(2), ggnew(2), gnorm(2), gg(2), alphacg(2), ene(2)
      real(t_p)  zero, pt5, one, resnrm, term
      save    zero, pt5, one
      data    zero/0.0_ti_p/, pt5/0.50_ti_p/, one/1.0_ti_p/
      external matvec
      real(t_p), allocatable :: dipfield(:,:,:),dipfieldbis(:,:,:)
      real(t_p) time0,time1,time2
c
c     MPI
c
      real(t_p), allocatable :: buffermpi1(:,:,:),buffermpi2(:,:,:)
      real(t_p), allocatable :: buffermpimu1(:,:,:),buffermpimu2(:,:,:)
      integer iglob, iipole, ierr, tag, proc
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
      buffermpi1 = 0_ti_p
      allocate (buffermpimu1(3,nrhs,max(npoleloc,1)))
      buffermpimu1 = 0_ti_p
      allocate (buffermpi2(3,nrhs,max(npolerecloc,1)))
      buffermpi2 = 0_ti_p
      allocate (buffermpimu2(3,nrhs,max(npolerecloc,1)))
      buffermpimu2 = 0_ti_p
c
      allocate (dipfield(3,nrhs,max(1,npoleloc)))
      dipfield = 0_ti_p
      allocate (dipfieldbis(3,nrhs,max(1,npolerecloc)))
      dipfieldbis = 0_ti_p
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
          diag(i) = polarity(iipole)
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
      res = 0_ti_p
      pp = 0_ti_p
      zr = 0_ti_p
      h = 0_ti_p

c
c     now, compute the initial direction
c
      ggold = 0_ti_p
c
      time0 = mpi_wtime()
      call tmatxbrecip(mu,murec,nrhs,dipfield,dipfieldbis)
      time1 = mpi_wtime()
      timerecdip = timerecdip + time1-time0
      call matvec(nrhs,.true.,mu,h)
      time2 = mpi_wtime()
      timerealdip = timerealdip + time2-time1
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

      term = (4.0_ti_p/3.0_ti_p) * aewald**3 / sqrtpi
      h(:,:,1:npoleloc) = h(:,:,1:npoleloc) + dipfield(:,:,1:npoleloc)
     $ - term*mu(:,:,1:npoleloc)
      res(:,:,1:npoleloc) = ef(:,:,1:npoleloc)-h(:,:,1:npoleloc)
cnew  do i = 1, 3
cnew     do j = 1,nrhs
cnew        do k = 1,npoleloc
cnew           res(i,j,k) = ef(k,i,j)-h(i,j,k)
cnew        enddo
cnew     enddo
cnew  enddo
      do k = 1, nrhs
        do j = 1, 3
          zr(j,k,1:npoleloc) = diag(1:npoleloc)*res(j,k,1:npoleloc)
        end do
      end do
      pp(:,:,1:npoleloc) = zr(:,:,1:npoleloc)
      do k = 1, nrhs
        ggold(k) = sum(res(:,k,:)*zr(:,k,:))
      end do
      call MPI_IALLREDUCE(MPI_IN_PLACE,ggold(1),nrhs,MPI_TPREC,MPI_SUM,
     $     COMM_TINKER,req1,ierr)
c
c     MPI : begin sending
c
      call commdirdir(nrhs,1,pp,reqrec,reqsend)
      call commrecdirdip(nrhs,1,murec,pp,buffermpimu1,
     $   buffermpimu2,req2rec,req2send)
      call commdirdir(nrhs,2,mu,reqrec,reqsend)
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
        time0 = mpi_wtime()
        call tmatxbrecip(pp,murec,nrhs,dipfield,dipfieldbis)
        time1 = mpi_wtime()
        call matvec(nrhs,.true.,pp,h)
        time2 = mpi_wtime()
        timerealdip = timerealdip + time2-time1
        timerecdip  = timerecdip  + time1-time0
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

        term = (4.0_ti_p/3.0_ti_p) * aewald**3 / sqrtpi
        h(:,:,1:npoleloc) = h(:,:,1:npoleloc) + dipfield(:,:,1:npoleloc)
     $   - term*pp(:,:,1:npoleloc)
        do k = 1, nrhs
          gg(k) = sum(pp(:,k,1:npoleloc)*h(:,k,1:npoleloc))
        end do
        call MPI_IALLREDUCE(MPI_IN_PLACE,gg(1),nrhs,MPI_TPREC,MPI_SUM,
     $    COMM_TINKER,req2,ierr)
        call MPI_WAIT(req2,status,ierr)
        do k = 1, nrhs
          if (gg(k).eq.zero) return
          alphacg(k)  = ggold(k)/gg(k)
          ggnew(k)  = zero
          ene(k)    = zero
        end do
        do k = 1, nrhs
          mu(:,k,1:npoleloc) = mu(:,k,1:npoleloc) + alphacg(k)*
     $     pp(:,k,1:npoleloc)
          res(:,k,1:npoleloc) = res(:,k,1:npoleloc) - alphacg(k)*
     $     h(:,k,1:npoleloc)
        end do
        do k = 1, nrhs
          do j = 1, 3
            zr(j,k,1:npoleloc) = diag(1:npoleloc)*res(j,k,1:npoleloc)
          end do
        end do
        do k = 1, nrhs
          ggnew(k) = sum(res(:,k,:)*zr(:,k,:))
          ene(k) = -pt5*sum(mu(1:3,k,1:npoleloc)*(res(1:3,k,1:npoleloc)+
     $      ef(:,k,1:npoleloc)))
cnew    do i = 1, 3
cnew    do j = 1,npoleloc
cnew       ggnew(k) = ggnew(k)+res(i,k,j)*zr(i,k,j)
cnew       ene(k) = ene(k) -pt5*mu(i,k,j)*(res(i,k,j)+ef(j,i,k))
cnew    end do
cnew    end do
        end do
        call MPI_IALLREDUCE(MPI_IN_PLACE,ggnew(1),nrhs,MPI_TPREC,
     $    MPI_SUM,COMM_TINKER,req3,ierr)
        call MPI_IALLREDUCE(MPI_IN_PLACE,ene(1),nrhs,MPI_TPREC,MPI_SUM,
     $    COMM_TINKER,req4,ierr)
        call MPI_WAIT(req3,status,ierr)
        call MPI_WAIT(req4,status,ierr)
        resnrm = zero
        do k = 1, nrhs
          gnorm(k) = sqrt(ggnew(k)/real(3*npolar,t_p))
          resnrm   = max(resnrm,gnorm(k))
        end do
        if (polprt.ge.2.and.rank.eq.0) write(iout,1010)
     $    it, (ene(k)*coulomb, gnorm(k), k = 1, nrhs)
        do k = 1, nrhs
          pp(:,k,1:npoleloc) = zr(:,k,1:npoleloc)+ggnew(k)/ggold(k)*
     $       pp(:,k,1:npoleloc)
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
      subroutine inducejac_pme2vec(matvec,nrhs,dodiis,ef,mu,murec)
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
      integer nrhs, info, proc
      real(t_p)  ef(3,nrhs,*), mu(3,nrhs,*), murec(3,nrhs,*)
cnew  real(t_p)  ef(max(1,npolebloc),3,nrhs), mu(3,nrhs,*), murec(3,nrhs,*)
      real(t_p), allocatable :: munew(:,:,:), h(:,:,:)
      real(t_p), allocatable :: xdiis(:,:), ediis(:,:), bmat(:,:),
     $   bloc(:,:), cex(:)
      integer i, j, k, it, ind, ndismx, nmat, lenb, tag
      logical dodiis
      parameter (ndismx=25)
      integer ipiv(ndismx+1)
      real(t_p)  zero, one, rnorm(2), rr, xx(1)
      real(t_p) term
      real(t_p) time0,time1,time2,time3
      save    zero, one, xx
      data    zero/0.0_ti_p/, one/1.0_ti_p/, xx/0.0_ti_p/
      external matvec
      real(t_p), allocatable :: dipfield(:,:,:),dipfieldbis(:,:,:)
c
c     MPI
c
      real(t_p), allocatable :: buffermpi1(:,:,:),buffermpi2(:,:,:)
      real(t_p), allocatable :: buffermpimu1(:,:,:),buffermpimu2(:,:,:)
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
      allocate (buffermpi1(3,nrhs,max(npoleloc,1)))
      buffermpi1 = 0_ti_p
      allocate (buffermpimu1(3,nrhs,max(npoleloc,1)))
      buffermpimu1 = 0_ti_p
      allocate (buffermpi2(3,nrhs,max(npolerecloc,1)))
      buffermpi2 = 0_ti_p
      allocate (buffermpimu2(3,nrhs,max(npolerecloc,1)))
      buffermpimu2 = 0_ti_p
c
      allocate (munew(3,nrhs,max(1,npolebloc)))
      munew = 0_ti_p
      allocate (dipfield(3,nrhs,max(1,npoleloc)))
      dipfield = 0_ti_p
      allocate (dipfieldbis(3,nrhs,max(1,npolerecloc)))
      dipfieldbis = 0_ti_p
      allocate (reqrecdirsend(nproc))
      allocate (reqrecdirrec(nproc))
      allocate (reqrec(nproc))
      allocate (reqsend(nproc))
      allocate (req2rec(nproc))
      allocate (req2send(nproc))
c
      allocate (h(3,nrhs,max(1,npolebloc)))
      h = 0_ti_p
      if (dodiis) then
        nmat = 1
        lenb = ndismx + 1
        allocate (xdiis(3*nrhs*max(npoleloc,1),ndismx))
        allocate (ediis(3*nrhs*max(npoleloc,1),ndismx))
        allocate (bmat(lenb,lenb))
        allocate (bloc(lenb,lenb))
        allocate (cex(lenb))
        bmat = 0_ti_p
      end if
c
c     main loop:
c
      do it = 1, politer
        h = 0_ti_p
        rnorm = 0_ti_p
c
        time0 = mpi_wtime()
        call tmatxbrecip(mu,murec,nrhs,dipfield,dipfieldbis)
        time1 = mpi_wtime()
        if (it.eq.1) timerecdip = timerecdip + time1-time0
        call matvec(nrhs,.false.,mu,h)
        time2 = mpi_wtime()
        if (it.eq.1) timerealdip = timerealdip + time2-time1
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
        time2 = mpi_wtime()
        h(1:3,1:nrhs,1:npoleloc) = h(1:3,1:nrhs,1:npoleloc) + 
     $   dipfield(1:3,1:nrhs,1:npoleloc)

        term = (4.0_ti_p/3.0_ti_p) * aewald**3 / sqrtpi
        h(1:3,1:nrhs,1:npoleloc) = h(1:3,1:nrhs,1:npoleloc) 
     $   - term*mu(1:3,1:nrhs,1:npoleloc)
c
        do i = 1, npoleloc
          iipole = poleglob(i)
          do k = 1, nrhs
            do j = 1, 3
              munew(j,k,i) = polarity(iipole)*(ef(j,k,i) - h(j,k,i))
cnew          munew(j,k,i) = polarity(iipole)*(ef(i,j,k) - h(j,k,i))
            end do
          end do
        end do
c
        do k = 1, nrhs
          rnorm(k) = sum((munew(1:3,k,1:npoleloc)-
     $      mu(1:3,k,1:npoleloc))**2)
        end do
c
        call MPI_IALLREDUCE(MPI_IN_PLACE,rnorm(1),nrhs,MPI_TPREC,
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
          cex = 0_ti_p
          cex(1) = one
          call dgesv(nmat,1,bloc,ndismx+1,ipiv,cex,nmat,info)
          munew = 0_ti_p
          call extrap(3*nrhs*npoleloc,nmat-1,xdiis,cex,munew)
        end if
        time2 = mpi_wtime()
c
        call commdirdir(nrhs,1,munew,reqrec,reqsend)
        call commrecdirdip(nrhs,1,murec,munew,buffermpimu1,
     $  buffermpimu2,req2rec,req2send)
        call commdirdir(nrhs,2,mu,reqrec,reqsend)
        call commrecdirdip(nrhs,2,murec,munew,buffermpimu1,
     $  buffermpimu2,req2rec,req2send)
        mu(:,:,1:npoleloc) = munew(:,:,1:npoleloc)

c
c     compute the norm of the increment.
c
        call MPI_WAIT(reqnorm,status,ierr)
        rr = zero
        do k = 1, nrhs
         rnorm(k) = sqrt(rnorm(k)/real(3*npolar,t_p))
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
      subroutine ulspredvec
      use mpole
      use tinheader
      use uprior
      implicit none
      integer i,j,k,m
      real(t_p) coeff,udk,upk
      real(t_p) amax,apmax
      real(t_p) b(maxualt)
      real(t_p) bp(maxualt)
      real(t_p) a(maxualt*(maxualt+1)/2)
      real(t_p) ap(maxualt*(maxualt+1)/2)
      real(t_p) c(maxualt,maxualt)
      real(t_p) cp(maxualt,maxualt)
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
            b(k) = 0.0_ti_p
            bp(k) = 0.0_ti_p
            do m = k, nualt
               c(k,m) = 0.0_ti_p
               cp(k,m) = 0.0_ti_p
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
         amax = 0.0_ti_p
         apmax = 0.0_ti_p
         do i = 1, k*(k+1)/2
            amax = max(amax,a(i))
            apmax = max(apmax,ap(i))
         end do
         if (amax .ne. 0.0_ti_p)  call cholesky (k,a,b)
         if (apmax .ne. 0.0_ti_p)  call cholesky (k,ap,bp)
c
c     transfer the final solution to the coefficient vector
c
         do k = 1, nualt-1
            bpred(k) = b(k)
            bpredp(k) = bp(k)
            bpreds(k) = b(k)
            bpredps(k) = bp(k)
         end do
         bpred(nualt) = 0.0_ti_p
         bpredp(nualt) = 0.0_ti_p
         bpreds(nualt) = 0.0_ti_p
         bpredps(nualt) = 0.0_ti_p
      end if
      return
      end
c
c     subroutine commfield : communicate some direct fields (Newton's third law)
c
      subroutine commfieldvec(nrhs,ef)
      use atoms
      use domdec
      use mpole
      use potent
      use mpi
      implicit none
      integer nrhs,i,j,k,tag,status(MPI_STATUS_SIZE)
      integer ierr,commloc
      integer, allocatable :: reqrec(:), reqsend(:)
      real(t_p) ef(3,nrhs,*)
cnew  real(t_p) ef(max(1,npolebloc),3,nrhs)
      real(t_p), allocatable :: buffer(:,:,:,:)
c
      if (use_pmecore) then
        commloc = comm_dir
      else
        commloc = COMM_TINKER
      end if
c
      allocate (buffer(3,nrhs,max(npoleloc,1),n_send1))
cnew  allocate (buffer(max(npoleloc,1),3,nrhs,n_send1))
      allocate (reqrec(nproc))
      allocate (reqsend(nproc))
c
      do i = 1, n_send1
        tag = nproc*rank + p_send1(i) + 1
        call MPI_IRECV(buffer(1,1,1,i),3*nrhs*npoleloc,MPI_TPREC,
     $   p_send1(i),tag,commloc,reqrec(i),ierr)
      end do
c
      do i = 1, n_recep1
        tag = nproc*p_recep1(i) + rank + 1
        call MPI_ISEND(ef(1,1,bufbegpole(p_recep1(i)+1)),
cnew    call MPI_ISEND(ef(bufbegpole(p_recep1(i)+1),1,1),
     $   3*nrhs*domlenpole(p_recep1(i)+1),MPI_TPREC,p_recep1(i),
     $  tag,commloc,reqsend(i),ierr)
      end do
c
      do i = 1, n_send1
        tag = nproc*rank + p_send1(i) + 1
        call MPI_WAIT(reqrec(i),status,ierr)
      end do
      do i = 1, n_recep1
        tag = nproc*p_recep1(i) + rank + 1
        call MPI_WAIT(reqsend(i),status,ierr)
      end do
c
      do i = 1, n_send1
        do j = 1, npoleloc
          do k = 1, nrhs
              ef(1,k,j) = ef(1,k,j) + buffer(1,k,j,i)
              ef(2,k,j) = ef(2,k,j) + buffer(2,k,j,i)
              ef(3,k,j) = ef(3,k,j) + buffer(3,k,j,i)
cnew        ef(j,1,k) = ef(j,1,k) + buffer(j,1,k,i)
cnew        ef(j,2,k) = ef(j,2,k) + buffer(j,2,k,i)
cnew        ef(j,3,k) = ef(j,3,k) + buffer(j,3,k,i)
          end do
        end do
      end do
c
      deallocate (buffer)
      deallocate (reqrec)
      deallocate (reqsend)
      return
      end
c
c     subroutine commfieldshort : communicate some short range direct fields (Newton's third law)
c
cnon  subroutine commfieldshort(nrhs,ef)
cnon  use atoms
cnon  use domdec
cnon  use mpole
cnon  use potent
cnon  use mpi
cnon  implicit none
cnon  integer nrhs,i,j,k,tag,status(MPI_STATUS_SIZE)
cnon  integer ierr,commloc
cnon  integer, allocatable :: reqrec(:), reqsend(:)
cnon  real(t_p) ef(3,nrhs,*)
cnon  real(t_p), allocatable :: buffer(:,:,:,:)
cnon
cnon  if (use_pmecore) then
cnon    commloc = comm_dir
cnon  else
cnon    commloc = COMM_TINKER
cnon  end if
cnon
cnon  allocate (buffer(3,nrhs,max(npoleloc,1),n_sendshort1))
cnon  allocate (reqrec(nproc))
cnon  allocate (reqsend(nproc))
cnon
cnon  do i = 1, n_sendshort1
cnon    tag = nproc*rank + p_sendshort1(i) + 1
cnon    call MPI_IRECV(buffer(1,1,1,i),3*nrhs*npoleloc,MPI_TPREC,
cnon $   p_sendshort1(i),tag,commloc,reqrec(i),ierr)
cnon  end do
cnon
cnon  do i = 1, n_recepshort1
cnon    tag = nproc*p_recepshort1(i) + rank + 1
cnon    call MPI_ISEND(ef(1,1,bufbegpole(p_recepshort1(i)+1)),
cnon $   3*nrhs*domlenpole(p_recepshort1(i)+1),MPI_TPREC,
cnon $   p_recepshort1(i),tag,commloc,reqsend(i),ierr)
cnon  end do
cnon
cnon  do i = 1, n_sendshort1
cnon    tag = nproc*rank + p_sendshort1(i) + 1
cnon    call MPI_WAIT(reqrec(i),status,ierr)
cnon  end do
cnon  do i = 1, n_recepshort1
cnon    tag = nproc*p_recepshort1(i) + rank + 1
cnon    call MPI_WAIT(reqsend(i),status,ierr)
cnon  end do
cnon
cnon  do i = 1, n_sendshort1
cnon    do j = 1, npoleloc
cnon      do k = 1, nrhs
cnon        ef(1,k,j) = ef(1,k,j) + buffer(1,k,j,i)
cnon        ef(2,k,j) = ef(2,k,j) + buffer(2,k,j,i)
cnon        ef(3,k,j) = ef(3,k,j) + buffer(3,k,j,i)
cnon      end do
cnon    end do
cnon  end do
cnon
cnon  deallocate (buffer)
cnon  deallocate (reqrec)
cnon  deallocate (reqsend)
cnon  return
cnon  end
c
!===================================================
!     sub diagvecvec
!===================================================
! Performs product of vector a with polarisabilities
!
      subroutine diagvecvec(nrhs, A, B)
      use atmlst
      use mpole
      use polar
      implicit none

      integer, intent(in) :: nrhs
      real(t_p), dimension(3,nrhs,npolebloc) :: A
      real(t_p), dimension(3,nrhs,npolebloc) :: B
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
