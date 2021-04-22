 
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
      subroutine newinduce_shortreal
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
      integer i, j, k, nrhs
c
c     MPI
c
      integer iipole
      integer, allocatable :: reqsend(:),reqrec(:)
      integer, allocatable :: req2send(:),req2rec(:)


c
      parameter (nrhs=2)
      real(t_p)  udsum, upsum

      real(t_p), allocatable :: ef(:,:,:), mu(:,:,:)
c
      external tmatxb_shortreal
c
      if (.not.use_polar) return
c
 1000 format(' illegal polalg in newinduce.')
c1010 format(' time for the ',a,F14.5)
c1020 format(' total elapsed time in newinduce: ',F14.5)
 
c
c     allocate some memory and clear the arrays:
c
      allocate (mu(3,nrhs,max(1,npolebloc)))
      mu = 0_ti_p
c
      allocate (ef(3,nrhs,max(1,npolebloc)))
      ef = 0_ti_p
c
      allocate (reqsend(nproc))
      allocate (reqrec(nproc))
      allocate (req2send(nproc))
      allocate (req2rec(nproc))
c
c     compute the electric fields:
c
ccc$    wtime0 = omp_get_wtime()
c
      call efld0_shortreal(nrhs,ef)
c
      call commfieldshort(nrhs,ef)
c
      call commdirdirshort(nrhs,0,mu,reqrec,reqsend)
c
ccc$    wtime1 = omp_get_wtime()
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
              udsum = udsum + bpred(k)*udshortalt(k,j,iipole)
              upsum = upsum + bpred(k)*upshortalt(k,j,iipole)
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
      call commdirdirshort(nrhs,1,mu,reqrec,reqsend)
      call commdirdirshort(nrhs,2,mu,reqrec,reqsend)
c
c     now, call the proper solver.
c
      if (polalg.eq.1) then
        call inducepcg_shortreal(tmatxb_shortreal,nrhs,.true.,ef,mu)
      else if (polalg.eq.2) then
        call inducejac_shortreal(tmatxb_shortreal,nrhs,.true.,ef,mu)
      else
         if (rank.eq.0) write(iout,1000) 
         call fatal
      end if
ccc$    wtime2 = omp_get_wtime()
c     if (polprt.ge.1.and.rank.eq.0) then
c       if (polprt.ge.2) then
c         write(iout,1010) 'fields:  ', wtime1-wtime0
c         write(iout,1010) 'dipoles: ', wtime2-wtime1
c       end if
c       write(iout,1020) wtime2 - wtime0
c     end if
c
         do i = 1, npolebloc
           iipole = poleglob(i)
           do j = 1, 3
             uind(j,iipole) = mu(j,1,i)
             uinp(j,iipole) = mu(j,2,i)
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
                     udshortalt(k,j,iipole) = udshortalt(k-1,j,iipole)
                     upshortalt(k,j,iipole) = upshortalt(k-1,j,iipole)
                  end do
                  udshortalt(1,j,iipole) = uind(j,iipole)
                  upshortalt(1,j,iipole) = uinp(j,iipole)
                end do
            end do
         end if
!      end if
      deallocate (reqrec)
      deallocate (reqsend)
      deallocate (req2rec)
      deallocate (req2send)
      deallocate (ef)
      deallocate (mu)
      return
      end
c
      subroutine inducepcg_shortreal(matvec,nrhs,precnd,ef,mu)
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
      real(t_p)  ef(3,nrhs,*), mu(3,nrhs,*)
      logical precnd
      real(t_p), allocatable :: res(:,:,:),h(:,:,:),pp(:,:,:),
     $          zr(:,:,:),diag(:)
      integer i, it, j, k
      real(t_p)  ggold(2), ggnew(2), gnorm(2), gg(2), alphacg(2),ene(2)
      real(t_p)  zero, pt5, one, resnrm
      save    zero, pt5, one
      data    zero/0.0_ti_p/, pt5/0.50_ti_p/, one/1.0_ti_p/
      external matvec
      real(t_p) time1,time2
c
c     MPI
c
      integer iglob, iipole, ierr
      integer req1, req2, req3, req4
      integer status(MPI_STATUS_SIZE)
      integer,allocatable :: reqrec(:),reqsend(:)
      integer,allocatable :: req2rec(:),req2send(:)
      integer,allocatable :: reqendrec(:),reqendsend(:)
      integer,allocatable :: req2endrec(:),req2endsend(:)
c
 1000 format(' cgiter shortreal converged after ',I3,' iterations.',/,
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
      time1 = mpi_wtime()
      call matvec(nrhs,.true.,mu,h)
      time2 = mpi_wtime()
      timerealdip = timerealdip + time2-time1
      call commfieldshort(nrhs,h)
c
      call commdirdirshort(nrhs,0,pp,reqrec,reqsend)
c
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
      call MPI_IALLREDUCE(MPI_IN_PLACE,ggold(1),nrhs,MPI_TPREC,MPI_SUM,
     $     COMM_TINKER,req1,ierr)
c
c     MPI : begin sending
c
      call commdirdirshort(nrhs,1,pp,reqrec,reqsend)
      call commdirdirshort(nrhs,2,mu,reqrec,reqsend)
      call MPI_WAIT(req1,status,ierr)
c
c     now, start the main loop:
c
      do it = 1, politer
        do k = 1, nrhs
          gg(k) = zero
        end do
c
        call matvec(nrhs,.true.,pp,h)
        call commfieldshort(nrhs,h)
c
c     MPI : begin reception
c
        call commdirdirshort(nrhs,0,pp,reqrec,reqsend)

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
          ene(k) = -pt5*sum(mu(:,k,1:npoleloc)*(res(:,k,1:npoleloc)+
     $      ef(:,k,1:npoleloc)))
        end do
        call MPI_IALLREDUCE(MPI_IN_PLACE,ggnew(1),nrhs,MPI_TPREC,
     $    MPI_SUM,COMM_TINKER,req3,ierr)
        call MPI_IALLREDUCE(MPI_IN_PLACE,ene(1),nrhs,MPI_TPREC,MPI_SUM,
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
        do k = 1, nrhs
          pp(:,k,1:npoleloc) = zr(:,k,1:npoleloc)+ggnew(k)/ggold(k)*
     $       pp(:,k,1:npoleloc)
        end do
c
        call commdirdirshort(nrhs,1,pp,reqrec,reqsend)
c
        call commdirdirshort(nrhs,2,pp,reqrec,reqsend)
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
      call commdirdirshort(nrhs,0,mu,reqendrec,reqendsend)
c
c     MPI : begin sending
c
      call commdirdirshort(nrhs,1,mu,reqendrec,reqendsend)
c
c
      call commdirdirshort(nrhs,2,mu,reqendrec,reqendsend)
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
      deallocate (reqsend)
      deallocate (reqrec)
      deallocate (req2send)
      deallocate (req2rec)
      deallocate (reqendsend)
      deallocate (reqendrec)
      deallocate (req2endsend)
      deallocate (req2endrec)
      deallocate (res)
      deallocate (h)
      deallocate (pp)
      deallocate (zr)
      deallocate (diag)
      return
      end
c
      subroutine inducejac_shortreal(matvec,nrhs,dodiis,ef,mu)
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
      real(t_p)  ef(3,nrhs,*), mu(3,nrhs,*)
      real(t_p), allocatable :: munew(:,:,:), h(:,:,:)
      real(t_p), allocatable :: xdiis(:,:), ediis(:,:), bmat(:,:),
     $   bloc(:,:), cex(:)
      integer i, j, k, it, ind, ndismx, nmat, lenb
      logical dodiis
      parameter (ndismx=25)
      integer ipiv(ndismx+1)
      real(t_p)  zero, one, rnorm(2), rr, xx(1)

      real(t_p) time1,time2
      save    zero, one, xx
      external matvec
c
c     MPI
c
      integer iglob, iipole, ierr
      integer reqnorm,reqdiis(2*ndismx+1)
      integer status(MPI_STATUS_SIZE)
      integer,allocatable :: reqrec(:),reqsend(:)
      integer,allocatable :: req2rec(:),req2send(:)
c
 1000 format(' itsolv shortreal converged after ',I3,' iterations.',/,
     $       ' final residual norm = ',3D14.7)
 1010 format(' residual norm at iteration ',I3,':',3D12.2)
 1020 format(' Jacobi/DIIS solver: induced dipoles',/,
     $  ' ipole       mux         muy         muz')
 1021 format(' Jacobi/DIIS solver: induced p-dipoles',/,
     $  ' ipole       mux         muy         muz')
 1030 format(i6,2x,f10.7,2x,f10.7,2x,f10.7)
c
c
      zero  = 0.0_ti_p
      one   = 1.0_ti_p
      xx(1) = 0.0_ti_p

      allocate (munew(3,nrhs,max(1,npolebloc)))
      munew = 0_ti_p
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
        time1 = mpi_wtime()
        call matvec(nrhs,.false.,mu,h)
        time2 = mpi_wtime()
        if (it.eq.1) timerealdip = timerealdip + time2-time1
        call commfieldshort(nrhs,h)
c
        call commdirdirshort(nrhs,0,mu,reqrec,reqsend)
c
c     jacobi step:
c
        time2 = mpi_wtime()
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
#ifdef _OPENACC
          call fatal_device("inducejac_shortreal")
#else
          call M_gesv(nmat,1,bloc,ndismx+1,ipiv,cex,nmat,info)
#endif
          munew = 0_ti_p
          call extrap(3*nrhs*npoleloc,nmat-1,xdiis,cex,munew)
        end if
        time2 = mpi_wtime()
c
        call commdirdirshort(nrhs,1,munew,reqrec,reqsend)
        call commdirdirshort(nrhs,2,mu,reqrec,reqsend)
        mu(:,:,1:npoleloc) = munew(:,:,1:npoleloc)

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
      deallocate (munew)
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
      subroutine efld0_shortreal(nrhs,ef)
c
c     Compute the direct space contribution to the permanent electric field.
c      Also compute the "p" field, which is used to
c     compute the energy according to the AMOEBA force field.
c
      use atmlst
      use atoms
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
      use shunt
      use tinheader ,only:ti_p,re_p
      implicit none
      integer i,iglob,kglob,kbis,nrhs
      real(t_p)  ef(3,nrhs,npolebloc)
      integer ii, j, k, kkk, iipole, kkpole, ipoleloc
      real(t_p)  dx, dy, dz, d, d2, d3, d5, pdi, pti, ck,
     $  dkx, dky, dkz, dkr, qkxx, qkyy, qkzz, qkxy, qkxz, qkyz,
     $  qkx, qky, qkz, qkr
      real(t_p) dix, diy, diz, dir, qixx, qiyy, qizz, qixy, qixz, qiyz,
     $  qix, qiy, qiz, qir, ci
      real(t_p) damp, pgamma, expdamp, scale3, scale5, scale7
      real(t_p) drr3,drr5,drr7
      real(t_p) prr3,prr5,prr7
      real(t_p) dsc3,dsc5,dsc7
      real(t_p) psc3,psc5,psc7
      real(t_p) zero, pt6, one, f50
      real(t_p) ralpha, alsq2, alsq2n, bfac, exp2a
      real(t_p) bn(0:3), fim(3), fid(3), fip(3)
      real(t_p) fkm(3), fkd(3), fkp(3)
      real(t_p) cutoff2
      real(t_p), allocatable :: dscale(:)
      real(t_p), allocatable :: pscale(:)

      save   zero, pt6, one, f50
      data   zero/0._ti_p/, pt6/0.6_ti_p/, one/1._ti_p/,f50/50._ti_p/
      character*10 mode
c
 1000 format(' Warning, system moved too much since last neighbor list
     $  update, try lowering nlupdate')
c
      mode = 'SHORTEWALD'
      call switch (mode)
      cutoff2 = cut2
c
      allocate (dscale(n))
      allocate (pscale(n))
      dscale = 1.0_ti_p
      pscale = 1.0_ti_p
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
        pdi = pdamp(iipole)
        pti = thole(iipole)
        ci   = rpole(1,iipole)
        dix  = rpole(2,iipole)
        diy  = rpole(3,iipole)
        diz  = rpole(4,iipole)
        qixx = rpole(5,iipole)
        qixy = rpole(6,iipole)
        qixz = rpole(7,iipole)
        qiyy = rpole(9,iipole)
        qiyz = rpole(10,iipole)
        qizz = rpole(13,iipole)
        do j = 1, n12(iglob)
           pscale(i12(j,iglob)) = p2scale
        end do
        do j = 1, n13(iglob)
           pscale(i13(j,iglob)) = p3scale
        end do
        do j = 1, n14(iglob)
           pscale(i14(j,iglob)) = p4scale
           do k = 1, np11(iglob)
              if (i14(j,iglob) .eq. ip11(k,iglob))
     &           pscale(i14(j,iglob)) = p4scale * p41scale
           end do
        end do
        do j = 1, n15(iglob)
           pscale(i15(j,iglob)) = p5scale
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
c
        do kkk = 1, nshortelst(ii)
          kkpole = shortelst(kkk,ii)
          kbis = poleloc(kkpole)
          kglob = ipole(kkpole)
          if ((kbis.eq.0).or.(kbis.gt.npolebloc)) then
            write(iout,1000)
            cycle
          end if
          dx = x(kglob) - x(iglob)
          dy = y(kglob) - y(iglob)
          dz = z(kglob) - z(iglob)
          call image(dx,dy,dz)
          d2 = dx*dx + dy*dy + dz*dz
          if (d2.le.off2) then
            d  = sqrt(d2)
c
c     calculate the error function damping terms
c
            ralpha = aewald * d
            bn(0) = erfc(ralpha) / d
            alsq2 = 2.0_ti_p * aewald**2
            alsq2n = 0.0_ti_p
            if (aewald .gt. 0.0_ti_p)
     &        alsq2n = 1.0_ti_p / (sqrtpi*aewald)
            exp2a = exp(-ralpha**2)
            do j = 1, 3
              bfac = real(j+j-1,t_p)
              alsq2n = alsq2 * alsq2n
              bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / d2
            end do
c
            d3   = d*d2
            d5   = d3*d2
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
            damp = pdi*pdamp(kkpole)
            scale3 = one
            scale5 = one
            scale7 = one
            if (damp.ne.zero) then
              pgamma = min(pti,thole(kkpole))
              damp = -pgamma*(d/damp)**3
              if (damp.gt.-F50) then
                expdamp = exp(damp)
                scale3 = one - expdamp
                scale5 = one - expdamp*(one-damp)
                scale7 = one - expdamp
     &                      *(one-damp + pt6*damp**2)
              end if
            end if
            dsc3 = scale3 * dscale(kglob)
            dsc5 = scale5 * dscale(kglob)
            dsc7 = scale7 * dscale(kglob)
            psc3 = scale3 * pscale(kglob)
            psc5 = scale5 * pscale(kglob)
            psc7 = scale7 * pscale(kglob)
            drr3 = (1.0_ti_p-dsc3) / (d*d2)
            drr5 = 3.0_ti_p * (1.0_ti_p-dsc5) / (d*d2*d2)
            drr7 = 15.0_ti_p * (1.0_ti_p-dsc7) / (d*d2*d2*d2)
            prr3 = (1.0_ti_p-psc3) / (d*d2)
            prr5 = 3.0_ti_p * (1.0_ti_p-psc5) / (d*d2*d2)
            prr7 = 15.0_ti_p * (1.0_ti_p-psc7) / (d*d2*d2*d2)
c
c     compute some intermediate quantities
c
            dir = dix*dx + diy*dy + diz*dz
            qix = qixx*dx + qixy*dy + qixz*dz
            qiy = qixy*dx + qiyy*dy + qiyz*dz
            qiz = qixz*dx + qiyz*dy + qizz*dz
            qir = qix*dx + qiy*dy + qiz*dz
            dkr = dkx*dx + dky*dy + dkz*dz
            qkx = qkxx*dx + qkxy*dy + qkxz*dz
            qky = qkxy*dx + qkyy*dy + qkyz*dz
            qkz = qkxz*dx + qkyz*dy + qkzz*dz
            qkr = qkx*dx + qky*dy + qkz*dz
c
            fim(1) = -dx*(bn(1)*ck-bn(2)*dkr+bn(3)*qkr)
     &                   - bn(1)*dkx + 2.0_ti_p*bn(2)*qkx
            fim(2) = -dy*(bn(1)*ck-bn(2)*dkr+bn(3)*qkr)
     &                   - bn(1)*dky + 2.0_ti_p*bn(2)*qky
            fim(3) = -dz*(bn(1)*ck-bn(2)*dkr+bn(3)*qkr)
     &                   - bn(1)*dkz + 2.0_ti_p*bn(2)*qkz
            fkm(1) = dx*(bn(1)*ci+bn(2)*dir+bn(3)*qir)
     &                  - bn(1)*dix - 2.0_ti_p*bn(2)*qix
            fkm(2) = dy*(bn(1)*ci+bn(2)*dir+bn(3)*qir)
     &                  - bn(1)*diy - 2.0_ti_p*bn(2)*qiy
            fkm(3) = dz*(bn(1)*ci+bn(2)*dir+bn(3)*qir)
     &                  - bn(1)*diz - 2.0_ti_p*bn(2)*qiz
            fid(1) = -dx*(drr3*ck-drr5*dkr+drr7*qkr)
     &                   - drr3*dkx + 2.0_ti_p*drr5*qkx
            fid(2) = -dy*(drr3*ck-drr5*dkr+drr7*qkr)
     &                  - drr3*dky + 2.0_ti_p*drr5*qky
            fid(3) = -dz*(drr3*ck-drr5*dkr+drr7*qkr)
     &                   - drr3*dkz + 2.0_ti_p*drr5*qkz
            fkd(1) = dx*(drr3*ci+drr5*dir+drr7*qir)
     &                  - drr3*dix - 2.0_ti_p*drr5*qix
            fkd(2) = dy*(drr3*ci+drr5*dir+drr7*qir)
     &                  - drr3*diy - 2.0_ti_p*drr5*qiy
            fkd(3) = dz*(drr3*ci+drr5*dir+drr7*qir)
     &                     - drr3*diz - 2.0_ti_p*drr5*qiz
            fip(1) = -dx*(prr3*ck-prr5*dkr+prr7*qkr)
     &                   - prr3*dkx + 2.0_ti_p*prr5*qkx
            fip(2) = -dy*(prr3*ck-prr5*dkr+prr7*qkr)
     &                   - prr3*dky + 2.0_ti_p*prr5*qky
            fip(3) = -dz*(prr3*ck-prr5*dkr+prr7*qkr)
     &                   - prr3*dkz + 2.0_ti_p*prr5*qkz
            fkp(1) = dx*(prr3*ci+prr5*dir+prr7*qir)
     &                  - prr3*dix - 2.0_ti_p*prr5*qix
            fkp(2) = dy*(prr3*ci+prr5*dir+prr7*qir)
     &                  - prr3*diy - 2.0_ti_p*prr5*qiy
            fkp(3) = dz*(prr3*ci+prr5*dir+prr7*qir)
     &                  - prr3*diz - 2.0_ti_p*prr5*qiz

            ef(1,1,ipoleloc) = ef(1,1,ipoleloc) + (fim(1) - fid(1))
            ef(2,1,ipoleloc) = ef(2,1,ipoleloc) + (fim(2) - fid(2))
            ef(3,1,ipoleloc) = ef(3,1,ipoleloc) + (fim(3) - fid(3))
            ef(1,2,ipoleloc) = ef(1,2,ipoleloc) + (fim(1) - fip(1))
            ef(2,2,ipoleloc) = ef(2,2,ipoleloc) + (fim(2) - fip(2))
            ef(3,2,ipoleloc) = ef(3,2,ipoleloc) + (fim(3) - fip(3))
c
            ef(1,1,kbis) = ef(1,1,kbis) + (fkm(1) - fkd(1))
            ef(2,1,kbis) = ef(2,1,kbis) + (fkm(2) - fkd(2))
            ef(3,1,kbis) = ef(3,1,kbis) + (fkm(3) - fkd(3))
            ef(1,2,kbis) = ef(1,2,kbis) + (fkm(1) - fkp(1))
            ef(2,2,kbis) = ef(2,2,kbis) + (fkm(2) - fkp(2))
            ef(3,2,kbis) = ef(3,2,kbis) + (fkm(3) - fkp(3))
          end if
        end do
c
c     reset interaction scaling coefficients for connected atoms
c
        do j = 1, n12(iglob)
           pscale(i12(j,iglob)) = 1.0_ti_p
        end do
        do j = 1, n13(iglob)
           pscale(i13(j,iglob)) = 1.0_ti_p
        end do
        do j = 1, n14(iglob)
           pscale(i14(j,iglob)) = 1.0_ti_p
        end do
        do j = 1, n15(iglob)
           pscale(i15(j,iglob)) = 1.0_ti_p
        end do
        do j = 1, np11(iglob)
           dscale(ip11(j,iglob)) = 1.0_ti_p
        end do
        do j = 1, np12(iglob)
           dscale(ip12(j,iglob)) = 1.0_ti_p
        end do
        do j = 1, np13(iglob)
           dscale(ip13(j,iglob)) = 1.0_ti_p
        end do
        do j = 1, np14(iglob)
           dscale(ip14(j,iglob)) = 1.0_ti_p
        end do
      end do
c
      deallocate (dscale)
      deallocate (pscale)
c
      return
      end
c
      subroutine tmatxb_shortreal(nrhs,dodiag,mu,efi)
c
c     Compute the direct space contribution to the electric field due to the current value
c     of the induced dipoles
c
      use atmlst
      use atoms
      use domdec
      use ewald
      use math
      use mpole
      use neigh
      use polar
      use polgrp
      use polpot
      use shunt
      use tinheader ,only:ti_p,re_p
      implicit none
      integer i,iipole,nrhs,iglob,kglob,kkpole,kkpoleloc
      integer ipoleloc
      real(t_p)  mu(3,nrhs,*), efi(3,nrhs,npolebloc)
      logical dodiag
      integer j, ii, kkk, irhs
      real(t_p)  dx, dy, dz, d, d2, damp, expdamp, pgamma,
     $  scale3, scale5, pdi, pti 
      real(t_p) ralpha, alsq2, alsq2n, bfac, exp2a
      real(t_p) rr3, rr5, dukx, duky, dukz, pukx, puky, pukz,
     $  puir, pukr, duir, dukr
      real(t_p) duix, duiy, duiz, puix, puiy, puiz
      real(t_p) bn(0:3), fid(3), fip(3), fimd(3), fimp(3)
      real(t_p) fkd(3), fkp(3), fkmd(3), fkmp(3)
      real(t_p)  zero, one, f50
      real(t_p), allocatable :: dscale(:)
      real(t_p)  cutoff2
      save    zero, one, f50
      data    zero/0._ti_p/, one/1._ti_p/, f50/50._ti_p/
      character*10 mode
c
c     initialize the result vector
c
      allocate (dscale(n))
      do i = 1, n
         dscale(i) = 1.0_ti_p
      end do
c
      do i = 1, npolebloc
        do irhs = 1, nrhs
          do j = 1, 3
            efi(j,irhs,i) = zero
          end do
        end do
      end do
c
c     gather some parameters, then set up the damping factors.
c
      mode = 'SHORTEWALD'
      call switch (mode)
      cutoff2 = cut2
c
      do ii = 1, npolelocnl
        iipole = poleglobnl(ii)
        iglob = ipole(iipole)
        i = loc(iglob)
        ipoleloc = poleloc(iipole)
        if (i.eq.0) cycle
        pdi = pdamp(iipole)
        pti = thole(iipole)
        duix = mu(1,1,ipoleloc)
        duiy = mu(2,1,ipoleloc)
        duiz = mu(3,1,ipoleloc)
        puix = mu(1,2,ipoleloc)
        puiy = mu(2,2,ipoleloc)
        puiz = mu(3,2,ipoleloc)
        do j = 1, np11(iglob)
           dscale(ip11(j,iglob)) = u1scale
        end do
        do j = 1, np12(iglob)
           dscale(ip12(j,iglob)) = u2scale
        end do
        do j = 1, np13(iglob)
           dscale(ip13(j,iglob)) = u3scale
        end do
        do j = 1, np14(iglob)
           dscale(ip14(j,iglob)) = u4scale
        end do
c
        do kkk = 1,nshortelst(ii)
          kkpole = shortelst(kkk,ii)
          kkpoleloc = poleloc(kkpole)
          kglob = ipole(kkpole)
          if (kkpoleloc.eq.0) cycle
          dx = x(kglob) - x(iglob)
          dy = y(kglob) - y(iglob)
          dz = z(kglob) - z(iglob)
          call image(dx,dy,dz)
          d2 = dx*dx + dy*dy + dz*dz
          if (d2.le.off2) then
c
c     compute the distances and the scaling factors according to
c     Thole's model.
c
            d  = sqrt(d2)
            dukx = mu(1,1,kkpoleloc)
            duky = mu(2,1,kkpoleloc)
            dukz = mu(3,1,kkpoleloc)
            pukx = mu(1,2,kkpoleloc)
            puky = mu(2,2,kkpoleloc)
            pukz = mu(3,2,kkpoleloc)
c
            ralpha = aewald * d
            bn(0) = erfc(ralpha) / d
            alsq2 = 2.0_ti_p * aewald**2
            alsq2n = 0.0_ti_p
            if (aewald .gt. 0.0_ti_p)
     &        alsq2n = 1.0_ti_p / (sqrtpi*aewald)
            exp2a = exp(-ralpha**2)
            do j = 1, 2
              bfac = dble(j+j-1)
              alsq2n = alsq2 * alsq2n
              bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / d2
            end do
c
            scale3 = dscale(kglob)
            scale5 = dscale(kglob)
            damp = pdi*pdamp(kkpole)
            if (damp.ne.zero) then
              pgamma = min(pti,thole(kkpole))
              damp = -pgamma*(d/damp)**3
              if (damp .gt. -f50) then
                expdamp = exp(damp)
                scale3 = scale3 * (one - expdamp)
                scale5 = scale5 * (one - expdamp*(one - damp))
              end if
            end if
c
c     compute the field.
c
            rr3 = (1.0_ti_p-scale3) / (d*d2)
            rr5 = 3.0_ti_p * (1.0_ti_p-scale5) / (d*d2*d2)
            duir = dx*duix + dy*duiy + dz*duiz
            dukr = dx*dukx + dy*duky + dz*dukz
            puir = dx*puix + dy*puiy + dz*puiz
            pukr = dx*pukx + dy*puky + dz*pukz
            fimd(1) = -bn(1)*dukx + bn(2)*dukr*dx
            fimd(2) = -bn(1)*duky + bn(2)*dukr*dy
            fimd(3) = -bn(1)*dukz + bn(2)*dukr*dz
            fkmd(1) = -bn(1)*duix + bn(2)*duir*dx
            fkmd(2) = -bn(1)*duiy + bn(2)*duir*dy
            fkmd(3) = -bn(1)*duiz + bn(2)*duir*dz
            fimp(1) = -bn(1)*pukx + bn(2)*pukr*dx
            fimp(2) = -bn(1)*puky + bn(2)*pukr*dy
            fimp(3) = -bn(1)*pukz + bn(2)*pukr*dz
            fkmp(1) = -bn(1)*puix + bn(2)*puir*dx
            fkmp(2) = -bn(1)*puiy + bn(2)*puir*dy
            fkmp(3) = -bn(1)*puiz + bn(2)*puir*dz
            fid(1) = -rr3*dukx + rr5*dukr*dx
            fid(2) = -rr3*duky + rr5*dukr*dy
            fid(3) = -rr3*dukz + rr5*dukr*dz
            fkd(1) = -rr3*duix + rr5*duir*dx
            fkd(2) = -rr3*duiy + rr5*duir*dy
            fkd(3) = -rr3*duiz + rr5*duir*dz
            fip(1) = -rr3*pukx + rr5*pukr*dx
            fip(2) = -rr3*puky + rr5*pukr*dy
            fip(3) = -rr3*pukz + rr5*pukr*dz
            fkp(1) = -rr3*puix + rr5*puir*dx
            fkp(2) = -rr3*puiy + rr5*puir*dy
            fkp(3) = -rr3*puiz + rr5*puir*dz
            do j = 1,3
              efi(j,1,ipoleloc) = efi(j,1,ipoleloc) - fimd(j) + fid(j)
              efi(j,2,ipoleloc) = efi(j,2,ipoleloc) - fimp(j) + fip(j)
              efi(j,1,kkpoleloc) = efi(j,1,kkpoleloc) - fkmd(j) + fkd(j)
              efi(j,2,kkpoleloc) = efi(j,2,kkpoleloc) - fkmp(j) + fkp(j)
            end do
          end if
        end do
c
c     reset interaction scaling coefficients for connected atoms
c
        do j = 1, np11(iglob)
           dscale(ip11(j,iglob)) = 1.0_ti_p
        end do
        do j = 1, np12(iglob)
           dscale(ip12(j,iglob)) = 1.0_ti_p
        end do
        do j = 1, np13(iglob)
           dscale(ip13(j,iglob)) = 1.0_ti_p
        end do
        do j = 1, np14(iglob)
           dscale(ip14(j,iglob)) = 1.0_ti_p
        end do
      end do
      if(dodiag) then
c
c     if dodiag is true, also compute the "self-induced" field,
c     i.e., the diagonal portion of the matrix/vector product.
c
        do i = 1, npoleloc
          iipole = poleglob(i)
          if (polarity(iipole) == 0_ti_p) then
             cycle
          else
            do irhs = 1, nrhs
              do j = 1, 3
                efi(j,irhs,i) = efi(j,irhs,i) +
     $             mu(j,irhs,i)/polarity(iipole)
              end do
            end do
          end if
        end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (dscale)
c
      return
      end

      subroutine dbletmatxb_shortreal(nrhs,dodiag,mu, mu2, efi, efi2)
c
c     Compute the direct space contribution to the electric field due to the current value
c     of the induced dipoles 
c     Now on two sets of 'dipoles' at once
c
      use atmlst
      use atoms
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
      use tinheader ,only:ti_p,re_p
      implicit none
      integer i,iipole,nrhs,iglob,kglob,kkpole,kkpoleloc
      integer ipoleloc
      real(t_p)  mu(3,nrhs,*), efi(3,nrhs,npolebloc)
      real(t_p)  mu2(3,nrhs,npolebloc), efi2(3,nrhs,npolebloc)
      logical dodiag
      integer j, ii, kkk, irhs
      real(t_p)  dx, dy, dz, d, d2, damp, expdamp, pgamma,
     $  scale3, scale5, pdi, pti 
      real(t_p) ralpha, alsq2, alsq2n, bfac, exp2a
      real(t_p) rr3, rr5, dukx, duky, dukz, pukx, puky, pukz,
     $  puir, pukr, duir, dukr
      real(t_p) :: dukx2, duky2, dukz2, pukx2, puky2, pukz2,
     $          puir2, pukr2,  duir2, dukr2
      real(t_p) duix, duiy, duiz, puix, puiy, puiz
      real(t_p) :: duix2, duiy2, duiz2, puix2, puiy2, puiz2
      real(t_p) bn(0:3)
      real(t_p), dimension(3,2) :: fid, fip, fimd, fimp, fkd, fkp,
     $                          fkmd, fkmp
      real(t_p)  zero, one, f50
      real(t_p), allocatable :: dscale(:)
      real(t_p)  cutoff2
      save    zero, one, f50
      data    zero/0._ti_p/, one/1._ti_p/, f50/50._ti_p/
      character*10 mode
c
c     initialize the result vector
c
      allocate (dscale(n))
      do i = 1, n
         dscale(i) = 1.0_ti_p
      end do
c
      do i = 1, npolebloc
        do irhs = 1, nrhs
          do j = 1, 3
            efi(j,irhs,i) = zero
            efi2(j,irhs,i) = zero
          end do
        end do
      end do
c
c     gather some parameters, then set up the damping factors.
c
      mode = 'SHORTEWALD'
      call switch (mode)
      cutoff2 = cut2
c
      do ii = 1, npolelocnl
        iipole = poleglobnl(ii)
        iglob = ipole(iipole)
        i = loc(iglob)
        ipoleloc = poleloc(iipole)
        if (i.eq.0) cycle
        pdi = pdamp(iipole)
        pti = thole(iipole)
        duix = mu(1,1,ipoleloc)
        duiy = mu(2,1,ipoleloc)
        duiz = mu(3,1,ipoleloc)
        puix = mu(1,2,ipoleloc)
        puiy = mu(2,2,ipoleloc)
        puiz = mu(3,2,ipoleloc)
        duix2 = mu2(1,1,ipoleloc)
        duiy2 = mu2(2,1,ipoleloc)
        duiz2 = mu2(3,1,ipoleloc)
        puix2 = mu2(1,2,ipoleloc)
        puiy2 = mu2(2,2,ipoleloc)
        puiz2 = mu2(3,2,ipoleloc)
        do j = 1, np11(iglob)
           dscale(ip11(j,iglob)) = u1scale
        end do
        do j = 1, np12(iglob)
           dscale(ip12(j,iglob)) = u2scale
        end do
        do j = 1, np13(iglob)
           dscale(ip13(j,iglob)) = u3scale
        end do
        do j = 1, np14(iglob)
           dscale(ip14(j,iglob)) = u4scale
        end do
c
        do kkk = 1,nshortelst(ii)
          kkpole = shortelst(kkk,ii)
          kkpoleloc = poleloc(kkpole)
          kglob = ipole(kkpole)
          if (kkpoleloc.eq.0) cycle
          dx = x(kglob) - x(iglob)
          dy = y(kglob) - y(iglob)
          dz = z(kglob) - z(iglob)
          call image(dx,dy,dz)
          d2 = dx*dx + dy*dy + dz*dz
          if (d2.le.off2) then
c
c     compute the distances and the scaling factors according to
c     Thole's model.
c
            d  = sqrt(d2)
            dukx = mu(1,1,kkpoleloc)
            duky = mu(2,1,kkpoleloc)
            dukz = mu(3,1,kkpoleloc)
            pukx = mu(1,2,kkpoleloc)
            puky = mu(2,2,kkpoleloc)
            pukz = mu(3,2,kkpoleloc)
            dukx2 = mu2(1,1,kkpoleloc)
            duky2 = mu2(2,1,kkpoleloc)
            dukz2 = mu2(3,1,kkpoleloc)
            pukx2 = mu2(1,2,kkpoleloc)
            puky2 = mu2(2,2,kkpoleloc)
            pukz2 = mu2(3,2,kkpoleloc)
c
            ralpha = aewald * d
            bn(0) = erfc(ralpha) / d
            alsq2 = 2.0_ti_p * aewald**2
            alsq2n = 0.0_ti_p
            if (aewald .gt. 0.0_ti_p)
     &        alsq2n = 1.0_ti_p / (sqrtpi*aewald)
            exp2a = exp(-ralpha**2)
            do j = 1, 2
              bfac = dble(j+j-1)
              alsq2n = alsq2 * alsq2n
              bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / d2
            end do
c
            scale3 = dscale(kglob)
            scale5 = dscale(kglob)
            damp = pdi*pdamp(kkpole)
            if (damp.ne.zero) then
              pgamma = min(pti,thole(kkpole))
              damp = -pgamma*(d/damp)**3
              if (damp .gt. -f50) then
                expdamp = exp(damp)
                scale3 = scale3 * (one - expdamp)
                scale5 = scale5 * (one - expdamp*(one - damp))
              end if
            end if
c
c     compute the field.
c
            rr3 = (1.0_ti_p-scale3) / (d*d2)
            rr5 = 3.0_ti_p * (1.0_ti_p-scale5) / (d*d2*d2)
            duir = dx*duix + dy*duiy + dz*duiz
            dukr = dx*dukx + dy*duky + dz*dukz
            puir = dx*puix + dy*puiy + dz*puiz
            pukr = dx*pukx + dy*puky + dz*pukz
            duir2 = dx*duix2 + dy*duiy2 + dz*duiz2
            dukr2 = dx*dukx2 + dy*duky2 + dz*dukz2
            puir2 = dx*puix2 + dy*puiy2 + dz*puiz2
            pukr2 = dx*pukx2 + dy*puky2 + dz*pukz2
            fimd(1,1) = -bn(1)*dukx + bn(2)*dukr*dx
            fimd(2,1) = -bn(1)*duky + bn(2)*dukr*dy
            fimd(3,1) = -bn(1)*dukz + bn(2)*dukr*dz
            fkmd(1,1) = -bn(1)*duix + bn(2)*duir*dx
            fkmd(2,1) = -bn(1)*duiy + bn(2)*duir*dy
            fkmd(3,1) = -bn(1)*duiz + bn(2)*duir*dz
            fimp(1,1) = -bn(1)*pukx + bn(2)*pukr*dx
            fimp(2,1) = -bn(1)*puky + bn(2)*pukr*dy
            fimp(3,1) = -bn(1)*pukz + bn(2)*pukr*dz
            fkmp(1,1) = -bn(1)*puix + bn(2)*puir*dx
            fkmp(2,1) = -bn(1)*puiy + bn(2)*puir*dy
            fkmp(3,1) = -bn(1)*puiz + bn(2)*puir*dz
            fimd(1,2) = -bn(1)*dukx2 + bn(2)*dukr2*dx
            fimd(2,2) = -bn(1)*duky2 + bn(2)*dukr2*dy
            fimd(3,2) = -bn(1)*dukz2 + bn(2)*dukr2*dz
            fkmd(1,2) = -bn(1)*duix2 + bn(2)*duir2*dx
            fkmd(2,2) = -bn(1)*duiy2 + bn(2)*duir2*dy
            fkmd(3,2) = -bn(1)*duiz2 + bn(2)*duir2*dz
            fimp(1,2) = -bn(1)*pukx2 + bn(2)*pukr2*dx
            fimp(2,2) = -bn(1)*puky2 + bn(2)*pukr2*dy
            fimp(3,2) = -bn(1)*pukz2 + bn(2)*pukr2*dz
            fkmp(1,2) = -bn(1)*puix2 + bn(2)*puir2*dx
            fkmp(2,2) = -bn(1)*puiy2 + bn(2)*puir2*dy
            fkmp(3,2) = -bn(1)*puiz2 + bn(2)*puir2*dz
            fid(1,1) = -rr3*dukx + rr5*dukr*dx
            fid(2,1) = -rr3*duky + rr5*dukr*dy
            fid(3,1) = -rr3*dukz + rr5*dukr*dz
            fkd(1,1) = -rr3*duix + rr5*duir*dx
            fkd(2,1) = -rr3*duiy + rr5*duir*dy
            fkd(3,1) = -rr3*duiz + rr5*duir*dz
            fip(1,1) = -rr3*pukx + rr5*pukr*dx
            fip(2,1) = -rr3*puky + rr5*pukr*dy
            fip(3,1) = -rr3*pukz + rr5*pukr*dz
            fkp(1,1) = -rr3*puix + rr5*puir*dx
            fkp(2,1) = -rr3*puiy + rr5*puir*dy
            fkp(3,1) = -rr3*puiz + rr5*puir*dz
            fid(1,2) = -rr3*dukx2 + rr5*dukr2*dx
            fid(2,2) = -rr3*duky2 + rr5*dukr2*dy
            fid(3,2) = -rr3*dukz2 + rr5*dukr2*dz
            fkd(1,2) = -rr3*duix2 + rr5*duir2*dx
            fkd(2,2) = -rr3*duiy2 + rr5*duir2*dy
            fkd(3,2) = -rr3*duiz2 + rr5*duir2*dz
            fip(1,2) = -rr3*pukx2 + rr5*pukr2*dx
            fip(2,2) = -rr3*puky2 + rr5*pukr2*dy
            fip(3,2) = -rr3*pukz2 + rr5*pukr2*dz
            fkp(1,2) = -rr3*puix2 + rr5*puir2*dx
            fkp(2,2) = -rr3*puiy2 + rr5*puir2*dy
            fkp(3,2) = -rr3*puiz2 + rr5*puir2*dz
            do j = 1,3
              efi(j,1,ipoleloc) =efi(j,1,ipoleloc) - fimd(j,1) +fid(j,1)
              efi(j,2,ipoleloc) =efi(j,2,ipoleloc) - fimp(j,1) +fip(j,1)
              efi(j,1,kkpoleloc)=efi(j,1,kkpoleloc)- fkmd(j,1) +fkd(j,1)
              efi(j,2,kkpoleloc)=efi(j,2,kkpoleloc)- fkmp(j,1) +fkp(j,1)

              efi2(j,1,ipoleloc)  = efi2(j,1,ipoleloc)  
     $                              - fimd(j,2) + fid(j,2)
              efi2(j,2,ipoleloc)  = efi2(j,2,ipoleloc)  
     $                              - fimp(j,2) + fip(j,2)
              efi2(j,1,kkpoleloc) = efi2(j,1,kkpoleloc) 
     $                              - fkmd(j,2) + fkd(j,2)
              efi2(j,2,kkpoleloc) = efi2(j,2,kkpoleloc) 
     $                              - fkmp(j,2) + fkp(j,2)
            end do
          end if
        end do
c
c     reset interaction scaling coefficients for connected atoms
c
        do j = 1, np11(iglob)
           dscale(ip11(j,iglob)) = 1.0_ti_p
        end do
        do j = 1, np12(iglob)
           dscale(ip12(j,iglob)) = 1.0_ti_p
        end do
        do j = 1, np13(iglob)
           dscale(ip13(j,iglob)) = 1.0_ti_p
        end do
        do j = 1, np14(iglob)
           dscale(ip14(j,iglob)) = 1.0_ti_p
        end do
      end do
      if(dodiag) then
c
c     if dodiag is true, also compute the "self-induced" field,
c     i.e., the diagonal portion of the matrix/vector product.
c
        do i = 1, npoleloc
          iipole = poleglob(i)
          if (polarity(iipole) == 0_ti_p) then
             cycle
          else
            do irhs = 1, nrhs
              do j = 1, 3
                efi(j,irhs,i) = efi(j,irhs,i) +
     $             mu(j,irhs,i)/polarity(iipole)
                efi2(j,irhs,i) = efi2(j,irhs,i) +
     $             mu2(j,irhs,i)/polarity(iipole)
              end do
            end do
          end if
        end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (dscale)
c
      return
      end
