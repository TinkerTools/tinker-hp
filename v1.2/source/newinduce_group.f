 
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c      
c     evaluate induced dipole moments and the polarization energy
c     of a group of atom without PBC
c     using a (preconditioned) conjugate gradient algorithm 
c     or a JI/DIIS solver
c     useful to define the polarization to be retrieved from the total
c     one in "group-scaled" interactions
c      
c     literature reference:
c     "Scalable Evaluation of Polarization Energy and Associated Forces
c     in Polarizable Molecular Dynamics: II. Toward Massively Parallel
c     Computations Using Smooth Particle Mesh Ewald",L. Lagardere et al.,
c     J. Chem. Theory Comput., 2015, 11 (6), pp 2589–2599
c
      subroutine newinduce_group
      use atmlst
      use domdec
      use group
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
      integer i, j, k, nrhs
      integer iglob,iglobgroup
c
c     MPI
c
      integer iipole
      integer, allocatable :: reqsend(:),reqrec(:)
      integer, allocatable :: req2send(:),req2rec(:)
c
      parameter (nrhs=2)

      real*8, allocatable :: ef(:,:,:), mu(:,:,:)
c
      external tmatxb_group
c
      if (.not.use_polar) return
c
 1000 format(' illegal polalg in newinduce.')
c
c
c     allocate some memory and clear the arrays:
c
      allocate (mu(3,nrhs,max(1,npolegroup)))
      mu = 0d0
c
      allocate (ef(3,nrhs,max(1,npolegroup)))
      ef = 0d0
c
      allocate (reqsend(nproc))
      allocate (reqrec(nproc))
      allocate (req2send(nproc))
      allocate (req2rec(nproc))
 
      call orderbuffergroup
c
c     compute the electric fields:
c
      call efld0_group(nrhs,ef)
c
      call commfieldfull(nrhs,ef)
c
c     guess the dipoles:
c
c     predicted values for always stable predictor-corrector method
c
      do i = 1, npolelocgroup
        iglobgroup = ipolegroup(poleglobgroup(i))
        iglob = globglobgroup(iglobgroup)
        iipole = pollist(iglob)
        do k = 1, nrhs
          do j = 1, 3
            mu(j,k,i) = polarity(iipole)*ef(j,k,i)
          end do
        end do
      end do
c
      call commdirdirfull(nrhs,0,mu,reqrec,reqsend)
      call commdirdirfull(nrhs,1,mu,reqrec,reqsend)
      call commdirdirfull(nrhs,2,mu,reqrec,reqsend)
c
      
c
c     now, call the proper solver.
c
      if (polalg.eq.1) then
        call inducepcg_group(tmatxb_group,nrhs,.true.,ef,mu)
      else if (polalg.eq.2) then
        call inducejac_group(tmatxb_group,nrhs,.true.,ef,mu)
      else
         if (rank.eq.0) write(iout,1000) 
         call fatal
      end if

c
c     move the computed dipoles in the module.
c
      do i = 1, npolegroup
        iglobgroup = poleglobgroup(i)
        do j = 1, 3
          uindgroup(j,iglobgroup) = mu(j,1,i)
          uinpgroup(j,iglobgroup) = mu(j,2,i)
        end do
      end do

      deallocate (reqrec)
      deallocate (reqsend)
      deallocate (req2rec)
      deallocate (req2send)
      deallocate (ef)
      deallocate (mu)
      return
      end
c
      subroutine inducepcg_group(matvec,nrhs,precnd,ef,mu)
      use atmlst
      use domdec
      use group
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
      real*8  ef(3,nrhs,*), mu(3,nrhs,*)
      logical precnd
      real*8, allocatable :: res(:,:,:), h(:,:,:), pp(:,:,:), zr(:,:,:),
     $  diag(:)
      integer i, it, j, k
      real*8  ggold(2), ggnew(2), gnorm(2), gg(2), alphacg(2), ene(2)
      real*8  zero, pt5, one, resnrm
      save    zero, pt5, one
      data    zero/0.0d0/, pt5/0.50d0/, one/1.0d0/
      external matvec
c
c     MPI
c
      integer iglob, iipole, ierr
      integer iglobgroup
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
      allocate (res(3,nrhs,max(1,npolelocgroup)))
      allocate (h(3,nrhs,max(1,npolegroup)))
      allocate (pp(3,nrhs,max(1,npolegroup)))
      allocate (zr(3,nrhs,max(1,npolelocgroup)))
      allocate (diag(npolelocgroup))
      if (precnd) then
        do i = 1, npolelocgroup
          iglobgroup = ipolegroup(poleglobgroup(i))
          iglob = globglobgroup(iglobgroup)
          iipole = pollist(iglob)
          diag(i) = polarity(iipole)
        end do
        if (polprt.ge.2.and.rank.eq.0) write (iout,1040)
      else
        do i = 1, npolelocgroup
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
      call matvec(nrhs,.true.,mu,h)
      call commfieldfull(nrhs,h)
c
      res(:,:,1:npolelocgroup) = ef(:,:,1:npolelocgroup)-
     $   h(:,:,1:npolelocgroup)
      do k = 1, nrhs
        do j = 1, 3
          zr(j,k,1:npolelocgroup) = diag(1:npolelocgroup)*
     $      res(j,k,1:npolelocgroup)
        end do
      end do
      pp(:,:,1:npolelocgroup) = zr(:,:,1:npolelocgroup)
      do k = 1, nrhs
        ggold(k) = sum(res(:,k,:)*zr(:,k,:))
      end do
      call MPI_IALLREDUCE(MPI_IN_PLACE,ggold(1),nrhs,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,req1,ierr)
c
      call commdirdirfull(nrhs,0,pp,reqrec,reqsend)
      call commdirdirfull(nrhs,1,pp,reqrec,reqsend)
      call commdirdirfull(nrhs,2,pp,reqrec,reqsend)
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
        call commfieldfull(nrhs,h)
c
        do k = 1, nrhs
          gg(k) = sum(pp(:,k,1:npolelocgroup)*h(:,k,1:npolelocgroup))
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
        do k = 1, nrhs
          mu(:,k,1:npolelocgroup) = mu(:,k,1:npolelocgroup)+ alphacg(k)*
     $     pp(:,k,1:npolelocgroup)
          res(:,k,1:npolelocgroup) =res(:,k,1:npolelocgroup)-alphacg(k)*
     $     h(:,k,1:npolelocgroup)
        end do
        do k = 1, nrhs
          do j = 1, 3
            zr(j,k,1:npolelocgroup) = diag(1:npolelocgroup)
     $        *res(j,k,1:npolelocgroup)
          end do
        end do
        do k = 1, nrhs
          ggnew(k) = sum(res(:,k,:)*zr(:,k,:))
          ene(k) = -pt5*sum(mu(:,k,1:npolelocgroup)*
     $      (res(:,k,1:npolelocgroup)+ ef(:,k,1:npolelocgroup)))
        end do
        call MPI_IALLREDUCE(MPI_IN_PLACE,ggnew(1),nrhs,MPI_REAL8,
     $    MPI_SUM,COMM_TINKER,req3,ierr)
        call MPI_IALLREDUCE(MPI_IN_PLACE,ene(1),nrhs,MPI_REAL8,MPI_SUM,
     $    COMM_TINKER,req4,ierr)
        call MPI_WAIT(req3,status,ierr)
        call MPI_WAIT(req4,status,ierr)
        resnrm = zero
        do k = 1, nrhs
          gnorm(k) = sqrt(ggnew(k)/dble(3*npolegroup))
          resnrm   = max(resnrm,gnorm(k))
        end do
        if (polprt.ge.2.and.rank.eq.0) write(iout,1010)
     $    it, (ene(k)*coulomb, gnorm(k), k = 1, nrhs)
        do k = 1, nrhs
          pp(:,k,1:npolelocgroup) = zr(:,k,1:npolelocgroup)
     $    +ggnew(k)/ggold(k)*pp(:,k,1:npolelocgroup)
        end do
c
        call commdirdirfull(nrhs,0,pp,reqrec,reqsend)
        call commdirdirfull(nrhs,1,pp,reqrec,reqsend)
        call commdirdirfull(nrhs,2,pp,reqrec,reqsend)
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
      call commdirdirfull(nrhs,0,mu,reqendrec,reqendsend)
      call commdirdirfull(nrhs,1,mu,reqendrec,reqendsend)
      call commdirdirfull(nrhs,2,mu,reqendrec,reqendsend)
c
c     finalize and debug printing:
c
      if (polprt.ge.3) then
        write(iout,1020)
        do i = 1, npolelocgroup
          iglobgroup = ipolegroup(poleglobgroup(i))
          iglob = globglobgroup(iglobgroup)
          write(iout,1030) iglob, (mu(j,1,i), j = 1, 3)
        end do
      end if
      if (polprt.ge.4) then
        write(iout,1021)
        do i = 1, npolelocgroup
          iglobgroup = ipolegroup(poleglobgroup(i))
          iglob = globglobgroup(iglobgroup)
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
      subroutine inducejac_group(matvec,nrhs,dodiis,ef,mu)
      use atmlst
      use domdec
      use group
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
      real*8  ef(3,nrhs,*), mu(3,nrhs,*)
      real*8, allocatable :: munew(:,:,:), h(:,:,:)
      real*8, allocatable :: xdiis(:,:), ediis(:,:), bmat(:,:),
     $   bloc(:,:), cex(:)
      integer i, j, k, it, ind, ndismx, nmat, lenb
      logical dodiis
      parameter (ndismx=25)
      integer ipiv(ndismx+1)
      real*8  zero, one, rnorm(2), rr, xx(1)

      save    zero, one, xx
      external matvec
c
c     MPI
c
      integer iglob, iipole, ierr
      integer iglobgroup
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
      zero  = 0.0d0
      one   = 1.0d0
      xx(1) = 0.0d0

      allocate (munew(3,nrhs,max(1,npolegroup)))
      munew = 0d0
      allocate (reqrec(nproc))
      allocate (reqsend(nproc))
      allocate (req2rec(nproc))
      allocate (req2send(nproc))
c
      allocate (h(3,nrhs,max(1,npolegroup)))
      h = 0d0
      if (dodiis) then
        nmat = 1
        lenb = ndismx + 1
        allocate (xdiis(3*nrhs*max(npolelocgroup,1),ndismx))
        allocate (ediis(3*nrhs*max(npolelocgroup,1),ndismx))
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
        call matvec(nrhs,.false.,mu,h)
        call commfieldfull(nrhs,h)
c
c     jacobi step:
c
c
        do i = 1, npolelocgroup
          iglobgroup = ipolegroup(poleglobgroup(i))
          iglob = globglobgroup(iglobgroup)
          iipole = pollist(iglob)
          do k = 1, nrhs
            do j = 1, 3
              munew(j,k,i) = polarity(iipole)*(ef(j,k,i) - h(j,k,i))
            end do
          end do
        end do
c
        do k = 1, nrhs
          rnorm(k) = sum((munew(1:3,k,1:npolelocgroup)-
     $      mu(1:3,k,1:npolelocgroup))**2)
        end do
c
        call MPI_IALLREDUCE(MPI_IN_PLACE,rnorm(1),nrhs,MPI_REAL8,
     $    MPI_SUM,COMM_TINKER,reqnorm,ierr)
        if (dodiis) then
          ind = 0
          do i = 1, npolelocgroup
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
          call diis(ndismx,3*nrhs*npolelocgroup,xdiis,ediis,
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
          call extrap(3*nrhs*npolelocgroup,nmat-1,xdiis,cex,munew)
        end if
c
        call commdirdirfull(nrhs,0,mu,reqrec,reqsend)
        call commdirdirfull(nrhs,1,munew,reqrec,reqsend)
        call commdirdirfull(nrhs,2,mu,reqrec,reqsend)
        mu(:,:,1:npolelocgroup) = munew(:,:,1:npolelocgroup)

c
c     compute the norm of the increment.
c
        call MPI_WAIT(reqnorm,status,ierr)
        rr = zero
        do k = 1, nrhs
         rnorm(k) = sqrt(rnorm(k)/dble(3*npolegroup))
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
        do i = 1, npolelocgroup
          iglobgroup = ipolegroup(poleglobgroup(i))
          iglob = globglobgroup(iglobgroup)
          write(iout,1030) iglob, (mu(j,1,i), j = 1, 3)
        end do
      end if
      if (polprt.ge.4) then
        write(iout,1021)
        do i = 1, npolelocgroup
          iglobgroup = ipolegroup(poleglobgroup(i))
          iglob = globglobgroup(iglobgroup)
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
c
      subroutine efld0_group(nrhs,ef)
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
      use group
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
      integer i,iglob,kglob,nrhs
      integer iglobgroup,kglobgroup
      real*8  ef(3,nrhs,npolegroup)
      integer ii, j, k, kk, kkk, iipole, kkpole
      real*8  dx, dy, dz, d, d2, d3, d5, pdi, pti, ck,
     $  dkx, dky, dkz, dkr, qkxx, qkyy, qkzz, qkxy, qkxz, qkyz,
     $  qkx, qky, qkz, qkr
      real*8 dix, diy, diz, dir, qixx, qiyy, qizz, qixy, qixz, qiyz,
     $  qix, qiy, qiz, qir, ci
      real*8 damp, pgamma, expdamp, scale3, scale5, scale7
      real*8 drr3,drr5,drr7
      real*8 prr3,prr5,prr7
      real*8 dsc3,dsc5,dsc7
      real*8 psc3,psc5,psc7
      real*8 zero, pt6, one, f50
      real*8 fid(3), fip(3)
      real*8 fkd(3), fkp(3)
      real*8 scaled,scalep
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: pscale(:)
      save   zero, pt6, one, f50
      data   zero/0.d0/, pt6/0.6d0/, one/1.d0/,f50/50.d0/
c
 1000 format(' Warning, system moved too much since last neighbor list
     $  update, try lowering nlupdate')
c
c
      allocate (dscale(n))
      allocate (pscale(n))
      do i = 1, n
         dscale(i) = 1.0d0
         pscale(i) = 1.0d0
      end do
c
      do ii = 1, npolelocgroup
        iglobgroup = ipolegroup(poleglobgroup(ii))
        iglob = globglobgroup(iglobgroup)
        iipole = pollist(iglob)
        pdi = pdamp(iipole)
        pti = thole(iipole)
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
        do kkk = iglobgroup+1, npolegroup
          kglobgroup = ipolegroup(kkk)
          kglob = globglobgroup(kglobgroup)
          kkpole = pollist(kglob)
          kk = polelocgroup(kkk)
          dx = x(kglob) - x(iglob)
          dy = y(kglob) - y(iglob)
          dz = z(kglob) - z(iglob)
          d2 = dx*dx + dy*dy + dz*dz
          d  = sqrt(d2)
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
     &                    *(one-damp + pt6*damp**2)
            end if
          end if
          scaled = dscale(kglob)
          scalep = pscale(kglob)
          dsc3 = scale3 * scaled
          dsc5 = scale5 * scaled
          dsc7 = scale7 * scaled
          psc3 = scale3 * scalep
          psc5 = scale5 * scalep
          psc7 = scale7 * scalep
          drr3 = dsc3 / (d*d2)
          drr5 = 3.0d0 * dsc5 / (d*d2*d2)
          drr7 = 15.0d0 * dsc7 / (d*d2*d2*d2)
          prr3 = psc3 / (d*d2)
          prr5 = 3.0d0 * psc5 / (d*d2*d2)
          prr7 = 15.0d0 * psc7 / (d*d2*d2*d2)
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
          fid(1) = -dx*(drr3*ck-drr5*dkr+drr7*qkr)
     &                 - drr3*dkx + 2.0d0*drr5*qkx
          fid(2) = -dy*(drr3*ck-drr5*dkr+drr7*qkr)
     &                - drr3*dky + 2.0d0*drr5*qky
          fid(3) = -dz*(drr3*ck-drr5*dkr+drr7*qkr)
     &                 - drr3*dkz + 2.0d0*drr5*qkz
          fkd(1) = dx*(drr3*ci+drr5*dir+drr7*qir)
     &                - drr3*dix - 2.0d0*drr5*qix
          fkd(2) = dy*(drr3*ci+drr5*dir+drr7*qir)
     &                - drr3*diy - 2.0d0*drr5*qiy
          fkd(3) = dz*(drr3*ci+drr5*dir+drr7*qir)
     &                   - drr3*diz - 2.0d0*drr5*qiz
          fip(1) = -dx*(prr3*ck-prr5*dkr+prr7*qkr)
     &                 - prr3*dkx + 2.0d0*prr5*qkx
          fip(2) = -dy*(prr3*ck-prr5*dkr+prr7*qkr)
     &                 - prr3*dky + 2.0d0*prr5*qky
          fip(3) = -dz*(prr3*ck-prr5*dkr+prr7*qkr)
     &                 - prr3*dkz + 2.0d0*prr5*qkz
          fkp(1) = dx*(prr3*ci+prr5*dir+prr7*qir)
     &                - prr3*dix - 2.0d0*prr5*qix
          fkp(2) = dy*(prr3*ci+prr5*dir+prr7*qir)
     &                - prr3*diy - 2.0d0*prr5*qiy
          fkp(3) = dz*(prr3*ci+prr5*dir+prr7*qir)
     &                - prr3*diz - 2.0d0*prr5*qiz

          ef(1,1,ii) = ef(1,1,ii) + fid(1)
          ef(2,1,ii) = ef(2,1,ii) + fid(2)
          ef(3,1,ii) = ef(3,1,ii) + fid(3)
          ef(1,2,ii) = ef(1,2,ii) + fip(1)
          ef(2,2,ii) = ef(2,2,ii) + fip(2)
          ef(3,2,ii) = ef(3,2,ii) + fip(3)
c
          ef(1,1,kk) = ef(1,1,kk) + fkd(1)
          ef(2,1,kk) = ef(2,1,kk) + fkd(2)
          ef(3,1,kk) = ef(3,1,kk) + fkd(3)
          ef(1,2,kk) = ef(1,2,kk) + fkp(1)
          ef(2,2,kk) = ef(2,2,kk) + fkp(2)
          ef(3,2,kk) = ef(3,2,kk) + fkp(3)

        end do
c
c     reset interaction scaling coefficients for connected atoms
c
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
      end do
c
      deallocate (dscale)
      deallocate (pscale)

      return
      end
c
      subroutine tmatxb_group(nrhs,dodiag,mu,efi)
c
c     Compute the direct space contribution to the electric field due to the current value
c     of the induced dipoles
c
      use sizes
      use atmlst
      use atoms
      use domdec
      use group
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
      integer i,nrhs,iglob,kglob,iipole,kkpole
      integer iglobgroup,kglobgroup
      real*8  mu(3,nrhs,npolegroup), efi(3,nrhs,npolegroup)
      logical dodiag
      integer j, ii, kk, kkk, irhs
      real*8  dx, dy, dz, d, d2, damp, expdamp, pgamma
      real*8  scale3, scale5, pdi, pti, pol
      real*8  rr3, rr5
      real*8  dukx, duky, dukz, pukx, puky, pukz
      real*8  puir, pukr,duir, dukr
      real*8 duix, duiy, duiz, puix, puiy, puiz
      real*8 fid(3), fip(3)
      real*8 fkd(3), fkp(3)
      real*8  zero, one, f50
      real*8, allocatable :: dscale(:)
      real*8  scale
      save    zero, one, f50

      zero = 0.0d0
      one  = 1.0d0
      f50  = 50.d0

c
c     initialize the result vector
c
      do i = 1, npolegroup
        do irhs = 1, nrhs
          do j = 1, 3
            efi(j,irhs,i) = zero
          end do
        end do
      end do
      allocate (dscale(n))
      do i = 1, n
         dscale(i) = 1.0d0
      end do
c
c     gather some parameters, then set up the damping factors.
c
c
      do ii = 1, npolelocgroup
        iglobgroup = ipolegroup(poleglobgroup(ii))
        iglob = globglobgroup(iglobgroup)
        iipole = pollist(iglob)
        pdi = pdamp(iipole)
        pti = thole(iipole)
        duix = mu(1,1,ii)
        duiy = mu(2,1,ii)
        duiz = mu(3,1,ii)
        puix = mu(1,2,ii)
        puiy = mu(2,2,ii)
        puiz = mu(3,2,ii)
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
        do kkk = iglobgroup+1, npolegroup
          kglobgroup = ipolegroup(kkk)
          kglob = globglobgroup(kglobgroup)
          kkpole = pollist(kglob)
          kk = polelocgroup(kkk)
          dx = x(kglob) - x(iglob)
          dy = y(kglob) - y(iglob)
          dz = z(kglob) - z(iglob)
          d2 = dx*dx + dy*dy + dz*dz
c
c     compute the distances and the scaling factors according to
c     Thole's model.
c
          d  = sqrt(d2)
          dukx = mu(1,1,kk)
          duky = mu(2,1,kk)
          dukz = mu(3,1,kk)
          pukx = mu(1,2,kk)
          puky = mu(2,2,kk)
          pukz = mu(3,2,kk)
c
          scale = dscale(kglob)
c
          scale3 = scale
          scale5 = scale
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
          rr3 = -scale3 / (d*d2)
          rr5 = -3.0d0 * scale5 / (d*d2*d2)
          duir = dx*duix + dy*duiy + dz*duiz
          dukr = dx*dukx + dy*duky + dz*dukz
          puir = dx*puix + dy*puiy + dz*puiz
          pukr = dx*pukx + dy*puky + dz*pukz

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

          efi(1,1,ii)  = efi(1,1,ii)  + fid(1)
          efi(1,2,ii)  = efi(1,2,ii)  + fip(1)
          efi(2,1,ii)  = efi(2,1,ii)  + fid(2)
          efi(2,2,ii)  = efi(2,2,ii)  + fip(2)
          efi(3,1,ii)  = efi(3,1,ii)  + fid(3)
          efi(3,2,ii)  = efi(3,2,ii)  + fip(3)
c
          efi(1,1,kk) = efi(1,1,kk)  + fkd(1)
          efi(1,2,kk) = efi(1,2,kk)  + fkp(1)
          efi(2,1,kk) = efi(2,1,kk)  + fkd(2)
          efi(2,2,kk) = efi(2,2,kk)  + fkp(2)
          efi(3,1,kk) = efi(3,1,kk)  + fkd(3)
          efi(3,2,kk) = efi(3,2,kk)  + fkp(3)
        end do
c
c     reset interaction scaling coefficients for connected atoms
c
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
      end do
      if(dodiag) then
c
c     if dodiag is true, also compute the "self-induced" field,
c     i.e., the diagonal portion of the matrix/vector product.
c
        do i = 1, npolelocgroup
          iglobgroup = ipolegroup(poleglobgroup(i))
          iglob = globglobgroup(iglobgroup)
          iipole = pollist(iglob)
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
      deallocate (dscale)
      return
      end

c
c    subroutine orderbuffergroup : get sizes and indexes of all the domains
c
      
      subroutine orderbuffergroup
      use atoms
      use domdec
      use group
      use mpole
      use potent
      use mpi
      implicit none
      integer i,iproc,tag,status(MPI_STATUS_SIZE)
      integer ierr,commloc
      integer count
      integer, allocatable :: reqrec(:), reqsend(:)
c
      if (use_pmecore) then
        commloc = comm_dir
      else
        commloc = COMM_TINKER
      end if
      allocate (reqrec(nproc))
      allocate (reqsend(nproc))
c
c     deal with atoms first
c
      domlengroup(rank+1) = nlocatgroup

      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*rank + iproc + 1
        call MPI_IRECV(domlengroup(iproc+1),1,MPI_INT,iproc,tag,
     $   COMM_TINKER,reqrec(iproc+1),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*iproc + rank + 1
        call MPI_ISEND(domlengroup(rank+1),1,MPI_INT,iproc,tag,
     $   COMM_TINKER,reqsend(iproc+1),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
          call MPI_WAIT(reqsend(iproc+1),status,ierr)
          call MPI_WAIT(reqrec(iproc+1),status,ierr)
        end if
      end do
c
      bufbeggroup(rank+1) = 1
      count = domlengroup(rank+1)
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
          if (domlengroup(iproc+1).eq.0) then
            bufbeggroup(iproc+1) = 1
          else
            bufbeggroup(iproc+1) = count + 1
          end if
          count = count + domlengroup(iproc+1)
        end if
      end do
c
c     get the indexes of the atoms in the group frame
c
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
          tag = nproc*rank + iproc + 1
          call MPI_IRECV(globgroup(bufbeggroup(iproc+1)),
     $     domlengroup(iproc+1),MPI_INT,iproc,tag,
     $     COMM_TINKER,reqrec(iproc+1),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
          tag = nproc*iproc + rank + 1
          call MPI_ISEND(globgroup,domlengroup(rank+1),MPI_INT,
     $     iproc,tag,COMM_TINKER,reqsend(iproc+1),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
          call MPI_WAIT(reqrec(iproc+1),status,ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
          call MPI_WAIT(reqsend(iproc+1),status,ierr)
        end if
      end do

      do iproc = 1, nproc
        if (iproc.ne.(rank+1)) then
          do i = 1, domlengroup(iproc)
            locgroup(globgroup(bufbeggroup(iproc)+i-1)) =
     $        bufbeggroup(iproc)+i-1
          end do
        end if
      end do
c
c     then deal with multipoles
c
      domlenpolegroup(rank+1) = npolelocgroup
c
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*rank + iproc + 1
        call MPI_IRECV(domlenpolegroup(iproc+1),1,MPI_INT,iproc,tag,
     $   COMM_TINKER,reqrec(iproc+1),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*iproc + rank + 1
        call MPI_ISEND(domlenpolegroup(rank+1),1,MPI_INT,iproc,tag,
     $   COMM_TINKER,reqsend(iproc+1),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
          call MPI_WAIT(reqsend(iproc+1),status,ierr)
          call MPI_WAIT(reqrec(iproc+1),status,ierr)
        end if
      end do
c
      bufbegpolegroup(rank+1) = 1
      count = domlenpolegroup(rank+1)
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
          if (domlenpolegroup(iproc+1).eq.0) then
            bufbegpolegroup(iproc+1) = 1
          else
            bufbegpolegroup(iproc+1) = count + 1
          end if
          count = count + domlenpolegroup(iproc+1)
        end if
      end do
c
c     get the indexes of the mulipoles in the group frame
c
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
          tag = nproc*rank + iproc + 1
          call MPI_IRECV(poleglobgroup(bufbegpolegroup(iproc+1)),
     $     domlenpolegroup(iproc+1),MPI_INT,iproc,tag,
     $     COMM_TINKER,reqrec(iproc+1),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
          tag = nproc*iproc + rank + 1
          call MPI_ISEND(poleglobgroup,domlenpolegroup(rank+1),MPI_INT,
     $     iproc,tag,COMM_TINKER,reqsend(iproc+1),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
          call MPI_WAIT(reqrec(iproc+1),status,ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
          call MPI_WAIT(reqsend(iproc+1),status,ierr)
        end if
      end do

      do iproc = 1, nproc
        if (iproc.ne.(rank+1)) then
          do i = 1, domlenpolegroup(iproc)
            polelocgroup(poleglobgroup(bufbegpolegroup(iproc)+i-1)) =
     $        bufbegpolegroup(iproc)+i-1
          end do
        end if
      end do
c
      deallocate (reqrec)
      deallocate (reqsend)
      return
      end
c
c     subroutine commfieldfull : communicate some direct fields (Newton's third law)
c
      subroutine commfieldfull(nrhs,ef)
      use atoms
      use domdec
      use group
      use mpole
      use potent
      use mpi
      implicit none
      integer nrhs,i,j,k,tag,status(MPI_STATUS_SIZE)
      integer ierr,commloc
      integer, allocatable :: reqrec(:), reqsend(:)
      real*8 ef(3,nrhs,*)
      real*8, allocatable :: buffer(:,:,:,:)
c
      if (use_pmecore) then
        commloc = comm_dir
      else
        commloc = COMM_TINKER
      end if
c
      allocate (buffer(3,nrhs,max(npolelocgroup,1),nproc))
      allocate (reqrec(nproc))
      allocate (reqsend(nproc))
c
      do i = 1, nproc
        if (i.ne.(rank+1)) then
          tag = nproc*rank + i
          call MPI_IRECV(buffer(1,1,1,i),3*nrhs*npolelocgroup,MPI_REAL8,
     $     i-1,tag,commloc,reqrec(i),ierr)
        end if
      end do
c
      do i = 1, nproc
        if (i.ne.(rank+1)) then
          tag = nproc*(i-1) + rank + 1
          call MPI_ISEND(ef(1,1,bufbegpolegroup(i)),
     $     3*nrhs*domlenpolegroup(i),MPI_REAL8,i-1,
     $    tag,commloc,reqsend(i),ierr)
        end if
      end do
c
      do i = 1, nproc
        if (i.ne.(rank+1)) then
          tag = nproc*rank + i
          call MPI_WAIT(reqrec(i),status,ierr)
        end if
      end do
      do i = 1, nproc
        if (i.ne.(rank+1)) then
          tag = nproc*(i-1) + rank + 1
          call MPI_WAIT(reqsend(i),status,ierr)
        end if
      end do
c
      do i = 1, nproc
        if (i.ne.(rank+1)) then
          do j = 1, npolelocgroup
            do k = 1, nrhs
              ef(1,k,j) = ef(1,k,j) + buffer(1,k,j,i)
              ef(2,k,j) = ef(2,k,j) + buffer(2,k,j,i)
              ef(3,k,j) = ef(3,k,j) + buffer(3,k,j,i)
            end do
          end do
        end if
      end do
c
      deallocate (buffer)
      deallocate (reqrec)
      deallocate (reqsend)
      return
      end
c
c    subroutine commdirdirfull : deal with communications of direct dipoles
c
c    rule determines what to do:
c        - 0: start reception of direct dipoles
c        - 1: start sending of direct dipoles 
c        - 2: wait for the communications to be done
c
      subroutine commdirdirfull(nrhs,rule,mu,reqrec,reqsend)
      use domdec
      use group
      use iounit
      use mpole
      use mpi
      implicit none
      integer nrhs,rule,ierr,status(MPI_STATUS_SIZE),tag,i
      integer :: reqrec(nproc),reqsend(nproc)
      real*8 mu(3,nrhs,npolegroup)
 1000 format(' illegal rule in commdirdirfull.')

      if (rule.eq.0) then
c
c     MPI : begin reception
c
        do i = 1, nproc
          if (i.ne.(rank+1)) then
            tag = nproc*rank + i
            call MPI_IRECV(mu(1,1,bufbegpolegroup(i)),
     $       3*nrhs*domlenpolegroup(i),MPI_REAL8,i-1,tag,
     $       COMM_TINKER,reqrec(i),ierr)
          end if
        end do

      else if (rule.eq.1) then
c
        do i = 1, nproc
          if (i.ne.(rank+1)) then
            tag = nproc*(i-1) + rank + 1
            call MPI_ISEND(mu,3*nrhs*npolelocgroup,
     $         MPI_REAL8,i-1,tag,COMM_TINKER,
     $         reqsend(i),ierr)
          end if
        end do
      else if (rule.eq.2) then
        do i = 1, nproc
          if (i.ne.(rank+1)) then
            call MPI_WAIT(reqrec(i),status,ierr)
          end if
        end do
        do i = 1, nproc
          if (i.ne.(rank+1)) then
            call MPI_WAIT(reqsend(i),status,ierr)
          end if
        end do
      else
        if (rank.eq.0) write(iout,1000) 
        call fatal 
      end if

      return
      end
c
c     subroutine commforcesgroup : deal with communications of "polarization group" forces
c
      subroutine commforcesgroup
      use sizes
      use atoms
      use deriv
      use domdec
      use group
      use potent
      use mpi
      implicit none
      integer i,tag,ierr
      integer, allocatable :: reqsend(:),reqrec(:)
      real*8, allocatable :: buffer(:,:,:)
      integer status(MPI_STATUS_SIZE)
c
c     allocate some arrays
c
      allocate (reqsend(nproc))
      allocate (reqrec(nproc))
      allocate (buffer(3,max(1,nlocatgroup),nproc))
c      buffer = 0d0
c      buffers = 0d0
c
c     MPI : begin reception in buffer
c
      do i = 1, nproc
        if (i.ne.(rank+1)) then
          tag = nproc*rank + i
          call MPI_IRECV(buffer(1,1,i),3*nlocatgroup,
     $     MPI_REAL8,i-1,tag,COMM_TINKER,reqrec(i),ierr)
        end if
      end do
c
c     MPI : begin sending
c
      do i = 1, nproc
        if (i.ne.(rank+1)) then
          tag = nproc*(i-1) + rank + 1
          call MPI_ISEND(depgroup(1,bufbeggroup(i)),
     $     3*domlengroup(i),
     $     MPI_REAL8,i-1,tag,COMM_TINKER,
     $     reqsend(i),ierr)
        end if
      end do
c
      do i = 1, nproc
        if (i.ne.(rank+1)) then
          call MPI_WAIT(reqrec(i),status,ierr)
        end if
      end do
      do i = 1, nproc
        if (i.ne.(rank+1)) then
          call MPI_WAIT(reqsend(i),status,ierr)
        end if
      end do
c
c     MPI : move in global arrays
c
      do i = 1, nproc
        if (i.ne.(rank+1)) then
          depgroup(:,1:nlocatgroup) = depgroup(:,1:nlocatgroup) + 
     $     buffer(:,1:nlocatgroup,i)
        end if
      end do
      deallocate (reqsend)
      deallocate (reqrec)
      deallocate (buffer)
      return
      end
