

!===================================================
!     sub diagvec2
!===================================================
! Performs product of vector a with polarisabilities
! BUt now on two vectors at once !

      !TODO 1.2 merge this file with the code precision and mpi comm
      subroutine diagvec2(nrhs, A,  c, b, d)
      use atmlst
      use mpole
      use polar
      implicit none

      integer, intent(in) :: nrhs
      real*8, dimension(3,nrhs,npolebloc) :: A, c
      real*8, dimension(3,nrhs,npolebloc) :: B, d
      integer :: i,iipole, irhs, j

      do i = 1, npolebloc
         iipole = poleglob(i)
         do irhs = 1, nrhs
            do j = 1,3
               B(j,irhs,i) = A(j,irhs,i)*polarity(iipole)
               D(j,irhs,i) = C(j,irhs,i)*polarity(iipole)
            end do
         end do
      end do

      return
      end

!===================================================
!     sub diagvec3
!===================================================
! Performs product of vector a with polarisabilities
! Why not three vecs
! Why not three vecs
! Why not three vecs

      subroutine diagvec3(nrhs, A1, a2, a3, b1, b2, b3)
      use atmlst
      use mpole
      use polar
      implicit none

      integer, intent(in) :: nrhs
      real*8, dimension(3,nrhs,npolebloc) :: A1, a2, a3
      real*8, dimension(3,nrhs,npolebloc) :: b1, b2, b3
      integer :: i,iipole, irhs, j

      do i = 1, npolebloc
         iipole = poleglob(i)
         do irhs = 1, nrhs
            do j = 1,3
               B1(j,irhs,i) = A1(j,irhs,i)*polarity(iipole)
               B2(j,irhs,i) = A2(j,irhs,i)*polarity(iipole)
               B3(j,irhs,i) = A3(j,irhs,i)*polarity(iipole)
            end do
         end do
      end do

      return
      end

      subroutine tmatxbrecipsave(mu,murec,nrhs,dipfield,dipfieldbis,
     $ fphivec)
c
c     Compute the reciprocal space contribution to the electric field due to the current  
c     value of the induced dipoles. Saves the potential and subsequent
!     derivatives in a fphi(20,n) vector.
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
      integer ierr,iglob,iloc
      integer status(MPI_STATUS_SIZE),tag
      integer nrhs,iipole

      integer i,j,k


      real*8 fuind(3),fuinp(3)
      real*8 term
      real*8 a(3,3)
      real*8 fdip_phi1(20), fdip_phi2(20)
      real*8 dipfield(3,nrhs,*),dipfieldbis(3,nrhs,*)
      real*8 mu(3,nrhs,*),murec(3,nrhs,*)
      integer, allocatable :: reqbcastrec(:),reqbcastsend(:)
      integer, allocatable :: reqrec(:),reqsend(:)
      real*8, allocatable :: qgridmpi(:,:,:,:,:)
      integer nprocloc,commloc,rankloc,proc
      real*8  :: fphivec(20, 2, max(1,npolerecloc))
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

      !(F.A.)

      allocate (qgridmpi(2,n1mpimax,n2mpimax,n3mpimax,nrec_recep))
      allocate (reqbcastrec(nprocloc))
      allocate (reqbcastsend(nprocloc))
      allocate (reqrec(nprocloc))
      allocate (reqsend(nprocloc))
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
c
c     zero out the PME charge grid
c
      qgrid2in_2d = 0d0
c     fill the pme grid, loop over the induced dipoles sites
      do j = 1, 3
        a(1,j) = dble(nfft1) * recip(j,1)
        a(2,j) = dble(nfft2) * recip(j,2)
        a(3,j) = dble(nfft3) * recip(j,3)
      end do
c
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
c     MPI : Begin reception
c
      do i = 1, nrec_recep
        tag = nprocloc*rankloc + prec_recep(i) + 1
        call MPI_IRECV(qgridmpi(1,1,1,1,i),2*n1mpimax*n2mpimax*n3mpimax,
     $   MPI_REAL8,prec_recep(i),tag,commloc,reqrec(i),ierr)
      end do
c
c     MPI : begin sending
c
      do i = 1, nrec_send
        proc = prec_send(i)
        tag = nprocloc*prec_send(i) + rankloc + 1
        call MPI_ISEND(qgrid2in_2d(1,1,1,1,i+1),
     $   2*n1mpimax*n2mpimax*n3mpimax,MPI_REAL8,proc,tag,commloc,
     $   reqsend(i),ierr)
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
        proc = prec_send(i)
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
     $   2*n1mpimax*n2mpimax*n3mpimax,
     $   MPI_REAL8,prec_recep(i),tag,commloc,reqbcastsend(i),ierr)
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
        call fphi_uind_big(iglob,i,fdip_phi1,fdip_phi2)

        !(F.A.)
        fphivec(:,1,i) = fdip_phi1
        fphivec(:,2,i) = fdip_phi2

        if (repart(iglob).ne.rank) then
c
c     convert the dipole fields from fractional to Cartesian
c
          do k = 1, 3
             dipfieldbis(k,1,i) = a(k,1)*fdip_phi1(2)
     &                           + a(k,2)*fdip_phi1(3)
     &                           + a(k,3)*fdip_phi1(4)
             dipfieldbis(k,2,i) = a(k,1)*fdip_phi2(2)
     &                           + a(k,2)*fdip_phi2(3)
     $                           + a(k,3)*fdip_phi2(4)
          end do
        else
          do k = 1, 3
             dipfield(k,1,iloc) = a(k,1)*fdip_phi1(2)
     &                           + a(k,2)*fdip_phi1(3)
     &                           + a(k,3)*fdip_phi1(4)
             dipfield(k,2,iloc) = a(k,1)*fdip_phi2(2)
     &                           + a(k,2)*fdip_phi2(3)
     $                           + a(k,3)*fdip_phi2(4)
          end do
        end if
      end do
!
      deallocate (qgridmpi)
      deallocate (reqrec)
      deallocate (reqsend)
      deallocate (reqbcastsend)
      deallocate (reqbcastrec)
      return
      end


      subroutine dbletmatxb_pme(nrhs,dodiag,mu, mu2, efi, efi2)
c
c     Compute the direct space contribution to the electric field due to the current value
c     of the induced dipoles, for two input vectors mu and mu2
c
      use atoms
      use atmlst
      use domdec
      use ewald
      use math
      use mpole
      use neigh
      use polar
      use polgrp
      use polpot
      use shunt
      use mpi
      implicit none
      integer i,iipole,nrhs,iglob,kglob,kkpole,kkpoleloc
      integer ipoleloc
      real*8  mu(3,nrhs,*), efi(3,nrhs,npolebloc)
      real*8  mu2(3,nrhs,*), efi2(3,nrhs,npolebloc)
      logical dodiag
      integer j, ii, kkk, irhs, inl
      real*8  dx, dy, dz, d, d2, damp, expdamp, pgamma,
     $  scale3, scale5, pdi, pti 
      real*8 ralpha, alsq2, alsq2n, bfac, exp2a
      real*8 dukx, duky, dukz, pukx, puky, pukz, puir, pukr,
     $  duir, dukr

      real*8 rr3, rr5, dukx2, duky2, dukz2, pukx2, puky2, pukz2, puir2,
     $   pukr2, duir2, dukr2
      real*8 duix, duiy, duiz, puix, puiy, puiz
      real*8 duix2, duiy2, duiz2, puix2, puiy2, puiz2
      real*8 bn(0:3), fid(3,2), fip(3,2), fimd(3,2), fimp(3,2)
      real*8 fkd(3,2), fkp(3,2), fkmd(3,2), fkmp(3,2)
      real*8  zero, one, f50
      real*8, allocatable :: dscale(:)
      real*8  cutoff2
      save    zero, one, f50
      data    zero/0.d0/, one/1.d0/, f50/50.d0/
      character*10 mode
c
c     initialize the result vector
c
      allocate (dscale(n))
      do i = 1, n
         dscale(i) = 1.0d0
      end do
c
      do i = 1, npolebloc
        do irhs = 1, nrhs
          do j = 1, 3
            efi(j,irhs,i) = zero
          end do
        end do
      end do
      efi2 = 0d0
c
c     gather some parameters, then set up the damping factors.
c
      mode = 'EWALD'
      call switch (mode)
      cutoff2 = cut2
c
      do ii = 1, npolelocnl
        iipole = poleglobnl(ii)
        iglob = ipole(iipole)
        i = loc(iglob)
        ipoleloc = poleloc(iipole)
        inl = polelocnl(iipole)
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
        do kkk = 1,nelst(inl)
          kkpole = elst(kkk,inl)
          kkpoleloc = poleloc(kkpole)
          kglob = ipole(kkpole)
          if (kkpoleloc.eq.0) cycle
          dx = x(kglob) - x(iglob)
          dy = y(kglob) - y(iglob)
          dz = z(kglob) - z(iglob)
          call image(dx,dy,dz)
          d2 = dx*dx + dy*dy + dz*dz
          if (d2.le.cutoff2) then
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
            alsq2 = 2.0d0 * aewald**2
            alsq2n = 0.0d0
            if (aewald .gt. 0.0d0)
     &        alsq2n = 1.0d0 / (sqrtpi*aewald)
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
            rr3 = (1.0d0-scale3) / (d*d2)
            rr5 = 3.0d0 * (1.0d0-scale5) / (d*d2*d2)
            duir = dx*duix + dy*duiy + dz*duiz
            dukr = dx*dukx + dy*duky + dz*dukz
            puir = dx*puix + dy*puiy + dz*puiz
            pukr = dx*pukx + dy*puky + dz*pukz
            duir2 = dx*duix2 + dy*duiy2 + dz*duiz2
            dukr2 = dx*dukx2 + dy*duky2 + dz*dukz2
            puir2 = dx*puix2 + dy*puiy2 + dz*puiz2
            pukr2 = dx*pukx2 + dy*puky2 + dz*pukz2


!            fimd = 0d0
!            fimp = 0d0
!            fkmd = 0d0
!            fkmp = 0d0
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
!            fid = 0d0
!            fkd = 0d0
!            fip = 0d0
!            fkp = 0d0
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
              efi(j,1,ipoleloc) = efi(j,1,ipoleloc) -fimd(j,1)+ fid(j,1)
              efi(j,2,ipoleloc) = efi(j,2,ipoleloc) -fimp(j,1)+ fip(j,1)
              efi(j,1,kkpoleloc)= efi(j,1,kkpoleloc)-fkmd(j,1)+ fkd(j,1)
              efi(j,2,kkpoleloc)= efi(j,2,kkpoleloc)-fkmp(j,1)+ fkp(j,1)
              efi2(j,1,ipoleloc) = efi2(j,1,ipoleloc) 
     $                                - fimd(j,2) +fid(j,2)
              efi2(j,2,ipoleloc) = efi2(j,2,ipoleloc) 
     $                                - fimp(j,2) +fip(j,2)
              efi2(j,1,kkpoleloc)= efi2(j,1,kkpoleloc)
     $                                - fkmd(j,2) +fkd(j,2)
              efi2(j,2,kkpoleloc)= efi2(j,2,kkpoleloc)
     $                                - fkmp(j,2) +fkp(j,2)
            end do
          end if
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
        do i = 1, npoleloc
          iipole = poleglob(i)
c
c     if no polarisability, take a negligeable value to allow convergence
c
          if (polarity(iipole) == 0d0) then
             do irhs = 1, nrhs
               do j = 1, 3
                 efi(j,irhs,i) = efi(j,irhs,i) +
     $              mu(j,irhs,i)*100000
                 efi2(j,irhs,i) = efi2(j,irhs,i) +
     $              mu2(j,irhs,i)*100000
               end do
             end do
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

!===================================================
!     sub fftthatplz2(vec,vecbis,fphivec)
! Quite self-explanatory, really. Puts vec on the grid, 
! returns the potential and successive derivatives in 
! fphivec.
!===================================================
      subroutine fftthatplz2(vec,vecbis,fphivec)
      use atmlst
      use bound
      use boxes
      use chgpot
      use domdec
      use ewald
      use fft
      use math
      use mpole
      use pme
      use potent
      use mpi
      implicit none

      real*8, intent(in), dimension(3, 2,npoleloc) :: vec
      real*8, intent(in), dimension(3, 2,npolerecloc) :: vecbis
      real*8, intent(out), dimension(20,2, npolerecloc) :: fphivec

      integer :: i, iipole, iglob,
     $ k, j , rankloc, ierr, tag, iloc
      integer :: status(MPI_STATUS_SIZE), commloc, nprocloc
      real*8 :: term

      real*8, allocatable, dimension(:,:) :: vec_frac1, vec_frac2
      real*8, allocatable, dimension(:,:,:,:) :: qgridout2 
      real*8, allocatable, dimension(:,:,:,:,:) :: qgridin2
      real*8, allocatable :: qgridmpi(:,:,:,:,:)
      integer, allocatable :: reqsend(:),reqrec(:)
      integer, allocatable :: req2send(:),req2rec(:)
      
      if (use_pmecore) then
        nprocloc = nrec
        commloc  = comm_rec
        rankloc  = rank_bis
      else
        nprocloc = nproc
        commloc  = COMM_TINKER
        rankloc  = rank
      end if

      allocate(vec_frac1(3,npolerecloc), vec_frac2(3,npolerecloc))
      allocate(qgridin2(2,n1mpimax,n2mpimax,n3mpimax,nrec_send+1))
!               qgridin3(2,n1mpimax,n2mpimax,n3mpimax,nrec_send+1))
      allocate(qgridout2(2,isize2(rankloc+1), jsize2(rankloc+1),
     $          ksize2(rankloc+1)))
      allocate (qgridmpi(2,n1mpimax,n2mpimax,n3mpimax,nrec_recep))
!               qgridout3(2,isize2(rankloc+1), jsize2(rankloc+1),
!     $          ksize2(rankloc+1)))
      allocate (reqsend(nproc))
      allocate (reqrec(nproc))
      allocate (req2send(nproc))
      allocate (req2rec(nproc))


      qgridin2 = 0d0
!      qgridin3 = 0d0
      fphivec = 0d0

      do i = 1 , npolerecloc
         iipole = polerecglob(i)
         iglob = ipole(iipole)
         iloc  = poleloc(iipole)

         if (repart(iglob).ne.rank) then
           call cart_to_frac_vec(vecbis(:,1,i), vec_frac1(:,i))
           call cart_to_frac_vec(vecbis(:,2,i), vec_frac2(:,i))
         else
           call cart_to_frac_vec(vec(:,1,iloc), vec_frac1(:,i))
           call cart_to_frac_vec(vec(:,2,iloc), vec_frac2(:,i))
         end if 
         call grid_uind_site(iglob,i, vec_frac1(:,i), vec_frac2(:,i), 
     $                         qgridin2)
      end do
c
c     MPI : begin reception
c
      do i = 1, nrec_recep
        tag = nprocloc*rankloc + prec_recep(i) + 1
        call MPI_IRECV(qgridmpi(1,1,1,1,i),2*n1mpimax*n2mpimax*
     $   n3mpimax,MPI_REAL8,prec_recep(i),tag,
     $   commloc,reqrec(i),ierr)
      end do
c
c     MPI : begin sending
c
      do i = 1, nrec_send
        tag = nprocloc*prec_send(i) + rankloc + 1
        call MPI_ISEND(qgridin2(1,1,1,1,i+1),
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
        qgridin2(:,:,:,:,1) = qgridin2(:,:,:,:,1)+
     $   qgridmpi(:,:,:,:,i)
      end do
      call fft2d_frontmpi(qgridin2,qgridout2,n1mpimax,n2mpimax,
     $ n3mpimax)
c      ntot = nfft1*nfft2*nfft3
c      if ((istart2(rankloc+1).eq.1).and.(jstart2(rankloc+1).eq.1).and.
c     $   (kstart2(rankloc+1).eq.1)) then
c           qfac_2d(1,1,1) = 0.0d0
c      end if
c      pterm = (pi/aewald)**2
c      volterm = pi * volbox
c      nf1 = (nfft1+1) / 2
c      nf2 = (nfft2+1) / 2
c      nf3 = (nfft3+1) / 2
c      do k3 = kstart2(rankloc+1),kend2(rankloc+1)
c         do k2 = jstart2(rankloc+1),jend2(rankloc+1)
c            do k1 = istart2(rankloc+1),iend2(rankloc+1)
c               m1 = k1 - 1
c               m2 = k2 - 1
c               m3 = k3 - 1
c               if (k1 .gt. nf1)  m1 = m1 - nfft1
c               if (k2 .gt. nf2)  m2 = m2 - nfft2
c               if (k3 .gt. nf3)  m3 = m3 - nfft3
c               if ((m1.eq.0).and.(m2.eq.0).and.(m3.eq.0)) goto 10
c               r1 = dble(m1)
c               r2 = dble(m2)
c               r3 = dble(m3)
c               h1 = recip(1,1)*r1 + recip(1,2)*r2 + recip(1,3)*r3
c               h2 = recip(2,1)*r1 + recip(2,2)*r2 + recip(2,3)*r3
c               h3 = recip(3,1)*r1 + recip(3,2)*r2 + recip(3,3)*r3
c               hsq = h1*h1 + h2*h2 + h3*h3
c               term = -pterm * hsq
c               expterm = 0.0d0
c               if ((term .gt. -50.0d0)) then
c                   denom = volterm*hsq*bsmod1(k1)*bsmod2(k2)*bsmod3(k3)
c                   expterm = exp(term) / denom
c                   if (.not. use_bounds) then
c                      expterm = expterm * (1.0d0-cos(pi*xbox*sqrt(hsq)))
c                   else if (octahedron) then
c                      if (mod(m1+m2+m3,2) .ne. 0)  expterm = 0.0d0
c                   end if
c                   qfac_2d(k1-istart2(rankloc+1)+1,k2-jstart2(rankloc+1)
c     $                     +1,k3-kstart2(rankloc+1)+1) = expterm
c               end if
c 10            continue
c            end do
c         end do
c      end do
      do k = 1, ksize2(rankloc+1)
         do j = 1, jsize2(rankloc+1)
           do i = 1, isize2(rankloc+1)
             term = qfac_2d(i,j,k)
             qgridout2(1,i,j,k) = term*qgridout2(1,i,j,k)
             qgridout2(2,i,j,k) = term*qgridout2(2,i,j,k)
           end do
         end do
      end do

      call fft2d_backmpi(qgridin2,qgridout2,n1mpimax,n2mpimax,
     $ n3mpimax)
c
c     MPI : Begin reception
c
      do i = 1, nrec_send
        tag = nprocloc*rankloc + prec_send(i) + 1
        call MPI_IRECV(qgridin2(1,1,1,1,i+1),
     $   2*n1mpimax*n2mpimax*n3mpimax,MPI_REAL8,
     $   prec_recep(i),tag,commloc,req2rec(i),ierr)
      end do
c
c     MPI : begin sending
c
      do i = 1, nrec_recep
        tag = nprocloc*prec_recep(i) + rankloc + 1
        call MPI_ISEND(qgridin2(1,1,1,1,1),
     $   2*n1mpimax*n2mpimax*n3mpimax,MPI_REAL8,
     $   prec_send(i),tag,commloc,req2send(i),ierr)
      end do
c
      do i = 1, nrec_send
        call MPI_WAIT(req2rec(i),status,ierr)
      end do
      do i = 1, nrec_recep
        call MPI_WAIT(req2send(i),status,ierr)
      end do

      do i = 1, npolerecloc
        iipole = polerecglob(i)
        iglob = ipole(iipole)
        call fphi_dervec_site2(iglob,i,qgridin2, fphivec(:,:,i))
      end do
      return
      end

!===================================================
!     sub fftthatplz(vec, fphivec)
! Quite self-explanatory, really. Puts vec on the grid, 
! returns the potential and successive derivatives in 
! fphivec.
!===================================================
      subroutine fftthatplz(vec,vecbis,fphivec)
      use atmlst
      use bound
      use boxes
      use chgpot
      use domdec
      use ewald
      use fft
      use math
      use mpole
      use pme
      use potent
      use mpi
      implicit none

      real*8, intent(in), dimension(3,npoleloc) :: vec
      real*8, intent(in), dimension(3,npolerecloc) :: vecbis
      real*8, intent(out), dimension(20, npolerecloc) :: fphivec

      integer :: i, iipole, iglob,
     $k, j , rankloc, ierr, tag, iloc
      integer :: status(MPI_STATUS_SIZE), commloc, nprocloc
      real*8 :: term

      real*8, allocatable, dimension(:,:) :: vec_frac
      real*8, allocatable, dimension(:,:,:,:) :: qgridout2 
      real*8, allocatable, dimension(:,:,:,:,:) :: qgridin2
      real*8, allocatable :: qgridmpi(:,:,:,:,:)
      integer, allocatable :: reqsend(:),reqrec(:)
      integer, allocatable :: req2send(:),req2rec(:)
      
      if (use_pmecore) then
        nprocloc = nrec
        commloc  = comm_rec
        rankloc  = rank_bis
      else
        nprocloc = nproc
        commloc  = COMM_TINKER
        rankloc  = rank
      end if

      allocate(vec_frac(3,npolerecloc))
      allocate(qgridin2(2,n1mpimax,n2mpimax,n3mpimax,nrec_send+1))
!               qgridin3(2,n1mpimax,n2mpimax,n3mpimax,nrec_send+1))
      allocate(qgridout2(2,isize2(rankloc+1), jsize2(rankloc+1),
     $          ksize2(rankloc+1)))
      allocate (qgridmpi(2,n1mpimax,n2mpimax,n3mpimax,nrec_recep))
!               qgridout3(2,isize2(rankloc+1), jsize2(rankloc+1),
!     $          ksize2(rankloc+1)))
      allocate (reqsend(nproc))
      allocate (reqrec(nproc))
      allocate (req2send(nproc))
      allocate (req2rec(nproc))


      qgridin2 = 0d0
!      qgridin3 = 0d0
      fphivec = 0d0

      do i = 1 , npolerecloc
         iipole = polerecglob(i)
         iglob = ipole(iipole)
         iloc  = poleloc(iipole)

         if (repart(iglob).ne.rank) then
           call cart_to_frac_vec(vecbis(:,i), vec_frac(:,i))
         else
           call cart_to_frac_vec(vec(:,iloc), vec_frac(:,i))
         end if
         call grid_uind_site(iglob,i, vec_frac(:,i), vec_frac(:,i), 
     $                       qgridin2)
      end do
c
c     MPI : begin reception
c
      do i = 1, nrec_recep
        tag = nprocloc*rankloc + prec_recep(i) + 1
        call MPI_IRECV(qgridmpi(1,1,1,1,i),2*n1mpimax*n2mpimax*
     $   n3mpimax,MPI_REAL8,prec_recep(i),tag,
     $   commloc,reqrec(i),ierr)
      end do
c
c     MPI : begin sending
c
      do i = 1, nrec_send
        tag = nprocloc*prec_send(i) + rankloc + 1
        call MPI_ISEND(qgridin2(1,1,1,1,i+1),
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
        qgridin2(:,:,:,:,1) = qgridin2(:,:,:,:,1)+
     $   qgridmpi(:,:,:,:,i)
      end do
      call fft2d_frontmpi(qgridin2,qgridout2,n1mpimax,n2mpimax,
     $ n3mpimax)
c      ntot = nfft1*nfft2*nfft3
c      if ((istart2(rankloc+1).eq.1).and.(jstart2(rankloc+1).eq.1).and.
c     $   (kstart2(rankloc+1).eq.1)) then
c           qfac_2d(1,1,1) = 0.0d0
c      end if
c      pterm = (pi/aewald)**2
c      volterm = pi * volbox
c      nf1 = (nfft1+1) / 2
c      nf2 = (nfft2+1) / 2
c      nf3 = (nfft3+1) / 2
c      do k3 = kstart2(rankloc+1),kend2(rankloc+1)
c         do k2 = jstart2(rankloc+1),jend2(rankloc+1)
c            do k1 = istart2(rankloc+1),iend2(rankloc+1)
c               m1 = k1 - 1
c               m2 = k2 - 1
c               m3 = k3 - 1
c               if (k1 .gt. nf1)  m1 = m1 - nfft1
c               if (k2 .gt. nf2)  m2 = m2 - nfft2
c               if (k3 .gt. nf3)  m3 = m3 - nfft3
c               if ((m1.eq.0).and.(m2.eq.0).and.(m3.eq.0)) goto 10
c               r1 = dble(m1)
c               r2 = dble(m2)
c               r3 = dble(m3)
c               h1 = recip(1,1)*r1 + recip(1,2)*r2 + recip(1,3)*r3
c               h2 = recip(2,1)*r1 + recip(2,2)*r2 + recip(2,3)*r3
c               h3 = recip(3,1)*r1 + recip(3,2)*r2 + recip(3,3)*r3
c               hsq = h1*h1 + h2*h2 + h3*h3
c               term = -pterm * hsq
c               expterm = 0.0d0
c               if ((term .gt. -50.0d0)) then
c                   denom = volterm*hsq*bsmod1(k1)*bsmod2(k2)*bsmod3(k3)
c                   expterm = exp(term) / denom
c                   if (.not. use_bounds) then
c                      expterm = expterm * (1.0d0-cos(pi*xbox*sqrt(hsq)))
c                   else if (octahedron) then
c                      if (mod(m1+m2+m3,2) .ne. 0)  expterm = 0.0d0
c                   end if
c                   qfac_2d(k1-istart2(rankloc+1)+1,k2-jstart2(rankloc+1)
c     $                     +1,k3-kstart2(rankloc+1)+1) = expterm
c               end if
c 10            continue
c            end do
c         end do
c      end do
      do k = 1, ksize2(rankloc+1)
         do j = 1, jsize2(rankloc+1)
           do i = 1, isize2(rankloc+1)
             term = qfac_2d(i,j,k)
             qgridout2(1,i,j,k) = term*qgridout2(1,i,j,k)
             qgridout2(2,i,j,k) = term*qgridout2(2,i,j,k)
           end do
         end do
      end do

      call fft2d_backmpi(qgridin2,qgridout2,n1mpimax,n2mpimax,
     $ n3mpimax)
c
c     MPI : Begin reception
c
      do i = 1, nrec_send
        tag = nprocloc*rankloc + prec_send(i) + 1
        call MPI_IRECV(qgridin2(1,1,1,1,i+1),
     $   2*n1mpimax*n2mpimax*n3mpimax,MPI_REAL8,
     $   prec_recep(i),tag,commloc,req2rec(i),ierr)
      end do
c
c     MPI : begin sending
c
      do i = 1, nrec_recep
        tag = nprocloc*prec_recep(i) + rankloc + 1
        call MPI_ISEND(qgridin2(1,1,1,1,1),
     $   2*n1mpimax*n2mpimax*n3mpimax,MPI_REAL8,
     $   prec_send(i),tag,commloc,req2send(i),ierr)
      end do
c
      do i = 1, nrec_send
        call MPI_WAIT(req2rec(i),status,ierr)
      end do
      do i = 1, nrec_recep
        call MPI_WAIT(req2send(i),status,ierr)
      end do

      do i = 1, npolerecloc
        iipole = polerecglob(i)
        iglob = ipole(iipole)
        call fphi_dervec_site(iglob,i,qgridin2, fphivec(:,i))
      end do
      return
      end

      ! sub cart_to_frac_vec
      ! given a vector a_in, returns the vector in fractional coords

      subroutine cart_to_frac_vec(a_in,a_out)
      implicit none
      integer j,k
      real*8 ctf(10,10)
      real*8,intent(in) :: a_in(3)
      real*8 ::  a_out(3)
c
c
c     find the matrix to convert Cartesian to fractional
c
      call cart_to_frac (ctf)
c
c     apply the transformation to get the fractional multipoles
c
      !a_out(1) = ctf(1,1) * a_in(1)
      do j = 1, 3
         a_out(j) = 0.0d0
         do k = 1, 3
            a_out(j) = a_out(j) + ctf(j+1,k+1)*a_in(k)
         end do
      end do

      end


!===================================================
!     sub fphi_uind_big
!===================================================
! Returns the pot and successive derivatives of it... 
!  ... up to 3rd order ! 
      subroutine fphi_uind_big(isite,impi,fdip_phi1,fdip_phi2)
      use domdec
      use fft
      use pme
      use potent
      use mpi
      implicit none
      integer istart,iend,jstart,jend,kstart,kend
      integer i,j,k,impi
      integer iproc,proc
      integer isite,iatm
      integer i0,j0,k0
      integer it1,it2,it3
      integer igrd0,jgrd0,kgrd0
      real*8 v0,v1,v2,v3
      real*8 u0,u1,u2,u3
      real*8 t0_1,t0_2,t1_1,t1_2, t3_1, t3_2
      real*8 t2_1,t2_2,tq_1,tq_2
      real*8 tu00_1,tu01_1,tu10_1
      real*8 tu00_2,tu01_2,tu10_2
      real*8 tu20_1,tu11_1,tu02_1
      real*8 tu20_2,tu11_2,tu02_2
      real*8 tu30_1,tu21_1,tu12_1,tu03_1
      real*8 tu30_2,tu21_2,tu12_2,tu03_2
      real*8 :: tuv000_1, tuv000_2
      real*8 tuv100_1,tuv010_1,tuv001_1
      real*8 tuv100_2,tuv010_2,tuv001_2
      real*8 tuv200_1,tuv020_1,tuv002_1
      real*8 tuv110_1,tuv101_1,tuv011_1
      real*8 tuv200_2,tuv020_2,tuv002_2
      real*8 tuv110_2,tuv101_2,tuv011_2
      real*8 tuv300_1,tuv030_1
      real*8 tuv003_1,tuv210_1,tuv201_1,tuv120_1
      real*8 tuv021_1,tuv102_1,tuv012_1,tuv111_1
      real*8 tuv300_2,tuv030_2
      real*8 tuv003_2,tuv210_2,tuv201_2,tuv120_2
      real*8 tuv021_2,tuv102_2,tuv012_2,tuv111_2
      real*8 fdip_phi1(20)
      real*8 fdip_phi2(20)
c
      iatm = isite
      igrd0 = igrid(1,iatm)
      jgrd0 = igrid(2,iatm)
      kgrd0 = igrid(3,iatm)
      tuv000_1 = 0.0d0
      tuv100_1 = 0.0d0
      tuv010_1 = 0.0d0
      tuv001_1 = 0.0d0
      tuv200_1 = 0.0d0
      tuv020_1 = 0.0d0
      tuv002_1 = 0.0d0
      tuv110_1 = 0.0d0
      tuv101_1 = 0.0d0
      tuv011_1 = 0.0d0
      tuv300_1 = 0.0d0
      tuv030_1 = 0.0d0
      tuv003_1 = 0.0d0
      tuv210_1 = 0.0d0
      tuv201_1 = 0.0d0
      tuv120_1 = 0.0d0
      tuv021_1 = 0.0d0
      tuv102_1 = 0.0d0
      tuv012_1 = 0.0d0
      tuv111_1 = 0.0d0
      tuv000_2 = 0.0d0
      tuv100_2 = 0.0d0
      tuv010_2 = 0.0d0
      tuv001_2 = 0.0d0
      tuv200_2 = 0.0d0
      tuv020_2 = 0.0d0
      tuv002_2 = 0.0d0
      tuv110_2 = 0.0d0
      tuv101_2 = 0.0d0
      tuv011_2 = 0.0d0
      tuv300_2 = 0.0d0
      tuv030_2 = 0.0d0
      tuv003_2 = 0.0d0
      tuv210_2 = 0.0d0
      tuv201_2 = 0.0d0
      tuv120_2 = 0.0d0
      tuv021_2 = 0.0d0
      tuv102_2 = 0.0d0
      tuv012_2 = 0.0d0
      tuv111_2 = 0.0d0
      k0 = kgrd0
      do it3 = 1, bsorder
         k0 = k0 + 1
         k = k0 + 1 + (nfft3-isign(nfft3,k0))/2
         v0 = thetai3(1,it3,impi)
         v1 = thetai3(2,it3,impi)
         v2 = thetai3(3,it3,impi)
         v3 = thetai3(4,it3,impi)
         tu00_1 = 0.0d0
         tu01_1 = 0.0d0
         tu10_1 = 0.0d0
         tu20_1 = 0.0d0
         tu11_1 = 0.0d0
         tu02_1 = 0.0d0
         tu30_1 = 0.0d0
         tu21_1 = 0.0d0
         tu12_1 = 0.0d0
         tu03_1 = 0.0d0
         tu00_2 = 0.0d0
         tu01_2 = 0.0d0
         tu10_2 = 0.0d0
         tu20_2 = 0.0d0
         tu11_2 = 0.0d0
         tu02_2 = 0.0d0
         tu30_2 = 0.0d0
         tu21_2 = 0.0d0
         tu12_2 = 0.0d0
         tu03_2 = 0.0d0
         j0 = jgrd0
         do it2 = 1, bsorder
            j0 = j0 + 1
            j = j0 + 1 + (nfft2-isign(nfft2,j0))/2
            u0 = thetai2(1,it2,impi)
            u1 = thetai2(2,it2,impi)
            u2 = thetai2(3,it2,impi)
            u3 = thetai2(4,it2,impi)
            t0_1 = 0.0d0
            t1_1 = 0.0d0
            t2_1 = 0.0d0
            t3_1 = 0.0d0
            t0_2 = 0.0d0
            t1_2 = 0.0d0
            t2_2 = 0.0d0
            t3_2 = 0.0d0
            i0 = igrd0
            do it1 = 1, bsorder
               i0 = i0 + 1
               i = i0 + 1 + (nfft1-isign(nfft1,i0))/2
c
               if (use_pmecore) then
                 kstart = kstart1(rank_bis+1)
                 kend = kend1(rank_bis+1)
                 jstart = jstart1(rank_bis+1)
                 jend = jend1(rank_bis+1)
                 istart = istart1(rank_bis+1)
                 iend = iend1(rank_bis+1)
               else
                 kstart = kstart1(rank+1)
                 kend = kend1(rank+1)
                 jstart = jstart1(rank+1)
                 jend = jend1(rank+1)
                 istart = istart1(rank+1)
                 iend = iend1(rank+1)
               end if
               if (((k.ge.kstart).and.(k.le.kend)).and.
     $           ((j.ge.jstart).and.(j.le.jend)).and.
     $           ((i.ge.istart).and.(i.le.iend))) then
                 tq_1  = qgrid2in_2d(1,i-istart+1,j-jstart+1,k-kstart+1,
     $             1)
                 tq_2  = qgrid2in_2d(2,i-istart+1,j-jstart+1,k-kstart+1,
     $             1)
                 goto 10
               end if
               do iproc = 1, nrec_send
                 proc = prec_send(iproc)
                 kstart = kstart1(proc+1)
                 kend = kend1(proc+1)
                 jstart = jstart1(proc+1)
                 jend = jend1(proc+1)
                 istart = istart1(proc+1)
                 iend = iend1(proc+1)
                 if (((k.ge.kstart).and.(k.le.kend)).and.
     $             ((j.ge.jstart).and.(j.le.jend)).and.
     $             ((i.ge.istart).and.(i.le.iend))) then
                   tq_1  = qgrid2in_2d(1,i-istart+1,j-jstart+1,
     $               k-kstart+1,iproc+1)
                   tq_2  = qgrid2in_2d(2,i-istart+1,j-jstart+1,
     $               k-kstart+1,iproc+1)
                   goto 10
                 end if
               end do
 10            continue
c
               t0_1 = t0_1 + tq_1*thetai1(1,it1,impi)
               t1_1 = t1_1 + tq_1*thetai1(2,it1,impi)
               t2_1 = t2_1 + tq_1*thetai1(3,it1,impi)
               t3_1 = t3_1 + tq_1*thetai1(4,it1,impi)
               t0_2 = t0_2 + tq_2*thetai1(1,it1,impi)
               t1_2 = t1_2 + tq_2*thetai1(2,it1,impi)
               t2_2 = t2_2 + tq_2*thetai1(3,it1,impi)
               t3_2 = t3_2 + tq_2*thetai1(4,it1,impi)
            end do
            tu00_1 = tu00_1 + t0_1*u0
            tu10_1 = tu10_1 + t1_1*u0
            tu01_1 = tu01_1 + t0_1*u1
            tu20_1 = tu20_1 + t2_1*u0
            tu11_1 = tu11_1 + t1_1*u1
            tu02_1 = tu02_1 + t0_1*u2
            tu30_1 = tu30_1 + t3_1*u0
            tu21_1 = tu21_1 + t2_1*u1
            tu12_1 = tu12_1 + t1_1*u2
            tu03_1 = tu03_1 + t0_1*u3
            tu00_2 = tu00_2 + t0_2*u0
            tu10_2 = tu10_2 + t1_2*u0
            tu01_2 = tu01_2 + t0_2*u1
            tu20_2 = tu20_2 + t2_2*u0
            tu11_2 = tu11_2 + t1_2*u1
            tu02_2 = tu02_2 + t0_2*u2
            tu30_2 = tu30_2 + t3_2*u0
            tu21_2 = tu21_2 + t2_2*u1
            tu12_2 = tu12_2 + t1_2*u2
            tu03_2 = tu03_2 + t0_2*u3
         end do
         tuv000_1 = tuv000_1 + tu00_1*v0
         tuv100_1 = tuv100_1 + tu10_1*v0
         tuv010_1 = tuv010_1 + tu01_1*v0
         tuv001_1 = tuv001_1 + tu00_1*v1
         tuv200_1 = tuv200_1 + tu20_1*v0
         tuv020_1 = tuv020_1 + tu02_1*v0
         tuv002_1 = tuv002_1 + tu00_1*v2
         tuv110_1 = tuv110_1 + tu11_1*v0
         tuv101_1 = tuv101_1 + tu10_1*v1
         tuv011_1 = tuv011_1 + tu01_1*v1
         tuv300_1 = tuv300_1 + tu30_1*v0
         tuv030_1 = tuv030_1 + tu03_1*v0
         tuv003_1 = tuv003_1 + tu00_1*v3
         tuv210_1 = tuv210_1 + tu21_1*v0
         tuv201_1 = tuv201_1 + tu20_1*v1
         tuv120_1 = tuv120_1 + tu12_1*v0
         tuv021_1 = tuv021_1 + tu02_1*v1
         tuv102_1 = tuv102_1 + tu10_1*v2
         tuv012_1 = tuv012_1 + tu01_1*v2
         tuv111_1 = tuv111_1 + tu11_1*v1
         tuv000_2 = tuv000_2 + tu00_2*v0
         tuv100_2 = tuv100_2 + tu10_2*v0
         tuv010_2 = tuv010_2 + tu01_2*v0
         tuv001_2 = tuv001_2 + tu00_2*v1
         tuv200_2 = tuv200_2 + tu20_2*v0
         tuv020_2 = tuv020_2 + tu02_2*v0
         tuv002_2 = tuv002_2 + tu00_2*v2
         tuv110_2 = tuv110_2 + tu11_2*v0
         tuv101_2 = tuv101_2 + tu10_2*v1
         tuv011_2 = tuv011_2 + tu01_2*v1
         tuv300_2 = tuv300_2 + tu30_2*v0
         tuv030_2 = tuv030_2 + tu03_2*v0
         tuv003_2 = tuv003_2 + tu00_2*v3
         tuv210_2 = tuv210_2 + tu21_2*v0
         tuv201_2 = tuv201_2 + tu20_2*v1
         tuv120_2 = tuv120_2 + tu12_2*v0
         tuv021_2 = tuv021_2 + tu02_2*v1
         tuv102_2 = tuv102_2 + tu10_2*v2
         tuv012_2 = tuv012_2 + tu01_2*v2
         tuv111_2 = tuv111_2 + tu11_2*v1
      end do
      fdip_phi1(1)  = tuv000_1
      fdip_phi1(2)  = tuv100_1
      fdip_phi1(3)  = tuv010_1
      fdip_phi1(4)  = tuv001_1
      fdip_phi1(5)  = tuv200_1
      fdip_phi1(6)  = tuv020_1
      fdip_phi1(7)  = tuv002_1
      fdip_phi1(8)  = tuv110_1
      fdip_phi1(9)  = tuv101_1
      fdip_phi1(10) = tuv011_1
      fdip_phi1(11) = tuv300_1
      fdip_phi1(12) = tuv030_1
      fdip_phi1(13) = tuv003_1
      fdip_phi1(14) = tuv210_1
      fdip_phi1(15) = tuv201_1
      fdip_phi1(16) = tuv120_1
      fdip_phi1(17) = tuv021_1
      fdip_phi1(18) = tuv102_1
      fdip_phi1(19) = tuv012_1
      fdip_phi1(20) = tuv111_1
      fdip_phi2(1)  = tuv000_2
      fdip_phi2(2)  = tuv100_2
      fdip_phi2(3)  = tuv010_2
      fdip_phi2(4)  = tuv001_2
      fdip_phi2(5)  = tuv200_2
      fdip_phi2(6)  = tuv020_2
      fdip_phi2(7)  = tuv002_2
      fdip_phi2(8)  = tuv110_2
      fdip_phi2(9)  = tuv101_2
      fdip_phi2(10) = tuv011_2
      fdip_phi2(11) = tuv300_2
      fdip_phi2(12) = tuv030_2
      fdip_phi2(13) = tuv003_2
      fdip_phi2(14) = tuv210_2
      fdip_phi2(15) = tuv201_2
      fdip_phi2(16) = tuv120_2
      fdip_phi2(17) = tuv021_2
      fdip_phi2(18) = tuv102_2
      fdip_phi2(19) = tuv012_2
      fdip_phi2(20) = tuv111_2
c
      return
      end


c
c       "fphi_dervec_site2" extracts the "equivalent" potential on the i-th site from
c       the given particle mesh Ewald grid for both scalings p and d
c
c
      subroutine fphi_dervec_site2(isite,impi, grid, fdervec)
      use domdec
      use fft
      use mpole
      use pme
      use potent
      use mpi
      implicit none
      integer istart,iend,jstart,jend,kstart,kend
      integer i,j,k,impi,rankloc
      integer iproc,proc
      integer isite,iatm
      integer i0,j0,k0
      integer it1,it2,it3
      integer igrd0,jgrd0,kgrd0
      real*8 v0,v1,v2,v3
      real*8 u0,u1,u2,u3
      real*8 t0,t1,t2,t3,tq
      real*8 t0_2,t1_2,t2_2,t3_2,tq_2
      real*8 tu00,tu10,tu01,tu20,tu11
      real*8 tu02,tu21,tu12,tu30,tu03
      real*8 tu00_2,tu10_2,tu01_2,tu20_2,tu11_2
      real*8 tu02_2,tu21_2,tu12_2,tu30_2,tu03_2
      real*8 tuv000,tuv100,tuv010,tuv001
      real*8 tuv200,tuv020,tuv002,tuv110
      real*8 tuv101,tuv011,tuv300,tuv030
      real*8 tuv003,tuv210,tuv201,tuv120
      real*8 tuv021,tuv102,tuv012,tuv111
      real*8 tuv000_2,tuv100_2,tuv010_2,tuv001_2
      real*8 tuv200_2,tuv020_2,tuv002_2,tuv110_2
      real*8 tuv101_2,tuv011_2,tuv300_2,tuv030_2
      real*8 tuv003_2,tuv210_2,tuv201_2,tuv120_2
      real*8 tuv021_2,tuv102_2,tuv012_2,tuv111_2
c      real*8 fphi(20)

      real*8, dimension(20,2) :: fdervec
      real*8, dimension(2,n1mpimax,n2mpimax,n3mpimax,nrec_send+1) :: 
     $        grid
c
      if (use_pmecore) then
        rankloc  = rank_bis
      else
        rankloc  = rank
      end if
c
      iatm = isite
      igrd0 = igrid(1,iatm)
      jgrd0 = igrid(2,iatm)
      kgrd0 = igrid(3,iatm)
      tuv000 = 0.0d0
      tuv001 = 0.0d0
      tuv010 = 0.0d0
      tuv100 = 0.0d0
      tuv200 = 0.0d0
      tuv020 = 0.0d0
      tuv002 = 0.0d0
      tuv110 = 0.0d0
      tuv101 = 0.0d0
      tuv011 = 0.0d0
      tuv300 = 0.0d0
      tuv030 = 0.0d0
      tuv003 = 0.0d0
      tuv210 = 0.0d0
      tuv201 = 0.0d0
      tuv120 = 0.0d0
      tuv021 = 0.0d0
      tuv102 = 0.0d0
      tuv012 = 0.0d0
      tuv111 = 0.0d0
      tuv000_2 = 0.0d0
      tuv001_2 = 0.0d0
      tuv010_2 = 0.0d0
      tuv100_2 = 0.0d0
      tuv200_2 = 0.0d0
      tuv020_2 = 0.0d0
      tuv002_2 = 0.0d0
      tuv110_2 = 0.0d0
      tuv101_2 = 0.0d0
      tuv011_2 = 0.0d0
      tuv300_2 = 0.0d0
      tuv030_2 = 0.0d0
      tuv003_2 = 0.0d0
      tuv210_2 = 0.0d0
      tuv201_2 = 0.0d0
      tuv120_2 = 0.0d0
      tuv021_2 = 0.0d0
      tuv102_2 = 0.0d0
      tuv012_2 = 0.0d0
      tuv111_2 = 0.0d0
      k0 = kgrd0
      do it3 = 1, bsorder
         k0 = k0 + 1
         k = k0 + 1 + (nfft3-isign(nfft3,k0))/2
         v0 = thetai3(1,it3,impi)
         v1 = thetai3(2,it3,impi)
         v2 = thetai3(3,it3,impi)
         v3 = thetai3(4,it3,impi)
         tu00 = 0.0d0
         tu10 = 0.0d0
         tu01 = 0.0d0
         tu20 = 0.0d0
         tu11 = 0.0d0
         tu02 = 0.0d0
         tu30 = 0.0d0
         tu21 = 0.0d0
         tu12 = 0.0d0
         tu03 = 0.0d0
         tu00_2 = 0.0d0
         tu10_2 = 0.0d0
         tu01_2 = 0.0d0
         tu20_2 = 0.0d0
         tu11_2 = 0.0d0
         tu02_2 = 0.0d0
         tu30_2 = 0.0d0
         tu21_2 = 0.0d0
         tu12_2 = 0.0d0
         tu03_2 = 0.0d0
         j0 = jgrd0
         do it2 = 1, bsorder
            j0 = j0 + 1
            j = j0 + 1 + (nfft2-isign(nfft2,j0))/2
            u0 = thetai2(1,it2,impi)
            u1 = thetai2(2,it2,impi)
            u2 = thetai2(3,it2,impi)
            u3 = thetai2(4,it2,impi)
            t0 = 0.0d0
            t1 = 0.0d0
            t2 = 0.0d0
            t3 = 0.0d0
            t0_2 = 0.0d0
            t1_2 = 0.0d0
            t2_2 = 0.0d0
            t3_2 = 0.0d0
            i0 = igrd0
            do it1 = 1, bsorder
               i0 = i0 + 1
               i = i0 + 1 + (nfft1-isign(nfft1,i0))/2
c
                tq = 0.0d0
                kstart = kstart1(rankloc+1)
                kend = kend1(rankloc+1)
                jstart = jstart1(rankloc+1)
                jend = jend1(rankloc+1)
                istart = istart1(rankloc+1)
                iend = iend1(rankloc+1)
               if (((k.ge.kstart).and.(k.le.kend)).and.
     $           ((j.ge.jstart).and.(j.le.jend)).and.
     $           ((i.ge.istart).and.(i.le.iend))) then
                 tq  = grid(1,i-istart+1,j-jstart+1,k-kstart+1,1)
                 tq_2  = grid(2,i-istart+1,j-jstart+1,k-kstart+1,1)
                 goto 10
               end if
               do iproc = 1, nrec_send
                 proc = prec_send(iproc)
                 kstart = kstart1(proc+1)
                 kend = kend1(proc+1)
                 jstart = jstart1(proc+1)
                 jend = jend1(proc+1)
                 istart = istart1(proc+1)
                 iend = iend1(proc+1)
                 if (((k.ge.kstart).and.(k.le.kend)).and.
     $             ((j.ge.jstart).and.(j.le.jend)).and.
     $             ((i.ge.istart).and.(i.le.iend))) then
                   tq  = grid(1,i-istart+1,j-jstart+1,k-kstart+1,
     $              iproc+1)
                   tq_2 =grid(2,i-istart+1,j-jstart+1,k-kstart+1,
     $              iproc+1)
                   goto 10
                 end if
               end do
 10            continue
c
               t0 = t0 + tq*thetai1(1,it1,impi)
               t1 = t1 + tq*thetai1(2,it1,impi)
               t2 = t2 + tq*thetai1(3,it1,impi)
               t3 = t3 + tq*thetai1(4,it1,impi)
               t0_2 = t0_2 + tq_2*thetai1(1,it1,impi)
               t1_2 = t1_2 + tq_2*thetai1(2,it1,impi)
               t2_2 = t2_2 + tq_2*thetai1(3,it1,impi)
               t3_2 = t3_2 + tq_2*thetai1(4,it1,impi)
            end do
            tu00 = tu00 + t0*u0
            tu10 = tu10 + t1*u0
            tu01 = tu01 + t0*u1
            tu20 = tu20 + t2*u0
            tu11 = tu11 + t1*u1
            tu02 = tu02 + t0*u2
            tu30 = tu30 + t3*u0
            tu21 = tu21 + t2*u1
            tu12 = tu12 + t1*u2
            tu03 = tu03 + t0*u3
            tu00_2 = tu00_2 + t0_2*u0
            tu10_2 = tu10_2 + t1_2*u0
            tu01_2 = tu01_2 + t0_2*u1
            tu20_2 = tu20_2 + t2_2*u0
            tu11_2 = tu11_2 + t1_2*u1
            tu02_2 = tu02_2 + t0_2*u2
            tu30_2 = tu30_2 + t3_2*u0
            tu21_2 = tu21_2 + t2_2*u1
            tu12_2 = tu12_2 + t1_2*u2
            tu03_2 = tu03_2 + t0_2*u3
         end do
         tuv000 = tuv000 + tu00*v0
         tuv100 = tuv100 + tu10*v0
         tuv010 = tuv010 + tu01*v0
         tuv001 = tuv001 + tu00*v1
         tuv200 = tuv200 + tu20*v0
         tuv020 = tuv020 + tu02*v0
         tuv002 = tuv002 + tu00*v2
         tuv110 = tuv110 + tu11*v0
         tuv101 = tuv101 + tu10*v1
         tuv011 = tuv011 + tu01*v1
         tuv300 = tuv300 + tu30*v0
         tuv030 = tuv030 + tu03*v0
         tuv003 = tuv003 + tu00*v3
         tuv210 = tuv210 + tu21*v0
         tuv201 = tuv201 + tu20*v1
         tuv120 = tuv120 + tu12*v0
         tuv021 = tuv021 + tu02*v1
         tuv102 = tuv102 + tu10*v2
         tuv012 = tuv012 + tu01*v2
         tuv111 = tuv111 + tu11*v1
         tuv000_2 = tuv000_2 + tu00_2*v0
         tuv100_2 = tuv100_2 + tu10_2*v0
         tuv010_2 = tuv010_2 + tu01_2*v0
         tuv001_2 = tuv001_2 + tu00_2*v1
         tuv200_2 = tuv200_2 + tu20_2*v0
         tuv020_2 = tuv020_2 + tu02_2*v0
         tuv002_2 = tuv002_2 + tu00_2*v2
         tuv110_2 = tuv110_2 + tu11_2*v0
         tuv101_2 = tuv101_2 + tu10_2*v1
         tuv011_2 = tuv011_2 + tu01_2*v1
         tuv300_2 = tuv300_2 + tu30_2*v0
         tuv030_2 = tuv030_2 + tu03_2*v0
         tuv003_2 = tuv003_2 + tu00_2*v3
         tuv210_2 = tuv210_2 + tu21_2*v0
         tuv201_2 = tuv201_2 + tu20_2*v1
         tuv120_2 = tuv120_2 + tu12_2*v0
         tuv021_2 = tuv021_2 + tu02_2*v1
         tuv102_2 = tuv102_2 + tu10_2*v2
         tuv012_2 = tuv012_2 + tu01_2*v2
         tuv111_2 = tuv111_2 + tu11_2*v1
      end do
      fdervec( 1,1) = tuv000
      fdervec( 2,1) = tuv100
      fdervec( 3,1) = tuv010
      fdervec( 4,1) = tuv001
      fdervec( 5,1) = tuv200
      fdervec( 6,1) = tuv020
      fdervec( 7,1) = tuv002
      fdervec( 8,1) = tuv110
      fdervec( 9,1) = tuv101
      fdervec(10,1) = tuv011
      fdervec(11,1) = tuv300
      fdervec(12,1) = tuv030
      fdervec(13,1) = tuv003
      fdervec(14,1) = tuv210
      fdervec(15,1) = tuv201
      fdervec(16,1) = tuv120
      fdervec(17,1) = tuv021
      fdervec(18,1) = tuv102
      fdervec(19,1) = tuv012
      fdervec(20,1) = tuv111
      fdervec( 1,2) = tuv000_2
      fdervec( 2,2) = tuv100_2
      fdervec( 3,2) = tuv010_2
      fdervec( 4,2) = tuv001_2
      fdervec( 5,2) = tuv200_2
      fdervec( 6,2) = tuv020_2
      fdervec( 7,2) = tuv002_2
      fdervec( 8,2) = tuv110_2
      fdervec( 9,2) = tuv101_2
      fdervec(10,2) = tuv011_2
      fdervec(11,2) = tuv300_2
      fdervec(12,2) = tuv030_2
      fdervec(13,2) = tuv003_2
      fdervec(14,2) = tuv210_2
      fdervec(15,2) = tuv201_2
      fdervec(16,2) = tuv120_2
      fdervec(17,2) = tuv021_2
      fdervec(18,2) = tuv102_2
      fdervec(19,2) = tuv012_2
      fdervec(20,2) = tuv111_2
      return
      end

c
c       "fphi_dervec_site" extracts the "equivalent" potential on the i-th site from
c       the given particle mesh Ewald grid
c
c
      subroutine fphi_dervec_site(isite,impi, grid, fdervec)
      use domdec
      use fft
      use mpole
      use pme
      use potent
      use mpi
      implicit none
      integer istart,iend,jstart,jend,kstart,kend
      integer i,j,k,impi,rankloc
      integer iproc,proc
      integer isite,iatm
      integer i0,j0,k0
      integer it1,it2,it3
      integer igrd0,jgrd0,kgrd0
      real*8 v0,v1,v2,v3
      real*8 u0,u1,u2,u3
      real*8 t0,t1,t2,t3,tq

      real*8 tu00,tu10,tu01,tu20,tu11
      real*8 tu02,tu21,tu12,tu30,tu03


      real*8 tuv000,tuv100,tuv010,tuv001
      real*8 tuv200,tuv020,tuv002,tuv110
      real*8 tuv101,tuv011,tuv300,tuv030
      real*8 tuv003,tuv210,tuv201,tuv120
      real*8 tuv021,tuv102,tuv012,tuv111





c      real*8 fphi(20)

      real*8, dimension(20) :: fdervec
      real*8, dimension(2,n1mpimax,n2mpimax,n3mpimax,nrec_send+1) :: 
     $        grid
c
      if (use_pmecore) then
        rankloc  = rank_bis
      else
        rankloc  = rank
      end if
c
      iatm = isite
      igrd0 = igrid(1,iatm)
      jgrd0 = igrid(2,iatm)
      kgrd0 = igrid(3,iatm)
      tuv000 = 0.0d0
      tuv001 = 0.0d0
      tuv010 = 0.0d0
      tuv100 = 0.0d0
      tuv200 = 0.0d0
      tuv020 = 0.0d0
      tuv002 = 0.0d0
      tuv110 = 0.0d0
      tuv101 = 0.0d0
      tuv011 = 0.0d0
      tuv300 = 0.0d0
      tuv030 = 0.0d0
      tuv003 = 0.0d0
      tuv210 = 0.0d0
      tuv201 = 0.0d0
      tuv120 = 0.0d0
      tuv021 = 0.0d0
      tuv102 = 0.0d0
      tuv012 = 0.0d0
      tuv111 = 0.0d0
      k0 = kgrd0
      do it3 = 1, bsorder
         k0 = k0 + 1
         k = k0 + 1 + (nfft3-isign(nfft3,k0))/2
         v0 = thetai3(1,it3,impi)
         v1 = thetai3(2,it3,impi)
         v2 = thetai3(3,it3,impi)
         v3 = thetai3(4,it3,impi)
         tu00 = 0.0d0
         tu10 = 0.0d0
         tu01 = 0.0d0
         tu20 = 0.0d0
         tu11 = 0.0d0
         tu02 = 0.0d0
         tu30 = 0.0d0
         tu21 = 0.0d0
         tu12 = 0.0d0
         tu03 = 0.0d0
         j0 = jgrd0
         do it2 = 1, bsorder
            j0 = j0 + 1
            j = j0 + 1 + (nfft2-isign(nfft2,j0))/2
            u0 = thetai2(1,it2,impi)
            u1 = thetai2(2,it2,impi)
            u2 = thetai2(3,it2,impi)
            u3 = thetai2(4,it2,impi)
            t0 = 0.0d0
            t1 = 0.0d0
            t2 = 0.0d0
            t3 = 0.0d0
            i0 = igrd0
            do it1 = 1, bsorder
               i0 = i0 + 1
               i = i0 + 1 + (nfft1-isign(nfft1,i0))/2
c
                tq = 0.0d0
                kstart = kstart1(rankloc+1)
                kend = kend1(rankloc+1)
                jstart = jstart1(rankloc+1)
                jend = jend1(rankloc+1)
                istart = istart1(rankloc+1)
                iend = iend1(rankloc+1)
               if (((k.ge.kstart).and.(k.le.kend)).and.
     $           ((j.ge.jstart).and.(j.le.jend)).and.
     $           ((i.ge.istart).and.(i.le.iend))) then
                 tq  = grid(1,i-istart+1,j-jstart+1,k-kstart+1,1)
                 goto 10
               end if
               do iproc = 1, nrec_send
                 proc = prec_send(iproc)
                 kstart = kstart1(proc+1)
                 kend = kend1(proc+1)
                 jstart = jstart1(proc+1)
                 jend = jend1(proc+1)
                 istart = istart1(proc+1)
                 iend = iend1(proc+1)
                 if (((k.ge.kstart).and.(k.le.kend)).and.
     $             ((j.ge.jstart).and.(j.le.jend)).and.
     $             ((i.ge.istart).and.(i.le.iend))) then
                   tq  = grid(1,i-istart+1,j-jstart+1,k-kstart+1,
     $              iproc+1)
                   goto 10
                 end if
               end do
 10            continue
c
               t0 = t0 + tq*thetai1(1,it1,impi)
               t1 = t1 + tq*thetai1(2,it1,impi)
               t2 = t2 + tq*thetai1(3,it1,impi)
               t3 = t3 + tq*thetai1(4,it1,impi)
            end do
            tu00 = tu00 + t0*u0
            tu10 = tu10 + t1*u0
            tu01 = tu01 + t0*u1
            tu20 = tu20 + t2*u0
            tu11 = tu11 + t1*u1
            tu02 = tu02 + t0*u2
            tu30 = tu30 + t3*u0
            tu21 = tu21 + t2*u1
            tu12 = tu12 + t1*u2
            tu03 = tu03 + t0*u3
         end do
         tuv000 = tuv000 + tu00*v0
         tuv100 = tuv100 + tu10*v0
         tuv010 = tuv010 + tu01*v0
         tuv001 = tuv001 + tu00*v1
         tuv200 = tuv200 + tu20*v0
         tuv020 = tuv020 + tu02*v0
         tuv002 = tuv002 + tu00*v2
         tuv110 = tuv110 + tu11*v0
         tuv101 = tuv101 + tu10*v1
         tuv011 = tuv011 + tu01*v1
         tuv300 = tuv300 + tu30*v0
         tuv030 = tuv030 + tu03*v0
         tuv003 = tuv003 + tu00*v3
         tuv210 = tuv210 + tu21*v0
         tuv201 = tuv201 + tu20*v1
         tuv120 = tuv120 + tu12*v0
         tuv021 = tuv021 + tu02*v1
         tuv102 = tuv102 + tu10*v2
         tuv012 = tuv012 + tu01*v2
         tuv111 = tuv111 + tu11*v1
      end do
      fdervec(1) = tuv000
      fdervec(2) = tuv100
      fdervec(3) = tuv010
      fdervec(4) = tuv001
      fdervec(5) = tuv200
      fdervec(6) = tuv020
      fdervec(7) = tuv002
      fdervec(8) = tuv110
      fdervec(9) = tuv101
      fdervec(10) = tuv011
      fdervec(11) = tuv300
      fdervec(12) = tuv030
      fdervec(13) = tuv003
      fdervec(14) = tuv210
      fdervec(15) = tuv201
      fdervec(16) = tuv120
      fdervec(17) = tuv021
      fdervec(18) = tuv102
      fdervec(19) = tuv012
      fdervec(20) = tuv111
      return
      end

! ================================================================
!   subroutine torque_prods
! ================================================================
! A small utility sub to generate the products (well especially to
! lighten the main sub)
!
! d(in)(3) : dipole (unrot)
! q2(in)(3,3) : quad (local, unrotated) to be derivated
! rmat(in)(3,3) : rotation matrix
! dri(in)(3,3,3) : deriv. of rotation matrix w.r. to i
! di(out) (3,3)  : derrotated d
! qi(out) (3,3,3): derrotated q2
!

      subroutine torque_prods(d, q2, rmat, dri,di, qi )
      implicit none

      integer :: l1, l2, l3, n1, n2
      real*8, dimension(3)   :: d
      real*8, dimension(3,3) :: di, q2, rmat
      real*8, dimension(3,3,3) :: qi, dri

      call aclear(9,di)
      call aclear(27,qi)
      do l1 = 1,3
         do l2 = 1,3
            do l3 = 1,3
               di(l3,l1) = di(l3,l1)
     $          + dri(l3,l1,l2)*d(l2)
            end do
            do n1 = 1,3
               do n2 = 1,3
                  do l3 = 1,3
                     qi(l3,l1,l2) = qi(l3,l1,l2)
     $                  + q2(n1,n2)*(dri(l3,l1,n1)
     $                  *rmat(l2,n2)
     $                  + dri(l3,l2,n2)*rmat(l1,n1))
                  end do
               end do
            end do
         end do
      end do
      return
      end subroutine
c
c     subroutine torque_dir: computes contribution of torques to derivatives of 
c      polarization energy with tcg (direct part)
c
      subroutine torquetcg_dir(torq_mu,torq_t,denedmu,denedt)
      use atmlst
      use domdec
      use mpole
      implicit none
      integer i,ii,iipole,iglob,iz,ix,iy
      integer jj,jjpole,jglob,jloc,jpoleloc
      integer jx,jy,jz,alphac,betac,jxloc,jyloc,jzloc
      logical doi,doix,doiy,doiz
      real*8 :: torq_mu(3,nbloc),torq_t(3,nbloc)
      real*8 :: denedmu(3,npolebloc),denedt(3,3,npolebloc)
      real*8, allocatable :: dtorquemu(:,:,:,:),dtorquetheta(:,:,:,:,:)
      real*8 sprod
      real*8, dimension(3) :: d
      real*8, dimension(3,3) :: q2, qr, di, rmat
      real*8, dimension(3,3,3) :: qi, dri, drix, driy, driz
c
      allocate (dtorquemu(3,4,3,nbloc))
      allocate (dtorquetheta(3,4,3,3,nbloc))
c
      torq_mu = 0d0
      torq_t = 0d0
      dtorquemu = 0d0
      dtorquetheta = 0d0

      do ii = 1, npolelocnl
         iipole = poleglobnl(ii)
         iglob = ipole(iipole)
         i = loc(iglob)
         iz = zaxis(iipole)
         ix = xaxis(iipole)
         iy = yaxis(iipole)
         do alphac = 1,3
            d(alphac) = pole(alphac + 1,iipole)
            q2(1,alphac) = pole(4+alphac,iipole)
            q2(2,alphac) = pole(7+alphac, iipole)
            q2(3,alphac) = pole(10+alphac, iipole)
            qr(1,alphac) = rpole(4+alphac,iipole)
            qr(2,alphac) = rpole(7+alphac, iipole)
            qr(3,alphac) = rpole(10+alphac, iipole)
         end do
        
         call derrot(iipole,.true.,iglob,iz,ix,iy,rmat,dri,
     $   driz,drix,driy)
         doi = .true.
         doix = .false.
         doiy = .false.
         doiz = .false.
         if (iz.gt.0) doiz = .true.
         if (ix.gt.0) doix = .true.
         if (iy.gt.0) doiy = .true.

         if (doi) then
            call torque_prods(d,q2,rmat,dri,di,qi)
            dtorquemu(:,1,:,i) = di
            dtorquetheta(:,1,:,:,i) = qi
             
         end if
         if (doix) then
            call torque_prods(d,q2,rmat,drix,di,qi)
            dtorquemu(:,2,:,i) = di
            dtorquetheta(:,2,:,:,i) = qi

         end if
         if (doiy) then
            call torque_prods(d,q2,rmat,driy,di,qi)
            dtorquemu(:,3,:,i) = di
            dtorquetheta(:,3,:,:,i) = qi

         end if
         if (doiz) then
            call torque_prods(d,q2,rmat,driz,di,qi)
            dtorquemu(:,4,:,i) = di
            dtorquetheta(:,4,:,:,i) = qi

         end if
      end do ! ii

      do jj = 1, npolelocnl
         jjpole  = poleglobnl(jj)
         jglob = ipole(jjpole)
         jloc  = loc(jglob)
         jpoleloc = poleloc(jjpole)
         jx = xaxis(jjpole)
         if (jx.gt.0) jxloc = loc(jx)
         jy = yaxis(jjpole)
         if (jy.gt.0) jyloc = loc(jy)
         jz = zaxis(jjpole)
         if (jz.gt.0) jzloc = loc(jz)
         do betac = 1,3
            torq_mu(betac,jloc) = torq_mu(betac,jloc)
     $      + sprod(3,dEnedmu(:,jpoleloc),dtorquemu(betac,1,:,jloc))

            torq_t(betac,jloc) = torq_t(betac,jloc)
     $      + sprod(9,dEnedt(:,:,jpoleloc),
     $      dtorquetheta(betac,1,:,:,jloc))

            if (jx .gt. 0) then
               torq_mu(betac,jxloc) = torq_mu(betac,jxloc)
     $       + sprod(3,dEnedmu(:,jpoleloc), dtorquemu(betac,2,:,jloc))
             torq_t(betac,jxloc) = torq_t(betac,jxloc)
     $       + sprod(9,dEnedt(:,:,jpoleloc), 
     $       dtorquetheta(betac,2,:,:,jloc))
            end if
            if (jy .gt. 0) then
               torq_mu(betac,jyloc) = torq_mu(betac,jyloc)
     $       + sprod(3,dEnedmu(:,jpoleloc), dtorquemu(betac,3,:,jloc))
             torq_t(betac,jyloc) = torq_t(betac,jyloc)
     $       + sprod(9,dEnedt(:,:,jpoleloc), 
     $       dtorquetheta(betac,3,:,:,jloc))
            end if
            if (jz .gt. 0) then
               torq_mu(betac,jzloc) = torq_mu(betac,jzloc)
     $       + sprod(3,dEnedmu(:,jpoleloc), dtorquemu(betac,4,:,jloc))
             torq_t(betac,jzloc) = torq_t(betac,jzloc)
     $       + sprod(9,dEnedt(:,:,jpoleloc), 
     $          dtorquetheta(betac,4,:,:,jloc))
            end if
         end do
      end do
c
      deallocate (dtorquemu,dtorquetheta)
      return
      end subroutine
c
c     subroutine torque_rec: computes contribution of torques to derivatives of 
c      polarization energy with tcg (reciprocal part)
c
      subroutine torquetcg_rec(torq_mu,torq_t,denedmu,denedt)
      use atmlst
      use domdec
      use mpole
      implicit none
      integer ii,iipole,iglob,iz,ix,iy
      integer jj,jjj,jjpole,jglob
      integer jx,jy,jz,alphac,betac,jxloc,jyloc,jzloc
      logical doi,doix,doiy,doiz
      real*8 :: torq_mu(3,nblocrec),torq_t(3,nblocrec)
      real*8 :: denedmu(3,nlocrec),denedt(3,3,nlocrec)
      real*8, allocatable :: dtorquemu(:,:,:,:),dtorquetheta(:,:,:,:,:)
      real*8 sprod
      real*8, dimension(3) :: d
      real*8, dimension(3,3) :: q2, qr, di, rmat
      real*8, dimension(3,3,3) :: qi, dri, drix, driy, driz
c
      allocate (dtorquemu(3,4,3,npolerecloc))
      allocate (dtorquetheta(3,4,3,3,npolerecloc))
c
      torq_mu = 0d0
      torq_t = 0d0
      dtorquemu = 0d0
      dtorquetheta = 0d0

      do ii = 1, npolerecloc
         iipole = polerecglob(ii)
         iglob = ipole(iipole)
         iz = zaxis(iipole)
         ix = xaxis(iipole)
         iy = yaxis(iipole)
         do alphac = 1,3
            d(alphac) = pole(alphac + 1,iipole)
            q2(1,alphac) = pole(4+alphac,iipole)
            q2(2,alphac) = pole(7+alphac, iipole)
            q2(3,alphac) = pole(10+alphac, iipole)
            qr(1,alphac) = rpole(4+alphac,iipole)
            qr(2,alphac) = rpole(7+alphac, iipole)
            qr(3,alphac) = rpole(10+alphac, iipole)
         end do
        
         call derrot(iipole,.true.,iglob,iz,ix,iy,rmat,dri,
     $   driz,drix,driy)
         doi = .true.
         doix = .false.
         doiy = .false.
         doiz = .false.
         if (iz.gt.0) doiz = .true.
         if (ix.gt.0) doix = .true.
         if (iy.gt.0) doiy = .true.

         if (doi) then
            call torque_prods(d,q2,rmat,dri,di,qi)
            dtorquemu(:,1,:,ii) = di
            dtorquetheta(:,1,:,:,ii) = qi
             
         end if
         if (doix) then
            call torque_prods(d,q2,rmat,drix,di,qi)
            dtorquemu(:,2,:,ii) = di
            dtorquetheta(:,2,:,:,ii) = qi

         end if
         if (doiy) then
            call torque_prods(d,q2,rmat,driy,di,qi)
            dtorquemu(:,3,:,ii) = di
            dtorquetheta(:,3,:,:,ii) = qi

         end if
         if (doiz) then
            call torque_prods(d,q2,rmat,driz,di,qi)
            dtorquemu(:,4,:,ii) = di
            dtorquetheta(:,4,:,:,ii) = qi

         end if
      end do ! ii

      do jj = 1, npolerecloc
         jjpole  = polerecglob(jj)
         jglob = ipole(jjpole)
         jjj = locrec1(jglob)
         jx = xaxis(jjpole)
         if (jx.gt.0) jxloc = locrec1(jx)
         jy = yaxis(jjpole)
         if (jy.gt.0) jyloc = locrec1(jy)
         jz = zaxis(jjpole)
         if (jz.gt.0) jzloc = locrec1(jz)
         do betac = 1,3
            torq_mu(betac,jjj) = torq_mu(betac,jjj)
     $      + sprod(3,dEnedmu(:,jj),dtorquemu(betac,1,:,jj))

            torq_t(betac,jjj) = torq_t(betac,jjj)
     $      + sprod(9,dEnedt(:,:,jj),
     $      dtorquetheta(betac,1,:,:,jj))

            if (jx .gt. 0) then
               torq_mu(betac,jxloc) = torq_mu(betac,jxloc)
     $       + sprod(3,dEnedmu(:,jj), dtorquemu(betac,2,:,jj))
             torq_t(betac,jxloc) = torq_t(betac,jxloc)
     $       + sprod(9,dEnedt(:,:,jj), 
     $       dtorquetheta(betac,2,:,:,jj))
            end if
            if (jy .gt. 0) then
               torq_mu(betac,jyloc) = torq_mu(betac,jyloc)
     $       + sprod(3,dEnedmu(:,jj), dtorquemu(betac,3,:,jj))
             torq_t(betac,jyloc) = torq_t(betac,jyloc)
     $       + sprod(9,dEnedt(:,:,jj), 
     $       dtorquetheta(betac,3,:,:,jj))
            end if
            if (jz .gt. 0) then
               torq_mu(betac,jzloc) = torq_mu(betac,jzloc)
     $       + sprod(3,dEnedmu(:,jj), dtorquemu(betac,4,:,jj))
             torq_t(betac,jzloc) = torq_t(betac,jzloc)
     $       + sprod(9,dEnedt(:,:,jj), 
     $          dtorquetheta(betac,4,:,:,jj))
            end if
         end do
      end do
c
      deallocate (dtorquemu,dtorquetheta)
      return
      end subroutine
