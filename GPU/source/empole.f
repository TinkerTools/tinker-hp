c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine empole3  --  mpole/polar energy & analysis  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "empole3" calculates the electrostatic energy due to
c     atomic multipole and dipole polarizability interactions
c
c
#include "tinker_macro.h"
      subroutine empole
      implicit none
c
c     choose the method for summing over multipole interactions
c
      call empolec
c
      return
      end
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine empole0c  --  Ewald multipole energy via list  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "empole0d" calculates the atomic multipole interaction energy
c     using particle mesh Ewald summation and a neighbor list
c
c
      subroutine empolec
      use sizes
      use atmlst
      use atoms
      use boxes
      use chgpot
      use domdec
      use energi
      use ewald
      use math
      use mpole
      use potent
      use tinheader ,only:ti_p,re_p
      use mpi
      implicit none
      integer i,ii,iglob,iipole,ierr
      real(t_p) e,f
      real(t_p) term,fterm
      real(t_p) cii,dii,qii
      real(t_p) xd,yd,zd
      real(t_p) ci,dix,diy,diz
      real(t_p) qixx,qixy,qixz
      real(t_p) qiyy,qiyz,qizz
c
c
c     zero out the total atomic multipole energy
c
      em = 0.0_re_p
      if (npole .eq. 0)  return
c
c     set Ewald coefficient
c
      aewald = aeewald
c
c     set the energy unit conversion factor
c
      f = electric / dielec
c
c     check the sign of multipole components at chiral sites
c
      call chkpole(.false.)
c
c     rotate the multipole components into the global frame
c
      call rotpole
c
c     compute the reciprocal space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.gt.ndir-1))
     $  then
        if (use_mrec) then
          call emrecip
        end if
      end if
c
c     compute the real space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.le.ndir-1))
     $   then
        if (use_mreal) then
           call emreal0d
        end if

        if (use_mself) then
c
c     compute the self-energy portion of the Ewald summation
c
        term = 2.0_ti_p * aewald * aewald
        fterm = -f * aewald / sqrtpi
        do ii = 1, npoleloc
           iipole = poleglob(ii)
           iglob = ipole(iipole)
           i = loc(iglob)
           ci = rpole(1,iipole)
           dix = rpole(2,iipole)
           diy = rpole(3,iipole)
           diz = rpole(4,iipole)
           qixx = rpole(5,iipole)
           qixy = rpole(6,iipole)
           qixz = rpole(7,iipole)
           qiyy = rpole(9,iipole)
           qiyz = rpole(10,iipole)
           qizz = rpole(13,iipole)
           cii = ci*ci
           dii = dix*dix + diy*diy + diz*diz
           qii = 2.0_ti_p*(qixy*qixy+qixz*qixz+qiyz*qiyz)
     &              + qixx*qixx + qiyy*qiyy + qizz*qizz
           e = fterm * (cii + term*(dii/3.0_ti_p+
     &                              2.0_ti_p*term*qii/5.0_ti_p))
           em = em + e
        end do
c
c       compute the cell dipole boundary correction term
c
        if (boundary .eq. 'VACUUM') then
           xd = 0.0_ti_p
           yd = 0.0_ti_p
           zd = 0.0_ti_p
           do ii = 1, npoleloc
              iipole = poleglob(ii)
              iglob = ipole(iipole)
              dix = rpole(2,iipole)
              diy = rpole(3,iipole)
              diz = rpole(4,iipole)
              xd = xd + dix + rpole(1,iipole)*x(iglob)
              yd = yd + diy + rpole(1,iipole)*y(iglob)
              zd = zd + diz + rpole(1,iipole)*z(iglob)
           end do
             call MPI_ALLREDUCE(MPI_IN_PLACE,xd,1,MPI_TPREC,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,yd,1,MPI_TPREC,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,zd,1,MPI_TPREC,MPI_SUM,
     $          COMM_TINKER,ierr)
             if (rank.eq.0) then
           term = (2.0_ti_p/3.0_ti_p) * f * (pi/volbox)
           e = term * (xd*xd+yd*yd+zd*zd)
           em = em + e
             end if
          end if
        end if
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine emreal0d  --  real space mpole energy via list  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "emreal0d" evaluates the real space portion of the Ewald sum
c     energy due to atomic multipoles using a neighbor list
c
c     literature reference:
c
c     W. Smith, "Point Multipoles in the Ewald Summation (Revisited)",
c     CCP5 Newsletter, 46, 18-30, 1998  (see http://www.ccp5.org/)
c
c
      subroutine emreal0d
      use sizes
      use atmlst
      use atoms
      use bound
      use chgpot
      use couple
      use cutoff
      use domdec
      use energi
      use ewald
      use group
      use math
      use mpole
      use mplpot
      use neigh
      use shunt
      use tinheader ,only:ti_p,re_p
      use mpi
      integer i,j,iglob,kglob,nnelst
      integer ii,kkk,iipole,kkpole
      real(t_p) e,f
      real(t_p) scalek
      real(t_p) xi,yi,zi
      real(t_p) xr,yr,zr
      real(t_p) r,r2,rr1,rr3
      real(t_p) rr5,rr7,rr9
      real(t_p) rr1i,rr3i,rr5i
      real(t_p) rr1k,rr3k,rr5k
      real(t_p) rr1ik,rr3ik,rr5ik
      real(t_p) rr7ik,rr9ik
      real(t_p) ci,dix,diy,diz
      real(t_p) qixx,qixy,qixz
      real(t_p) qiyy,qiyz,qizz
      real(t_p) ck,dkx,dky,dkz
      real(t_p) qkxx,qkxy,qkxz
      real(t_p) qkyy,qkyz,qkzz
      real(t_p) dir,dkr,dik,qik
      real(t_p) qix,qiy,qiz,qir
      real(t_p) qkx,qky,qkz,qkr
      real(t_p) diqk,dkqi,qiqk
      real(t_p) corei,corek
      real(t_p) vali,valk
      real(t_p) alphai,alphak
      real(t_p) term1,term2,term3
      real(t_p) term4,term5
      real(t_p) term1i,term2i,term3i
      real(t_p) term1k,term2k,term3k
      real(t_p) term1ik,term2ik,term3ik
      real(t_p) term4ik,term5ik
      real(t_p) dmpi(9),dmpk(9)
      real(t_p) dmpik(9),dmpe(9)
      real(t_p) fgrp
      real(t_p) s,ds,mpoleshortcut2
      real(t_p) facts
      logical testcut,shortrange,longrange,fullrange
      real(t_p), allocatable :: mscale(:)
      character*11 mode
      character*80 :: RoutineName


c     compute the short, long, or full real space part of the summation
      shortrange = use_mpoleshortreal
      longrange  = use_mpolelong
      fullrange  = .not.(shortrange.or.longrange)
      if (shortrange) then 
         RoutineName = 'emrealshort3d'
         mode        = 'SHORTEWALD'
      else if (longrange) then
         RoutineName = 'emreallong3d'
         mode        = 'EWALD'
      else
         RoutineName = 'emreal3d'
         mode        = 'EWALD'
      endif

c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
c
c     initialize connected atom exclusion coefficients
c
      mscale = 1.0
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      call switch (mode)
c
c     compute the real space portion of the Ewald summation
c
      do ii = 1, npolelocnl
         iipole = poleglobnl(ii)
         iglob  = ipole(iipole)
         i      = loc  (iglob)
         xi = x(iglob)
         yi = y(iglob)
         zi = z(iglob)
         ci = rpole(1,iipole)
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
         do j = 1, n12(iglob)
            mscale(i12(j,iglob)) = m2scale
         end do
         do j = 1, n13(iglob)
            mscale(i13(j,iglob)) = m3scale
         end do
         do j = 1, n14(iglob)
            mscale(i14(j,iglob)) = m4scale
         end do
         do j = 1, n15(iglob)
            mscale(i15(j,iglob)) = m5scale
         end do
c
c     evaluate all sites within the cutoff distance
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
            kglob = ipole(kkpole)
            if (use_group)  call groups (fgrp,iglob,kglob,0,0,0,0)
            xr = x(kglob) - xi
            yr = y(kglob) - yi
            zr = z(kglob) - zi
            if (use_bounds)  call image (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            testcut = merge(r2 .le. off2.and.r2.ge.mpoleshortcut2,
     &                      r2 .le. off2,
     &                      longrange
     &                     )
            if (testcut) then
               r = sqrt(r2)
               ck = rpole(1,kkpole)
               dkx = rpole(2,kkpole)
               dky = rpole(3,kkpole)
               dkz = rpole(4,kkpole)
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
               dik = dix*dkx + diy*dky + diz*dkz
               qik = qix*qkx + qiy*qky + qiz*qkz
               diqk = dix*qkx + diy*qky + diz*qkz
               dkqi = dkx*qix + dky*qiy + dkz*qiz
               qiqk = 2.0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                   + qixx*qkxx + qiyy*qkyy + qizz*qkzz
c
c     get reciprocal distance terms for this interaction
c
               rr1 = f / r
               rr3 = rr1 / r2
               rr5 = 3.0 * rr3 / r2
               rr7 = 5.0 * rr5 / r2
               rr9 = 7.0 * rr7 / r2
c
c     calculate real space Ewald error function damping
c
               call dampewald (9,r,r2,f,dmpe)
c
c     find damped multipole intermediates and energy value
c
               if (use_chgpen) then
                  corek = pcore(kkpole)
                  valk = pval(kkpole)
                  alphak = palpha(kkpole)
                  term1 = corei*corek
                  term1i = corek*vali
                  term2i = corek*dir
                  term3i = corek*qir
                  term1k = corei*valk
                  term2k = -corei*dkr
                  term3k = corei*qkr
                  term1ik = vali*valk
                  term2ik = valk*dir - vali*dkr + dik
                  term3ik = vali*qkr + valk*qir - dir*dkr
     &                         + 2.0*(dkqi-diqk+qiqk)
                  term4ik = dir*qkr - dkr*qir - 4.0*qik
                  term5ik = qir*qkr
                  call damppole (r,9,alphai,alphak,
     &                            dmpi,dmpk,dmpik)
                  scalek = mscale(kglob)
                  if (use_group)  scalek = scalek * fgrp
                  rr1i = dmpe(1) - (1.0-scalek*dmpi(1))*rr1
                  rr3i = dmpe(3) - (1.0-scalek*dmpi(3))*rr3
                  rr5i = dmpe(5) - (1.0-scalek*dmpi(5))*rr5
                  rr1k = dmpe(1) - (1.0-scalek*dmpk(1))*rr1
                  rr3k = dmpe(3) - (1.0-scalek*dmpk(3))*rr3
                  rr5k = dmpe(5) - (1.0-scalek*dmpk(5))*rr5
                  rr1ik = dmpe(1) - (1.0-scalek*dmpik(1))*rr1
                  rr3ik = dmpe(3) - (1.0-scalek*dmpik(3))*rr3
                  rr5ik = dmpe(5) - (1.0-scalek*dmpik(5))*rr5
                  rr7ik = dmpe(7) - (1.0-scalek*dmpik(7))*rr7
                  rr9ik = dmpe(9) - (1.0-scalek*dmpik(9))*rr9
                  rr1 = dmpe(1) - (1.0-scalek)*rr1
                  e = term1*rr1 + term4ik*rr7ik + term5ik*rr9ik
     &                   + term1i*rr1i + term1k*rr1k + term1ik*rr1ik
     &                   + term2i*rr3i + term2k*rr3k + term2ik*rr3ik
     &                   + term3i*rr5i + term3k*rr5k + term3ik*rr5ik
c
c     find standard multipole intermediates and energy value
c
               else
                  term1 = ci*ck
                  term2 = ck*dir - ci*dkr + dik
                  term3 = ci*qkr + ck*qir - dir*dkr
     &                       + 2.0*(dkqi-diqk+qiqk)
                  term4 = dir*qkr - dkr*qir - 4.0*qik
                  term5 = qir*qkr
                  scalek = 1.0 - mscale(kglob)
                  if (use_group)  then
                    scalek = 1.0 - mscale(kglob)*fgrp
                  end if
                  rr1 = dmpe(1) - scalek*rr1
                  rr3 = dmpe(3) - scalek*rr3
                  rr5 = dmpe(5) - scalek*rr5
                  rr7 = dmpe(7) - scalek*rr7
                  rr9 = dmpe(9) - scalek*rr9
                  e = term1*rr1 + term2*rr3 + term3*rr5
     &                   + term4*rr7 + term5*rr9
               end if

               if(shortrange .or. longrange)
     &            call switch_respa(r,ewaldshortcut,shortheal,s,ds)
c
c     fix the s factor, depending on the range
c
               if(shortrange) then
                  facts  =         s
               else if(longrange) then
                  facts  = 1.0 - s
               else
                  facts  = 1.0
               endif

               em = em + facts * e
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(iglob)
            mscale(i12(j,iglob)) = 1.0
         end do
         do j = 1, n13(iglob)
            mscale(i13(j,iglob)) = 1.0
         end do
         do j = 1, n14(iglob)
            mscale(i14(j,iglob)) = 1.0
         end do
         do j = 1, n15(iglob)
            mscale(i15(j,iglob)) = 1.0
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (mscale)
      return
      end
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine emrecip  --  PME recip space multipole energy  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "emrecip" evaluates the reciprocal space portion of the particle
c     mesh Ewald energy due to atomic multipole interactions
c
c     literature reference:
c
c     C. Sagui, L. G. Pedersen and T. A. Darden, "Towards an Accurate
c     Representation of Electrostatics in Classical Force Fields:
c     Efficient Implementation of Multipolar Interactions in
c     Biomolecular Simulations", Journal of Chemical Physics, 120,
c     73-87 (2004)
c
c     modifications for nonperiodic systems suggested by Tom Darden
c     during May 2007
c
c
      subroutine emrecip
      use atmlst
      use bound
      use boxes
      use chgpot
      use domdec
      use energi
      use ewald
      use fft
      use math
      use mpole
      use pme
      use potent
      use tinheader ,only:ti_p,re_p
      use mpi
      implicit none
      integer ierr,iipole,proc,iglob
      integer status(MPI_STATUS_SIZE),tag,commloc
      integer i,j,k
      integer k1,k2,k3
      integer m1,m2,m3
      integer nff
      integer nf1,nf2,nf3
      integer nprocloc,rankloc
      real(t_p) e,r1,r2,r3
      real(t_p) f,h1,h2,h3
      real(t_p) volterm,denom
      real(t_p) hsq,expterm
      real(t_p) term,pterm
      real(t_p) struc2
      real(t_p) cmp(10)
      real(t_p), allocatable ::  fmp(:,:)
      integer, allocatable :: reqsend(:),reqrec(:)
      integer, allocatable :: req2send(:),req2rec(:)
      real(t_p), allocatable :: qgridmpi(:,:,:,:,:)
c
      if (rank==0.and.tinkerdebug) write(*,*) ' emrecip'
      if (use_pmecore) then
        nprocloc = nrec
        rankloc = rank_bis
        commloc =  comm_rec
      else
        nprocloc = nproc
        rankloc = rank
        commloc = COMM_TINKER
      end if
c
c     dynamic allocation of global arrays
c
      if (allocated(fphirec)) deallocate(fphirec)
      allocate (fphirec(20,max(npolerecloc,1)))
      fphirec = 0_ti_p
c
c     dynamic allocation of local arrays
c
      allocate (fmp(20,npolerecloc))
      allocate (qgridmpi(2,n1mpimax,n2mpimax,n3mpimax,nrec_recep))
      allocate (reqsend(nproc))
      allocate (reqrec(nproc))
      allocate (req2send(nproc))
      allocate (req2rec(nproc))
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
      f = electric / dielec
c
c     copy the multipole moments into local storage areas
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
         cmp(8) = 2.0_ti_p * rpole(6,iipole)
         cmp(9) = 2.0_ti_p * rpole(7,iipole)
         cmp(10) = 2.0_ti_p * rpole(10,iipole)
         call cmp_to_fmp_site(cmp,fmp(1,i))
         call bspline_fill_site(iglob,i)
      end do
c
c     zero out the PME grid
c
      qgridin_2d = 0_ti_p
c
c     assign PME grid and perform 3-D FFT forward transform
c
      do i = 1, npolerecloc
        iipole = polerecglob(i)
        iglob = ipole(iipole)
        call grid_mpole_site(iglob,i,fmp(1,i))
      end do
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
c     MPI : begin sending
c
      do i = 1, nrec_send
        proc = prec_send(i)
        tag = nprocloc*prec_send(i) + rankloc + 1
        call MPI_ISEND(qgridin_2d(1,1,1,1,i+1),
     $   2*n1mpimax*n2mpimax*n3mpimax,MPI_TPREC,proc,tag,
     $   commloc,reqsend(i),ierr)
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

      do i = 1, nrec_recep
        qgridin_2d(:,:,:,:,1) = qgridin_2d(:,:,:,:,1) + 
     $    qgridmpi(:,:,:,:,i) 
      end do
c
c     Perform 3-D FFT forward transform
c
      call fft2d_frontmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
     $ n3mpimax)
c
c     make the scalar summation over reciprocal lattice
c
      pterm = (pi/aewald)**2
      volterm = pi * volbox
      nff = nfft1 * nfft2
      nf1 = (nfft1+1) / 2
      nf2 = (nfft2+1) / 2
      nf3 = (nfft3+1) / 2
c
      if ((istart2(rankloc+1).eq.1).and.(jstart2(rankloc+1).eq.1).and.
     $   (kstart2(rankloc+1).eq.1)) then
           qfac_2d(1,1,1) = 0.0_ti_p
      end if
c
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
            if (term .gt. -50.0_ti_p) then
               denom = volterm*hsq*bsmod1(k1)*bsmod2(k2)*bsmod3(k3)
               expterm = exp(term) / denom
               if (.not. use_bounds) then
                  expterm = expterm * (1.0_ti_p-cos(pi*xbox*sqrt(hsq)))
               else if (octahedron) then
                  if (mod(m1+m2+m3,2) .ne. 0)  expterm = 0.0_ti_p
               end if
            qfac_2d(k1-istart2(rankloc+1)+1,k2-jstart2(rankloc+1)+1,k3-
     $       kstart2(rankloc+1)+1) = expterm
            end if
 10         continue
          end do
        end do
      end do
c
c     account for zeroth grid point for nonperiodic system
c
      if ((istart2(rankloc+1).eq.1).and.(jstart2(rankloc+1).eq.1)
     $   .and.(kstart2(rankloc+1).eq.1)) then
        if (.not. use_bounds) then
           expterm = 0.5_ti_p * pi / xbox
           struc2 = qgrid2in_2d(1,1,1,1,1)**2 +
     $       qgrid2in_2d(2,1,1,1,1)**2
           e = f * expterm * struc2
           em = em + e
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
c     perform 3-D FFT backward transform and get potential
c
      call fft2d_backmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
     $ n3mpimax)
c
c     MPI : Begin reception
c
      do i = 1, nrec_send
        proc = prec_send(i)
        tag = nprocloc*rankloc + prec_send(i) + 1
        call MPI_IRECV(qgridin_2d(1,1,1,1,i+1),
     $   2*n1mpimax*n2mpimax*n3mpimax,MPI_TPREC,
     $   prec_send(i),tag,commloc,req2rec(i),ierr)
      end do
c
c     MPI : begin sending
c
      do i = 1, nrec_recep
        tag = nprocloc*prec_recep(i) + rankloc + 1
        call MPI_ISEND(qgridin_2d,
     $   2*n1mpimax*n2mpimax*n3mpimax,
     $   MPI_TPREC,prec_recep(i),tag,commloc,req2send(i),
     $   ierr)
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
        call fphi_mpole_site(iglob,i)
      end do
c
c     sum over multipoles and increment total multipole energy
c
      e = 0.0_ti_p
      do i = 1, npolerecloc
         do k = 1, 10
            e = e + fmp(k,i)*fphirec(k,i)
         end do
      end do
      e = 0.5_ti_p * f * e
      em = em + e
c
      deallocate (reqsend)
      deallocate (reqrec)
      deallocate (req2send)
      deallocate (req2rec)
      return
      end
