c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine echarge3  --  charge-charge energy & analysis  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "echarge3" calculates the charge-charge interaction energy
c     and partitions the energy among the atoms
c
c
#include "tinker_precision.h"
      subroutine echarge
      implicit none
c
c     choose the method for summing over pairwise interactions
c
      call echarge0c
c
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine echarge0c  --  Ewald charge analysis via list  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "echarge0c" calculates the charge-charge interaction energy
c     using a particle mesh Ewald summation
c
c
      subroutine echarge0c
      use atmlst
      use atoms
      use bound
      use boxes
      use charge
      use chgpot
      use couple
      use domdec
      use energi
      use ewald
      use group
      use inform
      use inter
      use iounit
      use math
      use molcul
      use neigh
      use potent
      use shunt
      use usage
      use mpi
      use tinheader
      implicit none
      integer i,j,k,iglob,iichg
      integer ii,kk,kkk,inl,kglob,kkchg
      integer in,kn
      real(t_p) e,efull
      real(t_p) f,fi,fik
      real(t_p) fs
      real(t_p) r,r2,rb,rew
      real(t_p) xi,yi,zi
      real(t_p) xr,yr,zr
      real(t_p) xd,yd,zd
      real(t_p) erfterm
      real(t_p) scale,scaleterm
      real(t_p), allocatable :: cscale(:)
      logical proceed,usei
      logical header,huge
      character*10 mode
c
c
c     zero out the Ewald summation energy and partitioning
c
      ec = 0.0_ti_p
c
      if (nion .eq. 0)  return
c
c     compute the reciprocal space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.gt.ndir-1))
     $  then
        call ecrecip
        if (use_pmecore) return
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (cscale(n))
c
c     initialize connected atom exclusion coefficients
c
      cscale = 1.0_ti_p
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'EWALD'
      call switch (mode)
c
c     compute the Ewald self-energy term over all the atoms
c
      fs = -f * aewald / sqrtpi
      do ii = 1, nionloc
         iichg = chgglob(ii)
         iglob = iion(iichg)
         i = loc(iglob)
         e = fs * pchg(iichg)**2
         ec = ec + e
      end do
c
c     compute the cell dipole boundary correction term
c
      if (boundary .eq. 'VACUUM') then
         xd = 0.0_ti_p
         yd = 0.0_ti_p
         zd = 0.0_ti_p
         do ii = 1, nionloc
            iichg = chgglob(ii)
            iglob = iion(iichg)
            i = loc(iglob)
            xd = xd + pchg(iichg)*x(iglob)
            yd = yd + pchg(iichg)*y(iglob)
            zd = zd + pchg(iichg)*z(iglob)
         end do
         e = (2.0_ti_p/3.0_ti_p) * f * (pi/volbox) * (xd*xd+yd*yd+zd*zd)
         ec = ec + e
      end if
c
c     compute the real space portion of the Ewald summation
c
      do ii = 1, nionlocnl
         iichg = chgglobnl(ii)
         iglob = iion(iichg)
         usei = use(iglob)
         xi = x(iglob)
         yi = y(iglob)
         zi = z(iglob)
         fi = f * pchg(iichg)
c
c     set exclusion coefficients for connected atoms
c
         do j = 1, n12(iglob)
            cscale(i12(j,iglob)) = c2scale
         end do
         do j = 1, n13(iglob)
            cscale(i13(j,iglob)) = c3scale
         end do
         do j = 1, n14(iglob)
            cscale(i14(j,iglob)) = c4scale
         end do
         do j = 1, n15(iglob)
            cscale(i15(j,iglob)) = c5scale
         end do
         do kkk = 1, nelst(ii)
            kkchg = elst(kkk,ii)
            kglob = iion(kkchg)
            if (kkchg.eq.0) cycle
c
c     compute the energy contribution for this interaction
c
               xr = xi - x(kglob)
               yr = yi - y(kglob)
               zr = zi - z(kglob)
c
c     find energy for interactions within real space cutoff
c
               call image (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  rb = r + ebuffer
                  fik = fi * pchg(kkchg)
                  rew = aewald * r
                  erfterm = erfc (rew)
                  scale = cscale(kglob)
                  scaleterm = scale - 1.0_ti_p
                  e = (fik/rb) * (erfterm+scaleterm)
                  ec = ec + e
c
c     increment the overall charge-charge energy component
c
                  efull = (fik/rb) * scale
               end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(iglob)
            cscale(i12(j,iglob)) = 1.0_ti_p
         end do
         do j = 1, n13(iglob)
            cscale(i13(j,iglob)) = 1.0_ti_p
         end do
         do j = 1, n14(iglob)
            cscale(i14(j,iglob)) = 1.0_ti_p
         end do
         do j = 1, n15(iglob)
            cscale(i15(j,iglob)) = 1.0_ti_p
         end do
      end do
c
c
c     perform deallocation of some local arrays
c
      deallocate (cscale)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine ecrecip  --  PME reciprocal space charge energy  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "ecrecip" evaluates the reciprocal space portion of the particle
c     mesh Ewald energy due to partial charges
c
c     literature reference:
c
c     U. Essmann, L. Perera, M. L Berkowitz, T. Darden, H. Lee and
c     L. G. Pedersen, "A Smooth Particle Mesh Ewald Method", Journal
c     of Chemical Physics, 103, 8577-8593 (1995)
c
c     modifications for nonperiodic systems suggested by Tom Darden
c     during May 2007
c
c
      subroutine ecrecip
      use atmlst
      use bound
      use boxes
      use charge
      use chgpot
      use domdec
      use energi
      use ewald
      use fft
      use math
      use pme
      use potent
      use mpi
      use tinheader
      implicit none
      integer i,j,k
      integer iichg,iglob
      integer k1,k2,k3
      integer m1,m2,m3
      integer nf1,nf2,nf3
      integer nff,npoint
      integer status(MPI_STATUS_SIZE),tag,ierr,proc
      integer nprocloc,rankloc,commloc
      real(t_p) e,f,denom
      real(t_p) term,expterm
      real(t_p) pterm,volterm
      real(t_p) hsq,struc2
      real(t_p) h1,h2,h3
      real(t_p) r1,r2,r3
      integer, allocatable :: req(:),reqbcast(:)
      real(t_p), allocatable :: qgridmpi(:,:,:,:,:)
c
      if (use_pmecore) then
        nprocloc = nrec
        rankloc = rank_bis
        commloc = comm_rec
      else
        nprocloc = nproc
        rankloc = rank
        commloc = MPI_COMM_WORLD
      end if
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
c
c     dynamic allocation of local arrays
c
      allocate (qgridmpi(2,n1mpimax,n2mpimax,n3mpimax,nrec_recep))
      allocate (req(nproc*nproc))
      allocate (reqbcast(nproc*nproc))
c
      do i = 1, nionrecloc
        iichg = chgrecglob(i)
        iglob = iion(iichg)
        call bspline_fill_site(iglob,i)
      end do
c
      qgridin_2d = 0_ti_p
c
c     MPI : Begin reception
c
      do i = 1, nrec_recep
        tag = nprocloc*rankloc + prec_recep(i) + 1
        call MPI_IRECV(qgridmpi(1,1,1,1,i),2*n1mpimax*n2mpimax*
     $   n3mpimax,MPI_TPREC,prec_recep(i),tag,commloc,req(tag),
     $   ierr)
      end do
      do i = 1, nionrecloc
        iichg = chgrecglob(i)
        iglob = iion(iichg)
        call grid_pchg_site(iglob,i,pchg(iichg))
      end do
c
c     MPI : begin sending
c
      do i = 1, nrec_send
        proc = prec_send(i)
        tag = nprocloc*prec_send(i) + rankloc + 1
        call MPI_ISEND(qgridin_2d(1,1,1,1,i+1),
     $   2*n1mpimax*n2mpimax*n3mpimax,MPI_TPREC,proc,tag,
     $   commloc,req(tag),ierr)
      end do
      do i = 1, nrec_recep
        tag = nprocloc*rankloc + prec_recep(i) + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
      do i = 1, nrec_send
        tag = nprocloc*prec_send(i) + rankloc + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
c
c     do the reduction 'by hand'
c
      do i = 1, nrec_recep
        qgridin_2d(:,:,:,:,1) = qgridin_2d(:,:,:,:,1) + 
     $    qgridmpi(:,:,:,:,i) 
      end do
      print*, sum(abs(qgridin_2d))
c
c     perform the 3-D FFT forward transformation
c
      call fft2d_frontmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
     $ n3mpimax)
c
c     use scalar sum to get the reciprocal space energy
c
      if ((istart2(rankloc+1).eq.1).and.(jstart2(rankloc+1).eq.1).and.
     $   (kstart2(rankloc+1).eq.1)) then
           qfac_2d(1,1,1) = 0.0_ti_p
      end if
      f = 0.5_ti_p * electric / dielec
      pterm = (pi/aewald)**2
      volterm = pi * volbox
      nff = nfft1 * nfft2
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
            if (term .gt. -50.0_ti_p) then
               denom = volterm*hsq*bsmod1(k1)*bsmod2(k2)*bsmod3(k3)
               expterm = exp(term) / denom
               if (.not. use_bounds) then
                  expterm = expterm * (1.0_ti_p-cos(pi*xbox*sqrt(hsq)))
               else if (octahedron) then
                  if (mod(m1+m2+m3,2) .ne. 0)  expterm = 0.0_ti_p
               end if
               struc2 = qgridout_2d(1,k1-istart2(rankloc+1)+1,
     $  k2-jstart2(rankloc+1)+1,k3-kstart2(rankloc+1)+1)**2 + 
     $  qgridout_2d(2,k1-istart2(rankloc+1)+1,k2-jstart2(rankloc+1)+1,
     $  k3-kstart2(rankloc+1)+1)**2
               e = f * expterm * struc2
               ec = ec + e
            end if
 10         continue
          end do
        end do
      end do
c
c     account for zeroth grid point for nonperiodic system
c
      if ((istart2(rankloc+1).eq.1).and.(jstart2(rankloc+1).eq.1).and.
     $   (kstart2(rankloc+1).eq.1)) then
        if (.not. use_bounds) then
           expterm = 0.5_ti_p * pi / xbox
           struc2 = qgridout_2d(1,1,1,1)**2 + qgridout_2d(2,1,1,1)**2
           e = f * expterm * struc2
           ec = ec + e
        end if
      end if
      deallocate (qgridmpi)
      deallocate (req)
      deallocate (reqbcast)
      return
      end
