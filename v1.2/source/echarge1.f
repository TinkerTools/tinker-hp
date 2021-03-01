c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine echarge1  --  charge-charge energy & derivs  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "echarge1" calculates the charge-charge interaction energy
c     and first derivatives with respect to Cartesian coordinates
c
c
      subroutine echarge1
      implicit none
c
c     choose the method for summing over pairwise interactions
c
      call echarge1c
c
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine echarge1c  --  Ewald charge derivs via list  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "echarge1c" calculates the charge-charge interaction energy
c     and first derivatives with respect to Cartesian coordinates
c     using a particle mesh Ewald summation and a neighbor list
c
c
      subroutine echarge1c
      use atmlst
      use atoms
      use bound
      use boxes
      use charge
      use chgpot
      use deriv
      use energi
      use ewald
      use domdec
      use iounit
      use inter
      use math
      use potent
      use timestat
      use usage
      use virial
      use mpi
      implicit none
      integer ii,i,iglob,iichg
      real*8 e,de,term
      real*8 f,fs
      real*8 xd,yd,zd
      real*8 dedx,dedy,dedz
      real*8 time0,time1
c
c
c     zero out the Ewald summation energy and derivatives
c
      ec = 0.0d0
      dec = 0.0d0
      if (nion .eq. 0)  return
c
c     compute the reciprocal space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.gt.ndir-1))
     $  then
        time0 = mpi_wtime()
        if (use_crec) then
          call ecrecip1
        end if
        time1 = mpi_wtime()
        timerec = timerec + time1 - time0
        if (use_pmecore) return
      end if
c
      if (use_cself) then
c
c     compute the Ewald self-energy term over all the atoms
c
        f = electric / dielec
        fs = -f * aewald / sqrtpi
        do ii = 1, nionloc
          iichg = chgglob(ii)
          e = fs * pchg(iichg)**2
          ec = ec + e
        end do
c
c     compute the cell dipole boundary correction term
c
        if (boundary .eq. 'VACUUM') then
           xd = 0.0d0
           yd = 0.0d0
           zd = 0.0d0
           do ii = 1, nionloc
             iichg = chgglob(ii)
             iglob = iion(iichg)
             i = loc(iglob)
             xd = xd + pchg(iichg)*x(iglob)
             yd = yd + pchg(iichg)*y(iglob)
             zd = zd + pchg(iichg)*z(iglob)
           end do
           term = (2.0d0/3.0d0) * f * (pi/volbox)
           e = term * (xd*xd+yd*yd+zd*zd)
           ec = ec + e
           do ii = 1, nionloc
              iichg = chgglob(ii)
              iglob = iion(iichg)
              i = loc(iglob)
              de = 2.0d0 * term * pchg(iichg)
              dedx = de * xd
              dedy = de * yd
              dedz = de * zd
              dec(1,i) = dec(1,i) + dedx
              dec(2,i) = dec(2,i) + dedy
              dec(3,i) = dec(3,i) + dedz
           end do
        end if
      end if
c
c     compute the real space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.le.ndir-1))
     $   then
        if (use_creal) then
          time0 = mpi_wtime()
          if (use_cshortreal) then
            call ecrealshort1d
          else if (use_clong) then
            call ecreallong1d
          else
            call ecreal1d
          end if
          time1 = mpi_wtime()
          timereal = timereal + time1 - time0
        end if
      end if
      return
      end


c
c     "ecreal1d" evaluates the real space portion of the Ewald sum
c     energy and forces due to atomic charge interactions, using a neighbor list
c
      subroutine ecreal1d
      use atmlst
      use atoms
      use bound
      use boxes
      use charge
      use chgpot
      use couple
      use deriv
      use domdec
      use energi
      use ewald
      use iounit
      use inter
      use math
      use molcul
      use neigh
      use potent
      use shunt
      use timestat
      use usage
      use virial
      use mpi
      implicit none
      integer i,j,k,iichg,iglob,kglob,kkchg
      integer ii,kkk

      real*8 e,de,efull
      real*8 f,fi,fik
      real*8 r,r2,rew
      real*8 rb,rb2

      real*8 xi,yi,zi
      real*8 xr,yr,zr

      real*8 erfc,erfterm
      real*8 scale,scaleterm
      real*8 dedx,dedy,dedz
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8, allocatable :: cscale(:)


      character*10 mode
      external erfc
 1000 format(' Warning, system moved too much since last neighbor list'
     $ ' update, try lowering nlupdate')
c
c     perform dynamic allocation of some local arrays
c
      allocate (cscale(n))
c
c     initialize connected atom scaling 
c
      cscale = 1.0d0
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'EWALD'
      call switch (mode)
c
c     compute the real space Ewald energy and first derivatives
c
      do ii = 1, nionlocnl
         iichg = chgglobnl(ii)
         iglob = iion(iichg)
         i = loc(iglob)
         if (i.eq.0) then
           write(iout,1000)
           cycle
         end if
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
            if (kkchg.eq.0) cycle
            kglob = iion(kkchg)
            k = loc(kglob)
            if (k.eq.0) then
              write(iout,1000)
              cycle
            end if
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
               rb2 = rb * rb
               fik = fi * pchg(kkchg)
               rew = aewald * r
               erfterm = erfc (rew)
               scale = cscale(kglob)
               scaleterm = scale - 1.0d0
               e = (fik/rb) * (erfterm+scaleterm)
               de = -fik * ((erfterm+scaleterm)/rb2
     &                 + (2.0d0*aewald/sqrtpi)*exp(-rew**2)/rb)
c
c     form the chain rule terms for derivative expressions
c
               de = de / r
               dedx = de * xr
               dedy = de * yr
               dedz = de * zr
c
c     increment the overall energy and derivative expressions
c
               ec = ec + e
               dec(1,i) = dec(1,i) + dedx
               dec(2,i) = dec(2,i) + dedy
               dec(3,i) = dec(3,i) + dedz
               dec(1,k) = dec(1,k) - dedx
               dec(2,k) = dec(2,k) - dedy
               dec(3,k) = dec(3,k) - dedz
c
c     increment the internal virial tensor components
c
               vxx = xr * dedx
               vyx = yr * dedx
               vzx = zr * dedx
               vyy = yr * dedy
               vzy = zr * dedy
               vzz = zr * dedz
               vir(1,1) = vir(1,1) + vxx
               vir(2,1) = vir(2,1) + vyx
               vir(3,1) = vir(3,1) + vzx
               vir(1,2) = vir(1,2) + vyx
               vir(2,2) = vir(2,2) + vyy
               vir(3,2) = vir(3,2) + vzy
               vir(1,3) = vir(1,3) + vzx
               vir(2,3) = vir(2,3) + vzy
               vir(3,3) = vir(3,3) + vzz
c
c     increment the total intramolecular energy
c
               if (molcule(iglob) .ne. molcule(kglob)) then
                  efull = (fik/rb) * scale
                  einter = einter + efull
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(iglob)
            cscale(i12(j,iglob)) = 1.0d0
         end do
         do j = 1, n13(iglob)
            cscale(i13(j,iglob)) = 1.0d0
         end do
         do j = 1, n14(iglob)
            cscale(i14(j,iglob)) = 1.0d0
         end do
         do j = 1, n15(iglob)
            cscale(i15(j,iglob)) = 1.0d0
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (cscale)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine ecrecip1  --  PME recip charge energy & derivs  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "ecrecip1" evaluates the reciprocal space portion of the particle
c     mesh Ewald summation energy and gradient due to partial charges
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
      subroutine ecrecip1
      use atmlst
      use atoms
      use bound
      use boxes
      use charge
      use chgpot
      use deriv
      use domdec
      use energi
      use ewald
      use fft
      use math
      use pme
      use potent
      use timestat
      use virial
      use mpi
      implicit none
      integer i,j,k
      integer iichg,iglob
      integer isite,iatm
      integer i0,j0,k0
      integer it1,it2,it3
      integer k1,k2,k3
      integer m1,m2,m3
      integer nf1,nf2,nf3
      integer nff
      integer status(MPI_STATUS_SIZE),tag,ierr
      integer nprocloc,rankloc,commloc,proc
      integer istart,iend,jstart,jend,kstart,kend
      integer iloc,iproc
      integer igrd0,jgrd0,kgrd0
      real*8 e,term,expterm
      real*8 vterm,pterm
      real*8 volterm
      real*8 f,fi,denom
      real*8 hsq,struc2
      real*8 de1,de2,de3
      real*8 dn1,dn2,dn3
      real*8 t1,t2,t3
      real*8 dt1,dt2,dt3
      real*8 h1,h2,h3
      real*8 r1,r2,r3
      integer, allocatable :: req(:),reqbcast(:)
      real*8, allocatable :: qgridmpi(:,:,:,:,:)
      real*8 time0,time1
c
      if (use_pmecore) then
        nprocloc = nrec
        commloc  = comm_rec
        rankloc = rank_bis
      else
        nprocloc = nproc
        commloc  = COMM_TINKER
        rankloc = rank
      end if
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
c
c     dynamic allocation of local arrays
c
      allocate (qgridmpi(2,n1mpimax,n2mpimax,n3mpimax,nrec_recep))
      allocate (req(nprocloc*nprocloc))
      allocate (reqbcast(nprocloc*nprocloc))
c
      do i = 1, nionrecloc
        iichg = chgrecglob(i)
        iglob = iion(iichg)
        call bspline_fill_site(iglob,i)
      end do
c
      qgridin_2d = 0d0
c
c     MPI : Begin reception
c
      do i = 1, nrec_recep
        tag = nprocloc*rankloc + prec_recep(i) + 1
        call MPI_IRECV(qgridmpi(1,1,1,1,i),2*n1mpimax*n2mpimax*
     $   n3mpimax,MPI_REAL8,prec_recep(i),tag,commloc,req(tag),
     $   ierr)
      end do
      time0 = mpi_wtime()
      do i = 1, nionrecloc
        iichg = chgrecglob(i)
        iglob = iion(iichg)
        call grid_pchg_site(iglob,i,pchg(iichg))
      end do
      time1 = mpi_wtime()
      timegrid1 = timegrid1 + time1-time0
c
c     MPI : begin sending
c
      time0 = mpi_wtime()
      do i = 1, nrec_send
        proc = prec_send(i)
        tag = nprocloc*prec_send(i) + rankloc + 1
        call MPI_ISEND(qgridin_2d(1,1,1,1,i+1),
     $   2*n1mpimax*n2mpimax*n3mpimax,MPI_REAL8,proc,tag,
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
      time1 = mpi_wtime()
      timerecreccomm = timerecreccomm + time1-time0
c
c     perform the 3-D FFT forward transformation
c
      time0 = mpi_wtime()
      call fft2d_frontmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
     $ n3mpimax)
      time1 = mpi_wtime()
      timeffts = timeffts + time1-time0
c
c     use scalar sum to get reciprocal space energy and virial
c
      time0 = mpi_wtime()
      if ((istart2(rankloc+1).eq.1).and.(jstart2(rankloc+1).eq.1).and.
     $   (kstart2(rankloc+1).eq.1)) then
           qfac_2d(1,1,1) = 0.0d0
      end if
      f = 0.5d0 * electric / dielec
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
            r1 = dble(m1)
            r2 = dble(m2)
            r3 = dble(m3)
            h1 = recip(1,1)*r1 + recip(1,2)*r2 + recip(1,3)*r3
            h2 = recip(2,1)*r1 + recip(2,2)*r2 + recip(2,3)*r3
            h3 = recip(3,1)*r1 + recip(3,2)*r2 + recip(3,3)*r3
            hsq = h1*h1 + h2*h2 + h3*h3
            term = -pterm * hsq
            expterm = 0.0d0
            if (term .gt. -50.0d0) then
               denom = volterm*hsq*bsmod1(k1)*bsmod2(k2)*bsmod3(k3)
               expterm = exp(term) / denom
               if (.not. use_bounds) then
                  expterm = expterm * (1.0d0-cos(pi*xbox*sqrt(hsq)))
               else if (octahedron) then
                  if (mod(m1+m2+m3,2) .ne. 0)  expterm = 0.0d0
               end if
               struc2 = qgridout_2d(1,k1-istart2(rankloc+1)+1,
     $  k2-jstart2(rankloc+1)+1,k3-kstart2(rankloc+1)+1)**2 + 
     $  qgridout_2d(2,k1-istart2(rankloc+1)+1,k2-jstart2(rankloc+1)+1,
     $  k3-kstart2(rankloc+1)+1)**2
               e = f * expterm * struc2
               ec = ec + e
               vterm = (2.0d0/hsq) * (1.0d0-term) * e
               vir(1,1) = vir(1,1) + h1*h1*vterm - e
               vir(2,1) = vir(2,1) + h1*h2*vterm
               vir(3,1) = vir(3,1) + h1*h3*vterm
               vir(1,2) = vir(1,2) + h2*h1*vterm
               vir(2,2) = vir(2,2) + h2*h2*vterm - e
               vir(3,2) = vir(3,2) + h2*h3*vterm
               vir(1,3) = vir(1,3) + h3*h1*vterm
               vir(2,3) = vir(2,3) + h3*h2*vterm
               vir(3,3) = vir(3,3) + h3*h3*vterm - e
            end if
            qgridout_2d(1,k1-istart2(rankloc+1)+1,
     $       k2-jstart2(rankloc+1)+1,k3-kstart2(rankloc+1)+1) = 
     $       expterm * qgridout_2d(1,k1-istart2(rankloc+1)+1,
     $       k2-jstart2(rankloc+1)+1,k3-kstart2(rankloc+1)+1)
            qgridout_2d(2,k1-istart2(rankloc+1)+1,
     $       k2-jstart2(rankloc+1)+1,k3-kstart2(rankloc+1)+1) = 
     $       expterm * qgridout_2d(2,k1-istart2(rankloc+1)+1,
     $       k2-jstart2(rankloc+1)+1,k3-kstart2(rankloc+1)+1)
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
           expterm = 0.5d0 * pi / xbox
           struc2 = qgridout_2d(1,1,1,1)**2 + qgridout_2d(2,1,1,1)**2
           e = f * expterm * struc2
           ec = ec + e
           qgridout_2d(1,1,1,1) = expterm * qgridout_2d(1,1,1,1)
           qgridout_2d(2,1,1,1) = expterm * qgridout_2d(2,1,1,1)
        end if
      end if
      time1 = mpi_wtime()
      timescalar = timescalar + time1-time0
c
c     perform the 3-D FFT backward transformation
c
      time0 = mpi_wtime()
      call fft2d_backmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
     $ n3mpimax)
      time1 = mpi_wtime()
      timeffts = timeffts + time1-time0
c
      time0 = mpi_wtime()
c     MPI : Begin reception
c
      do i = 1, nrec_send
        proc = prec_send(i)
        tag = nprocloc*rankloc + prec_send(i) + 1
        call MPI_IRECV(qgridin_2d(1,1,1,1,i+1),
     $   2*n1mpimax*n2mpimax*n3mpimax,MPI_REAL8,
     $   prec_send(i),tag,commloc,reqbcast(tag),ierr)
      end do
c
c     MPI : begin sending
c
      do i = 1, nrec_recep
        tag = nprocloc*prec_recep(i) + rankloc + 1
        call MPI_ISEND(qgridin_2d,
     $   2*n1mpimax*n2mpimax*n3mpimax,
     $   MPI_REAL8,prec_recep(i),tag,commloc,reqbcast(tag),
     $   ierr)
      end do
c
      do i = 1, nrec_send
        tag = nprocloc*rankloc + prec_send(i) + 1
        call MPI_WAIT(reqbcast(tag),status,ierr)
      end do
      do i = 1, nrec_recep
        tag = nprocloc*prec_recep(i) + rankloc + 1
        call MPI_WAIT(reqbcast(tag),status,ierr)
      end do
      time1 = mpi_wtime()
      timerecreccomm = timerecreccomm + time1-time0
c
c     get first derivatives of the reciprocal space energy
c
      time0 = mpi_wtime()
      f = electric / dielec
      dn1 = dble(nfft1)
      dn2 = dble(nfft2)
      dn3 = dble(nfft3)
      do isite = 1, nionrecloc
        iichg = chgrecglob(isite)
        iglob = iion(iichg)
        iloc = locrec1(iglob)
        iatm = iglob
        igrd0 = igrid(1,iatm)
        jgrd0 = igrid(2,iatm)
        kgrd0 = igrid(3,iatm)
        fi = f * pchg(iichg)
        de1 = 0.0d0
        de2 = 0.0d0
        de3 = 0.0d0
        k0 = kgrd0
        do it3 = 1, bsorder
           k0 = k0 + 1
           k = k0 + 1 + (nfft3-sign(nfft3,k0))/2
           t3 = thetai3(1,it3,isite)
           dt3 = dn3 * thetai3(2,it3,isite)
           j0 = jgrd0
           do it2 = 1, bsorder
              j0 = j0 + 1
              j = j0 + 1 + (nfft2-sign(nfft2,j0))/2
              t2 = thetai2(1,it2,isite)
              dt2 = dn2 * thetai2(2,it2,isite)
              i0 = igrd0
              do it1 = 1, bsorder
                 i0 = i0 + 1
                 i = i0 + 1 + (nfft1-sign(nfft1,i0))/2
                 t1 = thetai1(1,it1,isite)
                 dt1 = dn1 * thetai1(2,it1,isite)
c
                 kstart = kstart1(rankloc+1)
                 kend = kend1(rankloc+1)
                 jstart = jstart1(rankloc+1)
                 jend = jend1(rankloc+1)
                 istart = istart1(rankloc+1)
                 iend = iend1(rankloc+1)
                 if (((k.ge.kstart).and.(k.le.kend)).and.
     $             ((j.ge.jstart).and.(j.le.jend)).and.
     $             ((i.ge.istart).and.(i.le.iend))) then
                   term = qgridin_2d(1,i-istart+1,j-jstart+1,
     $                  k-kstart+1,1)
                   goto 100
                 end if
c
                 do iproc = 1, nrec_send
                   proc = prec_send(iproc)
                   kstart = kstart1(proc+1)
                   kend = kend1(proc+1)
                   jstart = jstart1(proc+1)
                   jend = jend1(proc+1)
                   istart = istart1(proc+1)
                   iend = iend1(proc+1)
                   if (((k.ge.kstart).and.(k.le.kend)).and.
     $               ((j.ge.jstart).and.(j.le.jend)).and.
     $               ((i.ge.istart).and.(i.le.iend))) then
                     term=qgridin_2d(1,i-istart+1,j-jstart+1,k-kstart+1,
     $                iproc+1)
                     goto 100
                   end if
                  end do
 100              continue
c
                 de1 = de1 + term*dt1*t2*t3
                 de2 = de2 + term*dt2*t1*t3
                 de3 = de3 + term*dt3*t1*t2
              end do
           end do
        end do
        decrec(1,iloc) =decrec(1,iloc)+fi*(recip(1,1)*de1+recip(1,2)*de2
     &                                     +recip(1,3)*de3)
        decrec(2,iloc) =decrec(2,iloc)+fi*(recip(2,1)*de1+recip(2,2)*de2
     &                                     +recip(2,3)*de3)
        decrec(3,iloc) =decrec(3,iloc)+fi*(recip(3,1)*de1+recip(3,2)*de2
     &                                     +recip(3,3)*de3)
      end do
      time1 = mpi_wtime()
      timegrid2 = timegrid2 + time1-time0
      deallocate (qgridmpi)
      deallocate (req)
      deallocate (reqbcast)
      return
      end
c
c     "ecrealshort1d" evaluates the short range real space portion of the Ewald sum
c     energy and forces due to atomic charge interactions, using a neighbor list
c
      subroutine ecrealshort1d
      use atmlst
      use atoms
      use bound
      use boxes
      use charge
      use chgpot
      use couple
      use cutoff
      use deriv
      use domdec
      use energi
      use ewald
      use iounit
      use inter
      use math
      use molcul
      use neigh
      use potent
      use shunt
      use timestat
      use usage
      use virial
      use mpi
      implicit none
      integer i,j,k,iichg,iglob,kglob,kkchg
      integer ii,kkk

      real*8 e,de,efull
      real*8 f,fi,fik
      real*8 r,r2,rew
      real*8 rb,rb2

      real*8 xi,yi,zi
      real*8 xr,yr,zr

      real*8 erfc,erfterm
      real*8 scale,scaleterm
      real*8 dedx,dedy,dedz
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8, allocatable :: cscale(:)

      real*8 s,ds

      character*10 mode
      external erfc
 1000 format(' Warning, system moved too much since last neighbor list'
     $ ' update, try lowering nlupdate')
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (cscale(n))
c
c     initialize connected atom scaling 
c
      cscale = 1.0d0
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'SHORTEWALD'
      call switch (mode)
c
c     compute the real space Ewald energy and first derivatives
c
      do ii = 1, nionlocnl
         iichg = chgglobnl(ii)
         iglob = iion(iichg)
         i = loc(iglob)
         if (i.eq.0) then
           write(iout,1000)
           cycle
         end if
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
         do kkk = 1, nshortelst(ii)
            kkchg = shortelst(kkk,ii)
            if (kkchg.eq.0) cycle
            kglob = iion(kkchg)
            k = loc(kglob)
            if (k.eq.0) then
              write(iout,1000)
              cycle
            end if
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
               rb2 = rb * rb
               fik = fi * pchg(kkchg)
               rew = aewald * r
               erfterm = erfc (rew)
               scale = cscale(kglob)
               scaleterm = scale - 1.0d0
               e = (fik/rb) * (erfterm+scaleterm)

               call switch_respa(r,off,shortheal,s,ds)
               
               de = -fik * ((erfterm+scaleterm)/rb2
     &                 + (2.0d0*aewald/sqrtpi)*exp(-rew**2)/rb)
c
c     form the chain rule terms for derivative expressions
c
               de = de / r
               dedx = de * xr *s -ds*xr*e/r
               dedy = de * yr *s -ds*yr*e/r
               dedz = de * zr *s -ds*zr*e/r
c
c     increment the overall energy and derivative expressions
c
               ec = ec + e*s

               dec(1,i) = dec(1,i) + dedx
               dec(2,i) = dec(2,i) + dedy
               dec(3,i) = dec(3,i) + dedz
               dec(1,k) = dec(1,k) - dedx
               dec(2,k) = dec(2,k) - dedy
               dec(3,k) = dec(3,k) - dedz
c
c     increment the internal virial tensor components
c
               vxx = xr * dedx
               vyx = yr * dedx
               vzx = zr * dedx
               vyy = yr * dedy
               vzy = zr * dedy
               vzz = zr * dedz
               vir(1,1) = vir(1,1) + vxx
               vir(2,1) = vir(2,1) + vyx
               vir(3,1) = vir(3,1) + vzx
               vir(1,2) = vir(1,2) + vyx
               vir(2,2) = vir(2,2) + vyy
               vir(3,2) = vir(3,2) + vzy
               vir(1,3) = vir(1,3) + vzx
               vir(2,3) = vir(2,3) + vzy
               vir(3,3) = vir(3,3) + vzz
c
c     increment the total intramolecular energy
c
               if (molcule(iglob) .ne. molcule(kglob)) then
                  efull = (fik/rb) * scale*s
                  einter = einter + efull
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(iglob)
            cscale(i12(j,iglob)) = 1.0d0
         end do
         do j = 1, n13(iglob)
            cscale(i13(j,iglob)) = 1.0d0
         end do
         do j = 1, n14(iglob)
            cscale(i14(j,iglob)) = 1.0d0
         end do
         do j = 1, n15(iglob)
            cscale(i15(j,iglob)) = 1.0d0
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (cscale)
      return
      end
c
c     "ecreallong1d" evaluates the long range real space portion of the Ewald sum
c     energy and forces due to atomic charge interactions, using a neighbor list
c
      subroutine ecreallong1d
      use atmlst
      use atoms
      use bound
      use boxes
      use charge
      use chgpot
      use couple
      use cutoff
      use deriv
      use domdec
      use energi
      use ewald
      use iounit
      use inter
      use math
      use molcul
      use neigh
      use potent
      use shunt
      use timestat
      use usage
      use virial
      use mpi
      implicit none
      integer i,j,k,iichg,iglob,kglob,kkchg
      integer ii,kkk

      real*8 e,de,efull
      real*8 f,fi,fik
      real*8 r,r2,rew
      real*8 rb,rb2

      real*8 xi,yi,zi
      real*8 xr,yr,zr

      real*8 erfc,erfterm
      real*8 scale,scaleterm
      real*8 dedx,dedy,dedz
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8, allocatable :: cscale(:)

      real*8 s,ds,cshortcut2

      character*10 mode
      external erfc
 1000 format(' Warning, system moved too much since last neighbor list'
     $ ' update, try lowering nlupdate')
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (cscale(n))
c
c     initialize connected atom scaling 
c
      cscale = 1.0d0
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'EWALD'
      call switch (mode)
      cshortcut2 = (chgshortcut-shortheal)**2
c
c     compute the real space Ewald energy and first derivatives
c
      do ii = 1, nionlocnl
         iichg = chgglobnl(ii)
         iglob = iion(iichg)
         i = loc(iglob)
         if (i.eq.0) then
           write(iout,1000)
           cycle
         end if
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
            if (kkchg.eq.0) cycle
            kglob = iion(kkchg)
            k = loc(kglob)
            if (k.eq.0) then
              write(iout,1000)
              cycle
            end if
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
            if ((r2 .le. off2).and.(r2.ge.cshortcut2)) then
               r = sqrt(r2)
               rb = r + ebuffer
               rb2 = rb * rb
               fik = fi * pchg(kkchg)
               rew = aewald * r
               erfterm = erfc (rew)
               scale = cscale(kglob)
               scaleterm = scale - 1.0d0
               e = (fik/rb) * (erfterm+scaleterm)
c
c     use energy switching if close the cutoff distance (at short range)
c
               call switch_respa(r,chgshortcut,shortheal,s,ds)
               de = -fik * ((erfterm+scaleterm)/rb2
     &                 + (2.0d0*aewald/sqrtpi)*exp(-rew**2)/rb)
               de = -e*ds+(1-s)*de

c
c     form the chain rule terms for derivative expressions
c
               de = de / r
               dedx = de * xr
               dedy = de * yr
               dedz = de * zr
c
c     increment the overall energy and derivative expressions
c
               ec = ec + (1-s)*e
               dec(1,i) = dec(1,i) + dedx
               dec(2,i) = dec(2,i) + dedy
               dec(3,i) = dec(3,i) + dedz
               dec(1,k) = dec(1,k) - dedx
               dec(2,k) = dec(2,k) - dedy
               dec(3,k) = dec(3,k) - dedz
c
c     increment the internal virial tensor components
c
               vxx = xr * dedx
               vyx = yr * dedx
               vzx = zr * dedx
               vyy = yr * dedy
               vzy = zr * dedy
               vzz = zr * dedz
               vir(1,1) = vir(1,1) + vxx
               vir(2,1) = vir(2,1) + vyx
               vir(3,1) = vir(3,1) + vzx
               vir(1,2) = vir(1,2) + vyx
               vir(2,2) = vir(2,2) + vyy
               vir(3,2) = vir(3,2) + vzy
               vir(1,3) = vir(1,3) + vzx
               vir(2,3) = vir(2,3) + vzy
               vir(3,3) = vir(3,3) + vzz
c
c     increment the total intramolecular energy
c
               if (molcule(iglob) .ne. molcule(kglob)) then
                  efull = (fik/rb) * scale
                  einter = einter + efull
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(iglob)
            cscale(i12(j,iglob)) = 1.0d0
         end do
         do j = 1, n13(iglob)
            cscale(i13(j,iglob)) = 1.0d0
         end do
         do j = 1, n14(iglob)
            cscale(i14(j,iglob)) = 1.0d0
         end do
         do j = 1, n15(iglob)
            cscale(i15(j,iglob)) = 1.0d0
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (cscale)
      return
      end
