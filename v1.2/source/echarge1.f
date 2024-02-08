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
      use potent
      implicit none
c
c     choose the method for summing over pairwise interactions
c
       call echarge1c
c
      return
      end
c
c     echarge1c: charge electrostatic interactions
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
      use mutant
      use pme
      use potent
      use timestat
      use usage
      use virial
      use mpi
      use sizes
      implicit none
      integer ii,i,iglob,iichg,ierr
      real*8 e,de,term,qtemp
      real*8 f,fs
      real*8 xd,yd,zd
      real*8 xdtemp,ydtemp,zdtemp
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
c     set Ewald coefficient
c
      aewald = aeewald
c
c     reset lambda elec forces for lambda dynamics
c
      if (use_lambdadyn) then
        delambdae = 0d0
      end if
c
c     compute the reciprocal space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.gt.ndir-1))
     $  then
        time0 = mpi_wtime()
        if (use_crec) then
          call ecrecip1
c
c     compute recip contribution to lambda derivative for lambda dynamics
c
          if (use_lambdadyn) then
            do i = 1, nionrecloc
              iichg = chgrecglob(i)
              iglob = iion(iichg)
              if (mut(iglob).and.elambda.gt.0) then
                delambdae = delambdae 
     $   +         (pchg(iichg)*cphirec(1,i))/elambda
              end if
            end do
          end if
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
          iglob = iion(iichg)
          e = fs * pchg(iichg)**2
          ec = ec + e
          if (use_lambdadyn.and.mut(iglob)) then
            qtemp =  pchg_orig(iichg)
            delambdae = delambdae + fs*2d0*elambda*qtemp**2
          end if
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
           call MPI_ALLREDUCE(MPI_IN_PLACE,xd,1,MPI_REAL8,MPI_SUM,
     $        COMM_TINKER,ierr)
           call MPI_ALLREDUCE(MPI_IN_PLACE,yd,1,MPI_REAL8,MPI_SUM,
     $        COMM_TINKER,ierr)
           call MPI_ALLREDUCE(MPI_IN_PLACE,zd,1,MPI_REAL8,MPI_SUM,
     $        COMM_TINKER,ierr)
           term = (2.0d0/3.0d0) * f * (pi/volbox)
           e = term * (xd*xd+yd*yd+zd*zd)
           if (rank.eq.0) then
             ec = ec + e
           end if
           delambdae = delambdae + term*(xdtemp**2+ydtemp**2+zdtemp**2)
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
          call ecreal1d
          time1 = mpi_wtime()
          timereal = timereal + time1 - time0
        end if
      end if
      return
      end
cc
cc     subroutine elambdacharge1c : charge electrostatic interactions during lambda dynamics
cc
c      subroutine elambdacharge1c
c      use atmlst
c      use atoms
c      use bound
c      use boxes
c      use charge
c      use chgpot
c      use deriv
c      use energi
c      use ewald
c      use domdec
c      use iounit
c      use inter
c      use math
c      use mutant
c      use potent
c      use timestat
c      use usage
c      use virial
c      use mpi
c      use potent
c      use sizes
c      use virial
c      implicit none
c      integer ii,i,iglob,iichg,ierr
c      real*8 e,de,term
c      real*8 f,fs
c      real*8 xd,yd,zd
c      real*8 xdtemp,ydtemp,zdtemp
c      real*8 dedx,dedy,dedz
c      real*8 time0,time1
c      real*8 elambdatemp
c      real*8, allocatable :: delambdarec0(:,:),delambdarec1(:,:)
c      real*8 :: elambdarec0,elambdarec1,qtemp
c      real*8 :: vir0(3,3),vir1(3,3),virtemp(3,3)
cc
c      allocate (delambdarec0(3,nlocrec2))
c      allocate (delambdarec1(3,nlocrec2))
c      elambdatemp = elambda  
cc
cc     zero out the Ewald summation energy and derivatives
cc
c      ec = 0.0d0
c      dec = 0.0d0
c      if (nion .eq. 0)  return
c      delambdae = 0d0
cc
cc     set Ewald coefficient
cc
c      aewald = aeewald
cc
cc     compute the reciprocal space part of the Ewald summation
cc
c      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.gt.ndir-1))
c     $  then
c        time0 = mpi_wtime()
c        if (use_crec) then
cc
cc         the reciprocal part is interpolated between 0 and 1
cc
c          elambda = 0d0
c          virtemp = vir
c          call MPI_BARRIER(hostcomm,ierr)
c          if (hostrank.eq.0) call altelec
c          call MPI_BARRIER(hostcomm,ierr)
c          ec = 0d0
c          decrec = 0d0
c          vir = 0d0
c          if (elambda.lt.1d0) then
c            call ecrecip1
c          end if
c          elambdarec0  = ec
c          delambdarec0 = decrec
c          vir0 = vir
c
c          elambda = 1d0
c          call MPI_BARRIER(hostcomm,ierr)
c          if (hostrank.eq.0) call altelec
c          call MPI_BARRIER(hostcomm,ierr)
c          ec = 0d0
c          decrec = 0d0
c          vir = 0d0
c          if (elambda.gt.0d0) then
c            call ecrecip1
c          end if
c          elambdarec1  = ec
c          delambdarec1 = decrec
c          vir1 = vir
c
c          elambda = elambdatemp 
c          ec = (1-elambda)*elambdarec0 + elambda*elambdarec1
c          decrec = (1-elambda)*delambdarec0+elambda*delambdarec1
c          delambdae = delambdae + elambdarec1-elambdarec0
c          vir = virtemp + (1-elambda)*vir0 + elambda*vir1
cc     
cc         reset lambda to initial value
cc
c          call MPI_BARRIER(hostcomm,ierr)
c          if (hostrank.eq.0) call altelec
c          call MPI_BARRIER(hostcomm,ierr)
c        end if
c        time1 = mpi_wtime()
c        timerec = timerec + time1 - time0
c        if (use_pmecore) return
c      end if
cc
c      if (use_cself) then
cc
cc     compute the Ewald self-energy term over all the atoms
cc
c        f = electric / dielec
c        fs = -f * aewald / sqrtpi
c        do ii = 1, nionloc
c          iichg = chgglob(ii)
c          iglob = iion(iichg)
c          e = fs * pchg(iichg)**2
c          ec = ec + e
c          if (mut(iglob)) then
c            qtemp =  pchg_orig(iichg)
c            delambdae = delambdae + fs*2d0*elambda*qtemp**2
c          end if
c        end do
cc
cc     compute the cell dipole boundary correction term
cc
c        if (boundary .eq. 'VACUUM') then
c           xd = 0.0d0
c           yd = 0.0d0
c           zd = 0.0d0
c           xdtemp = 0.0d0
c           ydtemp = 0.0d0
c           zdtemp = 0.0d0
c           do ii = 1, nionloc
c             iichg = chgglob(ii)
c             iglob = iion(iichg)
c             i = loc(iglob)
c             xd = xd + pchg(iichg)*x(iglob)
c             yd = yd + pchg(iichg)*y(iglob)
c             zd = zd + pchg(iichg)*z(iglob)
c             if (mut(iglob)) then
c               qtemp = pchg(iichg)/elambda
c               xdtemp = xdtemp + qtemp*x(iglob)
c               ydtemp = ydtemp + qtemp*y(iglob)
c               zdtemp = zdtemp + qtemp*z(iglob)
c             end if
c           end do
c           call MPI_ALLREDUCE(MPI_IN_PLACE,xd,1,MPI_REAL8,MPI_SUM,
c     $        COMM_TINKER,ierr)
c           call MPI_ALLREDUCE(MPI_IN_PLACE,yd,1,MPI_REAL8,MPI_SUM,
c     $        COMM_TINKER,ierr)
c           call MPI_ALLREDUCE(MPI_IN_PLACE,zd,1,MPI_REAL8,MPI_SUM,
c     $        COMM_TINKER,ierr)
c           term = (2.0d0/3.0d0) * f * (pi/volbox)
c           e = term * (xd*xd+yd*yd+zd*zd)
c           if (rank.eq.0) then
c             ec = ec + e
c           end if
c           delambdae = delambdae + term*(xdtemp**2+ydtemp**2+zdtemp**2)
c           do ii = 1, nionloc
c              iichg = chgglob(ii)
c              iglob = iion(iichg)
c              i = loc(iglob)
c              de = 2.0d0 * term * pchg(iichg)
c              dedx = de * xd
c              dedy = de * yd
c              dedz = de * zd
c              dec(1,i) = dec(1,i) + dedx
c              dec(2,i) = dec(2,i) + dedy
c              dec(3,i) = dec(3,i) + dedz
c           end do
c        end if
c      end if
cc
cc     compute the real space part of the Ewald summation
cc
c      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.le.ndir-1))
c     $   then
c        if (use_creal) then
c          time0 = mpi_wtime()
c          call ecreal1d
c          time1 = mpi_wtime()
c          timereal = timereal + time1 - time0
c        end if
c      end if
c      deallocate(delambdarec0,delambdarec1)
c      return
c      end
c
c
c     "ecreal1d" evaluates the real space portion of the Ewald sum
c     energy and first derivative due to atomic charge interactions
c     using a pairwise neighbor list
c
c     if longrange, calculates just the long range part
c     if shortrange, calculates just the short range part
c
      subroutine ecreal1d
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
      use group
      use inter
      use iounit
      use math
      use mutant
      use molcul
      use neigh
      use potent
      use shunt
      use timestat
      use usage
      use virial
      use mpi
      implicit none
      integer i,j,k,iichg,iglob,kglob,kkchg,nnchg
      integer ii,kkk
      real*8 e,de,efull
      real*8 f,fi,fik,fikbis
      real*8 r,r2,rew
      real*8 rb,rb2
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 qitemp,qktemp
      real*8 erfc,erfterm
      real*8 scale,scaleterm
      real*8 fgrp
      real*8 dedx,dedy,dedz
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8, allocatable :: cscale(:)
      real*8 s,ds,cshortcut2,facts,factds
      real*8 delambdaetemp
      logical testcut,shortrange,longrange,fullrange
      character*11 mode
      character*80 :: RoutineName
      external erfc
 1000 format(' Warning, system moved too much since last neighbor list'
     $ ' update, try lowering nlupdate')

c     compute the short, long, or full real space part of the Ewald summation
      shortrange = use_cshortreal
      longrange  = use_clong
      fullrange  = .not.(shortrange.or.longrange)

      if (shortrange) then 
         RoutineName = 'ecrealshort1d'
         mode        = 'SHORTEWALD'
      else if (longrange) then
         RoutineName = 'ecreallong1d'
         mode        = 'EWALD'
      else
         RoutineName = 'ecreal1d'
         mode        = 'EWALD'
      endif

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
      call switch (mode)
      cshortcut2 = (chgshortcut-shortheal)**2

c
c     compute the real space Ewald energy and first derivatives
c
      MAINLOOP:
     &do ii = 1, nionlocnl
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
         if (shortrange) then
           nnchg = nshortelst(ii)
         else
           nnchg = nelst(ii)
         end if
         do kkk = 1, nnchg
            if (shortrange) then
              kkchg = shortelst(kkk,ii)
            else
              kkchg = elst(kkk,ii)
            end if
            if (kkchg.eq.0) cycle
            kglob = iion(kkchg)
            if (use_group)  call groups (fgrp,iglob,kglob,0,0,0,0)
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
            testcut = merge(r2 .le. off2.and.r2.ge.cshortcut2,
     &                      r2 .le. off2,
     &                      longrange
     &                     )
            if (testcut) then
               r = sqrt(r2)
               rb = r + ebuffer
               rb2 = rb * rb
               fik = fi * pchg(kkchg)
               rew = aewald * r
               erfterm = erfc (rew)
               scale = cscale(kglob)
               if (use_group)  scale = scale * fgrp
               scaleterm = scale - 1.0d0
               e = (fik/rb) * (erfterm+scaleterm)
               de = -fik * ((erfterm+scaleterm)/rb2
     &                 + (2.0d0*aewald/sqrtpi)*exp(-rew**2)/rb)
c
c
               if(shortrange .or. longrange)
     &            call switch_respa(r,chgshortcut,shortheal,s,ds)

               if(shortrange) then
                  facts  =         s
                  factds =      + ds
               else if(longrange) then
                  facts  = 1.0d0 - s
                  factds =      - ds
               else
                  facts  = 1.0d0
                  factds = 0.0d0
               endif

               de = de * facts + e * factds
               e  =  e * facts
c
c              when two interacting atoms are mutated
c
               if (use_lambdadyn) then
                  delambdaetemp = 0d0
                  if (mut(iglob).and.mut(kglob)) then
                     qitemp = pchg_orig(iichg)
                     qktemp = pchg_orig(kkchg)
c                     delambdae = delambdae + 
                     delambdaetemp = 
     $    2d0*elambda*(f*qitemp*qktemp/rb)*(erfterm+scaleterm)*facts
                  else if ((mut(iglob).and..not.mut(kglob)).or.
     $                    (mut(kglob).and..not.mut(iglob))) then
                      fikbis = f*pchg_orig(iichg) * pchg_orig(kkchg)
c                      delambdae = delambdae + (fikbis/rb) * 
                      delambdaetemp = delambdaetemp + (fikbis/rb) * 
     $                   (erfterm+scaleterm)
                  end if
                  delambdaetemp = facts*delambdaetemp
                  if (shortrange) then
                    delambdaesave = delambdaesave + delambdaetemp
                  else
                    delambdae = delambdae + delambdaetemp
                  end if
               end if
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
      end do MAINLOOP
c
c     perform deallocation of some local arrays
c
      deallocate (cscale)
      return
      end
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
      integer i,j,k,ii
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
      real*8 h1,h2,h3,f1,f2,f3
      real*8 r1,r2,r3
      integer, allocatable :: req(:),reqbcast(:)
      real*8, allocatable :: qgridmpi(:,:,:,:,:)
      real*8 time0,time1
      time0 = mpi_wtime()
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
      do i = 1, nionrecloc
        iichg = chgrecglob(i)
        iglob = iion(iichg)
        call grid_pchg_site(iglob,i,pchg(iichg))
      end do
      time1 = mpi_wtime()
      timegrid = timegrid + time1-time0
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
      timerecreccomm = timerecreccomm + time1 - time0
c
c     perform the 3-D FFT forward transformation
c
      time0 = mpi_wtime()
      call fft2d_frontmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
     $ n3mpimax)
      time1 = mpi_wtime()
      timefft = timefft + time1-time0
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
      timefft = timefft + time1-time0
c
c     MPI : Begin reception
c
      time0 = mpi_wtime()
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
c     perform dynamic allocation of some global arrays
c
      if (allocated(cphirec)) deallocate (cphirec)
      allocate (cphirec(10,max(nionrecloc,1)))
      cphirec = 0d0
      if (allocated(fphirec)) deallocate(fphirec)
      allocate (fphirec(20,max(nionrecloc,1)))
      fphirec = 0d0
c
c     get first derivatives of the reciprocal space energy
c
      time0 = mpi_wtime()
      f = electric / dielec
      do i = 1, nionrecloc
        iichg = chgrecglob(i)
        iglob = iion(iichg)
        call fphi_chg_site(iglob,i)
      end do
      do i = 1, nionrecloc
         do j = 1, 4
            fphirec(j,i) = electric * fphirec(j,i)
         end do
        call fphi_to_cphi_site (fphirec(1,i),cphirec(1,i))
      end do
      do i = 1, nionrecloc
         iichg = chgrecglob(i)
         f1 = pchg(iichg)*fphirec(2,i)
         f2 = pchg(iichg)*fphirec(3,i)
         f3 = pchg(iichg)*fphirec(4,i)
         f1 = dble(nfft1) * f1
         f2 = dble(nfft2) * f2
         f3 = dble(nfft3) * f3
         h1 = recip(1,1)*f1 + recip(1,2)*f2 + recip(1,3)*f3
         h2 = recip(2,1)*f1 + recip(2,2)*f2 + recip(2,3)*f3
         h3 = recip(3,1)*f1 + recip(3,2)*f2 + recip(3,3)*f3
         iglob = iion(iichg)
         ii = locrec(iglob)
         decrec(1,ii) = decrec(1,ii) + h1
         decrec(2,ii) = decrec(2,ii) + h2
         decrec(3,ii) = decrec(3,ii) + h3
      end do
c      dn1 = dble(nfft1)
c      dn2 = dble(nfft2)
c      dn3 = dble(nfft3)
c      do isite = 1, nionrecloc
c        iichg = chgrecglob(isite)
c        iglob = iion(iichg)
c        iloc = locrec(iglob)
c        iatm = iglob
c        igrd0 = igrid(1,iatm)
c        jgrd0 = igrid(2,iatm)
c        kgrd0 = igrid(3,iatm)
c        fi = f * pchg(iichg)
c        de1 = 0.0d0
c        de2 = 0.0d0
c        de3 = 0.0d0
c        k0 = kgrd0
c        do it3 = 1, bsorder
c           k0 = k0 + 1
c           k = k0 + 1 + (nfft3-sign(nfft3,k0))/2
c           t3 = thetai3(1,it3,isite)
c           dt3 = dn3 * thetai3(2,it3,isite)
c           j0 = jgrd0
c           do it2 = 1, bsorder
c              j0 = j0 + 1
c              j = j0 + 1 + (nfft2-sign(nfft2,j0))/2
c              t2 = thetai2(1,it2,isite)
c              dt2 = dn2 * thetai2(2,it2,isite)
c              i0 = igrd0
c              do it1 = 1, bsorder
c                 i0 = i0 + 1
c                 i = i0 + 1 + (nfft1-sign(nfft1,i0))/2
c                 t1 = thetai1(1,it1,isite)
c                 dt1 = dn1 * thetai1(2,it1,isite)
cc
c                 kstart = kstart1(rankloc+1)
c                 kend = kend1(rankloc+1)
c                 jstart = jstart1(rankloc+1)
c                 jend = jend1(rankloc+1)
c                 istart = istart1(rankloc+1)
c                 iend = iend1(rankloc+1)
c                 if (((k.ge.kstart).and.(k.le.kend)).and.
c     $             ((j.ge.jstart).and.(j.le.jend)).and.
c     $             ((i.ge.istart).and.(i.le.iend))) then
c                   term = qgridin_2d(1,i-istart+1,j-jstart+1,
c     $                  k-kstart+1,1)
c                   goto 100
c                 end if
cc
c                 do iproc = 1, nrec_send
c                   proc = prec_send(iproc)
c                   kstart = kstart1(proc+1)
c                   kend = kend1(proc+1)
c                   jstart = jstart1(proc+1)
c                   jend = jend1(proc+1)
c                   istart = istart1(proc+1)
c                   iend = iend1(proc+1)
c                   if (((k.ge.kstart).and.(k.le.kend)).and.
c     $               ((j.ge.jstart).and.(j.le.jend)).and.
c     $               ((i.ge.istart).and.(i.le.iend))) then
c                     term=qgridin_2d(1,i-istart+1,j-jstart+1,k-kstart+1,
c     $                iproc+1)
c                     goto 100
c                   end if
c                  end do
c 100              continue
cc
c                 de1 = de1 + term*dt1*t2*t3
c                 de2 = de2 + term*dt2*t1*t3
c                 de3 = de3 + term*dt3*t1*t2
c              end do
c           end do
c        end do
c        decrec(1,iloc) =decrec(1,iloc)+fi*(recip(1,1)*de1+recip(1,2)*de2
c     &                                     +recip(1,3)*de3)
c        decrec(2,iloc) =decrec(2,iloc)+fi*(recip(2,1)*de1+recip(2,2)*de2
c     &                                     +recip(2,3)*de3)
c        decrec(3,iloc) =decrec(3,iloc)+fi*(recip(3,1)*de1+recip(3,2)*de2
c     &                                     +recip(3,3)*de3)
c      end do
      deallocate (qgridmpi)
      deallocate (req)
      deallocate (reqbcast)
      time1 = mpi_wtime()
      timegrid2 = timegrid2 + time1-time0
      return
      end
