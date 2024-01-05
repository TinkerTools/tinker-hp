c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine epolar1  --  polarization energy & derivs  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "epolar1" calculates the induced dipole polarization energy
c     and derivatives with respect to Cartesian coordinates
c
c
      subroutine epolar1
      use polpot
      use potent
      use mpi
      implicit none
c
c     choose the method for summing over polarization interactions
c
      if (use_lambdadyn) then
        call elambdapolar1c
      else
        if (use_polarshortreal) then
          if (polalgshort.eq.3) then
            call epolar1tcg
          else
            call epolar1c
          end if
        else
          if (polalg.eq.3) then
            call epolar1tcg
          else
            call epolar1c
          end if
        end if
      end if
      return
      end
c
c     ###################################################################
c     ##                                                               ##
c     ##  subroutine epolar1c  --  Ewald polarization derivs via list  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     "epolar1c" calculates the dipole polarization energy and
c     derivatives with respect to Cartesian coordinates using
c     particle mesh Ewald summation and a neighbor list
c
c
      subroutine epolar1c
      use atmlst
      use atoms
      use boxes
      use chgpot
      use deriv
      use domdec
      use energi
      use ewald
      use group
      use iounit
      use math
      use mpole
      use polar
      use polpot
      use potent
      use virial
      use mpi
      implicit none
      integer i,ii,iglob,iipole,ierr
      real*8 e,f,term,fterm
      real*8 dix,diy,diz
      real*8 uix,uiy,uiz,uii
      real*8 xd,yd,zd
      real*8 xq,yq,zq
      real*8 xu,yu,zu
      real*8 xup,yup,zup
      real*8 xv,yv,zv,vterm
      real*8 xufield,yufield
      real*8 zufield
      real*8 fix(3),fiy(3),fiz(3)
      real*8 trq(3)
c
c
c     zero out the polarization energy and derivatives
c
      ep = 0.0d0
      dep = 0d0
c
      if (npole .eq. 0)  return
      aewald = apewald
c
c     set the energy unit conversion factor
c
      f = electric / dielec
c
c     check the sign of multipole components at chiral sites
c
      if (.not. use_mpole)  call chkpole(.false.)
c
c     rotate the multipole components into the global frame
c
      if (.not. use_mpole)  call rotpole
c
c     compute the induced dipoles at each polarizable atom
c
      if (use_polarshortreal) then
        if (polalg.eq.5) then
          call dcinduce_shortreal
        else
          call newinduce_shortreal
        end if
      else if (use_pmecore) then
        if (polalg.eq.5) then
          call dcinduce_pme
        else
          call newinduce_pme
        end if
      else
        if (polalg.eq.5) then
          call dcinduce_pme2
        else
          call newinduce_pme2
        end if
      end if
c
c     compute the reciprocal space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.gt.ndir-1))
     $   then
        if (use_prec) then
          call eprecip1
        end if
      end if
c
c     compute the real space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.le.ndir-1))
     $   then
        if (use_preal) then
          call epreal1c
        end if

        if (use_pself) then
c
c     compute the Ewald self-energy term over all the atoms
c
          term = 2.0d0 * aewald * aewald
          fterm = -f * aewald / sqrtpi
          do ii = 1, npoleloc
             iipole = poleglob(ii)
             dix = rpole(2,iipole)
             diy = rpole(3,iipole)
             diz = rpole(4,iipole)
             uix = uind(1,iipole)
             uiy = uind(2,iipole)
             uiz = uind(3,iipole)
             uii = dix*uix + diy*uiy + diz*uiz
             e = fterm * term * uii / 3.0d0
             ep = ep + e
          end do
c
c       compute the self-energy torque term due to induced dipole
c
          term = (4.0d0/3.0d0) * f * aewald**3 / sqrtpi
          do ii = 1, npoleloc
             iipole = poleglob(ii)
             dix = rpole(2,iipole)
             diy = rpole(3,iipole)
             diz = rpole(4,iipole)
             uix = 0.5d0 * (uind(1,iipole)+uinp(1,iipole))
             uiy = 0.5d0 * (uind(2,iipole)+uinp(2,iipole))
             uiz = 0.5d0 * (uind(3,iipole)+uinp(3,iipole))
             trq(1) = term * (diy*uiz-diz*uiy)
             trq(2) = term * (diz*uix-dix*uiz)
             trq(3) = term * (dix*uiy-diy*uix)
             call torque (iipole,trq,fix,fiy,fiz,dep)
          end do
c
c       compute the cell dipole boundary correction term
c
          if (boundary .eq. 'VACUUM') then
             xd = 0.0d0
             yd = 0.0d0
             zd = 0.0d0
             xu = 0.0d0
             yu = 0.0d0
             zu = 0.0d0
             xup = 0.0d0
             yup = 0.0d0
             zup = 0.0d0
             do i = 1, npoleloc
                iipole = poleglob(i)
                iglob = ipole(iipole)
                xd = xd + rpole(2,iipole) + rpole(1,iipole)*x(iglob)
                yd = yd + rpole(3,iipole) + rpole(1,iipole)*y(iglob)
                zd = zd + rpole(4,iipole) + rpole(1,iipole)*z(iglob)
                xu = xu + uind(1,iipole)
                yu = yu + uind(2,iipole)
                zu = zu + uind(3,iipole)
                xup = xup + uinp(1,iipole)
                yup = yup + uinp(2,iipole)
                zup = zup + uinp(3,iipole)
             end do
             call MPI_ALLREDUCE(MPI_IN_PLACE,xd,1,MPI_REAL8,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,yd,1,MPI_REAL8,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,zd,1,MPI_REAL8,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,xu,1,MPI_REAL8,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,yu,1,MPI_REAL8,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,zu,1,MPI_REAL8,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,xup,1,MPI_REAL8,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,yup,1,MPI_REAL8,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,zup,1,MPI_REAL8,MPI_SUM,
     $          COMM_TINKER,ierr)
             term = (2.0d0/3.0d0) * f * (pi/volbox)
             if (rank.eq.0) then
               ep = ep + term*(xd*xu+yd*yu+zd*zu)
             end if
             do ii = 1, npoleloc
                iipole = poleglob(ii)
                iglob = ipole(iipole)
                i = loc(iglob)
                dep(1,i) = dep(1,i) + term*rpole(1,iipole)*(xu+xup)
                dep(2,i) = dep(2,i) + term*rpole(1,iipole)*(yu+yup)
                dep(3,i) = dep(3,i) + term*rpole(1,iipole)*(zu+zup)
             end do
             xufield = -term * (xu+xup)
             yufield = -term * (yu+yup)
             zufield = -term * (zu+zup)
             do i = 1, npoleloc
                iipole = poleglob(i)
              trq(1) = rpole(3,iipole)*zufield - rpole(4,iipole)*yufield
              trq(2) = rpole(4,iipole)*xufield - rpole(2,iipole)*zufield
              trq(3) = rpole(2,iipole)*yufield - rpole(3,iipole)*xufield
                call torque (iipole,trq,fix,fiy,fiz,dep)
             end do
c
c       boundary correction to virial due to overall cell dipole
c
             xd = 0.0d0
             yd = 0.0d0
             zd = 0.0d0
             xq = 0.0d0
             yq = 0.0d0
             zq = 0.0d0
             do i = 1, npoleloc
                iipole = poleglob(i)
                iglob = ipole(iipole)
                xd = xd + rpole(2,iipole)
                yd = yd + rpole(3,iipole)
                zd = zd + rpole(4,iipole)
                xq = xq + rpole(1,iipole)*x(iglob)
                yq = yq + rpole(1,iipole)*y(iglob)
                zq = zq + rpole(1,iipole)*z(iglob)
             end do
             call MPI_ALLREDUCE(MPI_IN_PLACE,xd,1,MPI_REAL8,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,yd,1,MPI_REAL8,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,zd,1,MPI_REAL8,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,xq,1,MPI_REAL8,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,yq,1,MPI_REAL8,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,zq,1,MPI_REAL8,MPI_SUM,
     $          COMM_TINKER,ierr)
             if (rank.eq.0) then
               xv = xq * (xu+xup)
               yv = yq * (yu+yup)
               zv = zq * (zu+zup)
               vterm = xv + yv + zv + xu*xup + yu*yup + zu*zup
     &                    + xd*(xu+xup) + yd*(yu+yup) + zd*(zu+zup)
               vterm = term * vterm
               vir(1,1) = vir(1,1) + term*xv + vterm
               vir(2,1) = vir(2,1) + term*xv
               vir(3,1) = vir(3,1) + term*xv
               vir(1,2) = vir(1,2) + term*yv
               vir(2,2) = vir(2,2) + term*yv + vterm
               vir(3,2) = vir(3,2) + term*yv
               vir(1,3) = vir(1,3) + term*zv
               vir(2,3) = vir(2,3) + term*zv
               vir(3,3) = vir(3,3) + term*zv + vterm
             end if
           end if
        end if
      end if
c
c     get group polarization if necessary
c
      if (use_group) call switch_group_grad

      return
      end
c
c
c
      subroutine elambdapolar1c
      use atmlst
      use deriv
      use domdec
      use energi
      use iounit
      use mpole
      use mutant
      use polar
      use polpot
      use potent
      use uprior
      use virial
      use mpi
      implicit none
      integer i,iipole,j,k,ierr
      real*8 elambdatemp,plambda
      real*8, allocatable :: delambdap0(:,:),delambdap1(:,:)
      real*8, allocatable :: delambdaprec0(:,:),delambdaprec1(:,:)
      real*8 :: elambdap0,elambdap1
      real*8 dplambdadelambdae,d2plambdad2elambdae
      real*8 :: vir0(3,3),vir1(3,3),virtemp(3,3)
c
      allocate (delambdaprec0(3,nlocrec2))
      allocate (delambdaprec1(3,nlocrec2))
      allocate (delambdap0(3,nbloc))
      allocate (delambdap1(3,nbloc))
      elambdatemp = elambda  
c
c     zero out the polarization energy and derivatives
c
      ep = 0.0d0
      dep = 0d0
      deprec = 0d0
c
      if (npole .eq. 0)  return
      if (.not.(use_mpole)) then
        delambdae = 0d0
      end if
      virtemp = vir

c
c     polarization is interpolated between elambda=1 and elambda=0, for lambda.gt.plambda,
c     otherwise the value taken is for elambda=0
c
      elambdap1 = 0d0
      delambdap1 = 0d0
      delambdaprec1 = 0d0
      vir1 = 0d0
      if (elambda.gt.bplambda) then
        elambda = 1d0
        call MPI_BARRIER(hostcomm,ierr)
        if (hostrank.eq.0) call altelec
        call MPI_BARRIER(hostcomm,ierr)
        call rotpole
        ep = 0d0
        dep = 0d0
        deprec = 0d0
        vir = 0d0
        call epolar1c
        elambdap1  = ep
        delambdap1 = dep
        delambdaprec1 = deprec
        vir1 = vir
      end if

      elambda = 0d0
      call MPI_BARRIER(hostcomm,ierr)
      if (hostrank.eq.0) call altelec
      call MPI_BARRIER(hostcomm,ierr)
      call rotpole
      ep = 0d0
      dep = 0d0
      deprec = 0d0
      vir = 0d0
      call epolar1c
      elambdap0  = ep
      delambdap0 = dep
      delambdaprec0 = deprec
      vir0 = vir
c
c     also store the dipoles to build ASPC guess
c
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
 
      elambda = elambdatemp 
c
c     interpolation of "plambda" between bplambda and 1 as a function of
c     elambda: 
c       plambda = 0 for elambda.le.bplambda
c       u = (elambda-bplambda)/(1-bplambda)
c       plambda = u**3 for elambda.gt.plambda
c       ep = (1-plambda)*ep0 +  plambda*ep1
c
      if (elambda.le.bplambda) then
        plambda = 0d0
        dplambdadelambdae = 0d0
        d2plambdad2elambdae = 0d0
      else
        plambda = ((elambda-bplambda)/(1d0-bplambda))**3
        dplambdadelambdae = 3d0*((elambda-bplambda)**2/(1-bplambda)**3)
      end if
      ep = plambda*elambdap1 + (1-plambda)*elambdap0
      deprec = (1-plambda)*delambdaprec0+plambda*delambdaprec1
      dep = (1-plambda)*delambdap0 + plambda*delambdap1
      delambdae = delambdae + (elambdap1-elambdap0)*dplambdadelambdae
      vir = virtemp + plambda*vir1 + (1-plambda)*vir0
c
c     reset lambda to initial value
c
      call MPI_BARRIER(hostcomm,ierr)
      if (hostrank.eq.0) call altelec
      call MPI_BARRIER(hostcomm,ierr)
      call rotpole
      return
      end
c
c
c     ###################################################################
c     ##                                                               ##
c     ##  subroutine eprecip1  --  PME recip polarize energy & derivs  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     "eprecip1" evaluates the reciprocal space portion of the particle
c     mesh Ewald summation energy and gradient due to dipole polarization
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
      subroutine eprecip1
      use atmlst
      use atoms
      use bound
      use boxes
      use chgpot
      use deriv
      use domdec
      use energi
      use ewald
      use fft
      use math
      use mpole
      use pme
      use polar
      use polpot
      use potent
      use virial
      use mpi
      implicit none
      integer status(MPI_STATUS_SIZE),tag,ierr
      integer i,j,k,ii,iipole,iglob,iloc,ilocrec
      integer j1,j2,j3
      integer k1,k2,k3
      integer m1,m2,m3
      integer ntot,nff
      integer nf1,nf2,nf3
      integer deriv1(10)
      integer deriv2(10)
      integer deriv3(10)
      real*8 e,eterm,f
      real*8 r1,r2,r3
      real*8 h1,h2,h3
      real*8 f1,f2,f3
      real*8 vxx,vyy,vzz
      real*8 vxy,vxz,vyz
      real*8 volterm,denom
      real*8 hsq,expterm
      real*8 term,pterm
      real*8 vterm,struc2
      real*8 trq(3),fix(3)
      real*8 fiy(3),fiz(3)
      real*8 cphim(4),cphid(4)
      real*8 cphip(4)
      real*8 a(3,3)
      real*8 xi,yi,zi,frcx,frcy,frcz
      real*8, allocatable :: cmp(:,:),fmp(:,:)
      real*8, allocatable :: fuind(:,:)
      real*8, allocatable :: fuinp(:,:)
      real*8, allocatable :: fphid(:,:)
      real*8, allocatable :: fphip(:,:)
      real*8, allocatable :: qgrip(:,:,:,:)
      real*8, allocatable :: qgridmpi(:,:,:,:,:)
      real*8, allocatable :: pot(:),potrec(:)
      real*8, allocatable :: decfx(:)
      real*8, allocatable :: decfy(:)
      real*8, allocatable :: decfz(:)
      integer, allocatable :: reqsend(:),reqrec(:)
      integer, allocatable :: req2send(:),req2rec(:)
      integer nprocloc,commloc,rankloc
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
c
c     set Ewald coefficient
c
      aewald = apewald
c
c     indices into the electrostatic field array
c
      data deriv1  / 2, 5,  8,  9, 11, 16, 18, 14, 15, 20 /
      data deriv2  / 3, 8,  6, 10, 14, 12, 19, 16, 20, 17 /
      data deriv3  / 4, 9, 10,  7, 15, 17, 13, 20, 18, 19 /
c
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
      f = electric / dielec
c
c     perform dynamic allocation of some global arrays
c
      allocate (fuind(3,npolerecloc))
      allocate (fuinp(3,npolerecloc))
      allocate (cmp(10,npolerecloc))
      allocate (fmp(10,npolerecloc))
      allocate (fphid(10,npolerecloc))
      allocate (fphip(10,npolerecloc))
      if (allocated(cphidprec)) deallocate (cphidprec)
      allocate (cphidprec(10,max(npolerecloc,1)))
      cphidprec = 0d0
      if (allocated(fphidprec)) deallocate(fphidprec)
      allocate (fphidprec(20,max(npolerecloc,1)))
      fphidprec = 0d0
      allocate (qgridmpi(2,n1mpimax,n2mpimax,n3mpimax,nrec_recep))
      allocate (reqsend(nproc))
      allocate (reqrec(nproc))
      allocate (req2send(nproc))
      allocate (req2rec(nproc))
c
c     zero out the temporary virial accumulation variables
c
      vxx = 0.0d0
      vxy = 0.0d0
      vxz = 0.0d0
      vyy = 0.0d0
      vyz = 0.0d0
      vzz = 0.0d0
cc
cc     get the fractional to Cartesian transformation matrix
cc
c      call frac_to_cart (ftc)
c
c     initialize variables required for the scalar summation
c
      ntot = nfft1 * nfft2 * nfft3
      pterm = (pi/aewald)**2
      volterm = pi * volbox
      nff = nfft1 * nfft2
      nf1 = (nfft1+1) / 2
      nf2 = (nfft2+1) / 2
      nf3 = (nfft3+1) / 2
c
c     remove scalar sum virial from prior multipole 3-D FFT
c
c      if (use_mpole) then
         vxx = -vmxx
         vxy = -vmxy
         vxz = -vmxz
         vyy = -vmyy
         vyz = -vmyz
         vzz = -vmzz
cc
cc     compute the arrays of B-spline coefficients
cc
c      else
c         call bspline_fill
c         call table_fill
cc
cc     assign only the permanent multipoles to the PME grid
cc     and perform the 3-D FFT forward transformation
cc
         do i = 1, npolerecloc
            iipole = polerecglob(i)
            cmp(1,i) = rpole(1,iipole)
            cmp(2,i) = rpole(2,iipole)
            cmp(3,i) = rpole(3,iipole)
            cmp(4,i) = rpole(4,iipole)
            cmp(5,i) = rpole(5,iipole)
            cmp(6,i) = rpole(9,iipole)
            cmp(7,i) = rpole(13,iipole)
            cmp(8,i) = 2.0d0 * rpole(6,iipole)
            cmp(9,i) = 2.0d0 * rpole(7,iipole)
            cmp(10,i) = 2.0d0 * rpole(10,iipole)
           call cmp_to_fmp_site(cmp(1,i),fmp(1,i))
         end do
c         call cmp_to_fmp (cmp,fmp)
c         call grid_mpole (fmp)
c         call fftfront
cc
cc     make the scalar summation over reciprocal lattice
cc
c         do i = 1, ntot-1
c            k3 = i/nff + 1
c            j = i - (k3-1)*nff
c            k2 = j/nfft1 + 1
c            k1 = j - (k2-1)*nfft1 + 1
c            m1 = k1 - 1
c            m2 = k2 - 1
c            m3 = k3 - 1
c            if (k1 .gt. nf1)  m1 = m1 - nfft1
c            if (k2 .gt. nf2)  m2 = m2 - nfft2
c            if (k3 .gt. nf3)  m3 = m3 - nfft3
c            r1 = dble(m1)
c            r2 = dble(m2)
c            r3 = dble(m3)
c            h1 = recip(1,1)*r1 + recip(1,2)*r2 + recip(1,3)*r3
c            h2 = recip(2,1)*r1 + recip(2,2)*r2 + recip(2,3)*r3
c            h3 = recip(3,1)*r1 + recip(3,2)*r2 + recip(3,3)*r3
c            hsq = h1*h1 + h2*h2 + h3*h3
c            term = -pterm * hsq
c            expterm = 0.0d0
c            if (term .gt. -50.0d0) then
c               denom = volterm*hsq*bsmod1(k1)*bsmod2(k2)*bsmod3(k3)
c               expterm = exp(term) / denom
c               if (.not. use_bounds) then
c                  expterm = expterm * (1.0d0-cos(pi*xbox*sqrt(hsq)))
c               else if (octahedron) then
c                  if (mod(m1+m2+m3,2) .ne. 0)  expterm = 0.0d0
c               end if
c               struc2 = qgrid(1,k1,k2,k3)**2 + qgrid(2,k1,k2,k3)**2
c               eterm = 0.5d0 * f * expterm * struc2
c               vterm = (2.0d0/hsq) * (1.0d0-term) * eterm
c               vxx = vxx - h1*h1*vterm + eterm
c               vxy = vxy - h1*h2*vterm
c               vxz = vxz - h1*h3*vterm
c               vyy = vyy - h2*h2*vterm + eterm
c               vyz = vyz - h2*h3*vterm
c               vzz = vzz - h3*h3*vterm + eterm
c            end if
c         end do
cc
cc     account for zeroth grid point for nonperiodic system
cc
c         qfac(1,1,1) = 0.0d0
c         if (.not. use_bounds) then
c            expterm = 0.5d0 * pi / xbox
c            struc2 = qgrid(1,1,1,1)**2 + qgrid(2,1,1,1)**2
c            e = f * expterm * struc2
c            qfac(1,1,1) = expterm
c         end if
cc
cc     complete the transformation of the PME grid
cc
c         do k = 1, nfft3
c            do j = 1, nfft2
c               do i = 1, nfft1
c                  term = qfac(i,j,k)
c                  qgrid(1,i,j,k) = term * qgrid(1,i,j,k)
c                  qgrid(2,i,j,k) = term * qgrid(2,i,j,k)
c               end do
c            end do
c         end do
cc
cc     perform 3-D FFT backward transform and get potential
cc
c         call fftback
c         call fphi_mpole (fphi)
c         do i = 1, npole
c            do j = 1, 20
c               fphi(j,i) = f * fphi(j,i)
c            end do
c         end do
c         call fphi_to_cphi (fphi,cphi)
c      end if
c
c     zero out the PME grid
c
      qgrid2in_2d = 0d0
c
c     convert Cartesian induced dipoles to fractional coordinates
c
      do i = 1, 3
         a(1,i) = dble(nfft1) * recip(i,1)
         a(2,i) = dble(nfft2) * recip(i,2)
         a(3,i) = dble(nfft3) * recip(i,3)
      end do
      do ii = 1, npolerecloc
         iipole = polerecglob(ii)
         iglob = ipole(iipole)
         do j = 1, 3
            fuind(j,ii) = a(j,1)*uind(1,iipole) + a(j,2)*uind(2,iipole)
     &                      + a(j,3)*uind(3,iipole)
            fuinp(j,ii) = a(j,1)*uinp(1,iipole) + a(j,2)*uinp(2,iipole)
     &                      + a(j,3)*uinp(3,iipole)
         end do
         call grid_uind_site(iglob,ii,fuind(1,ii),fuinp(1,ii),
     $    qgrid2in_2d)
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
        call MPI_ISEND(qgrid2in_2d(1,1,1,1,i+1),
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
        qgrid2in_2d(:,:,:,:,1) = qgrid2in_2d(:,:,:,:,1)+
     $   qgridmpi(:,:,:,:,i)
      end do
c
      call fft2d_frontmpi(qgrid2in_2d,qgrid2out_2d,n1mpimax,n2mpimax,
     $ n3mpimax)
c
c     account for zeroth grid point for nonperiodic system
c
      if (.not. use_bounds) then
         expterm = 0.5d0 * pi / xbox
         struc2 = qgrid2in_2d(1,1,1,1,1)**2 + qgrid2in_2d(2,1,1,1,1)**2
         e = f * expterm * struc2
         ep = ep + e
      end if
      if ((istart2(rankloc+1).eq.1).and.(jstart2(rankloc+1).eq.1)
     $   .and.(kstart2(rankloc+1).eq.1)) then
        if (.not. use_bounds) then
           expterm = 0.5d0 * pi / xbox
           struc2 = qgrid2in_2d(1,1,1,1,1)**2 +
     $       qgrid2in_2d(2,1,1,1,1)**2
           e = f * expterm * struc2
           ep = ep + e
        end if
      end if
c
c     complete the transformation of the PME grid
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
c     perform 3-D FFT backward transform and get potential
c
      call fft2d_backmpi(qgrid2in_2d,qgrid2out_2d,n1mpimax,n2mpimax,
     $ n3mpimax)
c
c     MPI : Begin reception
c
      do i = 1, nrec_send
        tag = nprocloc*rankloc + prec_send(i) + 1
        call MPI_IRECV(qgrid2in_2d(1,1,1,1,i+1),
     $   2*n1mpimax*n2mpimax*n3mpimax,MPI_REAL8,
     $   prec_recep(i),tag,commloc,req2rec(i),ierr)
      end do
c
c     MPI : begin sending
c
      do i = 1, nrec_recep
        tag = nprocloc*prec_recep(i) + rankloc + 1
        call MPI_ISEND(qgrid2in_2d(1,1,1,1,1),
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
       do ii = 1, npolerecloc
         iipole = polerecglob(ii)
         iglob = ipole(iipole)
         call fphi_uind_site(iglob,ii,fphid(1,ii),fphip(1,ii),
     $        fphidprec(1,ii))
         do j = 1, 10
            fphid(j,ii) = electric * fphid(j,ii)
            fphip(j,ii) = electric * fphip(j,ii)
         end do
         do j = 1, 20
            fphidprec(j,ii) = electric * fphidprec(j,ii)
         end do
         do j = 1, 20
            fphirec(j,ii) = electric * fphirec(j,ii)
         end do
       end do
c
c     increment the induced dipole energy and gradient
c
      e = 0.0d0
      do i = 1, npolerecloc
         f1 = 0.0d0
         f2 = 0.0d0
         f3 = 0.0d0
         do k = 1, 3
            j1 = deriv1(k+1)
            j2 = deriv2(k+1)
            j3 = deriv3(k+1)
            e = e + fuind(k,i)*fphirec(k+1,i)
            f1 = f1 + (fuind(k,i)+fuinp(k,i))*fphirec(j1,i)
     &              + fuind(k,i)*fphip(j1,i)
     &              + fuinp(k,i)*fphid(j1,i)
            f2 = f2 + (fuind(k,i)+fuinp(k,i))*fphirec(j2,i)
     &              + fuind(k,i)*fphip(j2,i)
     &              + fuinp(k,i)*fphid(j2,i)
            f3 = f3 + (fuind(k,i)+fuinp(k,i))*fphirec(j3,i)
     &              + fuind(k,i)*fphip(j3,i)
     &              + fuinp(k,i)*fphid(j3,i)
         end do
         do k = 1, 10
            f1 = f1 + fmp(k,i)*fphidprec(deriv1(k),i)
            f2 = f2 + fmp(k,i)*fphidprec(deriv2(k),i)
            f3 = f3 + fmp(k,i)*fphidprec(deriv3(k),i)
         end do
         f1 = 0.5d0 * dble(nfft1) * f1
         f2 = 0.5d0 * dble(nfft2) * f2
         f3 = 0.5d0 * dble(nfft3) * f3
         h1 = recip(1,1)*f1 + recip(1,2)*f2 + recip(1,3)*f3
         h2 = recip(2,1)*f1 + recip(2,2)*f2 + recip(2,3)*f3
         h3 = recip(3,1)*f1 + recip(3,2)*f2 + recip(3,3)*f3
         iipole = polerecglob(i)
         iglob = ipole(iipole)
         ii = locrec(iglob)
         deprec(1,ii) = deprec(1,ii) + h1
         deprec(2,ii) = deprec(2,ii) + h2
         deprec(3,ii) = deprec(3,ii) + h3
      end do
      e = 0.5d0 * e
      ep = ep + e
c
c     set the potential to be the induced dipole average
c
      do i = 1, npolerecloc
         do k = 1, 10
            fphidprec(k,i) = 0.5d0  * fphidprec(k,i)
         end do
         call fphi_to_cphi_site(fphidprec(1,i),cphidprec(1,i))
      end do
c
c     distribute torques into the induced dipole gradient
c
      do i = 1, npolerecloc
         iipole = polerecglob(i)
         trq(1) = cmp(4,i)*cphidprec(3,i) - cmp(3,i)*cphidprec(4,i)
     &               + 2.0d0*(cmp(7,i)-cmp(6,i))*cphidprec(10,i)
     &               + cmp(9,i)*cphidprec(8,i) +cmp(10,i)*cphidprec(6,i)
     &               - cmp(8,i)*cphidprec(9,i) -cmp(10,i)*cphidprec(7,i)
         trq(2) = cmp(2,i)*cphidprec(4,i) - cmp(4,i)*cphidprec(2,i)
     &               + 2.0d0*(cmp(5,i)-cmp(7,i))*cphidprec(9,i)
     &               + cmp(8,i)*cphidprec(10,i) +cmp(9,i)*cphidprec(7,i)
     &               - cmp(9,i)*cphidprec(5,i) -cmp(10,i)*cphidprec(8,i)
         trq(3) = cmp(3,i)*cphidprec(2,i) - cmp(2,i)*cphidprec(3,i)
     &               + 2.0d0*(cmp(6,i)-cmp(5,i))*cphidprec(8,i)
     &               + cmp(8,i)*cphidprec(5,i) +cmp(10,i)*cphidprec(9,i)
     &               - cmp(8,i)*cphidprec(6,i) -cmp(9,i)*cphidprec(10,i)
         call torque_rec (iipole,trq,fix,fiy,fiz,deprec)
      end do
c
c     induced dipole contribution to the internal virial
c
      do i = 1, npolerecloc
         iipole = polerecglob(i)
         do j = 2, 4
            cphim(j) = 0.0d0
            cphid(j) = 0.0d0
            cphip(j) = 0.0d0
            do k = 2, 4
               cphim(j) = cphim(j) + ftc(j,k)*fphirec(k,i)
               cphid(j) = cphid(j) + ftc(j,k)*fphid(k,i)
               cphip(j) = cphip(j) + ftc(j,k)*fphip(k,i)
            end do
         end do
         vxx = vxx - cphidprec(2,i)*cmp(2,i)
     &         - 0.5d0*(cphim(2)*(uind(1,iipole)+uinp(1,iipole))
     &         +cphid(2)*uinp(1,iipole)+cphip(2)*uind(1,iipole))
         vxy = vxy - 0.5d0*(cphidprec(2,i)*cmp(3,i)+cphidprec(3,i)*
     $         cmp(2,i))
     &         - 0.25d0*(cphim(2)*(uind(2,iipole)+uinp(2,iipole))
     &         +cphim(3)*(uind(1,iipole)+uinp(1,iipole))
     &         +cphid(2)*uinp(2,iipole)+cphip(2)*uind(2,iipole)
     &         +cphid(3)*uinp(1,iipole)+cphip(3)*uind(1,iipole))
         vxz = vxz - 0.5d0*(cphidprec(2,i)*cmp(4,i)+cphidprec(4,i)*
     $         cmp(2,i))
     &         - 0.25d0*(cphim(2)*(uind(3,iipole)+uinp(3,iipole))
     &         +cphim(4)*(uind(1,iipole)+uinp(1,iipole))
     &         +cphid(2)*uinp(3,iipole)+cphip(2)*uind(3,iipole)
     &         +cphid(4)*uinp(1,iipole)+cphip(4)*uind(1,iipole))
         vyy = vyy - cphidprec(3,i)*cmp(3,i)
     &         - 0.5d0*(cphim(3)*(uind(2,iipole)+uinp(2,iipole))
     &         +cphid(3)*uinp(2,iipole)+cphip(3)*uind(2,iipole))
         vyz = vyz - 0.5d0*(cphidprec(3,i)*cmp(4,i)+cphidprec(4,i)*
     $         cmp(3,i))
     &         - 0.25d0*(cphim(3)*(uind(3,iipole)+uinp(3,iipole))
     &         +cphim(4)*(uind(2,iipole)+uinp(2,iipole))
     &         +cphid(3)*uinp(3,iipole)+cphip(3)*uind(3,iipole)
     &         +cphid(4)*uinp(2,iipole)+cphip(4)*uind(2,iipole))
         vzz = vzz - cphidprec(4,i)*cmp(4,i)
     &         - 0.5d0*(cphim(4)*(uind(3,iipole)+uinp(3,iipole))
     &         +cphid(4)*uinp(3,iipole)+cphip(4)*uind(3,iipole))
         vxx = vxx - 2.0d0*cmp(5,i)*cphidprec(5,i) - cmp(8,i)*
     $         cphidprec(8,i)
     &         - cmp(9,i)*cphidprec(9,i)
         vxy = vxy - (cmp(5,i)+cmp(6,i))*cphidprec(8,i)
     &         - 0.5d0*(cmp(8,i)*(cphidprec(6,i)+cphidprec(5,i))
     &         +cmp(9,i)*cphidprec(10,i)+cmp(10,i)*cphidprec(9,i))
         vxz = vxz - (cmp(5,i)+cmp(7,i))*cphidprec(9,i)
     &         - 0.5d0*(cmp(9,i)*(cphidprec(5,i)+cphidprec(7,i))
     &          +cmp(8,i)*cphidprec(10,i)+cmp(10,i)*cphidprec(8,i))
         vyy = vyy - 2.0d0*cmp(6,i)*cphidprec(6,i) - cmp(8,i)*
     $         cphidprec(8,i)
     &         - cmp(10,i)*cphidprec(10,i)
         vyz = vyz - (cmp(6,i)+cmp(7,i))*cphidprec(10,i)
     &         - 0.5d0*(cmp(10,i)*(cphidprec(6,i)+cphidprec(7,i))
     &         +cmp(8,i)*cphidprec(9,i)+cmp(9,i)*cphidprec(8,i))
         vzz = vzz - 2.0d0*cmp(7,i)*cphidprec(7,i) -
     $             cmp(9,i)*cphidprec(9,i)
     &             - cmp(10,i)*cphidprec(10,i) 
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (fuind)
      deallocate (fuinp)
      deallocate (fphid)
      deallocate (fphip)
c
c     perform dynamic allocation of some local arrays
c
      allocate (qgrip(2,isize2(rankloc+1),jsize2(rankloc+1),
     $ ksize2(rankloc+1)))
c
c     assign permanent and induced multipoles to the PME grid
c     and perform the 3-D FFT forward transformation
c
c
c    zero out the grid
c
      qgridin_2d = 0d0
c
      do i = 1, npolerecloc
         iipole = polerecglob(i)
         iglob = ipole(iipole)
         do j = 2, 4
            cmp(j,i) = cmp(j,i) + uinp(j-1,iipole)
         end do
        call cmp_to_fmp_site (cmp(1,i),fmp(1,i))
        call grid_mpole_site(iglob,i,fmp(1,i))
      end do
c
c     MPI : Begin reception
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
        call MPI_ISEND(qgridin_2d(1,1,1,1,i+1),
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
        qgridin_2d(:,:,:,:,1) = qgridin_2d(:,:,:,:,1)+
     $   qgridmpi(:,:,:,:,i)
      end do
c
      call fft2d_frontmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
     $ n3mpimax)
      do k = 1, ksize2(rankloc+1)
         do j = 1, jsize2(rankloc+1)
            do i = 1, isize2(rankloc+1)
               qgrip(1,i,j,k) = qgridout_2d(1,i,j,k)
               qgrip(2,i,j,k) = qgridout_2d(2,i,j,k)
            end do
         end do
      end do
c
c     zero out the PME grid
c
      qgridin_2d = 0d0
      do i = 1, npolerecloc
         iipole = polerecglob(i)
         iglob = ipole(iipole)
         do j = 2, 4
            cmp(j,i) = cmp(j,i) + uind(j-1,iipole) - uinp(j-1,iipole)
         end do
         call cmp_to_fmp_site(cmp(1,i),fmp(1,i))
         call grid_mpole_site(iglob,i,fmp(1,i))
      end do
c
c     MPI : Begin reception
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
        call MPI_ISEND(qgridin_2d(1,1,1,1,i+1),
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
        qgridin_2d(:,:,:,:,1) = qgridin_2d(:,:,:,:,1)+
     $   qgridmpi(:,:,:,:,i)
      end do
c
      call fft2d_frontmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
     $ n3mpimax)
c
c     make the scalar summation over reciprocal lattice
c
      if ((istart2(rankloc+1).eq.1).and.(jstart2(rankloc+1).eq.1).and.
     $   (kstart2(rankloc+1).eq.1)) then
           qfac_2d(1,1,1) = 0.0d0
      end if
      do k3 = kstart2(rankloc+1),kend2(rankloc+1)
        do k2 = jstart2(rankloc+1),jend2(rankloc+1)
          do k1 = istart2(rankloc+1),iend2(rankloc+1)
           m1 = k1 - 1
           m2 = k2 - 1
           m3 = k3 - 1
           if (k1 .gt. nf1)  m1 = m1 - nfft1
           if (k2 .gt. nf2)  m2 = m2 - nfft2
           if (k3 .gt. nf3)  m3 = m3 - nfft3
           if ((m1.eq.0).and.(m2.eq.0).and.(m3.eq.0)) goto 20
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
     $         k2-jstart2(rankloc+1)+1,k3-kstart2(rankloc+1)+1)*
     $         qgrip(1,k1-istart2(rankloc+1)+1,
     $         k2-jstart2(rankloc+1)+1,
     $         k3-kstart2(rankloc+1)+1)
     &         + qgridout_2d(2,k1-istart2(rankloc+1)+1,k2-
     $         jstart2(rankloc+1)+1,k3-kstart2(rankloc+1)+1)*
     $         qgrip(2,k1-istart2(rankloc+1)+1,
     $         k2-jstart2(rankloc+1)+1,
     $         k3-kstart2(rankloc+1)+1)
              eterm = 0.5d0 * f * expterm * struc2
              vterm = (2.0d0/hsq) * (1.0d0-term) * eterm
              vxx = vxx + h1*h1*vterm - eterm
              vxy = vxy + h1*h2*vterm
              vxz = vxz + h1*h3*vterm
              vyy = vyy + h2*h2*vterm - eterm
              vyz = vyz + h2*h3*vterm
              vzz = vzz + h3*h3*vterm - eterm
           end if
           qfac_2d(k1-istart2(rankloc+1)+1,k2-jstart2(rankloc+1)+1,k3-
     $       kstart2(rankloc+1)+1) = expterm
 20         continue
         end do
       end do
      end do
c
c     perform dynamic allocation of some local arrays
c
      if (use_chgflx) then
         allocate (pot(nbloc))
         allocate (potrec(nlocrec))
         allocate (decfx(nbloc))
         allocate (decfy(nbloc))
         allocate (decfz(nbloc))
c
c     modify the gradient and virial for charge flux
c
         pot = 0d0
         potrec = 0d0
         do i = 1, npolerecloc
           iipole = polerecglob(i)
           iglob = ipole(iipole)
           ilocrec = locrec(iglob)
           potrec(ilocrec) = cphidprec(1,i)
         end do
c
c     communicate reciprocal potential to get local values
c
         call commpotrec(pot,potrec)
         call commpot(pot,1)
c
         call dcflux (pot,decfx,decfy,decfz)
         do i = 1, npolebloc
            iipole = poleglob(i)
            ii = ipole(iipole)
            iloc = loc(ii)
            xi = x(ii)
            yi = y(ii)
            zi = z(ii)
            frcx = decfx(iloc)
            frcy = decfy(iloc)
            frcz = decfz(iloc)
            dep(1,iloc) = dep(1,iloc) + frcx
            dep(2,iloc) = dep(2,iloc) + frcy
            dep(3,iloc) = dep(3,iloc) + frcz
            vxx = vxx + xi*frcx
            vxy = vxy + yi*frcx
            vxz = vxz + zi*frcx
            vyy = vyy + yi*frcy
            vyz = vyz + zi*frcy
            vzz = vzz + zi*frcz
         end do
c
c     perform deallocation of some local arrays
c
         deallocate (pot)
         deallocate (potrec)
         deallocate (decfx)
         deallocate (decfy)
         deallocate (decfz)
      end if
c
c     increment the internal virial tensor components
c
      vir(1,1) = vir(1,1) + vxx
      vir(2,1) = vir(2,1) + vxy
      vir(3,1) = vir(3,1) + vxz
      vir(1,2) = vir(1,2) + vxy
      vir(2,2) = vir(2,2) + vyy
      vir(3,2) = vir(3,2) + vyz
      vir(1,3) = vir(1,3) + vxz
      vir(2,3) = vir(2,3) + vyz
      vir(3,3) = vir(3,3) + vzz
c
c     perform deallocation of some local arrays
c
      deallocate (qgrip)
c      deallocate (fuind)
c      deallocate (fuinp)
      deallocate (cmp)
      deallocate (fmp)
c      deallocate (fphid)
c      deallocate (fphip)
c      deallocate (fphidp)
      deallocate (qgridmpi)
      deallocate (reqsend)
      deallocate (reqrec)
      deallocate (req2send)
      deallocate (req2rec)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine epreal1c  --  Ewald real space derivs via list  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "epreal1c" evaluates the real space portion of the Ewald
c     summation energy and gradient due to dipole polarization
c     via a neighbor list
c
c     if shortrange, calculates just the short range part
c
      subroutine epreal1c
      use atmlst
      use atoms
      use bound
      use chgpen
      use chgpot
      use couple
      use cutoff
      use deriv
      use domdec
      use energi
      use ewald
      use inter
      use iounit
      use math
      use molcul
      use mpole
      use neigh
      use polar
      use polgrp
      use polpot
      use potent
      use shunt
      use virial
      use mpi
      implicit none
      integer i,j,k,nnelst
      integer iglob,kglob,kbis
      integer ii,kkk
      integer it,kt
      integer iipole,kkpole
      integer ix,iy,iz
      real*8 e
      real*8 f,pgamma
      real*8 pdi,pti,ddi
      real*8 damp,expdamp
      real*8 temp3,temp5,temp7
      real*8 sc3,sc5,sc7
      real*8 psc3,psc5,psc7
      real*8 dsc3,dsc5,dsc7
      real*8 usc3,usc5
      real*8 psr3,psr5,psr7
      real*8 dsr3,dsr5,dsr7
      real*8 usr3,usr5
      real*8 rr3core,rr5core
      real*8 rr3i,rr5i
      real*8 rr7i,rr9i
      real*8 rr3k,rr5k
      real*8 rr7k,rr9k
      real*8 rr5ik,rr7ik
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 uix,uiy,uiz
      real*8 uixp,uiyp,uizp
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 ukx,uky,ukz
      real*8 ukxp,ukyp,ukzp
      real*8 dir,uir,uirp
      real*8 dkr,ukr,ukrp
      real*8 qix,qiy,qiz,qir
      real*8 qkx,qky,qkz,qkr
      real*8 corei,corek
      real*8 vali,valk
      real*8 alphai,alphak
      real*8 tuir,tukr
      real*8 tixx,tiyy,tizz
      real*8 tixy,tixz,tiyz
      real*8 tkxx,tkyy,tkzz
      real*8 tkxy,tkxz,tkyz
      real*8 tix3,tiy3,tiz3
      real*8 tix5,tiy5,tiz5
      real*8 tkx3,tky3,tkz3
      real*8 tkx5,tky5,tkz5
      real*8 term1,term2,term3
      real*8 term4,term5
      real*8 term6,term7
      real*8 term1core
      real*8 term1i,term2i,term3i
      real*8 term4i,term5i,term6i
      real*8 term7i,term8i
      real*8 term1k,term2k,term3k
      real*8 term4k,term5k,term6k
      real*8 term7k,term8k
      real*8 poti,potk
      real*8 depx,depy,depz
      real*8 frcx,frcy,frcz
      real*8 xix,yix,zix
      real*8 xiy,yiy,ziy
      real*8 xiz,yiz,ziz
      real*8 vxx,vyy,vzz
      real*8 vxy,vxz,vyz
      real*8 rc3(3),rc5(3),rc7(3)
      real*8 prc3(3),prc5(3),prc7(3)
      real*8 drc3(3),drc5(3),drc7(3)
      real*8 urc3(3),urc5(3),tep(3)
      real*8 fix(3),fiy(3),fiz(3)
      real*8 dmp3,dmp5,dmp7
      real*8 scalek
      real*8 dmpi(9),dmpk(9)
      real*8 dmpik(9),dmpe(9)
      real*8 fip(3),fkp(3)
      logical shortrange,longrange,fullrange
      real*8, allocatable :: pscale(:)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: uscale(:)
      real*8, allocatable :: wscale(:)
      real*8, allocatable :: ufld(:,:)
      real*8, allocatable :: dufld(:,:)
      real*8, allocatable :: pot(:)
      real*8, allocatable :: decfx(:)
      real*8, allocatable :: decfy(:)
      real*8, allocatable :: decfz(:)
      character*11 mode
      character*80 :: RoutineName
      external erfc
 1000 format(' Warning, system moved too much since last neighbor list'
     $   ' update, try lowering nlupdate')

c     compute the short, or full real space part of the summation
      shortrange = use_polarshortreal
      longrange  = .false.
      fullrange  = .not.(shortrange.or.longrange)

      if (shortrange) then 
         RoutineName = 'eprealshort1c'
         mode        = 'SHORTEWALD'
      else if (longrange) then
         RoutineName = 'epreallong1c'
         mode        = 'EWALD'
      else
         RoutineName = 'epreal1c'
         mode        = 'EWALD'
      endif
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (pscale(n))
      allocate (dscale(n))
      allocate (uscale(n))
      allocate (wscale(n))
      allocate (ufld(3,nbloc))
      allocate (dufld(6,nbloc))
      allocate (pot(nbloc))
      allocate (decfx(nbloc))
      allocate (decfy(nbloc))
      allocate (decfz(nbloc))
c
c     set exclusion coefficients and arrays to store fields
c
      pscale = 1.0d0
      dscale = 1.0d0
      uscale = 1.0d0
      wscale = 1.0d0
      ufld = 0.0d0
      dufld = 0.0d0
      pot = 0d0

c
c     set conversion factor, cutoff and switching coefficients
c
      f = 0.5d0 * electric / dielec
      call switch (mode)
c
c
c     compute the dipole polarization gradient components
c
      do ii = 1, npolelocnl
         iipole = poleglobnl(ii)
         iglob = ipole(iipole)
         i = loc(iglob)
         if ((i.eq.0).or.(i.gt.nbloc)) then
           write(iout,1000)
           cycle
         end if
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
         uix = uind(1,iipole)
         uiy = uind(2,iipole)
         uiz = uind(3,iipole)
         uixp = uinp(1,iipole)
         uiyp = uinp(2,iipole)
         uizp = uinp(3,iipole)
         if (use_thole) then
            pdi = pdamp(iipole)
            pti = thole(iipole)
            ddi = tholed(iipole)
         else if (use_chgpen) then
            corei = pcore(iipole)
            vali = pval(iipole)
            alphai = palpha(iipole)
         end if
c
c     set exclusion coefficients for connected atoms
c
         if (dpequal) then
            do j = 1, n12(iglob)
               pscale(i12(j,iglob)) = p2scale
               do k = 1, np11(iglob)
                  if (i12(j,iglob) .eq. ip11(k,iglob))
     &               pscale(i12(j,iglob)) = p2iscale
               end do
               dscale(i12(j,iglob)) = pscale(i12(j,iglob))
               wscale(i12(j,iglob)) = w2scale
            end do
            do j = 1, n13(iglob)
               pscale(i13(j,iglob)) = p3scale
               do k = 1, np11(iglob)
                  if (i13(j,iglob) .eq. ip11(k,iglob))
     &               pscale(i13(j,iglob)) = p3iscale
               end do
               dscale(i13(j,iglob)) = pscale(i13(j,iglob))
               wscale(i13(j,iglob)) = w3scale
            end do
            do j = 1, n14(iglob)
               pscale(i14(j,iglob)) = p4scale
               do k = 1, np11(iglob)
                   if (i14(j,iglob) .eq. ip11(k,iglob))
     &               pscale(i14(j,iglob)) = p4iscale
               end do
               dscale(i14(j,iglob)) = pscale(i14(j,iglob))
               wscale(i14(j,iglob)) = w4scale
            end do
            do j = 1, n15(iglob)
               pscale(i15(j,iglob)) = p5scale
               do k = 1, np11(iglob)
                  if (i15(j,iglob) .eq. ip11(k,iglob))
     &               pscale(i15(j,iglob)) = p5iscale
               end do
               dscale(i15(j,iglob)) = pscale(i15(j,iglob))
               wscale(i15(j,iglob)) = w5scale
            end do
            do j = 1, np11(iglob)
               uscale(ip11(j,iglob)) = u1scale
            end do
            do j = 1, np12(iglob)
               uscale(ip12(j,iglob)) = u2scale
            end do
            do j = 1, np13(iglob)
               uscale(ip13(j,iglob)) = u3scale
            end do
            do j = 1, np14(iglob)
               uscale(ip14(j,iglob)) = u4scale
            end do
         else
            do j = 1, n12(iglob)
               pscale(i12(j,iglob)) = p2scale
               do k = 1, np11(iglob)
                  if (i12(j,iglob) .eq. ip11(k,iglob))
     &               pscale(i12(j,iglob)) = p2iscale
               end do
               wscale(i12(j,iglob)) = w2scale
            end do
            do j = 1, n13(iglob)
               pscale(i13(j,iglob)) = p3scale
               do k = 1, np11(iglob)
                  if (i13(j,iglob) .eq. ip11(k,iglob))
     &               pscale(i13(j,iglob)) = p3iscale
               end do
               wscale(i13(j,iglob)) = w3scale
            end do
            do j = 1, n14(iglob)
               pscale(i14(j,iglob)) = p4scale
               do k = 1, np11(iglob)
                   if (i14(j,iglob) .eq. ip11(k,iglob))
     &               pscale(i14(j,iglob)) = p4iscale
               end do
               wscale(i14(j,iglob)) = w4scale
            end do
            do j = 1, n15(iglob)
               pscale(i15(j,iglob)) = p5scale
               do k = 1, np11(iglob)
                  if (i15(j,iglob) .eq. ip11(k,iglob))
     &               pscale(i15(j,iglob)) = p5iscale
               end do
               wscale(i15(j,iglob)) = w5scale
            end do
            do j = 1, np11(iglob)
               dscale(ip11(j,iglob)) = d1scale
               uscale(ip11(j,iglob)) = u1scale
            end do
            do j = 1, np12(iglob)
               dscale(ip12(j,iglob)) = d2scale
               uscale(ip12(j,iglob)) = u2scale
            end do
            do j = 1, np13(iglob)
               dscale(ip13(j,iglob)) = d3scale
               uscale(ip13(j,iglob)) = u3scale
            end do
            do j = 1, np14(iglob)
               dscale(ip14(j,iglob)) = d4scale
               uscale(ip14(j,iglob)) = u4scale
            end do
         end if
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
            kbis = loc(kglob)
            if (kbis.eq.0) then
              write(iout,1000)
              cycle
            end if
            xr = x(kglob) - xi
            yr = y(kglob) - yi
            zr = z(kglob) - zi
            if (use_bounds)  call image (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. off2) then
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
               ukx = uind(1,kkpole)
               uky = uind(2,kkpole)
               ukz = uind(3,kkpole)
               ukxp = uinp(1,kkpole)
               ukyp = uinp(2,kkpole)
               ukzp = uinp(3,kkpole)
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
               uir = uix*xr + uiy*yr + uiz*zr
               uirp = uixp*xr + uiyp*yr + uizp*zr
               ukr = ukx*xr + uky*yr + ukz*zr
               ukrp = ukxp*xr + ukyp*yr + ukzp*zr
c
c     get reciprocal distance terms for this interaction
c
               rr1 = f / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
c
c     calculate real space Ewald error function damping
c
               call dampewald (9,r,r2,f,dmpe)
c
c     apply Thole polarization damping to scale factors
c
               sc3 = 1.0d0
               sc5 = 1.0d0
               sc7 = 1.0d0
               do j = 1, 3
                  rc3(j) = 0.0d0
                  rc5(j) = 0.0d0
                  rc7(j) = 0.0d0
               end do
c
c     apply Thole polarization damping to scale factors
c
               if (use_thole) then
                  damp = pdi * pdamp(kkpole)
                  it = jpolar(iglob)
                  kt = jpolar(kglob)
                  if (use_tholed) then
                     pgamma = thdval(it,kt)
                     if (damp.ne.0.0d0 .and. pgamma.ne.0.0d0) then
                        damp = pgamma * (r/damp)**(1.5d0)
                        if (damp .lt. 50.0d0) then
                           expdamp = exp(-damp) 
                           sc3 = 1.0d0 - expdamp 
                           sc5 = 1.0d0 - expdamp*(1.0d0+0.5d0*damp)
                           sc7 = 1.0d0 - expdamp*(1.0d0+0.65d0*damp
     &                                      +0.15d0*damp**2)
                           temp3 = 1.5d0 * damp * expdamp / r2
                           temp5 = 0.5d0 * (1.0d0+damp)
                           temp7 = 0.7d0 + 0.15d0*damp**2/temp5
                           rc3(1) = xr * temp3
                           rc3(2) = yr * temp3
                           rc3(3) = zr * temp3
                           rc5(1) = rc3(1) * temp5
                           rc5(2) = rc3(2) * temp5
                           rc5(3) = rc3(3) * temp5
                           rc7(1) = rc5(1) * temp7
                           rc7(2) = rc5(2) * temp7
                           rc7(3) = rc5(3) * temp7
                        end if
                     end if
                  else
                     pgamma = thlval(it,kt)
                     if (damp.ne.0.0d0 .and. pgamma.ne.0.0d0) then
                        damp = pgamma * (r/damp)**3
                        if (damp .lt. 50.0d0) then
                           expdamp = exp(-damp)
                           sc3 = 1.0d0 - expdamp
                           sc5 = 1.0d0 - (1.0d0+damp)*expdamp
                           sc7 = 1.0d0 - (1.0d0+damp+0.6d0*damp**2)
     &                                          *expdamp
                           temp3 = 3.0d0 * damp * expdamp / r2
                           temp5 = damp
                           temp7 = -0.2d0 + 0.6d0*damp
                           rc3(1) = xr * temp3
                           rc3(2) = yr * temp3
                           rc3(3) = zr * temp3
                           rc5(1) = rc3(1) * temp5
                           rc5(2) = rc3(2) * temp5
                           rc5(3) = rc3(3) * temp5
                           rc7(1) = rc5(1) * temp7
                           rc7(2) = rc5(2) * temp7
                           rc7(3) = rc5(3) * temp7
                        end if
                     end if
                   end if
                   psc3 = 1.0d0 - sc3*pscale(kglob)
                   psc5 = 1.0d0 - sc5*pscale(kglob)
                   psc7 = 1.0d0 - sc7*pscale(kglob)
                   dsc3 = 1.0d0 - sc3*dscale(kglob)
                   dsc5 = 1.0d0 - sc5*dscale(kglob)
                   dsc7 = 1.0d0 - sc7*dscale(kglob)
                   usc3 = 1.0d0 - sc3*uscale(kglob)
                   usc5 = 1.0d0 - sc5*uscale(kglob)
                   psr3 = dmpe(3) - psc3*rr3
                   psr5 = dmpe(5) - psc5*rr5
                   psr7 = dmpe(7) - psc7*rr7
                   dsr3 = dmpe(3) - dsc3*rr3
                   dsr5 = dmpe(5) - dsc5*rr5
                   dsr7 = dmpe(7) - dsc7*rr7
                   usr3 = dmpe(3) - usc3*rr3
                   usr5 = dmpe(5) - usc5*rr5
                   do j = 1, 3
                      prc3(j) = rc3(j) * pscale(kglob)
                      prc5(j) = rc5(j) * pscale(kglob)
                      prc7(j) = rc7(j) * pscale(kglob)
                      drc3(j) = rc3(j) * dscale(kglob)
                      drc5(j) = rc5(j) * dscale(kglob)
                      drc7(j) = rc7(j) * dscale(kglob)
                      urc3(j) = rc3(j) * uscale(kglob)
                      urc5(j) = rc5(j) * uscale(kglob)
                   end do
c
c     apply charge penetration damping to scale factors
c
               else if (use_chgpen) then
                  corek = pcore(kkpole)
                  valk = pval(kkpole)
                  alphak = palpha(kkpole)
                  call damppole (r,9,alphai,alphak,dmpi,dmpk,dmpik)
                  rr3core = dmpe(3) - (1.0d0-dscale(kglob))*rr3
                  rr5core = dmpe(5) - (1.0d0-dscale(kglob))*rr5
                  rr3i = dmpe(3) - (1.0d0-dscale(kglob)*dmpi(3))*rr3
                  rr5i = dmpe(5) - (1.0d0-dscale(kglob)*dmpi(5))*rr5
                  rr7i = dmpe(7) - (1.0d0-dscale(kglob)*dmpi(7))*rr7
                  rr9i = dmpe(9) - (1.0d0-dscale(kglob)*dmpi(9))*rr9
                  rr3k = dmpe(3) - (1.0d0-dscale(kglob)*dmpk(3))*rr3
                  rr5k = dmpe(5) - (1.0d0-dscale(kglob)*dmpk(5))*rr5
                  rr7k = dmpe(7) - (1.0d0-dscale(kglob)*dmpk(7))*rr7
                  rr9k = dmpe(9) - (1.0d0-dscale(kglob)*dmpk(9))*rr9
                  rr5ik = dmpe(5) - (1.0d0-wscale(kglob)*dmpik(5))*rr5
                  rr7ik = dmpe(7) - (1.0d0-wscale(kglob)*dmpik(7))*rr7
               end if
c
c     store the potential at each site for use in charge flux
c
               if (use_chgflx) then
                  if (use_thole) then
                     poti = -ukr*psr3 - ukrp*dsr3
                     potk = uir*psr3 + uirp*dsr3
                  else if (use_chgpen) then
                     poti = -2.0d0 * ukr * rr3i
                     potk = 2.0d0 * uir * rr3k
                  end if
                  pot(i) = pot(i) + poti 
                  pot(kbis) = pot(kbis) + potk 
               end if
c
c     get the induced dipole field used for dipole torques
c
               if (use_thole) then
                  tix3 = psr3*ukx + dsr3*ukxp
                  tiy3 = psr3*uky + dsr3*ukyp
                  tiz3 = psr3*ukz + dsr3*ukzp
                  tkx3 = psr3*uix + dsr3*uixp
                  tky3 = psr3*uiy + dsr3*uiyp
                  tkz3 = psr3*uiz + dsr3*uizp
                  tuir = -psr5*ukr - dsr5*ukrp
                  tukr = -psr5*uir - dsr5*uirp
               else if (use_chgpen) then
                  tix3 = 2.0d0*rr3i*ukx
                  tiy3 = 2.0d0*rr3i*uky
                  tiz3 = 2.0d0*rr3i*ukz
                  tkx3 = 2.0d0*rr3k*uix
                  tky3 = 2.0d0*rr3k*uiy
                  tkz3 = 2.0d0*rr3k*uiz
                  tuir = -2.0d0*rr5i*ukr
                  tukr = -2.0d0*rr5k*uir
               end if
               ufld(1,i) = ufld(1,i) + tix3 + xr*tuir
               ufld(2,i) = ufld(2,i) + tiy3 + yr*tuir
               ufld(3,i) = ufld(3,i) + tiz3 + zr*tuir
               ufld(1,kbis) = ufld(1,kbis) + tkx3 + xr*tukr
               ufld(2,kbis) = ufld(2,kbis) + tky3 + yr*tukr
               ufld(3,kbis) = ufld(3,kbis) + tkz3 + zr*tukr
c
c     get permanent field to get polarization energy
c
               if (use_thole) then
                  call damptholed (iipole,kkpole,7,r,dmpik)
                  scalek = pscale(kglob)
                  dmp3 = dmpe(3) - (1.0d0-scalek*dmpik(3))*rr3
                  dmp5 = dmpe(5) - (1.0d0-scalek*dmpik(5))*rr5
                  dmp7 = dmpe(7) - (1.0d0-scalek*dmpik(7))*rr7
                  fip(1) = -xr*(dmp3*ck-dmp5*dkr+dmp7*qkr)
     &                        - dmp3*dkx + 2.0d0*dmp5*qkx
                  fip(2) = -yr*(dmp3*ck-dmp5*dkr+dmp7*qkr)
     &                        - dmp3*dky + 2.0d0*dmp5*qky
                  fip(3) = -zr*(dmp3*ck-dmp5*dkr+dmp7*qkr)
     &                        - dmp3*dkz + 2.0d0*dmp5*qkz
                  fkp(1) = xr*(dmp3*ci+dmp5*dir+dmp7*qir)
     &                        - dmp3*dix - 2.0d0*dmp5*qix
                  fkp(2) = yr*(dmp3*ci+dmp5*dir+dmp7*qir)
     &                        - dmp3*diy - 2.0d0*dmp5*qiy
                  fkp(3) = zr*(dmp3*ci+dmp5*dir+dmp7*qir)
     &                        - dmp3*diz - 2.0d0*dmp5*qiz
c
c     find the field components for charge penetration damping
c
               else if (use_chgpen) then
                  corek = pcore(kkpole)
                  valk = pval(kkpole)
                  alphak = palpha(kkpole)
                  call dampdir (r,alphai,alphak,dmpi,dmpk)
                  scalek = pscale(kglob)
                  rr3i = dmpe(3) - (1.0d0-scalek*dmpi(3))*rr3
                  rr5i = dmpe(5) - (1.0d0-scalek*dmpi(5))*rr5
                  rr7i = dmpe(7) - (1.0d0-scalek*dmpi(7))*rr7
                  rr3k = dmpe(3) - (1.0d0-scalek*dmpk(3))*rr3
                  rr5k = dmpe(5) - (1.0d0-scalek*dmpk(5))*rr5
                  rr7k = dmpe(7) - (1.0d0-scalek*dmpk(7))*rr7
                  rr3 = dmpe(3) - (1.0d0-scalek)*rr3
                  fip(1) = -xr*(rr3*corek + rr3k*valk
     &                        - rr5k*dkr + rr7k*qkr)
     &                        - rr3k*dkx + 2.0d0*rr5k*qkx
                  fip(2) = -yr*(rr3*corek + rr3k*valk
     &                        - rr5k*dkr + rr7k*qkr)
     &                        - rr3k*dky + 2.0d0*rr5k*qky
                  fip(3) = -zr*(rr3*corek + rr3k*valk
     &                        - rr5k*dkr + rr7k*qkr)
     &                        - rr3k*dkz + 2.0d0*rr5k*qkz
                  fkp(1) = xr*(rr3*corei + rr3i*vali
     &                        + rr5i*dir + rr7i*qir)
     &                        - rr3i*dix - 2.0d0*rr5i*qix
                  fkp(2) = yr*(rr3*corei + rr3i*vali
     &                        + rr5i*dir + rr7i*qir)
     &                        - rr3i*diy - 2.0d0*rr5i*qiy
                  fkp(3) = zr*(rr3*corei + rr3i*vali
     &                        + rr5i*dir + rr7i*qir)
     &                        - rr3i*diz - 2.0d0*rr5i*qiz
               end if

c
c     compute the energy contribution for this interaction
c
               e = 0d0
               do j = 1, 3
                 e = e - uind(j,iipole)*fip(j) - uind(j,kkpole)*fkp(j)
               end do
               ep = ep + e
c
c     get induced dipole field gradient used for quadrupole torques
c
               if (use_thole) then
                  tix5 = 2.0d0 * (psr5*ukx+dsr5*ukxp)
                  tiy5 = 2.0d0 * (psr5*uky+dsr5*ukyp)
                  tiz5 = 2.0d0 * (psr5*ukz+dsr5*ukzp)
                  tkx5 = 2.0d0 * (psr5*uix+dsr5*uixp)
                  tky5 = 2.0d0 * (psr5*uiy+dsr5*uiyp)
                  tkz5 = 2.0d0 * (psr5*uiz+dsr5*uizp)
                  tuir = -psr7*ukr - dsr7*ukrp
                  tukr = -psr7*uir - dsr7*uirp
               else if (use_chgpen) then
                  tix5 = 4.0d0 * (rr5i*ukx)
                  tiy5 = 4.0d0 * (rr5i*uky)
                  tiz5 = 4.0d0 * (rr5i*ukz)
                  tkx5 = 4.0d0 * (rr5k*uix)
                  tky5 = 4.0d0 * (rr5k*uiy)
                  tkz5 = 4.0d0 * (rr5k*uiz)
                  tuir = -2.0d0*rr7i*ukr 
                  tukr = -2.0d0*rr7k*uir 
               end if
               dufld(1,i) = dufld(1,i) + xr*tix5 + xr*xr*tuir
               dufld(2,i) = dufld(2,i) + xr*tiy5 + yr*tix5
     &                         + 2.0d0*xr*yr*tuir
               dufld(3,i) = dufld(3,i) + yr*tiy5 + yr*yr*tuir
               dufld(4,i) = dufld(4,i) + xr*tiz5 + zr*tix5
     &                         + 2.0d0*xr*zr*tuir
               dufld(5,i) = dufld(5,i) + yr*tiz5 + zr*tiy5
     &                         + 2.0d0*yr*zr*tuir
               dufld(6,i) = dufld(6,i) + zr*tiz5 + zr*zr*tuir
               dufld(1,kbis) = dufld(1,kbis) - xr*tkx5 - xr*xr*tukr
               dufld(2,kbis) = dufld(2,kbis) - xr*tky5 - yr*tkx5
     &                         - 2.0d0*xr*yr*tukr
               dufld(3,kbis) = dufld(3,kbis) - yr*tky5 - yr*yr*tukr
               dufld(4,kbis) = dufld(4,kbis) - xr*tkz5 - zr*tkx5
     &                         - 2.0d0*xr*zr*tukr
               dufld(5,kbis) = dufld(5,kbis) - yr*tkz5 - zr*tky5
     &                         - 2.0d0*yr*zr*tukr
               dufld(6,kbis) = dufld(6,kbis) - zr*tkz5 - zr*zr*tukr
c
c     get the dEd/dR terms used for direct polarization force
c
               if (use_thole) then
                  term1 = dmpe(5) - dsc3*rr5
                  term2 = dmpe(7) - dsc5*rr7
                  term3 = -dsr3 + term1*xr*xr - rr3*xr*drc3(1)
                  term4 = rr3*drc3(1) - term1*xr - dsr5*xr
                  term5 = term2*xr*xr - dsr5 - rr5*xr*drc5(1)
                  term6 = (dmpe(9)-dsc7*rr9)*xr*xr - dmpe(7)
     &                       - rr7*xr*drc7(1)
                  term7 = rr5*drc5(1) - 2.0d0*dmpe(7)*xr
     &                       + (dsc5+1.5d0*dsc7)*rr7*xr
                  tixx = ci*term3 + dix*term4 + dir*term5
     &                      + 2.0d0*dsr5*qixx + (qiy*yr+qiz*zr)*dsc7*rr7
     &                      + 2.0d0*qix*term7 + qir*term6
                  tkxx = ck*term3 - dkx*term4 - dkr*term5
     &                      + 2.0d0*dsr5*qkxx + (qky*yr+qkz*zr)*dsc7*rr7
     &                      + 2.0d0*qkx*term7 + qkr*term6
                  term3 = -dsr3 + term1*yr*yr - rr3*yr*drc3(2)
                  term4 = rr3*drc3(2) - term1*yr - dsr5*yr
                  term5 = term2*yr*yr - dsr5 - rr5*yr*drc5(2)
                  term6 = (dmpe(9)-dsc7*rr9)*yr*yr - dmpe(7)
     &                       - rr7*yr*drc7(2)
                  term7 = rr5*drc5(2) - 2.0d0*dmpe(7)*yr
     &                       + (dsc5+1.5d0*dsc7)*rr7*yr
                  tiyy = ci*term3 + diy*term4 + dir*term5
     &                      + 2.0d0*dsr5*qiyy + (qix*xr+qiz*zr)*dsc7*rr7
     &                      + 2.0d0*qiy*term7 + qir*term6
                  tkyy = ck*term3 - dky*term4 - dkr*term5
     &                      + 2.0d0*dsr5*qkyy + (qkx*xr+qkz*zr)*dsc7*rr7
     &                      + 2.0d0*qky*term7 + qkr*term6
                  term3 = -dsr3 + term1*zr*zr - rr3*zr*drc3(3)
                  term4 = rr3*drc3(3) - term1*zr - dsr5*zr
                  term5 = term2*zr*zr - dsr5 - rr5*zr*drc5(3)
                  term6 = (dmpe(9)-dsc7*rr9)*zr*zr - dmpe(7)
     &                       - rr7*zr*drc7(3)
                  term7 = rr5*drc5(3) - 2.0d0*dmpe(7)*zr
     &                       + (dsc5+1.5d0*dsc7)*rr7*zr
                  tizz = ci*term3 + diz*term4 + dir*term5
     &                      + 2.0d0*dsr5*qizz + (qix*xr+qiy*yr)*dsc7*rr7
     &                      + 2.0d0*qiz*term7 + qir*term6
                  tkzz = ck*term3 - dkz*term4 - dkr*term5
     &                      + 2.0d0*dsr5*qkzz + (qkx*xr+qky*yr)*dsc7*rr7
     &                      + 2.0d0*qkz*term7 + qkr*term6
                  term3 = term1*xr*yr - rr3*yr*drc3(1)
                  term4 = rr3*drc3(1) - term1*xr
                  term5 = term2*xr*yr - rr5*yr*drc5(1)
                  term6 = (dmpe(9)-dsc7*rr9)*xr*yr - rr7*yr*drc7(1)
                  term7 = rr5*drc5(1) - term2*xr
                  tixy = ci*term3 - dsr5*dix*yr + diy*term4 + dir*term5
     &                      + 2.0d0*dsr5*qixy - 2.0d0*dsr7*yr*qix
     &                      + 2.0d0*qiy*term7 + qir*term6
                  tkxy = ck*term3 + dsr5*dkx*yr - dky*term4 - dkr*term5
     &                      + 2.0d0*dsr5*qkxy - 2.0d0*dsr7*yr*qkx
     &                      + 2.0d0*qky*term7 + qkr*term6
                  term3 = term1*xr*zr - rr3*zr*drc3(1)
                  term5 = term2*xr*zr - rr5*zr*drc5(1)
                  term6 = (dmpe(9)-dsc7*rr9)*xr*zr - rr7*zr*drc7(1)
                  tixz = ci*term3 - dsr5*dix*zr + diz*term4 + dir*term5
     &                      + 2.0d0*dsr5*qixz - 2.0d0*dsr7*zr*qix
     &                      + 2.0d0*qiz*term7 + qir*term6
                  tkxz = ck*term3 + dsr5*dkx*zr - dkz*term4 - dkr*term5
     &                      + 2.0d0*dsr5*qkxz - 2.0d0*dsr7*zr*qkx
     &                      + 2.0d0*qkz*term7 + qkr*term6
                  term3 = term1*yr*zr - rr3*zr*drc3(2)
                  term4 = rr3*drc3(2) - term1*yr
                  term5 = term2*yr*zr - rr5*zr*drc5(2)
                  term6 = (dmpe(9)-dsc7*rr9)*yr*zr - rr7*zr*drc7(2)
                  term7 = rr5*drc5(2) - term2*yr
                  tiyz = ci*term3 - dsr5*diy*zr + diz*term4 + dir*term5
     &                      + 2.0d0*dsr5*qiyz - 2.0d0*dsr7*zr*qiy
     &                      + 2.0d0*qiz*term7 + qir*term6
                  tkyz = ck*term3 + dsr5*dky*zr - dkz*term4 - dkr*term5
     &                      + 2.0d0*dsr5*qkyz - 2.0d0*dsr7*zr*qky
     &                      + 2.0d0*qkz*term7 + qkr*term6
                  depx = tixx*ukxp + tixy*ukyp + tixz*ukzp
     &                      - tkxx*uixp - tkxy*uiyp - tkxz*uizp
                  depy = tixy*ukxp + tiyy*ukyp + tiyz*ukzp
     &                      - tkxy*uixp - tkyy*uiyp - tkyz*uizp
                  depz = tixz*ukxp + tiyz*ukyp + tizz*ukzp
     &                      - tkxz*uixp - tkyz*uiyp - tkzz*uizp
                  frcx = depx
                  frcy = depy
                  frcz = depz
c
c     get the dEp/dR terms used for direct polarization force
c
                  term1 = dmpe(5) - psc3*rr5
                  term2 = dmpe(7) - psc5*rr7
                  term3 = -psr3 + term1*xr*xr - rr3*xr*prc3(1)
                  term4 = rr3*prc3(1) - term1*xr - psr5*xr
                  term5 = term2*xr*xr - psr5 - rr5*xr*prc5(1)
                  term6 = (dmpe(9)-psc7*rr9)*xr*xr - dmpe(7)
     &                       - rr7*xr*prc7(1)
                  term7 = rr5*prc5(1) - 2.0d0*dmpe(7)*xr
     &                       + (psc5+1.5d0*psc7)*rr7*xr
                  tixx = ci*term3 + dix*term4 + dir*term5
     &                      + 2.0d0*psr5*qixx + (qiy*yr+qiz*zr)*psc7*rr7
     &                      + 2.0d0*qix*term7 + qir*term6
                  tkxx = ck*term3 - dkx*term4 - dkr*term5
     &                      + 2.0d0*psr5*qkxx + (qky*yr+qkz*zr)*psc7*rr7
     &                      + 2.0d0*qkx*term7 + qkr*term6
                  term3 = -psr3 + term1*yr*yr - rr3*yr*prc3(2)
                  term4 = rr3*prc3(2) - term1*yr - psr5*yr
                  term5 = term2*yr*yr - psr5 - rr5*yr*prc5(2)
                  term6 = (dmpe(9)-psc7*rr9)*yr*yr - dmpe(7)
     &                       - rr7*yr*prc7(2)
                  term7 = rr5*prc5(2) - 2.0d0*dmpe(7)*yr
     &                       + (psc5+1.5d0*psc7)*rr7*yr
                  tiyy = ci*term3 + diy*term4 + dir*term5
     &                      + 2.0d0*psr5*qiyy + (qix*xr+qiz*zr)*psc7*rr7
     &                      + 2.0d0*qiy*term7 + qir*term6
                  tkyy = ck*term3 - dky*term4 - dkr*term5
     &                      + 2.0d0*psr5*qkyy + (qkx*xr+qkz*zr)*psc7*rr7
     &                      + 2.0d0*qky*term7 + qkr*term6
                  term3 = -psr3 + term1*zr*zr - rr3*zr*prc3(3)
                  term4 = rr3*prc3(3) - term1*zr - psr5*zr
                  term5 = term2*zr*zr - psr5 - rr5*zr*prc5(3)
                  term6 = (dmpe(9)-psc7*rr9)*zr*zr - dmpe(7)
     &                       - rr7*zr*prc7(3)
                  term7 = rr5*prc5(3) - 2.0d0*dmpe(7)*zr
     &                       + (psc5+1.5d0*psc7)*rr7*zr
                  tizz = ci*term3 + diz*term4 + dir*term5
     &                      + 2.0d0*psr5*qizz + (qix*xr+qiy*yr)*psc7*rr7
     &                      + 2.0d0*qiz*term7 + qir*term6
                  tkzz = ck*term3 - dkz*term4 - dkr*term5
     &                      + 2.0d0*psr5*qkzz + (qkx*xr+qky*yr)*psc7*rr7
     &                      + 2.0d0*qkz*term7 + qkr*term6
                  term3 = term1*xr*yr - rr3*yr*prc3(1)
                  term4 = rr3*prc3(1) - term1*xr
                  term5 = term2*xr*yr - rr5*yr*prc5(1)
                  term6 = (dmpe(9)-psc7*rr9)*xr*yr - rr7*yr*prc7(1)
                  term7 = rr5*prc5(1) - term2*xr
                  tixy = ci*term3 - psr5*dix*yr + diy*term4 + dir*term5
     &                      + 2.0d0*psr5*qixy - 2.0d0*psr7*yr*qix
     &                      + 2.0d0*qiy*term7 + qir*term6
                  tkxy = ck*term3 + psr5*dkx*yr - dky*term4 - dkr*term5
     &                      + 2.0d0*psr5*qkxy - 2.0d0*psr7*yr*qkx
     &                      + 2.0d0*qky*term7 + qkr*term6
                  term3 = term1*xr*zr - rr3*zr*prc3(1)
                  term5 = term2*xr*zr - rr5*zr*prc5(1)
                  term6 = (dmpe(9)-psc7*rr9)*xr*zr - rr7*zr*prc7(1)
                  tixz = ci*term3 - psr5*dix*zr + diz*term4 + dir*term5
     &                      + 2.0d0*psr5*qixz - 2.0d0*psr7*zr*qix
     &                      + 2.0d0*qiz*term7 + qir*term6
                  tkxz = ck*term3 + psr5*dkx*zr - dkz*term4 - dkr*term5
     &                      + 2.0d0*psr5*qkxz - 2.0d0*psr7*zr*qkx
     &                      + 2.0d0*qkz*term7 + qkr*term6
                  term3 = term1*yr*zr - rr3*zr*prc3(2)
                  term4 = rr3*prc3(2) - term1*yr
                  term5 = term2*yr*zr - rr5*zr*prc5(2)
                  term6 = (dmpe(9)-psc7*rr9)*yr*zr - rr7*zr*prc7(2)
                  term7 = rr5*prc5(2) - term2*yr
                  tiyz = ci*term3 - psr5*diy*zr + diz*term4 + dir*term5
     &                      + 2.0d0*psr5*qiyz - 2.0d0*psr7*zr*qiy
     &                      + 2.0d0*qiz*term7 + qir*term6
                  tkyz = ck*term3 + psr5*dky*zr - dkz*term4 - dkr*term5
     &                      + 2.0d0*psr5*qkyz - 2.0d0*psr7*zr*qky
     &                      + 2.0d0*qkz*term7 + qkr*term6
                  depx = tixx*ukx + tixy*uky + tixz*ukz
     &                      - tkxx*uix - tkxy*uiy - tkxz*uiz
                  depy = tixy*ukx + tiyy*uky + tiyz*ukz
     &                      - tkxy*uix - tkyy*uiy - tkyz*uiz
                  depz = tixz*ukx + tiyz*uky + tizz*ukz
     &                      - tkxz*uix - tkyz*uiy - tkzz*uiz
                  frcx = frcx + depx
                  frcy = frcy + depy
                  frcz = frcz + depz
c
c     get the field gradient for direct polarization force
c
               else if (use_chgpen) then
                  term1i = rr3i - rr5i*xr*xr
                  term1core = rr3core - rr5core*xr*xr
                  term2i = 2.0d0*rr5i*xr 
                  term3i = rr7i*xr*xr - rr5i
                  term4i = 2.0d0*rr5i
                  term5i = 5.0d0*rr7i*xr
                  term6i = rr9i*xr*xr
                  term1k = rr3k - rr5k*xr*xr
                  term2k = 2.0d0*rr5k*xr
                  term3k = rr7k*xr*xr - rr5k
                  term4k = 2.0d0*rr5k
                  term5k = 5.0d0*rr7k*xr
                  term6k = rr9k*xr*xr
                  tixx = vali*term1i + corei*term1core  
     &                      + dix*term2i - dir*term3i
     &                      - qixx*term4i + qix*term5i - qir*term6i
     &                      + (qiy*yr+qiz*zr)*rr7i
                  tkxx = valk*term1k + corek*term1core
     &                      - dkx*term2k + dkr*term3k
     &                      - qkxx*term4k + qkx*term5k - qkr*term6k
     &                      + (qky*yr+qkz*zr)*rr7k
                  term1i = rr3i - rr5i*yr*yr
                  term1core = rr3core - rr5core*yr*yr
                  term2i = 2.0d0*rr5i*yr
                  term3i = rr7i*yr*yr - rr5i
                  term4i = 2.0d0*rr5i
                  term5i = 5.0d0*rr7i*yr
                  term6i = rr9i*yr*yr
                  term1k = rr3k - rr5k*yr*yr
                  term2k = 2.0d0*rr5k*yr
                  term3k = rr7k*yr*yr - rr5k
                  term4k = 2.0d0*rr5k
                  term5k = 5.0d0*rr7k*yr
                  term6k = rr9k*yr*yr
                  tiyy = vali*term1i + corei*term1core
     &                      + diy*term2i - dir*term3i
     &                      - qiyy*term4i + qiy*term5i - qir*term6i
     &                      + (qix*xr+qiz*zr)*rr7i
                  tkyy = valk*term1k + corek*term1core
     &                      - dky*term2k + dkr*term3k
     &                      - qkyy*term4k + qky*term5k - qkr*term6k
     &                      + (qkx*xr+qkz*zr)*rr7k
                  term1i = rr3i - rr5i*zr*zr
                  term1core = rr3core - rr5core*zr*zr
                  term2i = 2.0d0*rr5i*zr
                  term3i = rr7i*zr*zr - rr5i
                  term4i = 2.0d0*rr5i
                  term5i = 5.0d0*rr7i*zr
                  term6i = rr9i*zr*zr
                  term1k = rr3k - rr5k*zr*zr
                  term2k = 2.0d0*rr5k*zr
                  term3k = rr7k*zr*zr - rr5k
                  term4k = 2.0d0*rr5k
                  term5k = 5.0d0*rr7k*zr
                  term6k = rr9k*zr*zr
                  tizz = vali*term1i + corei*term1core
     &                      + diz*term2i - dir*term3i
     &                      - qizz*term4i + qiz*term5i - qir*term6i
     &                      + (qix*xr+qiy*yr)*rr7i
                  tkzz = valk*term1k + corek*term1core
     &                      - dkz*term2k + dkr*term3k
     &                      - qkzz*term4k + qkz*term5k - qkr*term6k
     &                      + (qkx*xr+qky*yr)*rr7k
                  term2i = rr5i*xr 
                  term1i = yr * term2i
                  term1core = rr5core*xr*yr
                  term3i = rr5i*yr
                  term4i = yr * (rr7i*xr)
                  term5i = 2.0d0*rr5i
                  term6i = 2.0d0*rr7i*xr
                  term7i = 2.0d0*rr7i*yr
                  term8i = yr*rr9i*xr
                  term2k = rr5k*xr
                  term1k = yr * term2k
                  term3k = rr5k*yr
                  term4k = yr * (rr7k*xr)
                  term5k = 2.0d0*rr5k
                  term6k = 2.0d0*rr7k*xr
                  term7k = 2.0d0*rr7k*yr
                  term8k = yr*rr9k*xr
                  tixy = -vali*term1i - corei*term1core 
     &                      + diy*term2i + dix*term3i
     &                      - dir*term4i - qixy*term5i + qiy*term6i
     &                      + qix*term7i - qir*term8i
                  tkxy = -valk*term1k - corek*term1core 
     &                      - dky*term2k - dkx*term3k
     &                      + dkr*term4k - qkxy*term5k + qky*term6k
     &                      + qkx*term7k - qkr*term8k
                  term2i = rr5i*xr
                  term1i = zr * term2i
                  term1core = rr5core*xr*zr
                  term3i = rr5i*zr
                  term4i = zr * (rr7i*xr)
                  term5i = 2.0d0*rr5i
                  term6i = 2.0d0*rr7i*xr
                  term7i = 2.0d0*rr7i*zr
                  term8i = zr*rr9i*xr
                  term2k = rr5k*xr
                  term1k = zr * term2k
                  term3k = rr5k*zr
                  term4k = zr * (rr7k*xr)
                  term5k = 2.0d0*rr5k
                  term6k = 2.0d0*rr7k*xr
                  term7k = 2.0d0*rr7k*zr
                  term8k = zr*rr9k*xr
                  tixz = -vali*term1i - corei*term1core
     &                      + diz*term2i + dix*term3i
     &                      - dir*term4i - qixz*term5i + qiz*term6i
     &                      + qix*term7i - qir*term8i
                  tkxz = -valk*term1k - corek*term1core
     &                      - dkz*term2k - dkx*term3k
     &                      + dkr*term4k - qkxz*term5k + qkz*term6k
     &                      + qkx*term7k - qkr*term8k
                  term2i = rr5i*yr
                  term1i = zr * term2i
                  term1core = rr5core*yr*zr
                  term3i = rr5i*zr
                  term4i = zr * (rr7i*yr)
                  term5i = 2.0d0*rr5i
                  term6i = 2.0d0*rr7i*yr
                  term7i = 2.0d0*rr7i*zr
                  term8i = zr*rr9i*yr
                  term2k = rr5k*yr
                  term1k = zr * term2k
                  term3k = rr5k*zr
                  term4k = zr * (rr7k*yr)
                  term5k = 2.0d0*rr5k
                  term6k = 2.0d0*rr7k*yr
                  term7k = 2.0d0*rr7k*zr
                  term8k = zr*rr9k*yr
                  tiyz = -vali*term1i - corei*term1core
     &                      + diz*term2i + diy*term3i
     &                      - dir*term4i - qiyz*term5i + qiz*term6i
     &                      + qiy*term7i - qir*term8i
                  tkyz = -valk*term1k - corek*term1core
     &                      - dkz*term2k - dky*term3k
     &                      + dkr*term4k - qkyz*term5k + qkz*term6k
     &                      + qky*term7k - qkr*term8k
                  depx = tixx*ukx + tixy*uky + tixz*ukz
     &                      - tkxx*uix - tkxy*uiy - tkxz*uiz
                  depy = tixy*ukx + tiyy*uky + tiyz*ukz
     &                      - tkxy*uix - tkyy*uiy - tkyz*uiz
                  depz = tixz*ukx + tiyz*uky + tizz*ukz
     &                      - tkxz*uix - tkyz*uiy - tkzz*uiz
                  frcx = -2.0d0 * depx
                  frcy = -2.0d0 * depy
                  frcz = -2.0d0 * depz
               end if
c
c     reset Thole values when alternate direct damping is used
c
               if (use_tholed) then
                  sc3 = 1.0d0
                  sc5 = 1.0d0
                  sc7 = 1.0d0
                  do j = 1, 3
                     rc3(j) = 0.0d0
                     rc5(j) = 0.0d0
                     rc7(j) = 0.0d0
                  end do
                  damp = pdi * pdamp(kkpole)
                  if (damp .ne. 0.0d0) then
                     pgamma = min(pti,thole(kkpole))
                     damp = pgamma * (r/damp)**3
                     if (damp .lt. 50.0d0) then
                        expdamp = exp(-damp)
                        sc3 = 1.0d0 - expdamp
                        sc5 = 1.0d0 - expdamp*(1.0d0+damp)
                        sc7 = 1.0d0 - expdamp*(1.0d0+damp
     &                                   +0.6d0*damp**2)
                        temp3 = 3.0d0 * damp * expdamp / r2
                        temp5 = damp
                        temp7 = -0.2d0 + 0.6d0*damp
                        rc3(1) = xr * temp3
                        rc3(2) = yr * temp3
                        rc3(3) = zr * temp3
                        rc5(1) = rc3(1) * temp5
                        rc5(2) = rc3(2) * temp5
                        rc5(3) = rc3(3) * temp5
                        rc7(1) = rc5(1) * temp7
                        rc7(2) = rc5(2) * temp7
                        rc7(3) = rc5(3) * temp7
                     end if
                  end if
                  usc3 = 1.0d0 - sc3*uscale(kglob)
                  usc5 = 1.0d0 - sc5*uscale(kglob)
                  usr3 = dmpe(3) - usc3*rr3
                  usr5 = dmpe(5) - usc5*rr5
                  do j = 1, 3
                     urc3(j) = rc3(j) * uscale(kglob)
                     urc5(j) = rc5(j) * uscale(kglob)
                  end do
               end if
c
c     get the dtau/dr terms used for mutual polarization force
c
               if (use_thole) then
                  term1 = dmpe(5) - usc3*rr5
                  term2 = dmpe(7) - usc5*rr7
                  term3 = usr5 + term1
                  term4 = rr3 * uscale(kglob)
                  term5 = -xr*term3 + rc3(1)*term4
                  term6 = -usr5 + xr*xr*term2 - rr5*xr*urc5(1)
                  tixx = uix*term5 + uir*term6
                  tkxx = ukx*term5 + ukr*term6
                  term5 = -yr*term3 + rc3(2)*term4
                  term6 = -usr5 + yr*yr*term2 - rr5*yr*urc5(2)
                  tiyy = uiy*term5 + uir*term6
                  tkyy = uky*term5 + ukr*term6
                  term5 = -zr*term3 + rc3(3)*term4
                  term6 = -usr5 + zr*zr*term2 - rr5*zr*urc5(3)
                  tizz = uiz*term5 + uir*term6
                  tkzz = ukz*term5 + ukr*term6
                  term4 = -usr5 * yr
                  term5 = -xr*term1 + rr3*urc3(1)
                  term6 = xr*yr*term2 - rr5*yr*urc5(1)
                  tixy = uix*term4 + uiy*term5 + uir*term6
                  tkxy = ukx*term4 + uky*term5 + ukr*term6
                  term4 = -usr5 * zr
                  term6 = xr*zr*term2 - rr5*zr*urc5(1)
                  tixz = uix*term4 + uiz*term5 + uir*term6
                  tkxz = ukx*term4 + ukz*term5 + ukr*term6
                  term5 = -yr*term1 + rr3*urc3(2)
                  term6 = yr*zr*term2 - rr5*zr*urc5(2)
                  tiyz = uiy*term4 + uiz*term5 + uir*term6
                  tkyz = uky*term4 + ukz*term5 + ukr*term6
                  depx = tixx*ukxp + tixy*ukyp + tixz*ukzp
     &                      + tkxx*uixp + tkxy*uiyp + tkxz*uizp
                  depy = tixy*ukxp + tiyy*ukyp + tiyz*ukzp
     &                      + tkxy*uixp + tkyy*uiyp + tkyz*uizp
                  depz = tixz*ukxp + tiyz*ukyp + tizz*ukzp
     &                      + tkxz*uixp + tkyz*uiyp + tkzz*uizp
                  frcx = frcx + depx
                  frcy = frcy + depy
                  frcz = frcz + depz
c
c     get the dtau/dr terms used for mutual polarization force
c
               else if (use_chgpen) then
                  term1 = 2.0d0 * rr5ik
                  term2 = term1*xr
                  term3 = rr5ik - rr7ik*xr*xr
                  tixx = uix*term2 + uir*term3
                  tkxx = ukx*term2 + ukr*term3
                  term2 = term1*yr
                  term3 = rr5ik - rr7ik*yr*yr
                  tiyy = uiy*term2 + uir*term3
                  tkyy = uky*term2 + ukr*term3
                  term2 = term1*zr
                  term3 = rr5ik - rr7ik*zr*zr
                  tizz = uiz*term2 + uir*term3
                  tkzz = ukz*term2 + ukr*term3
                  term1 = rr5ik*yr
                  term2 = rr5ik*xr
                  term3 = yr * (rr7ik*xr)
                  tixy = uix*term1 + uiy*term2 - uir*term3
                  tkxy = ukx*term1 + uky*term2 - ukr*term3
                  term1 = rr5ik * zr
                  term3 = zr * (rr7ik*xr)
                  tixz = uix*term1 + uiz*term2 - uir*term3
                  tkxz = ukx*term1 + ukz*term2 - ukr*term3
                  term2 = rr5ik*yr
                  term3 = zr * (rr7ik*yr)
                  tiyz = uiy*term1 + uiz*term2 - uir*term3
                  tkyz = uky*term1 + ukz*term2 - ukr*term3
                  depx = tixx*ukxp + tixy*ukyp + tixz*ukzp
     &                      + tkxx*uixp + tkxy*uiyp + tkxz*uizp
                  depy = tixy*ukxp + tiyy*ukyp + tiyz*ukzp
     &                      + tkxy*uixp + tkyy*uiyp + tkyz*uizp
                  depz = tixz*ukxp + tiyz*ukyp + tizz*ukzp
     &                      + tkxz*uixp + tkyz*uiyp + tkzz*uizp
                  frcx = frcx - depx
                  frcy = frcy - depy
                  frcz = frcz - depz
               end if
c
c     increment force-based gradient on the interaction sites
c
               dep(1,i) = dep(1,i) - frcx
               dep(2,i) = dep(2,i) - frcy
               dep(3,i) = dep(3,i) - frcz
               dep(1,kbis) = dep(1,kbis) + frcx
               dep(2,kbis) = dep(2,kbis) + frcy
               dep(3,kbis) = dep(3,kbis) + frcz
c
c     increment the virial due to pairwise Cartesian forces
c
               vxx = xr * frcx
               vxy = 0.5d0 * (yr*frcx+xr*frcy)
               vxz = 0.5d0 * (zr*frcx+xr*frcz)
               vyy = yr * frcy
               vyz = 0.5d0 * (zr*frcy+yr*frcz)
               vzz = zr * frcz
               vir(1,1) = vir(1,1) + vxx
               vir(2,1) = vir(2,1) + vxy
               vir(3,1) = vir(3,1) + vxz
               vir(1,2) = vir(1,2) + vxy
               vir(2,2) = vir(2,2) + vyy
               vir(3,2) = vir(3,2) + vyz
               vir(1,3) = vir(1,3) + vxz
               vir(2,3) = vir(2,3) + vyz
               vir(3,3) = vir(3,3) + vzz
               if (shortrange) then
                 virsave(1,1) = virsave(1,1) + vxx
                 virsave(2,1) = virsave(2,1) + vxy
                 virsave(3,1) = virsave(3,1) + vxz
                 virsave(1,2) = virsave(1,2) + vxy
                 virsave(2,2) = virsave(2,2) + vyy
                 virsave(3,2) = virsave(3,2) + vyz
                 virsave(1,3) = virsave(1,3) + vxz
                 virsave(2,3) = virsave(2,3) + vyz
                 virsave(3,3) = virsave(3,3) + vzz
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         if (dpequal) then
            do j = 1, n12(iglob)
               pscale(i12(j,iglob)) = 1.0d0
               dscale(i12(j,iglob)) = 1.0d0
               wscale(i12(j,iglob)) = 1.0d0
            end do
            do j = 1, n13(iglob)
               pscale(i13(j,iglob)) = 1.0d0
               dscale(i13(j,iglob)) = 1.0d0
               wscale(i13(j,iglob)) = 1.0d0
            end do
            do j = 1, n14(iglob)
               pscale(i14(j,iglob)) = 1.0d0
               dscale(i14(j,iglob)) = 1.0d0
               wscale(i14(j,iglob)) = 1.0d0
            end do
            do j = 1, n15(iglob)
               pscale(i15(j,iglob)) = 1.0d0
               dscale(i15(j,iglob)) = 1.0d0
               wscale(i15(j,iglob)) = 1.0d0
            end do
            do j = 1, np11(iglob)
               uscale(ip11(j,iglob)) = 1.0d0
            end do
            do j = 1, np12(iglob)
               uscale(ip12(j,iglob)) = 1.0d0
            end do
            do j = 1, np13(iglob)
               uscale(ip13(j,iglob)) = 1.0d0
            end do
            do j = 1, np14(iglob)
               uscale(ip14(j,iglob)) = 1.0d0
            end do
         else
            do j = 1, n12(iglob)
               pscale(i12(j,iglob)) = 1.0d0
               wscale(i12(j,iglob)) = 1.0d0
            end do
            do j = 1, n13(iglob)
               pscale(i13(j,iglob)) = 1.0d0
               wscale(i13(j,iglob)) = 1.0d0
            end do
            do j = 1, n14(iglob)
               pscale(i14(j,iglob)) = 1.0d0
               wscale(i14(j,iglob)) = 1.0d0
            end do
            do j = 1, n15(iglob)
               pscale(i15(j,iglob)) = 1.0d0
               wscale(i15(j,iglob)) = 1.0d0
            end do
            do j = 1, np11(iglob)
               dscale(ip11(j,iglob)) = 1.0d0
               uscale(ip11(j,iglob)) = 1.0d0
            end do
            do j = 1, np12(iglob)
               dscale(ip12(j,iglob)) = 1.0d0
               uscale(ip12(j,iglob)) = 1.0d0
            end do
            do j = 1, np13(iglob)
               dscale(ip13(j,iglob)) = 1.0d0
               uscale(ip13(j,iglob)) = 1.0d0
            end do
            do j = 1, np14(iglob)
               dscale(ip14(j,iglob)) = 1.0d0
               uscale(ip14(j,iglob)) = 1.0d0
            end do
         end if
      end do
c
c     torque is induced field and gradient cross permanent moments
c
      do ii = 1, npolelocnl
         iipole = poleglobnl(ii)
         iglob = ipole(iipole)
         i = loc(iglob)
         dix = rpole(2,iipole)
         diy = rpole(3,iipole)
         diz = rpole(4,iipole)
         qixx = rpole(5,iipole)
         qixy = rpole(6,iipole)
         qixz = rpole(7,iipole)
         qiyy = rpole(9,iipole)
         qiyz = rpole(10,iipole)
         qizz = rpole(13,iipole)
         tep(1) = diz*ufld(2,i) - diy*ufld(3,i)
     &               + qixz*dufld(2,i) - qixy*dufld(4,i)
     &               + 2.0d0*qiyz*(dufld(3,i)-dufld(6,i))
     &               + (qizz-qiyy)*dufld(5,i)
         tep(2) = dix*ufld(3,i) - diz*ufld(1,i)
     &               - qiyz*dufld(2,i) + qixy*dufld(5,i)
     &               + 2.0d0*qixz*(dufld(6,i)-dufld(1,i))
     &               + (qixx-qizz)*dufld(4,i)
         tep(3) = diy*ufld(1,i) - dix*ufld(2,i)
     &               + qiyz*dufld(4,i) - qixz*dufld(5,i)
     &               + 2.0d0*qixy*(dufld(1,i)-dufld(3,i))
     &               + (qiyy-qixx)*dufld(2,i)
         call torque (iipole,tep,fix,fiy,fiz,dep)
         iz = zaxis(iipole)
         ix = xaxis(iipole)
         iy = abs(yaxis(iipole))
         if (iz .eq. 0)  iz = iglob
         if (ix .eq. 0)  ix = iglob
         if (iy .eq. 0)  iy = iglob
         xiz = x(iz) - x(iglob)
         yiz = y(iz) - y(iglob)
         ziz = z(iz) - z(iglob)
         xix = x(ix) - x(iglob)
         yix = y(ix) - y(iglob)
         zix = z(ix) - z(iglob)
         xiy = x(iy) - x(iglob)
         yiy = y(iy) - y(iglob)
         ziy = z(iy) - z(iglob)
         vxx = xix*fix(1) + xiy*fiy(1) + xiz*fiz(1)
         vxy = 0.5d0 * (yix*fix(1) + yiy*fiy(1) + yiz*fiz(1)
     &                    + xix*fix(2) + xiy*fiy(2) + xiz*fiz(2))
         vxz = 0.5d0 * (zix*fix(1) + ziy*fiy(1) + ziz*fiz(1)
     &                    + xix*fix(3) + xiy*fiy(3) + xiz*fiz(3))
         vyy = yix*fix(2) + yiy*fiy(2) + yiz*fiz(2)
         vyz = 0.5d0 * (zix*fix(2) + ziy*fiy(2) + ziz*fiz(2)
     &                    + yix*fix(3) + yiy*fiy(3) + yiz*fiz(3))
         vzz = zix*fix(3) + ziy*fiy(3) + ziz*fiz(3)
         vir(1,1) = vir(1,1) + vxx
         vir(2,1) = vir(2,1) + vxy
         vir(3,1) = vir(3,1) + vxz
         vir(1,2) = vir(1,2) + vxy
         vir(2,2) = vir(2,2) + vyy
         vir(3,2) = vir(3,2) + vyz
         vir(1,3) = vir(1,3) + vxz
         vir(2,3) = vir(2,3) + vyz
         vir(3,3) = vir(3,3) + vzz
      end do
c
c
c     modify the gradient and virial for charge flux
c
c
c
      if (use_chgflx) then
c     first communicate potential
         call commpotsum(pot,0)        
         call commpot(pot,1)
         call dcflux (pot,decfx,decfy,decfz)
         do ii = 1, npolebloc
            iipole = poleglob(ii)
            iglob = ipole(iipole)
            i = loc(iglob)
            xi = x(iglob)
            yi = y(iglob)
            zi = z(iglob)
            frcx = decfx(i)
            frcy = decfy(i)
            frcz = decfz(i)
            dep(1,i) = dep(1,i) + frcx
            dep(2,i) = dep(2,i) + frcy
            dep(3,i) = dep(3,i) + frcz
            vxx = xi * frcx
            vxy = yi * frcx
            vxz = zi * frcx
            vyy = yi * frcy
            vyz = zi * frcy
            vzz = zi * frcz
            vir(1,1) = vir(1,1) + vxx
            vir(2,1) = vir(2,1) + vxy
            vir(3,1) = vir(3,1) + vxz
            vir(1,2) = vir(1,2) + vxy
            vir(2,2) = vir(2,2) + vyy
            vir(3,2) = vir(3,2) + vyz
            vir(1,3) = vir(1,3) + vxz
            vir(2,3) = vir(2,3) + vyz
            vir(3,3) = vir(3,3) + vzz
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (pscale)
      deallocate (dscale)
      deallocate (uscale)
      deallocate (wscale)
      deallocate (ufld)
      deallocate (dufld)
      deallocate (pot)
      deallocate (decfx)
      deallocate (decfy)
      deallocate (decfz)
      return
      end
