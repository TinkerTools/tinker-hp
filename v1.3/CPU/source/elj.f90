!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ##############################################################
!     ##                                                          ##
!     ##  subroutine elj  --  Lennard-Jones van der Waals energy  ##
!     ##                                                          ##
!     ##############################################################
!
!
!     "elj" calculates the Lennard-Jones 6-12 van der Waals energy
!
!
!> @brief 
!> calculates the Lennard-Jones 6-12 van der Waals energy
!> @param no params
subroutine elj
   use energi
   use inform
   use iounit
   use potent
   use vdwpot
   implicit none
   real*8 elrc
   character*11 mode
!
   if (deb_Path) write(iout,*), 'elj '
!
!     Sum over pairwise interactions
!
   call elj0c
!
!     apply long range van der Waals correction if desired
!
   if (use_vcorr) then
      mode = 'VDW'
      call evcorr (mode,elrc)
      ev = ev + elrc
   end if
   return
end
!
!     ###############################################################
!     ##                                                           ##
!     ##  subroutine elj0c  --  Lennard-Jones vdw energy via list  ##
!     ##                                                           ##
!     ###############################################################
!
!
!     "elj0c" calculates the Lennard-Jones 6-12 van der Waals energy
!     using a pairwise neighbor list
!
!     if shortrange, calculates just the long range part
!     if longrange , calculates just the short range part
!
!
!> @brief 
!> calculates the Lennard-Jones 6-12 van der Waals energy
!> using a pairwise neighbor list
!> @param no params
subroutine elj0c
   use atmlst
   use atoms
   use bound
   use couple
   use cutoff
   use domdec
   use energi
   use group
   use inform
   use inter
   use iounit
   use molcul
   use mutant
   use neigh
   use potent
   use shunt
   use usage
   use vdw
   use vdwpot
   use mpi
   implicit none
   integer i,j,iglob,iivdw,kglob,kbis,nnvlst
   integer ii,iv,it
   integer kk,kv,kt
   integer, allocatable :: iv14(:)
   real*8 e,eps,rdn,fgrp
   real*8 rv
   real*8 p6,p12
   real*8 xi,yi,zi
   real*8 xr,yr,zr
   real*8 rik,rik2,rik3
   real*8 rik4,rik5,taper
   real*8 rho,rhok,rvk,galpha,glamb1,gsc,evdw,lambdavt
   real*8 s,ds,vdwshortcut2
   real*8, allocatable :: xred(:)
   real*8, allocatable :: yred(:)
   real*8, allocatable :: zred(:)
   real*8, allocatable :: vscale(:)
   logical usei,proceed,muti,mutk,mutik
   logical testcut,shortrange,longrange,fullrange
   character*11 mode
   character*80 :: RoutineName
!
   if (deb_Path) write(iout,*), 'elj0c '
!
!
!     choose the method for summing over pairwise interactions
   shortrange = use_vdwshort
   longrange  = use_vdwlong
   fullrange  = .not.(shortrange.or.longrange)

   if (shortrange) then
      RoutineName = 'eljshort0c'
      mode        = 'SHORTVDW'
   else if (longrange) then
      RoutineName = 'eljlong0c'
      mode        = 'VDW'
   else
      RoutineName = 'elj0c'
      mode        = 'VDW'
   endif
!
!     zero out the van der Waals energy contribution
!
   ev = 0.0d0
!
!     perform dynamic allocation of some local arrays
!
   allocate (iv14(n))
   allocate (xred(nbloc))
   allocate (yred(nbloc))
   allocate (zred(nbloc))
   allocate (vscale(n))
!
!     set arrays needed to scale connected atom interactions
!
   vscale = 1.0d0
   iv14 = 0
!
!     set the coefficients for the switching function
!
   call switch (mode)
   vdwshortcut2 = (vdwshortcut-shortheal)**2
!
!     apply any reduction factor to the atomic coordinates
!
   do ii = 1, nvdwbloc
      iivdw = vdwglob(ii)
      iglob = ivdw(iivdw)
      i = loc(iglob)
      iv = ired(iglob)
      rdn = kred(iglob)
      xr = x(iglob) - x(iv)
      yr = y(iglob) - y(iv)
      zr = z(iglob) - z(iv)
      if (use_polymer) call image(xr,yr,zr)
      xred(i) = rdn*xr + x(iv)
      yred(i) = rdn*yr + y(iv)
      zred(i) = rdn*zr + z(iv)
   end do
!
!     find van der Waals energy via neighbor list search
!
   do ii = 1, nvdwlocnl
      iivdw = vdwglobnl(ii)
      iglob = ivdw(iivdw)
      i     = loc(iglob)
      muti  = mut(iglob)
      iv    = ired(iglob)
      it    = jvdw(iglob)
      xi    = xred(i)
      yi    = yred(i)
      zi    = zred(i)
      usei = (use(iglob) .or. use(iv))
!
!     set interaction scaling coefficients for connected atoms
!
      do j = 1, n12(iglob)
         vscale(i12(j,iglob)) = v2scale
      end do
      do j = 1, n13(iglob)
         vscale(i13(j,iglob)) = v3scale
      end do
      do j = 1, n14(iglob)
         vscale(i14(j,iglob)) = v4scale
         iv14(i14(j,iglob)) = iglob
      end do
      do j = 1, n15(iglob)
         vscale(i15(j,iglob)) = v5scale
      end do
!
!     decide whether to compute the current interaction
!

      if (shortrange) then
         nnvlst = nshortvlst(ii)
      else
         nnvlst = nvlst(ii)
      end if
      do kk = 1, nnvlst
         if (shortrange) then
            kglob = shortvlst(kk,ii)
         else
            kglob = vlst(kk,ii)
         end if
         if (use_group)  call groups(fgrp,iglob,kglob,0,0,0,0)
         kbis = loc(kglob)
         kv = ired(kglob)
         mutk = mut(kglob)
         proceed = (usei .or. use(kglob) .or. use(kv))
         if (.not.proceed) cycle
!
!     compute the energy contribution for this interaction
!

         kt = jvdw(kglob)
         xr = xi - xred(kbis)
         yr = yi - yred(kbis)
         zr = zi - zred(kbis)

         if (use_bounds) call image (xr,yr,zr)
         rik2 = xr*xr + yr*yr + zr*zr
!
!     check for an interaction distance less than the cutoff
!
         testcut = merge(rik2 .le. off2.and.rik2.ge.vdwshortcut2,&
         &rik2 .le. off2,&
         &longrange&
         &)
         if (testcut) then
            rik = sqrt(rik2)
            rv  = radmin(kt,it)
            eps = epsilon(kt,it)
            if (iv14(kglob) .eq. iglob) then
               rv = radmin4(kt,it)
               eps = epsilon4(kt,it)
            end if
            eps = eps * vscale(kglob)
!
!     set use of lambda scaling for decoupling or annihilation
!
            mutik = .false.
            if (muti .or. mutk) then
               if (vcouple .eq. 1) then
                  mutik = .true.
               else if (.not.muti .or. .not.mutk) then
                  mutik = .true.
               end if
            end if

            if (mutik) then
               rho = rik/rv
               rhok = rho**sck
               rvk  = rv**sck
               galpha = scalpha/(2d0**(sck/6d0))
               glamb1 = 1.0d0-vlambda
!                Softcore expression
               gsc = rv *(galpha * glamb1**scs + rhok )**(1d0/sck)

               evdw = eps*(rv**12/(gsc**12)-2d0*rv**6/(gsc**6))

               lambdavt =  vlambda ** sct
               e = lambdavt*evdw
            else
               p6 = rv**6 / rik2**3
               p12 = p6 * p6
               e = eps * (p12-2.0d0*p6)
            end if

!
!     use energy switching if near the cutoff distance at short range
!
            if(longrange.or.fullrange) then
               if (rik2 .gt. cut2) then
                  rik3  = rik2 * rik
                  rik4  = rik2 * rik2
                  rik5  = rik2 * rik3
                  taper =  c5*rik5 + c4*rik4 + c3*rik3&
                  &+ c2*rik2 + c1*rik + c0
                  e = e * taper
               end if
            endif
!
!     scale the interaction based on its group membership
!
            if (use_group)  e = e * fgrp

            if(shortrange .or. longrange)&
            &call switch_respa(rik,vdwshortcut,shortheal,s,ds)

!
!     increment the overall van der Waals energy components
!
            if(shortrange) then
               e  =   e * s
            else if(longrange) then
               e  =   e * (1.0d0 - s)
            endif

            if (e .ne. 0.0d0) then
               ev = ev + e
            end if
!
!     increment the total intermolecular energy
!
            if (molcule(iglob) .ne. molcule(kglob)) then
               einter = einter + e
            end if
!
         end if
      end do
!
!     reset interaction scaling coefficients for connected atoms
!
      do j = 1, n12(iglob)
         vscale(i12(j,iglob)) = 1.0d0
      end do
      do j = 1, n13(iglob)
         vscale(i13(j,iglob)) = 1.0d0
      end do
      do j = 1, n14(iglob)
         vscale(i14(j,iglob)) = 1.0d0
      end do
      do j = 1, n15(iglob)
         vscale(i15(j,iglob)) = 1.0d0
      end do
   end do
!
!     perform deallocation of some local arrays
!
   deallocate (iv14)
   deallocate (xred)
   deallocate (yred)
   deallocate (zred)
   deallocate (vscale)

   return
end
