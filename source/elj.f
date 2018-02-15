c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine elj  --  Lennard-Jones van der Waals energy  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "elj" calculates the Lennard-Jones 6-12 van der Waals energy
c
c
      subroutine elj
      implicit none
      real*8 elrc
      include 'sizes.i'
      include 'energi.i'
      include 'vdwpot.i'
c
c
c     choose the method for summing over pairwise interactions
c
      call elj0c
c
c     apply long range van der Waals correction if desired
c
      if (use_vcorr) then
         call evcorr (elrc)
         ev = ev + elrc
      end if
      return
      end
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine elj0c  --  Lennard-Jones vdw energy via list  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "elj0c" calculates the Lennard-Jones 6-12 van der Waals energy
c     using a pairwise neighbor list
c
c
      subroutine elj0c
      implicit none
      include 'sizes.i'
      include 'atmlst.i'
      include 'atoms.i'
      include 'bound.i'
      include 'couple.i'
      include 'energi.i'
      include 'group.i'
      include 'neigh.i'
      include 'shunt.i'
      include 'usage.i'
      include 'vdw.i'
      include 'vdwpot.i'
      include 'openmp.i'
      integer i,j,k,iglob,iivdw,kglob,kbis,inl
      integer ii,iv,it
      integer kk,kv,kt
      integer, allocatable :: iv14(:)
      real*8 e,evt
      real*8 p6,p12,eps
      real*8 rv,rdn,fgrp
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 rik,rik2,rik3
      real*8 rik4,rik5,taper
      real*8, allocatable :: xred(:)
      real*8, allocatable :: yred(:)
      real*8, allocatable :: zred(:)
      real*8, allocatable :: vscale(:)
      logical proceed,usei
      character*6 mode
c
c
c     zero out the van der Waals energy contribution
c
      ev = 0.0d0
c
c     perform dynamic allocation of some local arrays
c
      allocate (iv14(n))
      allocate (xred(nbloc))
      allocate (yred(nbloc))
      allocate (zred(nbloc))
      allocate (vscale(n))
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         vscale(i) = 1.0d0
         iv14(i) = 0
      end do
c
c     set the coefficients for the switching function
c
      mode = 'VDW'
      call switch (mode)
c
c     apply any reduction factor to the atomic coordinates
c
      do ii = 1, nvdwbloc
         iivdw = vdwglob(ii)
         iglob = ivdw(iivdw)
         i = loc(iglob)
         iv = ired(iglob)
         rdn = kred(iglob)
         xred(i) = rdn*(x(iglob)-x(iv)) + x(iv)
         yred(i) = rdn*(y(iglob)-y(iv)) + y(iv)
         zred(i) = rdn*(z(iglob)-z(iv)) + z(iv)
      end do
c
c     transfer global to local copies for OpenMP calculation
c
      evt = ev
c
c     set OpenMP directives for the major loop structure
c
c!$OMP PARALLEL default(private) shared(nvdw,ivdw,ired,kred,
c!$OMP& jvdw,xred,yred,zred,use,nvlst,vlst,n12,n13,n14,n15,
c!$OMP& i12,i13,i14,i15,v2scale,v3scale,v4scale,v5scale,
c!$OMP& use_group,off2,radmin,epsilon,radmin4,epsilon4,
c!$OMP& cut2,c0,c1,c2,c3,c4,c5) firstprivate(vscale,iv14)
c!$OMP& shared(evt)
c!$OMP DO reduction(+:evt) schedule(dynamic)
c
c     find the van der Waals energy via neighbor list search
c
      do ii = 1, nvdwlocnl
         iivdw = vdwglobnl(ii)
         iglob = ivdw(iivdw)
         i = loc(iglob)
         inl = vdwlocnl(iivdw)
         iv = ired(iglob)
         it = jvdw(iglob)
         xi = xred(i)
         yi = yred(i)
         zi = zred(i)
         usei = (use(iglob) .or. use(iv))
c
c     set interaction scaling coefficients for connected atoms
c
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
c
c     decide whether to compute the current interaction
c
         do kk = 1, nvlst(inl)
            kglob = vlst(kk,inl)
            kbis = loc(kglob)
            kv = ired(kglob)
            proceed = .true.
            if (use_group) call groups(proceed,fgrp,iglob,kglob,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(kglob) .or. use(kv))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               kt = jvdw(kglob)
               xr = xi - xred(kbis)
               yr = yi - yred(kbis)
               zr = zi - zred(kbis)
               if (use_bounds) call image (xr,yr,zr)
               rik2 = xr*xr + yr*yr + zr*zr
c
c     check for an interaction distance less than the cutoff
c
               if (rik2 .le. off2) then
                  rv = radmin(kt,it)
                  eps = epsilon(kt,it)
                  if (iv14(kglob) .eq. iglob) then
                     rv = radmin4(kt,it)
                     eps = epsilon4(kt,it)
                  end if
                  eps = eps * vscale(kglob)
                  p6 = rv**6 / rik2**3
                  p12 = p6 * p6
                  e = eps * (p12 - 2.0d0*p6)
c
c     use energy switching if near the cutoff distance
c
                  if (rik2 .gt. cut2) then
                     rik = sqrt(rik2)
                     rik3 = rik2 * rik
                     rik4 = rik2 * rik2
                     rik5 = rik2 * rik3
                     taper = c5*rik5 + c4*rik4 + c3*rik3
     &                          + c2*rik2 + c1*rik + c0
                     e = e * taper
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group)  e = e * fgrp
c
c     increment the overall van der Waals energy components
c
                  evt = evt + e
               end if
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
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
c
c     end OpenMP directives for the major loop structure
c
c!$OMP END DO
c!$OMP END PARALLEL
c
c     transfer local to global copies for OpenMP calculation
c
      ev = evt
c
c     perform deallocation of some local arrays
c
      deallocate (iv14)
      deallocate (xred)
      deallocate (yred)
      deallocate (zred)
      deallocate (vscale)
      return
      end
