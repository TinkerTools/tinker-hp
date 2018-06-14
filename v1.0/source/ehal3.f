c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine ehal3  --  buffered 14-7 vdw energy & analysis  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "ehal3" calculates the buffered 14-7 van der Waals energy
c     and partitions the energy among the atoms
c
c
      subroutine ehal3
      implicit none
      integer i
      real*8 elrc,aelrc
      real*8 time0,time1
      include 'sizes.i'
      include 'analyz.i'
      include 'atoms.i'
      include 'cutoff.i'
      include 'energi.i'
      include 'inform.i'
      include 'iounit.i'
      include 'vdwpot.i'
      include 'openmp.i'
      include 'mpif.h'
c
c
c     choose the method for summing over pairwise interactions
c
      time0 = mpi_wtime()
      call ehal3c
c      call ehal3cvec
      time1 = mpi_wtime()
c      write(*,*) 'time vdw = ',time1-time0
c
c     apply long range van der Waals correction if desired
c
      if (use_vcorr.and.rank.eq.0) then
         call evcorr (elrc)
         ev = ev + elrc
         aelrc = elrc / dble(n)
         do i = 1, nbloc
            aev(i) = aev(i) + aelrc
         end do
         if (verbose .and. elrc.ne.0.0d0) then
            write (iout,10)  elrc
   10       format (/,' Long Range vdw Correction :',9x,f12.4)
         end if
      end if
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine ehal3c  --  buffered 14-7 analysis via list  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "ehal3c" calculates the buffered 14-7 van der Waals energy
c     and also partitions the energy among the atoms using a
c     pairwise neighbor list
c
c
      subroutine ehal3c
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'analyz.i'
      include 'atmtyp.i'
      include 'atmlst.i'
      include 'atoms.i'
      include 'bound.i'
      include 'couple.i'
      include 'energi.i'
      include 'group.i'
      include 'inform.i'
      include 'inter.i'
      include 'iounit.i'
      include 'molcul.i'
      include 'mutant.i'
      include 'neigh.i'
      include 'shunt.i'
      include 'usage.i'
      include 'vdw.i'
      include 'vdwpot.i'
      include 'openmp.i'
      integer i,j,k,iglob,kglob,kbis
      integer ii,iv,it,iivdw,inl
      integer kk,kv,kt
      integer nevt
      integer, allocatable :: iv14(:)
      real*8 e,eps,rdn
      real*8 fgrp,rv,rv7
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 rho,tau,taper
      real*8 scal,t1,t2
      real*8 rik,rik2,rik3
      real*8 rik4,rik5,rik7
      real*8 evt,eintert
      real*8, allocatable :: xred(:)
      real*8, allocatable :: yred(:)
      real*8, allocatable :: zred(:)
      real*8, allocatable :: vscale(:)
      real*8, allocatable :: aevt(:)
      logical proceed,usei
      logical muti,mutk
      logical header,huge
      character*6 mode
c
c
c     zero out the van der Waals energy and partitioning terms
c
      nev = 0
      ev = 0.0d0
      aev = 0.0d0
      header = .true.
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
      vscale = 1.0d0
      iv14 = 0
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
         muti = mut(iglob)
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
            mutk = mut(kglob)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
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
                  rik = sqrt(rik2)
                  rv  = radmin(kt,it)
                  eps = epsilon(kt,it)
                  if (iv14(kglob) .eq. iglob) then
                     rv = radmin4(kt,it)
                     eps = epsilon4(kt,it)
                  end if
                  eps = eps * vscale(kglob)
c
c     get the interaction energy, via soft core if necessary
c
                  if ((muti .and. .not.mutk) .or.
     &                (mutk .and. .not.muti)) then
                     rho = rik / rv
                     eps = eps * vlambda**scexp
                     scal = scalpha * (1.0d0-vlambda)**2
                     t1 = (1.0d0+dhal)**7 / (scal+(rho+dhal)**7)
                     t2 = (1.0d0+ghal) / (scal+rho**7+ghal)
                     e = eps * t1 * (t2-2.0d0)
                  else
                     rv7 = rv**7
                     rik7 = rik**7
                     rho = rik7 + ghal*rv7
                     tau = (dhal+1.0d0) / (rik + dhal*rv)
                     e = eps * rv7 * tau**7
     &                      * ((ghal+1.0d0)*rv7/rho-2.0d0)
                  end if
c
c     use energy switching if near the cutoff distance
c
                  if (rik2 .gt. cut2) then
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
                  if (e .ne. 0.0d0) then
                     nev = nev + 1
                     ev = ev + e
                     aev(i) = aev(i) + 0.5d0*e
                     aev(kbis) = aev(kbis) + 0.5d0*e
                  end if
c
c     increment the total intermolecular energy
c
                  if (molcule(iglob) .ne. molcule(kglob)) then
                     einter = einter + e
                  end if
c
c     print a message if the energy of this interaction is large
c
                  huge = (e .gt. 10.0d0)
                  if ((debug.and.e.ne.0.0d0)
     &                  .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,10)
   10                   format (/,' Individual van der Waals',
     &                             ' Interactions :',
     &                          //,' Type',14x,'Atom Names',
     &                             20x,'Minimum',4x,'Actual',
     &                             6x,'Energy',/)
                     end if
                     write (iout,20)  iglob,name(iglob),kglob,
     $                     name(kglob),
     &                                rv,sqrt(rik2),e
   20                format (' VDW-Hal',3x,2(i7,'-',a3),
     &                          13x,2f10.4,f12.4)
                  end if
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
c     perform deallocation of some local arrays
c
      deallocate (iv14)
      deallocate (xred)
      deallocate (yred)
      deallocate (zred)
      deallocate (vscale)
      return
      end
c
c     "ehal3cvec" calculates the buffered 14-7 van der Waals energy
c     and also partitions the energy among the atoms using a
c     pairwise neighbor list
c
c
      subroutine ehal3cvec
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'analyz.i'
      include 'atmlst.i'
      include 'atoms.i'
      include 'couple.i'
      include 'energi.i'
      include 'mutant.i'
      include 'neigh.i'
      include 'shunt.i'
      include 'usage.i'
      include 'vdw.i'
      include 'vdwpot.i'
      include 'openmp.i'
      include 'mpif.h'
      integer i,j,k,iglob,kglob,kbis
      integer ii,iv,it,iivdw,inl
      integer kk,kv,kt
      integer nevt
      integer nlvdw
      integer, allocatable :: iv14(:)
      real*8 e,eps,rdn,etemp
      real*8 fgrp,rv,rv7
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 rho,tau,taper
      real*8 scal,t1,t2
      real*8 rik,rik2,rik3
      real*8 rik4,rik5,rik7
      real*8 evt,eintert
      real*8, allocatable :: xred(:)
      real*8, allocatable :: yred(:)
      real*8, allocatable :: zred(:)
      real*8, allocatable :: vscale(:)
      real*8, allocatable :: aevt(:)
      real*8, allocatable :: rik2vec(:), posvec(:,:),epsvec(:),rvvec(:)
      real*8, allocatable :: rik2vec1(:)
      real*8, allocatable :: rv7vec(:),rikvec(:),rhovec(:)
      real*8, allocatable :: tauvec(:),rik7vec(:)
      real*8, allocatable :: rik3vec(:),rik4vec(:),rik5vec(:)
      real*8, allocatable :: tapervec(:)
      integer, allocatable :: kbisvec(:)
      real*8 :: time0,time1,time2,time3,time4,time5,time6
      real*8 :: timeglob0,timeglob1,timeglob2,timeglob3,timeglob4
      real*8 :: timeglob5
      logical proceed,usei
      logical muti,mutk
      logical header,huge
      character*6 mode
c
      timeglob0 = 0.0d0
      timeglob1 = 0.0d0
      timeglob2 = 0.0d0
      timeglob3 = 0.0d0
      timeglob4 = 0.0d0
      timeglob5 = 0.0d0
c
c     zero out the van der Waals energy and partitioning terms
c
      nev = 0
      ev = 0.0d0
      aev = 0.0d0
      header = .true.
c
c     perform dynamic allocation of some local arrays
c
      allocate (iv14(n))
      allocate (xred(nbloc))
      allocate (yred(nbloc))
      allocate (zred(nbloc))
      allocate (vscale(n))
c
      allocate (rik2vec(maxvlst),posvec(3,maxvlst),epsvec(maxvlst))
      allocate (rik2vec1(maxvlst))
      allocate (rvvec(maxvlst),rv7vec(maxvlst),rikvec(maxvlst))
      allocate (rhovec(maxvlst),tauvec(maxvlst),rik7vec(maxvlst))
      allocate (rik3vec(maxvlst),rik4vec(maxvlst),rik5vec(maxvlst))
      allocate (tapervec(maxvlst),kbisvec(maxvlst))
c
c     set arrays needed to scale connected atom interactions
c
      vscale = 1
      iv14 = 0.0d0
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
         muti = mut(iglob)
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
c     fill vectors to compute interactions involving this atom
c
         time0 = mpi_wtime()
         do kk = 1, nvlst(inl)
           kglob = vlst(kk,inl)
           kbis = loc(kglob)
           kv = ired(kglob)
           posvec(1,kk) = xi - xred(kbis)
           posvec(2,kk) = yi - yred(kbis)
           posvec(3,kk) = zi - zred(kbis)
         end do
         time1 = mpi_wtime()
         call imagevec(posvec,nvlst(inl))
         time2 = mpi_wtime()
         timeglob1 = timeglob1 + time2-time0
         rik2vec1 = posvec(1,1:nvlst(inl))*posvec(1,1:nvlst(inl))
     $     + posvec(2,1:nvlst(inl))*posvec(2,1:nvlst(inl))
     $     + posvec(3,1:nvlst(inl))*posvec(3,1:nvlst(inl))
         nlvdw = 0
         time3 = mpi_wtime()
         timeglob2 = timeglob2 + time3-time2
         do kk = 1, nvlst(inl)
           kglob = vlst(kk,inl)
           rik2 = rik2vec1(kk)
           kt = jvdw(kglob)
           if (rik2.le.off2) then
             nlvdw = nlvdw+1
             kbisvec(nlvdw) = loc(kglob)
             rik2vec(nlvdw) = rik2vec1(kk)
             rvvec(nlvdw) = radmin(kt,it)
             epsvec(nlvdw) = epsilon(kt,it)*vscale(kglob)
           end if
         end do
         time4 = mpi_wtime()
         timeglob3 = timeglob3 + time4-time3
         rv7vec(1:nlvdw) = rvvec(1:nlvdw)**7
         rikvec(1:nlvdw) = sqrt(rik2vec(1:nlvdw))
         rik3vec(1:nlvdw) = rik2vec(1:nlvdw)*rikvec(1:nlvdw)
         rik4vec(1:nlvdw) = rik2vec(1:nlvdw)*rik2vec(1:nlvdw)
         rik5vec(1:nlvdw) = rik2vec(1:nlvdw)*rik3vec(1:nlvdw)
         rik7vec(1:nlvdw) = rik2vec(1:nlvdw)*rik5vec(1:nlvdw)
         rhovec(1:nlvdw) = rik7vec(1:nlvdw) + ghal*rv7vec(1:nlvdw)
         tauvec(1:nlvdw) = ((dhal+1.0d0)/(rikvec(1:nlvdw)
     $         + dhal*rvvec(1:nlvdw)))**7
         time5 = mpi_wtime()
         timeglob4 = timeglob4 + time5-time4
c
         tapervec(1:nlvdw) = c5*rik5vec(1:nlvdw)+c4*rik4vec(1:nlvdw)
     $    + c3*rik3vec(1:nlvdw)+c2*rik2vec(1:nlvdw)+c1*rikvec(1:nlvdw)
     $    + c0
c
         nlvdw = 0
         do kk = 1, nvlst(inl)
           kglob = vlst(kk,inl)
           rik2 = rik2vec1(kk)
           if (rik2.le.off2) then
             nlvdw = nlvdw + 1
             if (rik2.le.cut2) then
             tapervec(nlvdw) = 1.0d0
             end if
           end if
         end do
c
         e = sum(tapervec(1:nlvdw)*(epsvec(1:nlvdw)*rv7vec(1:nlvdw)*
     $   tauvec(1:nlvdw)*((ghal+1.0d0)*rv7vec(1:nlvdw)/rhovec(1:nlvdw)
     $   -2.0d0)))
         time6 = mpi_wtime()
         timeglob5 = timeglob5 + time6-time5
c
         nev = nev + nlvdw
         ev = ev + e
         aev(i) = aev(i) + e
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
c     perform deallocation of some local arrays
c
      deallocate (iv14)
      deallocate (xred)
      deallocate (yred)
      deallocate (zred)
      deallocate (vscale)
c      write(*,*) 'timeglob1 = ',timeglob1
c      write(*,*) 'timeglob2 = ',timeglob2
c      write(*,*) 'timeglob3 = ',timeglob3
c      write(*,*) 'timeglob4 = ',timeglob4
c      write(*,*) 'timeglob5 = ',timeglob5
c
      deallocate (rik2vec,posvec,epsvec,rvvec)
      return
      end
