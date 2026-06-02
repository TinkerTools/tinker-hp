!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!
!     ##################################################################
!     ##                                                              ##
!     ##  subroutine switch_group_ene  --  group polarization energy  ##
!     ##                                                              ##
!     ##################################################################
!
!
!     "switch_group_ene" calculates the group induced dipole polarization energy:
!      if A and B are two groups and their interactions are scaled by alpha,
!      then the Epol(A-B) = alpha*(Epol(A+B) - Epol(A) - Epol(B)) st:
!      Epol(A+B) = alpha*Epol(A+B) + (1-alpha)*(Epol(A) + Epol(B))
!
!
subroutine switch_group_ene
   use atoms
   use domdec
   use energi
   use group
   use inform
   use iounit
   use mpole
   use mpi
   implicit none
   integer i,igroup,jgroup,iglob,iipole
   real*8 epoltot
   real*8, allocatable :: epol_group(:)
   logical scaledinter
10 format('interactions between group',2x,i5,2x,'and group',2x,i5,&
      &2x,'scaled by ',f16.5)
!
!     return if more than 2 groups
!
   if (ngrp.gt.2) then
      if (rank.eq.0) then
         write(iout,*) 'polarization group scaling only possible&
            &        with 2 groups'
      end if
      return
   end if
!
!     return if no interactions between groups is scaled
!
   scaledinter = .false.
   do igroup = 1, ngrp
      do jgroup = 1, ngrp
         if (wgrp(igroup,jgroup).ne.1.0d0) scaledinter = .true.
      end do
   end do
   if (.not.(scaledinter)) return
!
   allocate (epol_group(ngrp))
!
!     save total unscaled polarization energy
!
   epoltot = ep
!
!     loop on the groups and add their energy
!
   do igroup = 1, ngrp
      natgroup = igrp(2,igroup) - igrp(1,igroup) + 1
      if (allocated(globglobgroup)) deallocate (globglobgroup)
      allocate (globglobgroup(natgroup))
      if (allocated(loclocgroup)) deallocate (loclocgroup)
      allocate (loclocgroup(n))
      if (allocated(globgroup)) deallocate (globgroup)
      allocate (globgroup(natgroup))
      if (allocated(locgroup)) deallocate (locgroup)
      allocate (locgroup(natgroup))
      if (allocated(domlengroup)) deallocate (domlengroup)
      allocate (domlengroup(nproc))
      if (allocated(bufbeggroup)) deallocate (bufbeggroup)
      allocate (bufbeggroup(nproc))
      if (allocated(domlenpolegroup)) deallocate (domlenpolegroup)
      allocate (domlenpolegroup(nproc))
      if (allocated(bufbegpolegroup)) deallocate (bufbegpolegroup)
      allocate (bufbegpolegroup(nproc))
      if (allocated(poleglobgroup)) deallocate (poleglobgroup)
      allocate (poleglobgroup(natgroup))
      if (allocated(polelocgroup)) deallocate (polelocgroup)
      allocate (polelocgroup(natgroup))
      if (allocated(ipolegroup)) deallocate (ipolegroup)
      allocate (ipolegroup(natgroup))
      if (allocated(pollistgroup)) deallocate (pollistgroup)
      allocate (pollistgroup(n))
      do i = 1, natgroup
         iglob = kgrp(igrp(1,igroup)+i-1)
         globglobgroup(i) = iglob
         loclocgroup(iglob) = i
      end do

      npolegroup = 0
      do i = 1, natgroup
         iglob = globglobgroup(i)
         iipole = pollist(iglob)
         if (iipole.eq.0) cycle
         npolegroup = npolegroup + 1
         pollistgroup(i) = npolegroup
         ipolegroup(npolegroup) = i
      end do


      nlocatgroup = 0
      do i = 1, natgroup
         iglob = globglobgroup(i)
         if (repart(iglob).eq.rank) then
            nlocatgroup = nlocatgroup + 1
            globgroup(nlocatgroup) = i
            locgroup(i) = nlocatgroup
         end if
      end do

      npolelocgroup = 0
      do i = 1, npolegroup
         iglob = ipolegroup(i)
         if (repart(iglob).eq.rank) then
            npolelocgroup = npolelocgroup + 1
            poleglobgroup(npolelocgroup) = i
            polelocgroup(i) = npolelocgroup
         end if
      end do
      call epolar3_group
      epol_group(igroup) = epgroup
   end do
!
!     add the energies
!
   ep = wgrp(1,2)*epoltot
   ep = ep + (wgrp(1,1)-wgrp(1,2))*epol_group(1)
   ep = ep + (wgrp(2,2)-wgrp(1,2))*epol_group(2)
   do igroup = 1, ngrp
      do jgroup = 1, ngrp
         if ((rank.eq.0).and.(verbose)) then
            write(iout,10) igroup,jgroup,wgrp(igroup,jgroup)
         end if
      end do
   end do

   deallocate (epol_group)
   return
end
!
!     #############################################################################
!     ##                                                                         ##
!     ##  subroutine switch_group_grad  --  group polarization energy and forces ##
!     ##                                                                         ##
!     #############################################################################
!
!
!     "switch_group_grad" calculates the group induced dipole polarization energy and forces:
!      if A and B are two groups and their interactions are scaled by alpha,
!      then the Epol(A-B) = alpha*(Epol(A+B) - Epol(A) - Epol(B)) st:
!      Epol(A+B) = alpha*Epol(A+B) + (1-alpha)*(Epol(A) + Epol(B))
!
!
subroutine switch_group_grad
   use atoms
   use deriv
   use domdec
   use energi
   use group
   use inform
   use iounit
   use mpole
   use mpi
   implicit none
   integer i,igroup,jgroup,iglob,iipole,iloc,iglobgroup
   real*8 epoltot
   real*8, allocatable :: epol_group(:)
   real*8, allocatable :: depol_group(:,:,:)
   logical scaledinter
10 format('interactions between group',2x,i5,2x,'and group',2x,i5,&
      &2x,'scaled by ',f16.5)
!
!     return if more than 2 groups
!
   if (ngrp.gt.2) then
      if (rank.eq.0) then
         write(iout,*) 'polarization group scaling only possible&
            &        with 2 groups'
      end if
      return
   end if
!
!     return if no interactions between groups is scaled
!
   scaledinter = .false.
   do igroup = 1, ngrp
      do jgroup = 1, ngrp
         if (wgrp(igroup,jgroup).ne.1.0d0) scaledinter = .true.
      end do
   end do
   if (.not.(scaledinter)) return
!
   allocate (epol_group(ngrp))
   allocate (depol_group(3,nloc,ngrp))
!
!     save total unscaled polarization energy
!
   epoltot = ep
   depol_group = 0d0
!
!     loop on the groups and add their energy
!
   do igroup = 1, ngrp
      natgroup = igrp(2,igroup) - igrp(1,igroup) + 1
      if (allocated(globglobgroup)) deallocate (globglobgroup)
      allocate (globglobgroup(natgroup))
      if (allocated(loclocgroup)) deallocate (loclocgroup)
      allocate (loclocgroup(n))
      if (allocated(globgroup)) deallocate (globgroup)
      allocate (globgroup(natgroup))
      if (allocated(locgroup)) deallocate (locgroup)
      allocate (locgroup(natgroup))
      if (allocated(domlengroup)) deallocate (domlengroup)
      allocate (domlengroup(nproc))
      if (allocated(bufbeggroup)) deallocate (bufbeggroup)
      allocate (bufbeggroup(nproc))
      if (allocated(domlenpolegroup)) deallocate (domlenpolegroup)
      allocate (domlenpolegroup(nproc))
      if (allocated(bufbegpolegroup)) deallocate (bufbegpolegroup)
      allocate (bufbegpolegroup(nproc))
      if (allocated(poleglobgroup)) deallocate (poleglobgroup)
      allocate (poleglobgroup(natgroup))
      if (allocated(polelocgroup)) deallocate (polelocgroup)
      allocate (polelocgroup(natgroup))
      if (allocated(ipolegroup)) deallocate (ipolegroup)
      allocate (ipolegroup(natgroup))
      if (allocated(pollistgroup)) deallocate (pollistgroup)
      allocate (pollistgroup(n))
      do i = 1, natgroup
         iglob = kgrp(igrp(1,igroup)+i-1)
         globglobgroup(i) = iglob
         loclocgroup(iglob) = i
      end do

      npolegroup = 0
      do i = 1, natgroup
         iglob = globglobgroup(i)
         iipole = pollist(iglob)
         if (iipole.eq.0) cycle
         npolegroup = npolegroup + 1
         pollistgroup(i) = npolegroup
         ipolegroup(npolegroup) = i
      end do


      nlocatgroup = 0
      do i = 1, natgroup
         iglob = globglobgroup(i)
         if (repart(iglob).eq.rank) then
            nlocatgroup = nlocatgroup + 1
            globgroup(nlocatgroup) = i
            locgroup(i) = nlocatgroup
         end if
      end do

      npolelocgroup = 0
      do i = 1, npolegroup
         iglob = ipolegroup(i)
         if (repart(iglob).eq.rank) then
            npolelocgroup = npolelocgroup + 1
            poleglobgroup(npolelocgroup) = i
            polelocgroup(i) = npolelocgroup
         end if
      end do
      call epolar1_group
      epol_group(igroup) = epgroup
      do i = 1, nlocatgroup
         iglobgroup = globgroup(i)
         iglob = globglobgroup(iglobgroup)
         iloc = loc(iglob)
         depol_group(:,iloc,igroup) = depgroup(:,i)
      end do
   end do
!
!     add the energies and the forces
!
   ep = wgrp(1,2)*epoltot
   dep = wgrp(1,2)*dep
   deprec = wgrp(1,2)*deprec
   dep(:,1:nloc) = dep(:,1:nloc) + (wgrp(1,1)-wgrp(1,2))*&
      &depol_group(:,1:nloc,1)
   ep = ep + (wgrp(1,1)-wgrp(1,2))*epol_group(1)
   dep(:,1:nloc) = dep(:,1:nloc) + (wgrp(2,2)-wgrp(1,2))*&
      &depol_group(:,1:nloc,2)
   ep = ep + (wgrp(2,2)-wgrp(1,2))*epol_group(2)
   do igroup = 1, ngrp
      do jgroup = 1, ngrp
         if ((rank.eq.0).and.(verbose)) then
            write(iout,10) igroup,jgroup,wgrp(igroup,jgroup)
         end if
      end do
   end do

   deallocate (epol_group)
   deallocate (depol_group)
   return
end
