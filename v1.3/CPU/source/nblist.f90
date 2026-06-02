!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###############################################################
!     ##                                                           ##
!     ##  subroutine nblist  --  maintain pairwise neighbor lists  ##
!     ##                                                           ##
!     ###############################################################
!
!
!     "nblist" constructs and maintains nonbonded pair neighbor lists
!     for vdw and electrostatic interactions
!
!
!> @brief 
!> constructs and maintains nonbonded pair neighbor lists
!> for vdw and electrostatic interactions
!> @param[in] istep: index of timestep
subroutine nblist(istep)
   use sizes
   use cutoff
   use domdec
   use inform
   use iounit
   use neigh
   use potent
   use mpi
   implicit none
   integer istep,modnl
   real*8 time0,time1
!
   if (deb_Path) write(iout,*), 'nblist '
!
!
!
!     check number of steps between nl updates
!
   modnl = mod(istep,ineigup)
   nblocloop = merge( nbloc,&
   &(int(nbloc/16)+1)*16,&
   &(mod(nbloc,16).eq.0))
   if (modnl.ne.0) return

   if ((use_clist).or.(use_mlist)) then
      if (allocated(nelst)) deallocate (nelst)
      allocate (nelst(nlocnl))
      if (allocated(elst)) deallocate (elst)
      allocate (elst(maxelst,nlocnl))
      nelst = 0
      elst = 0
   end if
   if ((use_shortclist).or.(use_shortmlist)) then
      if (allocated(nshortelst)) deallocate (nshortelst)
      allocate (nshortelst(nlocnl))
      if (allocated(shortelst)) deallocate (shortelst)
      allocate (shortelst(maxelst,nlocnl))
      nshortelst = 0
      shortelst = 0
   end if
   if ((use_vlist).or.(use_dlist)) then
      if (allocated(nvlst)) deallocate (nvlst)
      allocate (nvlst(nlocnl))
      if (allocated(vlst)) deallocate (vlst)
      allocate (vlst(maxvlst,nlocnl))
      nvlst = 0
      vlst = 0
   end if
   if ((use_shortvlist).or.(use_shortdlist)) then
      if (allocated(nshortvlst)) deallocate (nshortvlst)
      allocate (nshortvlst(nlocnl))
      if (allocated(shortvlst)) deallocate (shortvlst)
      allocate (shortvlst(maxvlst,nlocnl))
      nshortvlst = 0
      shortvlst = 0
   end if
   if ((use_pmecore).and.(rank.gt.ndir-1)) return
!
!
!     build the cells at the beginning and assign the particules to them
!
!
   if (use_shortclist) then
      call clistcell2
   else if (use_clist) then
      call clistcell
   end if
   if (use_shortvlist) then
      call vlistcell2
   else if (use_vlist) then
!        time0 = mpi_wtime()
        call vlistcell
!        time1 = mpi_wtime()
      time0 = mpi_wtime()
!      call vlistcellvec
      time1 = mpi_wtime()
   end if
   if (use_shortdlist) then
      call dlistcell2
   else if (use_dlist) then
      call dlistcell
   end if
   if (use_shortmlist) then
      call mlistcell2
   else if (use_mlist) then
      time0 = mpi_wtime()
      call mlistcell
      time1 = mpi_wtime()
   end if
!
   return
end
!
!    "mlistcell" performs a complete rebuild of the
!     electrostatic neighbor lists for multipoles using linked cells method
!
!> @brief 
!> "mlistcell" performs a complete rebuild of the
!>  electrostatic neighbor lists for multipoles using linked cells method
!> @param no params
subroutine mlistcell
   use sizes
   use atmlst
   use atoms
   use boxes
   use cutoff
   use domdec
   use inform
   use iounit
   use mpole
   use neigh
   use mpi
   implicit none
   integer iglob
   integer i,icell,j,k,nneigloc
   integer ineig,iipole,kkpole
   integer kcell,kglob
   integer ncell_loc
   integer, allocatable :: index(:),indcell_loc(:)
   real*8 xi,yi,zi,xk,yk,zk,r2
   real*8, allocatable :: pos(:,:),r2vec(:)
   real*8 boxedge2
   logical docompute

!
!     check size of the box and cutoff for minimum image convention
!
1000 format('Error in multipole neigbor list: max cutoff ',&
   &'buffer should be less than half one edge of the box')
1010 format('Multipole cutoff = ',F14.3)
1020 format('List buffer      = ',F14.3)
!
   if (deb_Path) write(iout,*), 'mlistcell '
!
   boxedge2 = min(xbox2,ybox2,zbox2)
   if (mbuf2.gt.boxedge2*boxedge2) then
      if (rank.eq.0) then
         write(iout,1000)
         write(iout,1010) mpolecut
         write(iout,1020) lbuffer
      end if
      call fatal
   end if
!
   allocate (index(nbloc))
   allocate (indcell_loc(nbloc))
   allocate(pos(3,nbloc))
   allocate(r2vec(nbloc))
!
!     perform a complete list build
!
   do i = 1, npolelocnl
      iipole = poleglobnl(i)
      iglob  = ipole(iipole)
      icell = repartcell(iglob)
!
!       align data of the local cell and the neighboring ones
!
      ncell_loc = cell_len(icell)
      indcell_loc(1:ncell_loc) =&
      &indcell(bufbegcell(icell):(bufbegcell(icell)+cell_len(icell)-1))
      do ineig = 1, numneigcell(icell)
         kcell = neigcell(ineig,icell)
         indcell_loc(ncell_loc+1:(ncell_loc+cell_len(kcell))) =&
         &indcell(bufbegcell(kcell):(bufbegcell(kcell)+cell_len(kcell)-1))
         ncell_loc = ncell_loc + cell_len(kcell)
      end do
!
!       do the neighbor search
!
      nneigloc = 0
      xi = x(iglob)
      yi = y(iglob)
      zi = z(iglob)
      do k = 1, ncell_loc
         kglob = indcell_loc(k)
         kkpole = pollist(kglob)
!
!   skip atom if it is not in the multipole list
!
         if (kkpole.eq.0) cycle
         if (kglob.le.iglob) cycle
         xk = x(kglob)
         yk = y(kglob)
         zk = z(kglob)
         pos(1,nneigloc+1) = xi - xk
         pos(2,nneigloc+1) = yi - yk
         pos(3,nneigloc+1) = zi - zk
         call midpointimagetriclinic(xi,yi,zi,xk,yk,zk,pos(1,nneigloc+1),&
         &pos(2,nneigloc+1),pos(3,nneigloc+1),docompute)
         if (docompute) then
            nneigloc = nneigloc + 1
            index(nneigloc) = kglob
         end if
      end do
!
!       compute the distances and build the list accordingly
!
      r2vec(1:nneigloc) = pos(1,1:nneigloc)*pos(1,1:nneigloc) +&
      &pos(2,1:nneigloc)*pos(2,1:nneigloc) +&
      &pos(3,1:nneigloc)*pos(3,1:nneigloc)

      j = 0
      do k = 1, nneigloc
         r2 = r2vec(k)
         kglob = index(k)
         if (r2 .le. mbuf2) then
            j = j + 1
            kkpole = pollist(kglob)
            elst(j,i) = kkpole
         end if
      end do
      nelst(i) = j
!
!     check to see if the neighbor list is too long
!
      if (nelst(i) .ge. maxelst) then
         if (rank.eq.0) then
            write (iout,10)
10          format (/,' MBUILD  --  Too many Neighbors;',&
            &' Increase MAXELST')
            call fatal
         end if
      end if
   end do
!
   deallocate (pos)
   deallocate (index)
   deallocate (indcell_loc)
   deallocate (r2vec)
   return
end
!
!    "mlistcell2" performs a complete rebuild of the
!     short range and regular electrostatic neighbor lists for
!     multipoles using linked cells method
!
!> @brief 
!> performs a complete rebuild of the
!> short range and regular electrostatic neighbor lists for
!> multipoles using linked cells method
!> @param no params
subroutine mlistcell2
   use sizes
   use atmlst
   use atoms
   use boxes
   use cutoff
   use domdec
   use inform
   use iounit
   use mpole
   use neigh
   use mpi
   implicit none
   integer iglob
   integer i,icell,j,j1,k,nneigloc
   integer ineig,iipole,kkpole
   integer kcell,kglob
   integer ncell_loc
   integer, allocatable :: index(:),indcell_loc(:)
   real*8 xi,yi,zi,xk,yk,zk,r2
   real*8 mbufbeg2
   real*8 boxedge2
   real*8, allocatable :: pos(:,:),r2vec(:)
   logical docompute
!
!     check size of the box and cutoff for minimum image convention
!
1000 format('Error in multipole neigbor list: max cutoff + ',&
   &'buffer should be less than half one edge of the box')
1010 format('Multipole cutoff = ',F14.3)
1020 format('List buffer      = ',F14.3)
!
   if (deb_Path) write(iout,*), 'mlistcell2 '
!
   boxedge2 = min(xbox2,ybox2,zbox2)
   if (mbuf2.gt.boxedge2*boxedge2) then
      if (rank.eq.0) then
         write(iout,1000)
         write(iout,1010) mpolecut
         write(iout,1020) lbuffer
      end if
      call fatal
   end if
!
   allocate (index(nbloc))
   allocate (indcell_loc(nbloc))
   allocate(pos(3,nbloc))
   allocate(r2vec(nbloc))
!
!     starting distances for long range real space interactions
!
   mbufbeg2 = (mpoleshortcut-lbuffer-shortheal)**2
!
!     perform a complete list build
!
   do i = 1, npolelocnl
      iipole = poleglobnl(i)
      iglob  = ipole(iipole)
      icell = repartcell(iglob)
!
!       align data of the local cell and the neighboring ones
!
      ncell_loc = cell_len(icell)
      indcell_loc(1:ncell_loc) =&
      &indcell(bufbegcell(icell):(bufbegcell(icell)+cell_len(icell)-1))
      do ineig = 1, numneigcell(icell)
         kcell = neigcell(ineig,icell)
         indcell_loc(ncell_loc+1:(ncell_loc+cell_len(kcell))) =&
         &indcell(bufbegcell(kcell):(bufbegcell(kcell)+cell_len(kcell)-1))
         ncell_loc = ncell_loc + cell_len(kcell)
      end do
!
!       do the neighbor search
!
      nneigloc = 0
      xi = x(iglob)
      yi = y(iglob)
      zi = z(iglob)
      do k = 1, ncell_loc
         kglob = indcell_loc(k)
         kkpole = pollist(kglob)
!
!   skip atom if it is not in the multipole list
!
         if (kkpole.eq.0) cycle
         if (kglob.le.iglob) cycle
         xk = x(kglob)
         yk = y(kglob)
         zk = z(kglob)
         pos(1,nneigloc+1) = xi - xk
         pos(2,nneigloc+1) = yi - yk
         pos(3,nneigloc+1) = zi - zk
         call midpointimagetriclinic(xi,yi,zi,xk,yk,zk,pos(1,nneigloc+1),&
         &pos(2,nneigloc+1),pos(3,nneigloc+1),docompute)
         if (docompute) then
            nneigloc = nneigloc + 1
            index(nneigloc) = kglob
         end if
      end do
!
!       compute the distances and build the list accordingly
!
      r2vec(1:nneigloc) = pos(1,1:nneigloc)*pos(1,1:nneigloc) +&
      &pos(2,1:nneigloc)*pos(2,1:nneigloc) +&
      &pos(3,1:nneigloc)*pos(3,1:nneigloc)

      j = 0
      j1 = 0
      do k = 1, nneigloc
         r2 = r2vec(k)
         kglob = index(k)
         if (r2 .le. mshortbuf2) then
            j1 = j1 + 1
            kkpole = pollist(kglob)
            shortelst(j1,i) = kkpole
         end if
         if (r2.le.mbuf2) then
            j = j + 1
            kkpole = pollist(kglob)
            elst(j,i) = kkpole
         end if
      end do
      nelst(i) = j
      nshortelst(i) = j1
!
!     check to see if the neighbor list is too long
!
      if (nelst(i) .ge. maxelst) then
         if (rank.eq.0) then
            write (iout,10)
10          format (/,' MBUILD  --  Too many Neighbors;',&
            &' Increase MAXELST')
            call fatal
         end if
      end if
   end do
!
   deallocate (pos)
   deallocate (index)
   deallocate (indcell_loc)
   deallocate (r2vec)
   return
end
!
!    subroutine initmpipme : build the arrays to communicate direct and reciprocal fields
!    during the calculation of the induced dipoles
!
!
!> @brief 
!> build the arrays to communicate direct and reciprocal fields
!> during the calculation of the induced dipoles
!> @param no params
subroutine initmpipme
   use atmlst
   use domdec
   use inform
   use iounit
   use mpole
   use pme
   use mpi
   implicit none
   integer ierr,iipole
   integer i,iproc,tag,iglob
   integer count1
   integer status(MPI_STATUS_SIZE)
   integer, allocatable :: req(:),req2(:),count(:)
   allocate (req(nproc*nproc))
   allocate (req2(nproc*nproc))
   allocate (count(nproc))
!
   if (deb_Path) write(iout,*), 'initmpipme '
!
!
   count = 0
!
!     deal with Direct-Recip communications
!
   if (allocated(buf1)) deallocate (buf1)
   allocate (buf1(nblocrecdir))
!      buf1 = 0
   if (allocated(buf2)) deallocate (buf2)
   allocate (buf2(nblocrecdir))
!      buf2 = 0
   if (allocated(buflen1)) deallocate (buflen1)
   allocate (buflen1(nproc))
   buflen1 = 0
   if (allocated(buflen2)) deallocate (buflen2)
   allocate (buflen2(nproc))
   buflen2 = 0
   if (allocated(bufbeg1)) deallocate (bufbeg1)
   allocate (bufbeg1(nproc))
   bufbeg1 = 0
   if (allocated(bufbeg2)) deallocate (bufbeg2)
   allocate (bufbeg2(nproc))
   bufbeg2 = 0
!
   do i = 1, npolerecloc
      iipole = polerecglob(i)
      iglob = ipole(iipole)
      if (repart(iglob).ne.rank) then
         buflen2(repart(iglob)+1) = buflen2(repart(iglob)+1)+1
      end if
   end do
   count1 = 0
   do iproc = 1, nrecdir_recep1
      if (precdir_recep1(iproc).ne.rank) then
         if (buflen2(precdir_recep1(iproc)+1).ne.0) then
            bufbeg2(precdir_recep1(iproc)+1) = count1 + 1
         else
            bufbeg2(precdir_recep1(iproc)+1) = 1
         end if
         count1 = count1 + buflen2(precdir_recep1(iproc)+1)
      end if
   end do
!
   do i = 1, npolerecloc
      iipole = polerecglob(i)
      iglob = ipole(iipole)
      if (repart(iglob).ne.rank) then
         buf2(bufbeg2(repart(iglob)+1)+count(repart(iglob)+1))=&
         &iipole
         count(repart(iglob)+1) = count(repart(iglob)+1) + 1
      end if
   end do
!
!     send and receive sizes of the buffers
!
   do i = 1, nrecdir_send1
      if (precdir_send1(i).ne.rank) then
         tag = nproc*rank + precdir_send1(i) + 1
         call MPI_IRECV(buflen1(precdir_send1(i)+1),1,MPI_INT,&
         &precdir_send1(i),tag,COMM_TINKER,req(tag),ierr)
      end if
   end do
   do i = 1, nrecdir_recep1
      if (precdir_recep1(i).ne.rank) then
         tag = nproc*precdir_recep1(i) + rank + 1
         call MPI_ISEND(buflen2(precdir_recep1(i)+1),1,MPI_INT,&
         &precdir_recep1(i),tag,COMM_TINKER,req(tag),ierr)
      end if
   end do
!
   do i = 1, nrecdir_send1
      if (precdir_send1(i).ne.rank) then
         tag = nproc*rank + precdir_send1(i) + 1
         call MPI_WAIT(req(tag),status,ierr)
      end if
   end do
   do i = 1, nrecdir_recep1
      if (precdir_recep1(i).ne.rank) then
         tag = nproc*precdir_recep1(i) + rank + 1
         call MPI_WAIT(req(tag),status,ierr)
      end if
   end do
   count1 = 0
   do iproc = 1, nrecdir_send1
      if (precdir_send1(iproc).ne.rank) then
         if (buflen1(precdir_send1(iproc)+1).ne.0) then
            bufbeg1(precdir_send1(iproc)+1) = count1 + 1
         else
            bufbeg1(precdir_send1(iproc)+1) = 1
         end if
         count1 = count1 + buflen1(precdir_send1(iproc)+1)
      end if
   end do
!
!     send and receive list of corresponding indexes
!
   do i = 1, nrecdir_send1
      if (precdir_send1(i).ne.rank) then
         tag = nproc*rank + precdir_send1(i) + 1
         call MPI_IRECV(buf1(bufbeg1(precdir_send1(i)+1)),&
         &buflen1(precdir_send1(i)+1),&
         &MPI_INT,precdir_send1(i),tag,COMM_TINKER,req2(tag),ierr)
      end if
   end do
   do i = 1, nrecdir_recep1
      if (precdir_recep1(i).ne.rank) then
         tag = nproc*precdir_recep1(i) + rank + 1
         call MPI_ISEND(buf2(bufbeg2(precdir_recep1(i)+1)),&
         &buflen2(precdir_recep1(i)+1),MPI_INT,precdir_recep1(i),tag,&
         &COMM_TINKER,req2(tag),ierr)
      end if
   end do
!
   do i = 1, nrecdir_send1
      if (precdir_send1(i).ne.rank) then
         tag = nproc*rank + precdir_send1(i) + 1
         call MPI_WAIT(req2(tag),status,ierr)
      end if
   end do
   do i = 1, nrecdir_recep1
      if (precdir_recep1(i).ne.rank) then
         tag = nproc*precdir_recep1(i) + rank + 1
         call MPI_WAIT(req2(tag),status,ierr)
      end if
   end do
!
   deallocate (req)
   deallocate (req2)
   deallocate (count)
   return
end
!
!     subroutine reinitnl : get the number of particules whose nl has to be computed
!     and the associated indexes
!
!> @brief 
!> get the number of particules whose nl has to be computed
!> and the associated indexes
!> @param[in] istep: index of timestep
subroutine reinitnl(istep)
   use atoms
   use cutoff
   use domdec
   use inform
   use iounit
   use neigh
   implicit none
   real*8 mbuf,vbuf,bigbuf
   integer modnl
   integer istep
!
   if (deb_Path) write(iout,*), 'reinitnl '
!
!
   mbuf = sqrt(mbuf2)+2.0d0
   vbuf = sqrt(vbuf2) + 2.0d0
   bigbuf = max(mbuf,vbuf,ddcut)
!
   modnl = mod(istep,ineigup)
   if (modnl.ne.0) return
!
   if (.not.allocated(ineignl)) allocate (ineignl(n))
   ineignl = 0
!
   call ddbuild(bigbuf,nbig_recep,pbig_recep,modnl.eq.0)
   return
end
!
!    "clistcell" performs a complete rebuild of the
!     electrostatic neighbor lists for charges using linked cells method
!
!> @brief 
!> performs a complete rebuild of the
!>  electrostatic neighbor lists for charges using linked cells method
subroutine clistcell
   use sizes
   use atmlst
   use atoms
   use boxes
   use charge
   use cutoff
   use domdec
   use inform
   use iounit
   use neigh
   use mpi
   implicit none
   integer iglob
   integer i,icell,j,k,nneigloc
   integer ineig,iichg,kkchg
   integer kcell,kglob
   integer ncell_loc
   integer, allocatable :: index(:),indcell_loc(:)
   real*8 xi,yi,zi,xk,yk,zk,r2
   real*8 boxedge2
   real*8, allocatable :: pos(:,:),r2vec(:)
   logical docompute
!
!     check size of the box and cutoff for minimum image convention
!
1000 format('Error in chargee neigbor list: max cutoff ',&
   &'buffer should be less than half one edge of the box')
1010 format('Charge cutoff = ',F14.3)
1020 format('List buffer      = ',F14.3)
!
   if (deb_Path) write(iout,*), 'clistcell '
!
   boxedge2 = min(xbox2,ybox2,zbox2)
   if (cbuf2.gt.boxedge2*boxedge2) then
      if (rank.eq.0) then
         write(iout,1000)
         write(iout,1010) chgcut
         write(iout,1020) lbuffer
      end if
      call fatal
   end if
!
   allocate (index(nbloc))
   allocate (indcell_loc(nbloc))
   allocate(pos(3,nbloc))
   allocate(r2vec(nbloc))
!
!     perform a complete list build
!
   do i = 1, nionlocnl
      iichg = chgglobnl(i)
      iglob  = iion(iichg)
      icell = repartcell(iglob)
!
!       align data of the local cell and the neighboring ones
!
      ncell_loc = cell_len(icell)
      indcell_loc(1:ncell_loc) =&
      &indcell(bufbegcell(icell):(bufbegcell(icell)+cell_len(icell)-1))
      do ineig = 1, numneigcell(icell)
         kcell = neigcell(ineig,icell)
         indcell_loc(ncell_loc+1:(ncell_loc+cell_len(kcell))) =&
         &indcell(bufbegcell(kcell):(bufbegcell(kcell)+cell_len(kcell)-1))
         ncell_loc = ncell_loc + cell_len(kcell)
      end do
!
!       do the neighbor search
!
      nneigloc = 0
      xi = x(iglob)
      yi = y(iglob)
      zi = z(iglob)
      do k = 1, ncell_loc
         kglob = indcell_loc(k)
         kkchg = chglist(kglob)
         if (kkchg.eq.0) cycle
         if (kglob.le.iglob) cycle
         xk = x(kglob)
         yk = y(kglob)
         zk = z(kglob)
         pos(1,nneigloc+1) = xi - xk
         pos(2,nneigloc+1) = yi - yk
         pos(3,nneigloc+1) = zi - zk
         call midpointimagetriclinic(xi,yi,zi,xk,yk,zk,pos(1,nneigloc+1),&
         &pos(2,nneigloc+1),pos(3,nneigloc+1),docompute)
         if (docompute) then
            nneigloc = nneigloc + 1
            index(nneigloc) = kglob
         end if
      end do
!
!       compute the distances and build the list accordingly
!
      r2vec(1:nneigloc) = pos(1,1:nneigloc)*pos(1,1:nneigloc) +&
      &pos(2,1:nneigloc)*pos(2,1:nneigloc) +&
      &pos(3,1:nneigloc)*pos(3,1:nneigloc)

      j = 0
      do k = 1, nneigloc
         r2 = r2vec(k)
         kglob = index(k)
         if (r2 .le. cbuf2) then
            j = j + 1
            kkchg = chglist(kglob)
            elst(j,i) = kkchg
         end if
      end do
      nelst(i) = j
!
!     check to see if the neighbor list is too long
!
      if (nelst(i) .ge. maxelst) then
         if (rank.eq.0) then
            write (iout,10)
10          format (/,' MBUILD  --  Too many Neighbors;',&
            &' Increase MAXELST')
            call fatal
         end if
      end if
   end do
!
   deallocate (pos)
   deallocate (index)
   deallocate (indcell_loc)
   deallocate (r2vec)
   return
end
!
!    "clistcell2" performs a complete rebuild of the
!     electrostatic short range and regular neighbor lists for charges
!     using linked cells method
!
!> @brief 
!> performs a complete rebuild of the
!>  electrostatic short range and regular neighbor lists for charges
!>  using linked cells method
!> @param[in] istep: index of timestep
subroutine clistcell2
   use sizes
   use atmlst
   use atoms
   use boxes
   use charge
   use cutoff
   use domdec
   use inform
   use iounit
   use neigh
   use mpi
   implicit none
   integer iglob
   integer i,icell,j,j1,k,nneigloc
   integer ineig,iichg,kkchg
   integer kcell,kglob
   integer ncell_loc
   integer, allocatable :: index(:),indcell_loc(:)
   real*8 xi,yi,zi,xk,yk,zk,r2
   real*8 cbufbeg2,boxedge2
   real*8, allocatable :: pos(:,:),r2vec(:)
   logical docompute
!
!     check size of the box and cutoff for minimum image convention
!
1000 format('Error in charge neigbor list: max cutoff ',&
   &'buffer should be less than half one edge of the box')
1010 format('Multipole cutoff = ',F14.3)
1020 format('List buffer      = ',F14.3)
!
   if (deb_Path) write(iout,*), 'clistcell2 '
!
   boxedge2 = min(xbox2,ybox2,zbox2)
   if (cbuf2.gt.boxedge2*boxedge2) then
      if (rank.eq.0) then
         write(iout,1000)
         write(iout,1010) chgcut
         write(iout,1020) lbuffer
      end if
      call fatal
   end if
!
   allocate (index(nbloc))
   allocate (indcell_loc(nbloc))
   allocate(pos(3,nbloc))
   allocate(r2vec(nbloc))
!
!     starting distances for long range real space interactions
!
   cbufbeg2 = (chgshortcut-lbuffer-shortheal)**2
!
!     perform a complete list build
!
   do i = 1, nionlocnl
      iichg = chgglobnl(i)
      iglob  = iion(iichg)
      icell = repartcell(iglob)
!
!       align data of the local cell and the neighboring ones
!
      ncell_loc = cell_len(icell)
      indcell_loc(1:ncell_loc) =&
      &indcell(bufbegcell(icell):(bufbegcell(icell)+cell_len(icell)-1))
      do ineig = 1, numneigcell(icell)
         kcell = neigcell(ineig,icell)
         indcell_loc(ncell_loc+1:(ncell_loc+cell_len(kcell))) =&
         &indcell(bufbegcell(kcell):(bufbegcell(kcell)+cell_len(kcell)-1))
         ncell_loc = ncell_loc + cell_len(kcell)
      end do
!
!       do the neighbor search
!
      nneigloc = 0
      xi = x(iglob)
      yi = y(iglob)
      zi = z(iglob)
      do k = 1, ncell_loc
         kglob = indcell_loc(k)
         kkchg = chglist(kglob)
         if (kkchg.eq.0) cycle
         if (kglob.le.iglob) cycle
         xk = x(kglob)
         yk = y(kglob)
         zk = z(kglob)
         pos(1,nneigloc+1) = xi - xk
         pos(2,nneigloc+1) = yi - yk
         pos(3,nneigloc+1) = zi - zk
         call midpointimagetriclinic(xi,yi,zi,xk,yk,zk,pos(1,nneigloc+1),&
         &pos(2,nneigloc+1),pos(3,nneigloc+1),docompute)
         if (docompute) then
            nneigloc = nneigloc + 1
            index(nneigloc) = kglob
         end if
      end do
!
!       compute the distances and build the list accordingly
!
      r2vec(1:nneigloc) = pos(1,1:nneigloc)*pos(1,1:nneigloc) +&
      &pos(2,1:nneigloc)*pos(2,1:nneigloc) +&
      &pos(3,1:nneigloc)*pos(3,1:nneigloc)

      j = 0
      j1 = 0
      do k = 1, nneigloc
         r2 = r2vec(k)
         kglob = index(k)
         if (r2 .le. cshortbuf2) then
            j1 = j1 + 1
            kkchg = chglist(kglob)
            shortelst(j1,i) = kkchg
         end if
         if ((r2 .le. cbuf2).and.(r2.ge.cbufbeg2)) then
            j = j + 1
            kkchg = chglist(kglob)
            elst(j,i) = kkchg
         end if
      end do
      nelst(i) = j
      nshortelst(i) = j1
!
!     check to see if the neighbor list is too long
!
      if ((nelst(i) .ge. maxelst).or.(nshortelst(i).ge.maxelst)) then
         if (rank.eq.0) then
            write (iout,10)
10          format (/,' MBUILD  --  Too many Neighbors;',&
            &' Increase MAXELST')
            call fatal
         end if
      end if
   end do
!
   deallocate (pos)
   deallocate (index)
   deallocate (indcell_loc)
   deallocate (r2vec)
   return
end

!
!    "vlistcell" performs a complete rebuild of the
!     vdw neighbor lists for charges using linked cells method
!
!> @brief 
!> "vlistcell" performs a complete rebuild of the
!>  vdw neighbor lists for charges using linked cells method
subroutine vlistcell
   use atmlst
   use atoms
   use boxes
   use bound
   use cutoff
   use domdec
   use inform
   use iounit
   use kvdws
   use neigh
   use vdw
   use mpi
   implicit none
   integer iglob,iloc
   integer i,ii,icell,j,k,nneigloc
   integer ineig,iivdw,iv
   integer kcell,kglob,kbis
   integer ncell_loc
   integer, allocatable :: index(:),indcell_loc(:)
   real*8 xi,yi,zi,xk,yk,zk,r2,rdn
   real*8 xr,yr,zr
   real*8, allocatable :: pos(:,:),r2vec(:)
   real*8, allocatable :: xred(:)
   real*8, allocatable :: yred(:)
   real*8, allocatable :: zred(:)
   real*8 boxedge2
   logical docompute
!
!     check size of the box and cutoff for minimum image convention
!
1000 format('Error in VDW neigbor list: max cutoff ',&
   &'buffer should be less than half one edge of the box')
1010 format('VDW cutoff = ',F14.3)
1020 format('List buffer      = ',F14.3)
!
   if (deb_Path) write(iout,*), 'vlistcell '
!
   boxedge2 = min(xbox2,ybox2,zbox2)
   if (vbuf2.gt.boxedge2*boxedge2) then
      if (rank.eq.0) then
         write(iout,1000)
         write(iout,1010) vdwcut
         write(iout,1020) lbuffer
      end if
      call fatal
   end if
!
   allocate (xred(nbloc))
   allocate (yred(nbloc))
   allocate (zred(nbloc))
!
   allocate (index(nbloc))
   allocate (indcell_loc(nbloc))
   allocate(pos(3,nbloc))
   allocate(r2vec(nbloc))
!
!     apply reduction factors to find coordinates for each site
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
!     perform a complete list build
!
   do i = 1, nvdwlocnl
      iivdw = vdwglobnl(i)
      iglob  = ivdw(iivdw)
      icell = repartcell(iglob)
      iloc = loc(iglob)
!
!       align data of the local cell and the neighboring ones
!
      ncell_loc = cell_len(icell)
      indcell_loc(1:ncell_loc) =&
      &indcell(bufbegcell(icell):(bufbegcell(icell)+cell_len(icell)-1))
      do ineig = 1, numneigcell(icell)
         kcell = neigcell(ineig,icell)
         indcell_loc(ncell_loc+1:(ncell_loc+cell_len(kcell))) =&
         &indcell(bufbegcell(kcell):(bufbegcell(kcell)+cell_len(kcell)-1))
         ncell_loc = ncell_loc + cell_len(kcell)
      end do
!
!       do the neighbor search
!
      nneigloc = 0
      xi = xred(iloc)
      yi = yred(iloc)
      zi = zred(iloc)
      do k = 1, ncell_loc
         kglob = indcell_loc(k)
         if (kglob.le.iglob) cycle
         if (rad(jvdw(kglob)).eq.0) cycle
         kbis = loc(kglob)
         xk = xred(kbis)
         yk = yred(kbis)
         zk = zred(kbis)
         call midpointimagetriclinic(xi,yi,zi,xk,yk,zk,pos(1,nneigloc+1),&
         &pos(2,nneigloc+1),pos(3,nneigloc+1),docompute)
         if (docompute) then
            nneigloc = nneigloc + 1
            index(nneigloc) = kglob
         end if
      end do
!
!       compute the distances and build the list accordingly
!
      r2vec(1:nneigloc) = pos(1,1:nneigloc)*pos(1,1:nneigloc) +&
      &pos(2,1:nneigloc)*pos(2,1:nneigloc) +&
      &pos(3,1:nneigloc)*pos(3,1:nneigloc)

      j = 0
      do k = 1, nneigloc
         r2 = r2vec(k)
         kglob = index(k)
         if (r2 .le. vbuf2) then
            j = j + 1
            vlst(j,i) = kglob
         end if
      end do
      nvlst(i) = j
!
!     check to see if the neighbor list is too long
!
      if (nvlst(i) .ge. maxvlst) then
         if (rank.eq.0) then
            write (iout,10)
10          format (/,' VBUILD  --  Too many Neighbors;',&
            &' Increase MAXVLST')
            call fatal
         end if
      end if
   end do
!
   deallocate (xred)
   deallocate (yred)
   deallocate (zred)
!
   deallocate (pos)
   deallocate (index)
   deallocate (indcell_loc)
   deallocate (r2vec)
   return
end
!
!    "vlistcell2" performs a complete rebuild of the
!     short range and regular vdw neighbor lists for charges
!      using linked cells method
!
!> @brief 
!> "vlistcell2" performs a complete rebuild of the
!>  short range and regular vdw neighbor lists for charges
!>   using linked cells method
subroutine vlistcell2
   use atmlst
   use atoms
   use boxes
   use bound
   use cutoff
   use domdec
   use inform
   use iounit
   use kvdws
   use neigh
   use vdw
   use mpi
   implicit none
   integer iglob,iloc
   integer i,ii,icell,j,j1,k,nneigloc
   integer ineig,iivdw,iv
   integer kcell,kglob,kbis
   integer ncell_loc
   integer, allocatable :: index(:),indcell_loc(:)
   real*8 xi,yi,zi,xk,yk,zk,r2,rdn
   real*8 xr,yr,zr
   real*8 boxedge2
   real*8, allocatable :: pos(:,:),r2vec(:)
   real*8, allocatable :: xred(:)
   real*8, allocatable :: yred(:)
   real*8, allocatable :: zred(:)
   logical docompute
!
!     check size of the box and cutoff for minimum image convention
!
1000 format('Error in VDW neigbor list: max cutoff ',&
   &'buffer should be less than half one edge of the box')
1010 format('VDW cutoff = ',F14.3)
1020 format('List buffer      = ',F14.3)
!
   if (deb_Path) write(iout,*), 'vlistcell2 '
!
   boxedge2 = min(xbox2,ybox2,zbox2)
   if (vbuf2.gt.boxedge2*boxedge2) then
      if (rank.eq.0) then
         write(iout,1000)
         write(iout,1010) vdwcut
         write(iout,1020) lbuffer
      end if
      call fatal
   end if
!
   allocate (xred(nbloc))
   allocate (yred(nbloc))
   allocate (zred(nbloc))
!
   allocate (index(nbloc))
   allocate (indcell_loc(nbloc))
   allocate(pos(3,nbloc))
   allocate(r2vec(nbloc))
!
!     starting distances for long range real space interactions
!
!
!     apply reduction factors to find coordinates for each site
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
!     perform a complete list build
!
   do i = 1, nvdwlocnl
      iivdw = vdwglobnl(i)
      iglob  = ivdw(iivdw)
      icell = repartcell(iglob)
      iloc = loc(iglob)
!
!       align data of the local cell and the neighboring ones
!
      ncell_loc = cell_len(icell)
      indcell_loc(1:ncell_loc) =&
      &indcell(bufbegcell(icell):(bufbegcell(icell)+cell_len(icell)-1))
      do ineig = 1, numneigcell(icell)
         kcell = neigcell(ineig,icell)
         indcell_loc(ncell_loc+1:(ncell_loc+cell_len(kcell))) =&
         &indcell(bufbegcell(kcell):(bufbegcell(kcell)+cell_len(kcell)-1))
         ncell_loc = ncell_loc + cell_len(kcell)
      end do
!
!       do the neighbor search
!
      nneigloc = 0
      xi = xred(iloc)
      yi = yred(iloc)
      zi = zred(iloc)
      do k = 1, ncell_loc
         kglob = indcell_loc(k)
         if (kglob.le.iglob) cycle
         if (rad(jvdw(kglob)).eq.0) cycle
         kbis = loc(kglob)
         xk = xred(kbis)
         yk = yred(kbis)
         zk = zred(kbis)
         pos(1,nneigloc+1) = xi - xk
         pos(2,nneigloc+1) = yi - yk
         pos(3,nneigloc+1) = zi - zk
         call midpointimagetriclinic(xi,yi,zi,xk,yk,zk,pos(1,nneigloc+1),&
         &pos(2,nneigloc+1),pos(3,nneigloc+1),docompute)
         if (docompute) then
            nneigloc = nneigloc + 1
            index(nneigloc) = kglob
         end if
      end do
!
!       compute the distances and build the list accordingly
!
      r2vec(1:nneigloc) = pos(1,1:nneigloc)*pos(1,1:nneigloc) +&
      &pos(2,1:nneigloc)*pos(2,1:nneigloc) +&
      &pos(3,1:nneigloc)*pos(3,1:nneigloc)

      j = 0
      j1 = 0
      do k = 1, nneigloc
         r2 = r2vec(k)
         kglob = index(k)
         if (r2 .le. vshortbuf2) then
            j1 = j1 + 1
            shortvlst(j1,i) = kglob
         end if
         if (r2.le.vbuf2) then
            j = j + 1
            vlst(j,i) = kglob
         end if
      end do
      nvlst(i) = j
      nshortvlst(i) = j1
!
!     check to see if the neighbor list is too long
!
      if ((nvlst(i).ge.maxvlst).or.(nshortvlst(i).ge.maxvlst)) then
         if (rank.eq.0) then
            write (iout,10)
10          format (/,' VBUILD  --  Too many Neighbors;',&
            &' Increase MAXVLST')
            call fatal
         end if
      end if
   end do
!
   deallocate (xred)
   deallocate (yred)
   deallocate (zred)
!
   deallocate (pos)
   deallocate (index)
   deallocate (indcell_loc)
   deallocate (r2vec)
   return
end
!
!    "dlistcell" performs a complete rebuild of the
!     dispersion neighbor lists using linked cells method
!
!> @brief 
!> "dlistcell" performs a complete rebuild of the
!>  dispersion neighbor lists using linked cells method
subroutine dlistcell
   use atmlst
   use atoms
   use boxes
   use cutoff
   use disp
   use domdec
   use inform
   use iounit
   use kvdws
   use neigh
   use vdw
   use mpi
   implicit none
   integer iglob,iloc
   integer i,icell,j,k,nneigloc
   integer ineig,iidisp
   integer kcell,kglob
   integer ncell_loc
   integer, allocatable :: index(:),indcell_loc(:)
   real*8 xi,yi,zi,xk,yk,zk,r2
   real*8, allocatable :: pos(:,:),r2vec(:)
   logical docompute
   real*8 boxedge2
!
!     check size of the box and cutoff for minimum image convention
!
1000 format('Error in dispersion neigbor list: max cutoff ',&
   &'buffer should be less than half one edge of the box')
1010 format('Dispersion cutoff = ',F14.3)
1020 format('List buffer      = ',F14.3)
!
   if (deb_Path) write(iout,*), 'dlistcell '
!
   boxedge2 = min(xbox2,ybox2,zbox2)
   if (dbuf2.gt.boxedge2*boxedge2) then
      if (rank.eq.0) then
         write(iout,1000)
         write(iout,1010) dispcut
         write(iout,1020) lbuffer
      end if
      call fatal
   end if
!
   allocate (index(nbloc))
   allocate (indcell_loc(nbloc))
   allocate(pos(3,nbloc))
   allocate(r2vec(nbloc))
!
!     perform a complete list build
!
!   write(*,*) 'ndispbloc = ',ndispbloc,'ndisplocnl = ',ndisplocnl
   do i = 1, ndisplocnl
      iidisp = dispglobnl(i)
      iglob  = idisp(iidisp)
      icell = repartcell(iglob)
      iloc = loc(iglob)
      xi = x(iglob)
      yi = y(iglob)
      zi = z(iglob)
!
!       align data of the local cell and the neighboring ones
!
      ncell_loc = cell_len(icell)
      indcell_loc(1:ncell_loc) =&
      &indcell(bufbegcell(icell):(bufbegcell(icell)+cell_len(icell)-1))
      do ineig = 1, numneigcell(icell)
         kcell = neigcell(ineig,icell)
         indcell_loc(ncell_loc+1:(ncell_loc+cell_len(kcell))) =&
         &indcell(bufbegcell(kcell):(bufbegcell(kcell)+cell_len(kcell)-1))
         ncell_loc = ncell_loc + cell_len(kcell)
      end do
!
!       do the neighbor search
!
      nneigloc = 0
      do k = 1, ncell_loc
         kglob = indcell_loc(k)
         if (kglob.le.iglob) cycle
         xk = x(kglob)
         yk = y(kglob)
         zk = z(kglob)
         pos(1,nneigloc+1) = xi - xk
         pos(2,nneigloc+1) = yi - yk
         pos(3,nneigloc+1) = zi - zk
         call midpointimagetriclinic(xi,yi,zi,xk,yk,zk,pos(1,nneigloc+1),&
         &pos(2,nneigloc+1),pos(3,nneigloc+1),docompute)
         if (docompute) then
            nneigloc = nneigloc + 1
            index(nneigloc) = kglob
         end if
      end do
!
!       compute the distances and build the list accordingly
!
      r2vec(1:nneigloc) = pos(1,1:nneigloc)*pos(1,1:nneigloc) +&
      &pos(2,1:nneigloc)*pos(2,1:nneigloc) +&
      &pos(3,1:nneigloc)*pos(3,1:nneigloc)

      j = 0
      do k = 1, nneigloc
         r2 = r2vec(k)
         kglob = index(k)
         if (r2 .le. dbuf2) then
            j = j + 1
            vlst(j,i) = displist(kglob)
         end if
      end do
      nvlst(i) = j
!
!     check to see if the neighbor list is too long
!
      if (nvlst(i) .ge. maxvlst) then
         if (rank.eq.0) then
            write (iout,10)
10          format (/,' DLISTCELL  --  Too many Neighbors;',&
            &' Increase MAXVLST')
            call fatal
         end if
      end if
   end do
!
   deallocate (pos)
   deallocate (index)
   deallocate (indcell_loc)
   deallocate (r2vec)
   return
end
!
!    "dlistcell2" performs a complete rebuild of the
!     dispersion neighbor lists using linked cells method
!
!> @brief 
!> "dlistcell2" performs a complete rebuild of the
!>  dispersion neighbor lists using linked cells method
subroutine dlistcell2
   use atmlst
   use atoms
   use boxes
   use cutoff
   use disp
   use domdec
   use inform
   use iounit
   use kvdws
   use neigh
   use vdw
   use mpi
   implicit none
   integer iglob,iloc
   integer i,icell,j,j1,k,nneigloc
   integer ineig,iidisp
   integer kcell,kglob
   integer ncell_loc
   integer, allocatable :: index(:),indcell_loc(:)
   real*8 xi,yi,zi,xk,yk,zk,r2
   real*8, allocatable :: pos(:,:),r2vec(:)
   logical docompute
   real*8 boxedge2
!
!     check size of the box and cutoff for minimum image convention
!
1000 format('Error in dispersion neigbor list: max cutoff ',&
   &'buffer should be less than half one edge of the box')
1010 format('Dispersion cutoff = ',F14.3)
1020 format('List buffer      = ',F14.3)
!
   if (deb_Path) write(iout,*), 'dlistcell2 '
!
   boxedge2 = min(xbox2,ybox2,zbox2)
   if (dbuf2.gt.boxedge2*boxedge2) then
      if (rank.eq.0) then
         write(iout,1000)
         write(iout,1010) dispcut
         write(iout,1020) lbuffer
      end if
      call fatal
   end if
!
   allocate (index(nbloc))
   allocate (indcell_loc(nbloc))
   allocate(pos(3,nbloc))
   allocate(r2vec(nbloc))
!
!     perform a complete list build
!
   do i = 1, ndisplocnl
      iidisp = dispglobnl(i)
      iglob  = idisp(iidisp)
      icell = repartcell(iglob)
      iloc = loc(iglob)
      xi = x(iglob)
      yi = y(iglob)
      zi = z(iglob)
!
!       align data of the local cell and the neighboring ones
!
      ncell_loc = cell_len(icell)
      indcell_loc(1:ncell_loc) =&
      &indcell(bufbegcell(icell):(bufbegcell(icell)+cell_len(icell)-1))
      do ineig = 1, numneigcell(icell)
         kcell = neigcell(ineig,icell)
         indcell_loc(ncell_loc+1:(ncell_loc+cell_len(kcell))) =&
         &indcell(bufbegcell(kcell):(bufbegcell(kcell)+cell_len(kcell)-1))
         ncell_loc = ncell_loc + cell_len(kcell)
      end do
!
!       do the neighbor search
!
      nneigloc = 0
      do k = 1, ncell_loc
         kglob = indcell_loc(k)
         if (kglob.le.iglob) cycle
         xk = x(kglob)
         yk = y(kglob)
         zk = z(kglob)
         pos(1,nneigloc+1) = xi - xk
         pos(2,nneigloc+1) = yi - yk
         pos(3,nneigloc+1) = zi - zk
         call midpointimagetriclinic(xi,yi,zi,xk,yk,zk,pos(1,nneigloc+1),&
         &pos(2,nneigloc+1),pos(3,nneigloc+1),docompute)
         if (docompute) then
            nneigloc = nneigloc + 1
            index(nneigloc) = kglob
         end if
      end do
!
!       compute the distances and build the list accordingly
!
      r2vec(1:nneigloc) = pos(1,1:nneigloc)*pos(1,1:nneigloc) +&
      &pos(2,1:nneigloc)*pos(2,1:nneigloc) +&
      &pos(3,1:nneigloc)*pos(3,1:nneigloc)

      j = 0
      j1 = 0
      do k = 1, nneigloc
         r2 = r2vec(k)
         kglob = index(k)
         if (r2 .le. dshortbuf2) then
            j1 = j1 + 1
            shortvlst(j1,i) = displist(kglob)
         end if
         if (r2 .le. dbuf2) then
            j = j + 1
            vlst(j,i) = displist(kglob)
         end if
      end do
      nvlst(i) = j
      nshortvlst(i) = j1
!
!     check to see if the neighbor list is too long
!
      if (nvlst(i) .ge. maxvlst) then
         if (rank.eq.0) then
            write (iout,10)
10          format (/,' DLISTCELL  --  Too many Neighbors;',&
            &' Increase MAXVLST')
            call fatal
         end if
      end if
   end do
!
   deallocate (pos)
   deallocate (index)
   deallocate (indcell_loc)
   deallocate (r2vec)
   return
end

!
!     ####################################################################
!     ##                                                                ##
!     ##  function imagevec3  --  compute the minimum image distance    ##
!     ##                                                                ##
!     ####################################################################
!
!
!     "imagevec3" takes the components of pairwise distances between
!     two points in a periodic box and converts to the components
!     of the minimum image distances. Indice i designs x, y or z
!     direction
!
!     xcell    length of the a-axis of the complete replicated cell
!     ycell    length of the b-axis of the complete replicated cell
!     zcell    length of the c-axis of the complete replicated cell
!     xcell2   half the length of the a-axis of the replicated cell
!     ycell2   half the length of the b-axis of the replicated cell
!     zcell2   half the length of the c-axis of the replicated cell

!     n REALLY SHOULD be a multiple of 16

subroutine imagevec3(imageout,pos,n,i)
   use cell
   implicit none
   integer, intent (in) ::i,n
   real*8,  intent (in) ::  pos(n)
   real*8,  intent (out) ::  imageout(n)
   integer k
   real*8 coordcell,coordcell2
!
   SELECT CASE (i)
    CASE (1)
      coordcell  = xcell
      coordcell2 = xcell2
    CASE (2)
      coordcell  = ycell
      coordcell2 = ycell2
    CASE (3)
      coordcell  = zcell
      coordcell2 = zcell2
    CASE DEFAULT
      coordcell  = xcell
      coordcell2 = xcell2
   END SELECT
   do k = 1, n
      imageout(k) =  pos(k)&
      &- int( (abs(pos(k)) - coordcell2) / coordcell&
      &+ 1.0d0&
      &) * sign (coordcell,pos(k))
   enddo
end subroutine imagevec3
