c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine nblist  --  maintain pairwise neighbor lists  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "nblist" constructs and maintains nonbonded pair neighbor lists
c     for vdw and electrostatic interactions
c
c
      subroutine nblist(istep)
      implicit none
      include 'sizes.i'
      include 'cutoff.i'
      include 'potent.i'
      include 'openmp.i'
      include 'mpif.h'
      include 'timestat.i'
      include 'neigh.i'
      integer istep,modnl
      real*8  time0,time1
c
c
c     check number of steps between nl updates
c
      modnl = mod(istep,ineigup)
      if (modnl.ne.0) return
      if (associated(nelst)) deallocate (nelst)
      allocate (nelst(nlocnl))
      if (associated(elst)) deallocate (elst)
      allocate (elst(maxelst,nlocnl))
      nelst = 0
      elst = 0
      if (associated(nvlst)) deallocate (nvlst)
      allocate (nvlst(nlocnl))
      if (associated(vlst)) deallocate (vlst)
      allocate (vlst(maxvlst,nlocnl))
      nvlst = 0
      vlst = 0
      if ((use_pmecore).and.(rank.gt.ndir-1)) return
c
      time0 = mpi_wtime()
c
c     build the cells at the beginning and assign the particules to them
c
      call build_cell_list(istep)
c
      if (use_mlist) call mlistcell
      if (use_vlist) call vlistcell
      time1 = mpi_wtime()
      timenl = timenl + time1 - time0
c
      return
      end
c
c    "mlistcell" performs a complete rebuild of the
c     electrostatic neighbor lists for atomic multipoles using linked cells method
c
      subroutine mlistcell
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'atmlst.i'
      include 'iounit.i'
      include 'mpole.i'
      include 'neigh.i'
      include 'openmp.i'
      include 'mpif.h'
      integer iglob
      integer i,icell,j,k
      integer ineig,iipole
      integer kcell,kloc,kglob
      real*8 xr,yr,zr,xi,yi,zi,xk,yk,zk,r2
      logical docompute
c
c
c     perform a complete list build
c
      do i = 1, npolelocnl
        j = 0
        iipole = poleglobnl(i)
        iglob  = ipole(iipole)
        icell = repartcell(iglob)
        xi = x(iglob)
        yi = y(iglob)
        zi = z(iglob)
c
c      search in the same cell
c
        do k = 1, cell_len(icell)
          kloc = bufbegcell(icell) + k - 1
          kglob = indcell(kloc)
          if (kglob.le.iglob) cycle
          xk = x(kglob)
          yk = y(kglob)
          zk = z(kglob)
          call midpoint(xi,yi,zi,xk,yk,zk,docompute)
          if (docompute) then
            xr = xi - xk
            yr = yi - yk
            zr = zi - zk
            call image (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. mbuf2) then
               j = j + 1
               elst(j,i) = kglob
            end if
          end if
        end do
c
c      search in the neighboring cells
c
        do ineig = 1, numneigcell(icell)
          kcell = neigcell(ineig,icell)
          do k = 1, cell_len(kcell)
            kloc = bufbegcell(kcell) + k - 1
            kglob = indcell(kloc)
            if (kglob.le.iglob) cycle
            xk = x(kglob)
            yk = y(kglob)
            zk = z(kglob)
            call midpoint(xi,yi,zi,xk,yk,zk,docompute)
            if (docompute) then
              xr = xi - xk
              yr = yi - yk
              zr = zi - zk
              call image (xr,yr,zr)
              r2 = xr*xr + yr*yr + zr*zr
              if (r2 .le. mbuf2) then
                 j = j + 1
                 elst(j,i) = kglob
              end if
            end if
          end do
        end do
        nelst(i) = j
c
c     check to see if the neighbor list is too long
c
        if (nelst(i) .ge. maxelst) then
           if (rank.eq.0) then
             write (iout,10)
   10        format (/,' MBUILD  --  Too many Neighbors;',
     &                  ' Increase MAXELST')
             call fatal
           end if
        end if
      end do
      return
      end
c
c    "vlistcell" performs a complete rebuild of the
c     vdw neighbor lists using linked cells method
c
      subroutine vlistcell
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'atmlst.i'
      include 'iounit.i'
      include 'vdw.i'
      include 'neigh.i'
      include 'openmp.i'
      include 'mpif.h'
      integer iglob,iloc
      integer i,icell,j,k
      integer ineig
      integer kcell,kloc,kglob,kbis
      integer ii,iivdw,iv
      real*8 xr,yr,zr,xi,yi,zi,xk,yk,zk,r2,rdn
      real*8, allocatable :: xred(:)
      real*8, allocatable :: yred(:)
      real*8, allocatable :: zred(:)
      logical docompute
c
      allocate (xred(nbloc))
      allocate (yred(nbloc))
      allocate (zred(nbloc))
c
c     apply reduction factors to find coordinates for each site
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
c     perform a complete list build
c
      do i = 1, nvdwlocnl
        j = 0
        iivdw = vdwglobnl(i)
        iglob  = ivdw(iivdw)
        icell = repartcell(iglob)
        iloc = loc(iglob)
        xi = xred(iloc)
        yi = yred(iloc)
        zi = zred(iloc)
c
c      search in the same cell
c
        do k = 1, cell_len(icell)
          kloc = bufbegcell(icell) + k - 1
          kglob = indcell(kloc)
          if (kglob.le.iglob) cycle
          kbis = loc(kglob)
          xk = xred(kbis)
          yk = yred(kbis)
          zk = zred(kbis)
          call midpoint(xi,yi,zi,xk,yk,zk,docompute)
          if (docompute) then
            xr = xi - xk
            yr = yi - yk
            zr = zi - zk
            call image (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. vbuf2) then
               j = j + 1
               vlst(j,i) = kglob
            end if
          end if
        end do
c
c      search in the neighboring cells
c
        do ineig = 1, numneigcell(icell)
          kcell = neigcell(ineig,icell)
          do k = 1, cell_len(kcell)
            kloc = bufbegcell(kcell) + k - 1
            kglob = indcell(kloc)
            if (kglob.le.iglob) cycle
            kbis = loc(kglob)
            xk = xred(kbis)
            yk = yred(kbis)
            zk = zred(kbis)
            call midpoint(xi,yi,zi,xk,yk,zk,docompute)
            if (docompute) then
              xr = xi - xk
              yr = yi - yk
              zr = zi - zk
              call image (xr,yr,zr)
              r2 = xr*xr + yr*yr + zr*zr
              if (r2 .le. vbuf2) then
                 j = j + 1
                 vlst(j,i) = kglob
              end if
            end if
          end do
        end do
        nvlst(i) = j
c
c     check to see if the neighbor list is too long
c
        if (nvlst(i) .ge. maxvlst) then
           if (rank.eq.0) then
             write (iout,10)
   10        format (/,' VBUILD  --  Too many Neighbors;',
     &                  ' Increase MAXVLST')
             call fatal
           end if
        end if
      end do
      deallocate (xred)
      deallocate (yred)
      deallocate (zred)
      return
      end
c
c    subroutine initmpipme : build the arrays to communicate direct and reciprocal fields
c    during the calculation of the induced dipoles
c
c
      subroutine initmpipme
      implicit none
      include 'sizes.i'
      include 'atmlst.i'
      include 'openmp.i'
      include 'mpole.i'
      include 'pme.i'
      include 'mpif.h'
      integer ierr,iipole
      integer i,iproc,tag,iglob
      integer count1
      integer status(MPI_STATUS_SIZE)
      integer, allocatable :: req(:),req2(:),count(:)
      allocate (req(nproc*nproc))
      allocate (req2(nproc*nproc))
      allocate (count(nproc))
c
      count = 0 
c
c     deal with Direct-Recip communications
c
      if (associated(buf1)) deallocate (buf1)
      allocate (buf1(nblocrecdir))
c      buf1 = 0
      if (associated(buf2)) deallocate (buf2)
      allocate (buf2(nblocrecdir))
c      buf2 = 0
      if (associated(buflen1)) deallocate (buflen1)
      allocate (buflen1(nproc))
      buflen1 = 0
      if (associated(buflen2)) deallocate (buflen2)
      allocate (buflen2(nproc))
      buflen2 = 0
      if (associated(bufbeg1)) deallocate (bufbeg1)
      allocate (bufbeg1(nproc))
      bufbeg1 = 0
      if (associated(bufbeg2)) deallocate (bufbeg2)
      allocate (bufbeg2(nproc))
      bufbeg2 = 0
c
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
c
      do i = 1, npolerecloc
        iipole = polerecglob(i)
        iglob = ipole(iipole)
        if (repart(iglob).ne.rank) then
          buf2(bufbeg2(repart(iglob)+1)+count(repart(iglob)+1))=
     $      iipole
          count(repart(iglob)+1) = count(repart(iglob)+1) + 1
        end if
      end do
c
c     send and receive sizes of the buffers
c
       do i = 1, nrecdir_send1
         if (precdir_send1(i).ne.rank) then
          tag = nproc*rank + precdir_send1(i) + 1
          call MPI_IRECV(buflen1(precdir_send1(i)+1),1,MPI_INT,
     $   precdir_send1(i),tag,MPI_COMM_WORLD,req(tag),ierr)
        end if
      end do
      do i = 1, nrecdir_recep1
        if (precdir_recep1(i).ne.rank) then
          tag = nproc*precdir_recep1(i) + rank + 1
          call MPI_ISEND(buflen2(precdir_recep1(i)+1),1,MPI_INT,
     $     precdir_recep1(i),tag,MPI_COMM_WORLD,req(tag),ierr)
        end if
      end do
c
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
c
c     send and receive list of corresponding indexes
c
      do i = 1, nrecdir_send1
        if (precdir_send1(i).ne.rank) then
          tag = nproc*rank + precdir_send1(i) + 1
          call MPI_IRECV(buf1(bufbeg1(precdir_send1(i)+1)),
     $     buflen1(precdir_send1(i)+1),
     $     MPI_INT,precdir_send1(i),tag,MPI_COMM_WORLD,req2(tag),ierr)
        end if
      end do
      do i = 1, nrecdir_recep1
        if (precdir_recep1(i).ne.rank) then
          tag = nproc*precdir_recep1(i) + rank + 1
          call MPI_ISEND(buf2(bufbeg2(precdir_recep1(i)+1)),
     $     buflen2(precdir_recep1(i)+1),MPI_INT,precdir_recep1(i),tag,
     $     MPI_COMM_WORLD,req2(tag),ierr)
        end if
      end do
c
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
c
      if (associated(thetai1)) deallocate (thetai1)
      if (associated(thetai2)) deallocate (thetai2)
      if (associated(thetai3)) deallocate (thetai3)
      allocate (thetai1(4,bsorder,nlocrec))
      allocate (thetai2(4,bsorder,nlocrec))
      allocate (thetai3(4,bsorder,nlocrec))
      deallocate (req)
      deallocate (req2)
      deallocate (count)
      return
      end
c
c     subroutine reinitnl : get the number of particules whose nl has to be computed
c     and the associated indexes
c
      subroutine reinitnl(istep)
      implicit none
      include 'atoms.i'
      include 'openmp.i'
      include 'neigh.i'
      real*8 d,mbuf,vbuf,torquebuf,bigbuf
      integer iproc,i,iglob,modnl
      integer iloc,istep
c
      mbuf = sqrt(mbuf2)
      vbuf = sqrt(vbuf2) + 2.0d0
      torquebuf = mbuf + lbuffer
      if (torquebuf.gt.(vbuf)) then
        bigbuf = torquebuf
      else
        bigbuf = vbuf
      end if
c
      modnl = mod(istep,ineigup)
      if (modnl.ne.0) return
c
      if (.not.associated(ineignl)) allocate (ineignl(n))
      ineignl = 0
c
      nlocnl = nloc
      ineignl(1:nloc) = glob(1:nloc)
c
      do iproc = 1, nbig_recep
        do i = 1, domlen(pbig_recep(iproc)+1)
          iloc = bufbeg(pbig_recep(iproc)+1)+i-1
          iglob = glob(iloc)
          call distprocpart(iglob,rank,d,.true.)
          if (d.le.(bigbuf/2)) then
            nlocnl = nlocnl + 1
            ineignl(nlocnl) = iglob
c            locnl(iglob) = nlocnl
          end if
        end do
      end do
      return
      end
c
c
c     subroutine build_cell_list : build the cells in order to build the non bonded neighbor
c     lists with the cell-list method
c
      subroutine build_cell_list(istep)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'neigh.i'
      include 'openmp.i'
      include 'mpif.h'
      integer i,proc,icell,j,k,p,q,r,istep,iglob
      integer count,iloc
      integer temp_x,temp_y,temp_z
      integer nx_cell,ny_cell,nz_cell
      integer temp,numneig,tempcell
      real*8 xmin,xmax,ymin,ymax,zmin,zmax
      real*8 lenx,leny,lenz
      real*8 mbuf,vbuf,bigbuf
      real*8 lenx_cell,leny_cell,lenz_cell
      real*8 xr,yr,zr
      real*8 eps1,eps2
      real*8, allocatable :: xbegcelltemp(:),ybegcelltemp(:)
      real*8, allocatable :: zbegcelltemp(:)
      real*8, allocatable :: xendcelltemp(:),yendcelltemp(:)
      real*8, allocatable :: zendcelltemp(:)
      integer, allocatable :: filledcell(:),indcelltemp(:)
      logical docompute
c
      eps1 = 1.0d-10
      eps2 = 1.0d-8
      mbuf = sqrt(mbuf2)
      vbuf = sqrt(vbuf2)+2.0
      if ((mbuf+lbuffer).gt.vbuf) then
        bigbuf = mbuf+lbuffer
      else
        bigbuf = vbuf
      end if
c
c
c     divide the searching domain in cells of size the multipole cutoff
c
      xmin = xbegproc(rank+1)
      xmax = xendproc(rank+1)
      ymin = ybegproc(rank+1)
      ymax = yendproc(rank+1)
      zmin = zbegproc(rank+1)
      zmax = zendproc(rank+1)
      do i = 1, nbig_recep
        proc = pbig_recep(i)
        if (xbegproc(proc+1).le.xmin) xmin = xbegproc(proc+1)
        if (xendproc(proc+1).ge.xmax) xmax = xendproc(proc+1)
        if (ybegproc(proc+1).le.ymin) ymin = ybegproc(proc+1)
        if (yendproc(proc+1).ge.ymax) ymax = yendproc(proc+1)
        if (zbegproc(proc+1).le.zmin) zmin = zbegproc(proc+1)
        if (zendproc(proc+1).ge.zmax) zmax = zendproc(proc+1)
      end do
c
      lenx = abs(xmax-xmin)
      nx_cell = max(1,int(lenx/(bigbuf)))
      lenx_cell = lenx/nx_cell
      leny = abs(ymax-ymin)
      ny_cell = max(1,int(leny/(bigbuf)))
      leny_cell = leny/ny_cell
      lenz = abs(zmax-zmin)
      nz_cell = max(1,int(lenz/(bigbuf)))
      lenz_cell = lenz/nz_cell
      ncell_tot = nx_cell*ny_cell*nz_cell
c
      allocate (xbegcelltemp(nx_cell))
      allocate (xendcelltemp(nx_cell))
      allocate (ybegcelltemp(ny_cell))
      allocate (yendcelltemp(ny_cell))
      allocate (zbegcelltemp(nz_cell))
      allocate (zendcelltemp(nz_cell))
      if (associated(xbegcell)) deallocate (xbegcell)
      allocate (xbegcell(ncell_tot))
      if (associated(ybegcell)) deallocate (ybegcell)
      allocate (ybegcell(ncell_tot))
      if (associated(zbegcell)) deallocate (zbegcell)
      allocate (zbegcell(ncell_tot))
      if (associated(xendcell)) deallocate (xendcell)
      allocate (xendcell(ncell_tot))
      if (associated(yendcell)) deallocate (yendcell)
      allocate (yendcell(ncell_tot))
      if (associated(zendcell)) deallocate (zendcell)
      allocate (zendcell(ncell_tot))
      if (associated(neigcell)) deallocate (neigcell)
      allocate (neigcell(26,ncell_tot))
      if (associated(numneigcell)) deallocate (numneigcell)
      allocate (numneigcell(ncell_tot))
      allocate (filledcell(ncell_tot))
c
      do i = 0, nx_cell-1
        xbegcelltemp(i+1) = xmin + i*lenx_cell
        xendcelltemp(i+1) = xmin + (i+1)*lenx_cell
      end do
      do i = 0, ny_cell-1
        ybegcelltemp(i+1) = ymin + i*leny_cell
        yendcelltemp(i+1) = ymin + (i+1)*leny_cell
      end do
      do i = 0, nz_cell-1
        zbegcelltemp(i+1) = zmin + i*lenz_cell
        zendcelltemp(i+1) = zmin + (i+1)*lenz_cell
      end do
c
c     assign cell
c
      do k = 1, nz_cell
        do j = 1, ny_cell
          do i = 1, nx_cell
              icell = (k-1)*ny_cell*nx_cell+(j-1)*nx_cell+i
              xbegcell(icell) = xbegcelltemp(i)
              xendcell(icell) = xendcelltemp(i)
              ybegcell(icell) = ybegcelltemp(j)
              yendcell(icell) = yendcelltemp(j)
              zbegcell(icell) = zbegcelltemp(k)
              zendcell(icell) = zendcelltemp(k)
              numneig = 0
              filledcell = 0
              filledcell(icell) = 1
c
              do p = -1,1
                do q = -1,1
                  do r = -1,1
                    if ((p.eq.0).and.(q.eq.0).and.(r.eq.0)) goto 10
c
                    temp_x = p+i
                    temp_y = q+j-1
                    temp_z = r+k-1
                    if ((i.eq.1).and.(p.eq.-1)) temp_x = nx_cell
                    if ((i.eq.nx_cell).and.(p.eq.1)) temp_x = 1
                    if ((j.eq.1).and.(q.eq.-1)) temp_y = ny_cell-1
                    if ((j.eq.ny_cell).and.(q.eq.1)) temp_y = 0
                    if ((k.eq.1).and.(r.eq.-1)) temp_z = nz_cell-1
                    if ((k.eq.nz_cell).and.(r.eq.1)) temp_z = 0
                    tempcell = temp_z*ny_cell*nx_cell+temp_y*nx_cell+
     $                temp_x
                    if (filledcell(tempcell).eq.1) goto 10
                    filledcell(tempcell) = 1
                    numneig = numneig+1
                    neigcell(numneig,icell) = tempcell
 10               continue
                  end do
                end do
              end do
              numneigcell(icell) = numneig
          end do
        end do
      end do
      deallocate (filledcell)
      deallocate (xbegcelltemp)
      deallocate (xendcelltemp)
      deallocate (ybegcelltemp)
      deallocate (yendcelltemp)
      deallocate (zbegcelltemp)
      deallocate (zendcelltemp)
c
c     assign the atoms to the cells
c
      if (associated(cell_len)) deallocate (cell_len)
      allocate (cell_len(ncell_tot))
      if (associated(indcell)) deallocate (indcell)
      allocate (indcell(n))
      if (associated(bufbegcell)) deallocate (bufbegcell)
      allocate (bufbegcell(ncell_tot))
      if (associated(repartcell)) deallocate (repartcell)
      allocate (repartcell(n))
      allocate (indcelltemp(n))
      cell_len = 0
c
      do i = 1, nlocnl
        iglob = ineignl(i)
        xr = x(iglob)
        yr = y(iglob)
        zr = z(iglob)
        if (use_bounds) call image(xr,yr,zr)
        if (abs(xr-xmax).lt.eps1) xr = xr-eps2
        if (abs(yr-ymax).lt.eps1) yr = yr-eps2
        if (abs(zr-zmax).lt.eps1) zr = zr-eps2
        do icell = 1, ncell_tot
          if ((zr.ge.zbegcell(icell)).and.
     $     (zr.lt.zendcell(icell)).and.(yr.ge.ybegcell(icell))
     $    .and.(yr.lt.yendcell(icell)).and.(xr.ge.xbegcell(icell))
     $    .and.(xr.lt.xendcell(icell))) then
            repartcell(iglob) = icell
            cell_len(icell) = cell_len(icell) + 1
            indcelltemp(iglob) = cell_len(icell)
          end if
        end do
      end do
c
      bufbegcell(1) = 1
      count = cell_len(1)
      do icell = 2, ncell_tot
        if (cell_len(icell).ne.0) then
          bufbegcell(icell) = count + 1
        else
          bufbegcell(icell) = 1
        end if
        count = count + cell_len(icell)
      end do
c
      do i = 1, nlocnl
        iglob = ineignl(i)
        icell = repartcell(iglob)
        iloc  = bufbegcell(icell) + indcelltemp(iglob) - 1
        indcell(iloc) = iglob
      end do
      deallocate (indcelltemp)
      return
      end
