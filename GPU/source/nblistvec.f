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
#include "tinker_precision.h"
      subroutine nblistvec(istep)
      use sizes
      use domdec
      use cutoff
      use neigh
      use potent
      use timestat
      use mpi
      implicit none
      integer istep,modnl,j,k
      real(t_p)  time0,time1
c
c
c     check number of steps between nl updates
c
      modnl = mod(istep,ineigup)
      if (modnl.ne.0) return
      if (allocated(nelst)) deallocate (nelst)
      allocate (nelst(nlocnl))
      if (allocated(elst)) deallocate (elst)
      allocate (elst(maxelst,nlocnl))
!DIR$ ASSUME (mod(maxelst,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
      do j = 1,maxelst
!DIR$ ASSUME (mod(nlocnl,16).eq.0)
         do k = 1, nlocnl
            elst(j,k) = 0
            nelst(k) = 0
         enddo
      enddo
      if (allocated(nvlst)) deallocate (nvlst)
      allocate (nvlst(nlocnl))
      if (allocated(vlst)) deallocate (vlst)
      allocate (vlst(maxvlst,nlocnl))
!DIR$ ASSUME (mod(maxvlst,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
      do j = 1,maxvlst
!DIR$ ASSUME (mod(nlocnl,16).eq.0)
         do k = 1, nlocnl
            vlst(j,k) = 0
            nvlst(k) = 0
         enddo
      enddo
      if ((use_pmecore).and.(rank.gt.ndir-1)) return
c
      time0 = mpi_wtime()
c
c     build the cells at the beginning and assign the particules to them
c
      call build_cell_list(istep)
c
      call timer_enter( timer_nl )   
      if (use_clist) call clistcellvec
      if (use_vlist) call vlistcellvec
      if (use_mlist) call mlistcellvec
      call timer_exit( timer_nl )
      timenl = timenl + timer_get_last( timer_nl )
c
      return
      end
c
c    "mlistcellvec" performs a complete rebuild of the
c     electrostatic neighbor lists for multipoles using linked cells method
c
      subroutine mlistcellvec
      use atmlst
      use atoms
      use domdec
      use iounit
      use mpole
      use neigh
      use mpi
      use vec_list
      use utilvec

      implicit none
      integer iglob
      integer iipole
      integer icell,kcell
      integer lencell,lenloop16
      integer nneigcell,nneigloc,nneloop16
      integer ncell_loc,nceloop16
      integer i,j,k,kk
      real(t_p) xi,yi,zi
 
!DIR$ ATTRIBUTES ALIGN:64:: indcell_locvec
      integer indcell_locvec(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64:: kglobvec
      integer kglobvec(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64:: kglobvec1
      integer kglobvec1(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64:: kkpolevec
      integer kkpolevec(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64:: xposvec
      real(t_p) xposvec(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64:: yposvec
      real(t_p) yposvec(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64:: zposvec
      real(t_p) zposvec(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64:: xkvec
      real(t_p) xkvec(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64:: ykvec
      real(t_p) ykvec(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64:: zkvec
      real(t_p) zkvec(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64:: mask
      logical mask(nblocloop)
c
      if(rank.eq.0.and.tinkerdebug) write(*,*) 'mlistcellvec'

c
c     perform a complete list build
c
      do i = 1, npolelocnl
         iipole = poleglobnl(i)
         iglob  = ipole(iipole)
         icell  = repartcell(iglob)
c
c       align data of the local cell and the neighboring ones
c
         ncell_loc = cell_len(icell)
         nceloop16 = (int(ncell_loc / 16) + 1) * 16! First multiple of 16
         nceloop16 = merge(ncell_loc,nceloop16, mod(ncell_loc,16).eq.0)

         nneigcell = numneigcell(icell)
         nneloop16 = (int(nneigcell / 16) + 1) * 16! First multiple of 16
         nneloop16 = merge(nneigcell,nneloop16, mod(nneigcell,16).eq.0)
!DIR$ ASSUME (mod(nceloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, nceloop16
            if(k.le.ncell_loc) then
              indcell_locvec(k) = indcell(bufbegcell(icell)+k-1)
            endif
         enddo
!DIR$ ASSUME (mod(nneloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
         do j = 1, nneloop16
            if(j.le.nneigcell) then
               kcell     = neigcell(j,icell)
               lencell   = cell_len(kcell)
               lenloop16 = (int(lencell / 16) + 1) * 16! First multiple of 16
!DIR$ ASSUME (mod(lenloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
               do k = 1 , lenloop16
                  if(k.le.lencell) then
                     indcell_locvec(ncell_loc+k) =
     &                                indcell(bufbegcell(kcell)+k-1)
                  endif
               enddo
               ncell_loc = ncell_loc + lencell
            endif
         end do
         nceloop16 = (int(ncell_loc / 16) + 1) * 16! First multiple of 16
         nceloop16 = merge(ncell_loc,nceloop16, mod(ncell_loc,16).eq.0)
c
c       do the neighbor search
c
         xi = x(iglob)
         yi = y(iglob)
         zi = z(iglob)
!DIR$ ASSUME (mod(nceloop16,16).eq.0)
!DIR$ SIMD
         do k = 1, nceloop16
            if(k.le.ncell_loc) then
               kglobvec(k) = indcell_locvec(k)
            else
               kglobvec(k) = iglob ! exclude value by default
            endif
         enddo
c
c   keep atom if it is in the multipole list
c
         kk = 0
!DIR$ ASSUME (mod(nceloop16,16).eq.0)
!DIR$ VECTOR NOMASK_READWRITE
         do k = 1, nceloop16
            if(pollist(kglobvec(k)).ne.0 .and.kglobvec(k).gt.iglob) then
              kk = kk + 1
              kglobvec1(kk)  = kglobvec(k)
            endif
         enddo
         ncell_loc = kk

         nceloop16 = (int(ncell_loc / 16) + 1) * 16! First multiple of 16
         nceloop16 = merge(ncell_loc,nceloop16, mod(ncell_loc,16).eq.0)

!DIR$ ASSUME (mod(nceloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
         do k = 1, nceloop16
            if(k.le.ncell_loc) then
               xkvec(k) = x(kglobvec1(k))
               ykvec(k) = y(kglobvec1(k))
               zkvec(k) = z(kglobvec1(k))
            endif
               xposvec(k) = xi - xkvec(k)
               yposvec(k) = yi - ykvec(k)
               zposvec(k) = zi - zkvec(k)
         enddo
!DIR$ ASSUME (mod(nceloop16,16).eq.0)
cold     call  imagevec(xposvec,nceloop16,1)
cold     call  imagevec(yposvec,nceloop16,2)
cold     call  imagevec(zposvec,nceloop16,3)
         call  image3dvec(xposvec,yposvec,zposvec,nceloop16)
!DIR$ ASSUME (mod(nceloop16,16).eq.0)
c
c   select  the atoms with the midpoint method
c
         call midpointimagevec(ncell_loc,
     &                          nceloop16,
     &                          xkvec,
     &                          ykvec,
     &                          zkvec,
     &                          xposvec,
     &                          yposvec,
     &                          zposvec,
     &                          mask)
c
c    combine the mask with the distances cutoff and build the list accordingly
c
         kk = 0
!DIR$ ASSUME (mod(nceloop16,16).eq.0)
         do k = 1, nceloop16
            mask(k) = (     mask(k)
     &                 .and.(  xposvec(k) ** 2 + yposvec(k) ** 2
     &                       + zposvec(k) ** 2 .le.mbuf2 ))
            if (mask(k)) then
               kk = kk + 1
               kkpolevec(kk) = pollist(kglobvec1(k))
            end if
         end do
         nelst(i) = kk
         nneloop16 = (int(nelst(i) / 16) + 1) * 16! First multiple of 16
         nneloop16 = merge(nelst(i),nneloop16, mod(nelst(i),16).eq.0)
!DIR$ ASSUME (mod(nneloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
        do k = 1,  nneloop16
             elst(k,i) = kkpolevec(k)
        enddo
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
!$acc data present(elst,nelst)
!$acc update device(elst(:,:),nelst(:))
!$acc end data
      return
      end subroutine mlistcellvec
c
c    subroutine initmpipme : build the arrays to communicate direct and reciprocal fields
c    during the calculation of the induced dipoles
c
c
      subroutine initmpipmevec
      use atmlst
      use domdec
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
c
!DIR$ VECTOR ALIGNED
      count = 0 
c
c     deal with Direct-Recip communications
c
      if (allocated(buf1)) deallocate (buf1)
      allocate (buf1(nblocrecdir))
c      buf1 = 0
      if (allocated(buf2)) deallocate (buf2)
      allocate (buf2(nblocrecdir))
c      buf2 = 0
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
c
!DIR$ VECTOR ALIGNED
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
!$acc update device(buf2,bufbeg2)
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
!$acc update device(bufbeg1)
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
!$acc update device(buf1)
      if (allocated(thetai1)) deallocate (thetai1)
      if (allocated(thetai2)) deallocate (thetai2)
      if (allocated(thetai3)) deallocate (thetai3)
      allocate (thetai1(4,bsorder,nlocrec))
      allocate (thetai2(4,bsorder,nlocrec))
      allocate (thetai3(4,bsorder,nlocrec))
!$acc enter data create(thetai1,thetai2,thetai3)
      deallocate (req)
      deallocate (req2)
      deallocate (count)
      return
      end
c
c     subroutine reinitnl : get the number of particules whose nl has to be computed
c     and the associated indexes
c
      subroutine reinitnlvec(istep)
      use atoms
      use domdec
      use neigh
      use tinheader ,only:ti_p,re_p
      implicit none
      real(t_p) d,mbuf,vbuf,torquebuf,bigbuf
      integer iproc,i,iglob,modnl
      integer iloc,istep
c
      mbuf = sqrt(mbuf2)
      vbuf = sqrt(vbuf2) + 2.0_ti_p
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
      if (.not.allocated(ineignl)) allocate (ineignl(n))
      ineignl = 0
c
      nlocnl = nloc
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
      do i = 1,nloc
c     ineignl(1:nloc) = glob(1:nloc)
         ineignl(i) = glob(i)
      enddo
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
      subroutine build_cell_listvec(istep)
      use atoms
      use bound
      use domdec
      use neigh
      use mpi
      implicit none
      integer i,proc,icell,j,k,p,q,r,istep,iglob
      integer count,iloc
      integer temp_x,temp_y,temp_z
      integer temp,numneig,tempcell
      real(t_p) xmin,xmax,ymin,ymax,zmin,zmax
      real(t_p) lenx,leny,lenz
      real(t_p) mbuf,vbuf,bigbuf
      real(t_p) lenx_cell,leny_cell,lenz_cell
      real(t_p) xr,yr,zr
      real(t_p) eps1,eps2
      real(t_p), allocatable :: xbegcelltemp(:),ybegcelltemp(:)
      real(t_p), allocatable :: zbegcelltemp(:)
      real(t_p), allocatable :: xendcelltemp(:),yendcelltemp(:)
      real(t_p), allocatable :: zendcelltemp(:)
      integer, allocatable :: filledcell(:),indcelltemp(:)
      logical docompute
c
      eps1 = 1.0d-10
      eps2 = 1.0d-8
      mbuf = sqrt(mbuf2)
      vbuf = sqrt(vbuf2)+2.0
      if ((mbuf+lbuffer).gt.vbuf) then
        bigbuf = (mbuf+lbuffer)/2
      else
        bigbuf = vbuf/2
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
!DIR$ VECTOR ALIGNED
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
      if (allocated(xbegcell)) deallocate (xbegcell)
      allocate (xbegcell(ncell_tot))
      if (allocated(ybegcell)) deallocate (ybegcell)
      allocate (ybegcell(ncell_tot))
      if (allocated(zbegcell)) deallocate (zbegcell)
      allocate (zbegcell(ncell_tot))
      if (allocated(xendcell)) deallocate (xendcell)
      allocate (xendcell(ncell_tot))
      if (allocated(yendcell)) deallocate (yendcell)
      allocate (yendcell(ncell_tot))
      if (allocated(zendcell)) deallocate (zendcell)
      allocate (zendcell(ncell_tot))
      if (allocated(neigcell)) deallocate (neigcell)
      allocate (neigcell(124,ncell_tot))
      if (allocated(numneigcell)) deallocate (numneigcell)
      allocate (numneigcell(ncell_tot))
      allocate (filledcell(ncell_tot))
c
!DIR$ VECTOR ALIGNED
      do i = 0, nx_cell-1
        xbegcelltemp(i+1) = xmin + i*lenx_cell
        xendcelltemp(i+1) = xmin + (i+1)*lenx_cell
      end do
!DIR$ VECTOR ALIGNED
      do i = 0, ny_cell-1
        ybegcelltemp(i+1) = ymin + i*leny_cell
        yendcelltemp(i+1) = ymin + (i+1)*leny_cell
      end do
!DIR$ VECTOR ALIGNED
      do i = 0, nz_cell-1
        zbegcelltemp(i+1) = zmin + i*lenz_cell
        zendcelltemp(i+1) = zmin + (i+1)*lenz_cell
      end do
c
c     assign cell
c
!DIR$ VECTOR ALIGNED
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
!DIR$ VECTOR ALIGNED
              filledcell = 0
              filledcell(icell) = 1
c
c              do p = -1,1
c                do q = -1,1
c                  do r = -1,1
              do p = -2,2
                do q = -2,2
                  do r = -2,2
                    if ((p.eq.0).and.(q.eq.0).and.(r.eq.0)) goto 10
c
                    temp_x = p+i
                    temp_y = q+j-1
                    temp_z = r+k-1
c                    if ((i.eq.1).and.(p.eq.-1)) temp_x = nx_cell
c                    if ((i.eq.nx_cell).and.(p.eq.1)) temp_x = 1
c                    if ((j.eq.1).and.(q.eq.-1)) temp_y = ny_cell-1
c                    if ((j.eq.ny_cell).and.(q.eq.1)) temp_y = 0
c                    if ((k.eq.1).and.(r.eq.-1)) temp_z = nz_cell-1
c                    if ((k.eq.nz_cell).and.(r.eq.1)) temp_z = 0

                    if ((i.eq.1).and.(p.eq.-2)) temp_x = nx_cell-1
                    if ((i.eq.1).and.(p.eq.-1)) temp_x = nx_cell
                    if ((i.eq.2).and.(p.eq.-2)) temp_x = nx_cell
                    if ((i.eq.nx_cell).and.(p.eq.1)) temp_x = 1
                    if ((i.eq.nx_cell).and.(p.eq.2)) temp_x = 2
                    if ((i.eq.nx_cell-1).and.(p.eq.2)) temp_x = 1
                    if ((j.eq.1).and.(q.eq.-2)) temp_y = ny_cell-2
                    if ((j.eq.1).and.(q.eq.-1)) temp_y = ny_cell-1
                    if ((j.eq.2).and.(q.eq.-2)) temp_y = ny_cell-1
                    if ((j.eq.ny_cell).and.(q.eq.1)) temp_y = 0
                    if ((j.eq.ny_cell).and.(q.eq.2)) temp_y = 1
                    if ((j.eq.ny_cell-1).and.(q.eq.2)) temp_y = 0
                    if ((k.eq.1).and.(r.eq.-2)) temp_z = nz_cell-2
                    if ((k.eq.1).and.(r.eq.-1)) temp_z = nz_cell-1
                    if ((k.eq.2).and.(r.eq.-2)) temp_z = nz_cell-1
                    if ((k.eq.nz_cell).and.(r.eq.1)) temp_z = 0
                    if ((k.eq.nz_cell).and.(r.eq.2)) temp_z = 1
                    if ((k.eq.nz_cell-1).and.(r.eq.2)) temp_z = 0
                    tempcell = temp_z*ny_cell*nx_cell+temp_y*nx_cell+
     $                temp_x
c                    write(*,*) 'tempcell = ',tempcell,i,j,k
c                    write(*,*) '******'
c                    write(*,*) 'p,q,r = ',p,q,r
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
      if (allocated(cell_len)) deallocate (cell_len)
      allocate (cell_len(ncell_tot))
      if (allocated(indcell)) deallocate (indcell)
      allocate (indcell(n))
      if (allocated(bufbegcell)) deallocate (bufbegcell)
      allocate (bufbegcell(ncell_tot))
      if (allocated(repartcell)) deallocate (repartcell)
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
!DIR$ VECTOR ALIGNED
        do icell = 1, ncell_tot
          if ((zr.ge.zbegcell(icell)).and.
     $     (zr.lt.zendcell(icell)).and.(yr.ge.ybegcell(icell))
     $    .and.(yr.lt.yendcell(icell)).and.(xr.ge.xbegcell(icell))
     $    .and.(xr.lt.xendcell(icell))) then
!DIR£ VECTOR ALIGNED
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
!DIR$ VECTOR ALIGNED
      do i = 1, nlocnl
        iglob = ineignl(i)
        icell = repartcell(iglob)
        iloc  = bufbegcell(icell) + indcelltemp(iglob) - 1
        indcell(iloc) = iglob
c       indcell(bufbegcell(icell) + indcelltemp(iglob) - 1) = iglob
      end do
      deallocate (indcelltemp)
      return
      end
c
c    "clistcellvec" performs a complete rebuild of the
c     electrostatic neighbor lists for charges using linked cells method
c
      subroutine clistcellvec
      use sizes
      use atmlst
      use atoms
      use domdec
      use iounit
      use charge
      use neigh
      use mpi
      use vec_list
      use utilvec

      implicit none
      integer iglob
      integer iichg
      integer icell,kcell
      integer lencell,lenloop16
      integer nneigcell,nneigloc,nneloop16
      integer ncell_loc,nceloop16
      integer i,j,k, kk

      real(t_p) xi,yi,zi
!DIR$ ATTRIBUTES ALIGN:64:: indcell_locvec
      integer indcell_locvec(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64:: kglobvec
      integer kglobvec(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64:: kglobvec1
      integer kglobvec1(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64:: kkchgvec
      integer kkchgvec(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64:: xposvec
      real(t_p) xposvec(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64:: yposvec
      real(t_p) yposvec(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64:: zposvec
      real(t_p) zposvec(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64:: xkvec
      real(t_p) xkvec(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64:: ykvec
      real(t_p) ykvec(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64:: zkvec
      real(t_p) zkvec(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64:: mask
      logical mask(nblocloop)
c
      if(rank.eq.0.and.tinkerdebug) write(*,*) 'clistcellvec'

c
c     perform a complete list build
c
      do i = 1, nionlocnl
         iichg  = chgglobnl(i)
         iglob  = iion(iichg)
         icell = repartcell(iglob)
c
c       align data of the local cell and the neighboring ones
c
         ncell_loc = cell_len(icell)
         nceloop16 = (int(ncell_loc / 16) + 1) * 16! First multiple of 16
         nceloop16 = merge(ncell_loc,nceloop16, mod(ncell_loc,16).eq.0)

         nneigcell = numneigcell(icell)
         nneloop16 = (int(nneigcell / 16) + 1) * 16! First multiple of 16
         nneloop16 = merge(nneigcell,nneloop16, mod(nneigcell,16).eq.0)
!DIR$ ASSUME (mod(nceloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, nceloop16
            if(k.le.ncell_loc) then
              indcell_locvec(k) = indcell(bufbegcell(icell)+k-1)
            endif
         enddo
!DIR$ ASSUME (mod(nneloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
         do j = 1, nneloop16
            if(j.le.nneigcell) then
               kcell     = neigcell(j,icell)
               lencell   = cell_len(kcell)
               lenloop16 = (int(lencell / 16) + 1) * 16! First multiple of 16
!DIR$ ASSUME (mod(lenloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
               do k = 1 , lenloop16
                  if(k.le.lencell) then
                     indcell_locvec(ncell_loc+k) =
     &                                   indcell(bufbegcell(kcell)+k-1)
                  endif
               enddo
               ncell_loc = ncell_loc + lencell
            endif
         end do
         nceloop16 = (int(ncell_loc / 16) + 1) * 16! First multiple of 16
         nceloop16 = merge(ncell_loc,nceloop16, mod(ncell_loc,16).eq.0)
c
c       do the neighbor search
c
         xi = x(iglob)
         yi = y(iglob)
         zi = z(iglob)
!DIR$ ASSUME (mod(nceloop16,16).eq.0)
!DIR$ SIMD
         do k = 1, nceloop16
            if(k.le.ncell_loc) then
               kglobvec(k) = indcell_locvec(k)
            else
               kglobvec(k) = iglob ! exclude value by default
            endif
         enddo
c
c   keep atom if it is in the charge list
c
         kk = 0
!DIR$ ASSUME (mod(nceloop16,16).eq.0)
!DIR$ VECTOR NOMASK_READWRITE
         do k = 1, nceloop16
            if(     chglist(kglobvec(k)).ne.0
     &         .and.kglobvec(k).gt.iglob) then
              kk = kk + 1
              kglobvec1(kk) = kglobvec(k)
            endif
         enddo
         ncell_loc = kk

         nceloop16 = (int(ncell_loc / 16) + 1) * 16! First multiple of 16
         nceloop16 = merge(ncell_loc,nceloop16, mod(ncell_loc,16).eq.0)

!DIR$ ASSUME (mod(nceloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
         do k = 1, nceloop16
            if(k.le.ncell_loc) then
               xkvec(k) = x(kglobvec1(k))
               ykvec(k) = y(kglobvec1(k))
               zkvec(k) = z(kglobvec1(k))
            endif
            xposvec(k) = xi - xkvec(k)
            yposvec(k) = yi - ykvec(k)
            zposvec(k) = zi - zkvec(k)
         enddo
!DIR$ ASSUME (mod(nceloop16,16).eq.0)
cold     call  imagevec(xposvec,nceloop16,1)
cold     call  imagevec(yposvec,nceloop16,2)
cold     call  imagevec(zposvec,nceloop16,3)
         call  image3dvec(xposvec,yposvec,zposvec,nceloop16)
c
c   select  the atoms with the midpoint method
c
         call midpointimagevec(ncell_loc,
     &                          nceloop16,
     &                          xkvec,
     &                          ykvec,
     &                          zkvec,
     &                          xposvec,
     &                          yposvec,
     &                          zposvec,
     &                          mask)
c
c    combine the mask with the distances cutoff and build the list accordingly
c
         kk = 0
!DIR$ ASSUME (mod(nceloop16,16).eq.0)
         do k = 1, nceloop16
            mask(k) = (     mask(k)
     &                 .and.(  xposvec(k) ** 2 + yposvec(k) ** 2
     &                       + zposvec(k) ** 2 .le.cbuf2 ))
            if (mask(k)) then
               kk = kk + 1
               kkchgvec(kk) = chglist(kglobvec1(k))
            end if
         end do
         nelst(i) = kk
         nneloop16 = (int(nelst(i) / 16) + 1) * 16! First multiple of 16
         nneloop16 = merge(nelst(i),nneloop16, mod(nelst(i),16).eq.0)
!DIR$ ASSUME (mod(nneloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
        do k = 1,  nneloop16
             elst(k,i) = kkchgvec(k)
        enddo
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
!$acc data present(elst,nelst)
!$acc update device(elst(:,:),nelst(:)) 
!$acc end data
      return
      end subroutine clistcellvec
c
c    "vlistcell" performs a complete rebuild of the
c     vdw neighbor lists for charges using linked cells method
c
      subroutine vlistcellvec
      use atmlst
      use atoms
      use domdec
      use iounit
      use kvdws
      use neigh
      use vdw
      use mpi
      use vec_list
      use utilvec

      implicit none
      integer iglob
      integer iglobdefault,kglobdefault
      integer iloc
      integer icell,kcell
      integer lencell,lenloop16
      integer nneigcell,nneigloc,nneloop16
      integer ncell_loc,nceloop16
      integer i,j,k, kk
      integer iivdw

      real(t_p) xi,yi,zi
!DIR$ ATTRIBUTES ALIGN:64:: indcell_locvec
      integer indcell_locvec(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64:: iglobvec
      integer iglobvec(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64:: ivec
      integer ivec(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64:: ivvec
      integer ivvec(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64:: kglobvec
      integer kglobvec(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64:: kglobvec1
      integer kglobvec1(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64:: kbisvec
      integer kbisvec(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64:: rdnvec
      real(t_p) rdnvec(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64:: xposvec
      real(t_p) xposvec(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64:: yposvec
      real(t_p) yposvec(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64:: zposvec
      real(t_p) zposvec(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64:: xkvec
      real(t_p) xkvec(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64:: ykvec
      real(t_p) ykvec(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64:: zkvec
      real(t_p) zkvec(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64:: xred
      real(t_p) xred(0:nblocloop)
!DIR$ ATTRIBUTES ALIGN:64:: yred
      real(t_p) yred(0:nblocloop)
!DIR$ ATTRIBUTES ALIGN:64:: zred
      real(t_p) zred(0:nblocloop)
!DIR$ ATTRIBUTES ALIGN:64:: mask
      logical mask(nblocloop)
      if(rank.eq.0.and.tinkerdebug) write(*,*) 'vlistcellvec'

c     set default values to point to for various indices
      iglobdefault = ivdw(vdwglob (nvdwbloc))

c
c     apply reduction factors to find coordinates for each site
c
!DIR$ ASSUME (mod(nvdwblocloop,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ VECTOR G2S 
!DIR$ SIMD
      do k = 1,nvdwblocloop
         if (k.le.nvdwbloc) then
            iglobvec (k) = ivdw (vdwglob (k))
         else
            iglobvec(k) = iglobdefault ! safe value by default
         endif
      enddo

!DIR$ ASSUME (mod(nvdwblocloop,16).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1 , nvdwblocloop
         ivec    (k) = loc  (iglobvec(k))
         ivvec   (k) = ired (iglobvec(k))
         rdnvec  (k) = kred (iglobvec(k))

         xred(ivec(k)) =  rdnvec(k) * x(iglobvec(k))
     &                  + (1.0_ti_p - rdnvec(k)) * x(ivvec(k))
         yred(ivec(k)) =  rdnvec(k) * y(iglobvec(k))
     &                  + (1.0_ti_p - rdnvec(k)) * y(ivvec(k))
         zred(ivec(k)) =  rdnvec(k) * z(iglobvec(k))
     &                  + (1.0_ti_p - rdnvec(k)) * z(ivvec(k))
      enddo

c
c     perform a complete list build
c
      do i = 1, nvdwlocnl
         iivdw = vdwglobnl(i)
         iglob = ivdw(iivdw)
         icell = repartcell(iglob)
         iloc  = loc(iglob)
         kglobdefault = iglob
c
c       align data of the local cell and the neighboring ones
c
         ncell_loc = cell_len(icell)
         nceloop16 = (int(ncell_loc / 16) + 1) * 16! First multiple of 16
         nceloop16 = merge(ncell_loc,nceloop16, mod(ncell_loc,16).eq.0)

         nneigcell = numneigcell(icell)
         nneloop16 = (int(nneigcell / 16) + 1) * 16! First multiple of 16
         nneloop16 = merge(nneigcell,nneloop16, mod(nneigcell,16).eq.0)
!DIR$ ASSUME (mod(nceloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, nceloop16
            if(k.le.ncell_loc) then
              indcell_locvec(k) = indcell(bufbegcell(icell)+k-1)
            endif
         enddo
!DIR$ ASSUME (mod(nneloop16,16).eq.0)
!DIR$ NOFUSION
!DIR$ VECTOR ALIGNED
         do j = 1, nneloop16
            if(j.le.nneigcell) then
               kcell     = neigcell(j,icell)
               lencell   = cell_len(kcell)
               lenloop16 = (int(lencell / 16) + 1) * 16! First multiple of 16
!DIR$ ASSUME (mod(lenloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
               do k = 1 , lenloop16
                  if(k.le.lencell) then
                     indcell_locvec(ncell_loc+k) =
     &                                   indcell(bufbegcell(kcell)+k-1)
                  endif
               enddo
               ncell_loc = ncell_loc + lencell
            endif
         end do

         nceloop16 = (int(ncell_loc / 16) + 1) * 16! First multiple of 16
         nceloop16 = merge(ncell_loc,nceloop16, mod(ncell_loc,16).eq.0)
c
c       do the neighbor search
c
         xi = xred(iloc)
         yi = yred(iloc)
         zi = zred(iloc)
!DIR$ ASSUME (mod(nceloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
         do k = 1, nceloop16
            if(k.le.ncell_loc) then
               kglobvec(k) = indcell_locvec(k)
            else
               kglobvec(k) = kglobdefault ! exclude value by default
            endif
         enddo
c   keep atom if it is in the vdw neighbor charge list
c
         kk = 0
!DIR$ ASSUME (mod(nceloop16,16).eq.0)
!DIR$ VECTOR NOMASK_READWRITE
!DIR$ VECTOR ALIGNED
         do k = 1, nceloop16
            if(     rad(jvdw(kglobvec(k))).ne.0
     &         .and.kglobvec(k).gt.iglob) then
              kk = kk + 1
              kglobvec1(kk) = kglobvec(k)
              kbisvec  (kk) = loc(kglobvec(k))
            endif
         enddo
         ncell_loc = kk
  
         nceloop16 = (int(ncell_loc / 16) + 1) * 16! First multiple of 16
         nceloop16 = merge(ncell_loc,nceloop16, mod(ncell_loc,16).eq.0)

!DIR$ ASSUME (mod(nceloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, nceloop16
            if(k.le.ncell_loc) then
               xkvec(k) = xred(kbisvec(k))
               ykvec(k) = yred(kbisvec(k))
               zkvec(k) = zred(kbisvec(k))
            endif
            xposvec(k) = xi - xkvec(k)
            yposvec(k) = yi - ykvec(k)
            zposvec(k) = zi - zkvec(k)
         enddo
!DIR$ ASSUME (mod(nceloop16,16).eq.0)
c        call imagevec(xposvec,nceloop16,1)
c        call imagevec(yposvec,nceloop16,2)
c        call imagevec(zposvec,nceloop16,3)
         call image3dvec(xposvec,yposvec,zposvec,nceloop16)
c
c   select  the atoms with the midpoint method
c
         call midpointimagevec(ncell_loc,
     &                          nceloop16,
     &                          xkvec,
     &                          ykvec,
     &                          zkvec,
     &                          xposvec,
     &                          yposvec,
     &                          zposvec,
     &                          mask)
c
c    combine the mask with the distances cutoff and build the list accordingly
c
         kk = 0
!DIR$ ASSUME (mod(nceloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, nceloop16
            mask(k) = (     mask(k)
     &                 .and.(  xposvec(k) ** 2 + yposvec(k) ** 2
     &                       + zposvec(k) ** 2 .le.vbuf2 ))
            if (mask(k)) then
               kk = kk + 1
               kglobvec(kk) = kglobvec1(k)
            end if
         end do
         nvlst(i) = kk
         nneloop16 = (int(nvlst(i) / 16) + 1) * 16! First multiple of 16
         nneloop16 = merge(nvlst(i),nneloop16, mod(nvlst(i),16).eq.0)
!DIR$ ASSUME (mod(nneloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
        do k = 1,  nneloop16
             vlst(k,i) = kglobvec(k)
        enddo
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
!$acc data present(vlst,nvlst)
!$acc update device(vlst(:,:),nvlst(:))
!$acc end data
      return
      end subroutine vlistcellvec
