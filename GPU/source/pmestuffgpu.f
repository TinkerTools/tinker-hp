c
c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c
c     "bspline_fill_site" finds B-spline coefficients and derivatives
c     for PME i-th atomic sites along the fractional coordinate axes
c
c
#include "tinker_precision.h"
      subroutine bspline_fill_sitegpu(config)
      use atmlst
      use atoms
      use boxes
      use charge   ,only:iion,nion,nionrecloc
      use domdec   ,only:nproc,nrec
      use inform   ,only:deb_Path
      use interfaces,only:grid_pchg_site_p,grid_pchg_sitecu
     &              ,grid_mpole_site_p,grid_mpole_sitecu
      use mpole
      use neigh    ,only:celle_pole,celle_chg
      use pme
      use tinheader,only:ti_p,prec1_eps
      use utilgpu  ,only:rec_queue,ngangs_rec
      use potent
      implicit none
      integer,intent(in),optional::config
      integer c_mpole,c_charge,c_n
      integer i,ifr,k,cfg,iipole,iichg
      integer,pointer,save:: glob_p(:),type_p(:)
      real(t_p) xi,yi,zi
      real(t_p) w,fr
      parameter(c_mpole=0,c_charge=1)
!$acc routine(bsplgen) seq
!$acc routine(bsplgen_chg) seq

      if (present(config)) then
         cfg = config
      else
         cfg = c_mpole
      end if

      if (cfg.eq.c_mpole) then
         glob_p => ipole
         type_p => polerecglob
         c_n    =  npolerecloc
#ifdef _OPENACC
         if (associated(grid_mpole_site_p,grid_mpole_sitecu)) then
            type_p => celle_pole
         end if
#endif
      else if (cfg.eq.c_charge) then
         glob_p => iion
         type_p => chgrecglob
         c_n    =  nionrecloc
#ifdef _OPENACC
         if (associated(grid_pchg_site_p,grid_pchg_sitecu)) then
            type_p => celle_chg
         end if
#endif
      end if

      if (cfg.eq.c_mpole) then

!$acc parallel loop num_gangs(4*ngangs_rec)
!$acc&         present(x,y,z,recip,igrid,thetai1,thetai2,
!$acc&  thetai3,type_p,glob_p)
!$acc&         async(rec_queue)
      do k=1,c_n
         iipole = type_p(k)
         i      = glob_p(iipole)
c
c     get the b-spline coefficients for the i-th atomic site
c
         xi   = x(i)
         yi   = y(i)
         zi   = z(i)
         w    = xi*recip(1,1) + yi*recip(2,1) + zi*recip(3,1)
         fr   = real(nfft1,t_p) * (w-anint(w)+0.5_ti_p)
         ifr  = int(fr-pme_eps)
         w    = fr - real(ifr,t_p)
         igrid(1,i) = ifr - bsorder
         call bsplgen (w,k,thetai1(1,1,1))
         w    = xi*recip(1,2) + yi*recip(2,2) + zi*recip(3,2)
         fr   = real(nfft2,t_p) * (w-anint(w)+0.5_ti_p)
         ifr  = int(fr-pme_eps)
         w    = fr - real(ifr,t_p)
         igrid(2,i) = ifr - bsorder
         call bsplgen (w,k,thetai2(1,1,1))
         w    = xi*recip(1,3) + yi*recip(2,3) + zi*recip(3,3)
         fr   = real(nfft3,t_p) * (w-anint(w)+0.5_ti_p)
         ifr  = int(fr-pme_eps)
         w    = fr - real(ifr,t_p)
         igrid(3,i) = ifr - bsorder
         call bsplgen (w,k,thetai3(1,1,1))
      end do

      else if (cfg.eq.c_charge) then

!$acc parallel loop num_gangs(4*ngangs_rec)
!$acc&         present(x,y,z,recip,igrid,thetai1,thetai2,
!$acc&  thetai3,type_p,glob_p)
!$acc&         async(rec_queue)
      do k=1,c_n
         iichg  = type_p(k)
         i      = glob_p(iichg)
c
c     get the b-spline coefficients for the i-th atomic site
c
         xi   = x(i)
         yi   = y(i)
         zi   = z(i)
         w    = xi*recip(1,1) + yi*recip(2,1) + zi*recip(3,1)
         fr   = real(nfft1,t_p) * (w-anint(w)+0.5_ti_p)
         ifr  = int(fr-pme_eps)
         w    = fr - real(ifr,t_p)
         igrid(1,i) = ifr - bsorder
         call bsplgen_chg (w,k,thetai1(1,1,1))
         w    = xi*recip(1,2) + yi*recip(2,2) + zi*recip(3,2)
         fr   = real(nfft2,t_p) * (w-anint(w)+0.5_ti_p)
         ifr  = int(fr-pme_eps)
         w    = fr - real(ifr,t_p)
         igrid(2,i) = ifr - bsorder
         call bsplgen_chg (w,k,thetai2(1,1,1))
         w    = xi*recip(1,3) + yi*recip(2,3) + zi*recip(3,3)
         fr   = real(nfft3,t_p) * (w-anint(w)+0.5_ti_p)
         ifr  = int(fr-pme_eps)
         w    = fr - real(ifr,t_p)
         igrid(3,i) = ifr - bsorder
         call bsplgen_chg (w,k,thetai3(1,1,1))
      end do

      end if
      end
c
c       "grid_mpole_site" places the i-th fractional atomic multipole onto
c       the particle mesh Ewald grid
c
c
      subroutine grid_mpole_sitegpu(fmpvec)
      use atmlst
      use chunks
      use domdec
      use fft
      use mpole
      use pme
      use potent
      use sizes
      use utils
      use utilgpu,only:rec_queue,ngangs_rec
      implicit none
      integer istart,iend,jstart,jend,kstart,kend
      integer istat,ied,jstat,jed,kstat,ked
      integer iipole
      integer i,j,k,m,impi,rankloc
      integer mk,mj
      integer iproc,proc
      integer ii,jj,kk
      integer isite,iatm
      integer offsetx,offsety
      integer offsetz
      integer location,twonlpts_1,twonlpts_12
      integer nlptsit
      real(t_p) v0,u0,t0
      real(t_p) v1,u1,t1
      real(t_p) v2,u2,t2
      real(t_p) vut(6)
      real(t_p) term,term0,term1,term2
      real(t_p) fmp(10)
      real(t_p) fmpvec(10,max(npolerecloc,1))
c
      if (use_pmecore) then
        rankloc  = rank_bis
      else
        rankloc  = rank
      end if
      kstat = kstart1(rankloc+1)
      ked   = kend1  (rankloc+1)
      jstat = jstart1(rankloc+1)
      jed   = jend1  (rankloc+1)
      istat = istart1(rankloc+1)
      ied   = iend1  (rankloc+1)
      twonlpts_1  = 2*nlpts+1
      twonlpts_12 = twonlpts_1**2
      nlptsit     = (2*nlpts+1)**3 - 1

c
c     put the permanent multipole moments onto the grid
c
!$acc parallel loop gang worker vector_length(32)
!$acc&         present(polerecglob,igrid,ipole,thetai1,thetai2)
!$acc&         present(thetai3,fmpvec,qgridin_2d)
!$acc&         async(rec_queue)
      do impi = 1,npolerecloc
        isite   = ipole(polerecglob(impi))
        offsetx = 1 - (igrid(1,isite) + grdoff - nlpts)
        offsety = 1 - (igrid(2,isite) + grdoff - nlpts)
        offsetz = 1 - (igrid(3,isite) + grdoff - nlpts)
c
c       put the induced dipole moments onto the grid
c
!$acc loop vector private(vut)
        do kk = 0, nlptsit
           k  = igrid(3,isite) + grdoff + kk/twonlpts_12 - nlpts
           j  = igrid(2,isite) + grdoff 
     &        + mod(kk/twonlpts_1,twonlpts_1) - nlpts
           i  = igrid(1,isite) + grdoff 
     &        + mod(kk,twonlpts_1)-nlpts
           mk = k + offsetz
           if (k .lt. 1) k = k + nfft3
           mj = j + offsety
           if (j .lt. 1) j = j + nfft2
           m  = i + offsetx
           if (i .lt. 1) i = i + nfft1
c
           vut(1)= thetai3(1,mk,impi)
           vut(2)= thetai3(2,mk,impi)
           vut(3)= thetai3(3,mk,impi)
           vut(4)= thetai2(1,mj,impi)
           vut(5)= thetai2(2,mj,impi)
           vut(6)= thetai2(3,mj,impi)
           term0 = fmpvec(1,impi)*vut(4)*vut(1)
     &           + fmpvec(3,impi)*vut(5)*vut(1)
     &           + fmpvec(4,impi)*vut(4)*vut(2)
     &           + fmpvec(6,impi)*vut(6)*vut(1)
     &           + fmpvec(7,impi)*vut(4)*vut(3)
     &           + fmpvec(10,impi)*vut(5)*vut(2)
           term1 = fmpvec(2,impi)*vut(4)*vut(1)
     &           + fmpvec(8,impi)*vut(5)*vut(1)
     &           + fmpvec(9,impi)*vut(4)*vut(2)
           term2 = fmpvec(5,impi)*vut(4)*vut(1)
           vut(1)= thetai1(1,m,impi)
           vut(2)= thetai1(2,m,impi)
           vut(3)= thetai1(3,m,impi)
           term  = term0*vut(1) + term1*vut(2) + term2*vut(3)

           if (((k.ge.kstat).and.(k.le.ked)).and.
     &         ((j.ge.jstat).and.(j.le.jed)).and.
     &         ((i.ge.istat).and.(i.le.ied))) then
c            location = 1+(i-istat)*2+(j-jstat)*2*n1mpimax
c    &                   +(k-kstat)*2*n1mpimax*n2mpimax
c            call atomic_adds(location,term,qgridin_p(1))
!$acc atomic update
             qgridin_2d(1,i-istat+1,j-jstat+1,k-kstat+1,1) =
     &       qgridin_2d(1,i-istat+1,j-jstat+1,k-kstat+1,1)
     &       + term
             cycle
           end if

!$acc loop seq
           do iproc = 1, nrec_send
              proc   = prec_send(iproc)
              kstart = kstart1(proc+1)
              kend   = kend1  (proc+1)
              jstart = jstart1(proc+1)
              jend   = jend1  (proc+1)
              istart = istart1(proc+1)
              iend   = iend1  (proc+1)
              if (((k.ge.kstart).and.(k.le.kend)).and.
     &            ((j.ge.jstart).and.(j.le.jend)).and.
     &            ((i.ge.istart).and.(i.le.iend))) then
c             location = 1+(i-istart)*2+(j-jstart)*2*n1mpimax
c    &                    +(k-kstart)*2*n1mpimax*n2mpimax
c    &                    +iproc*2*n1mpimax*n2mpimax*n3mpimax
c             call atomic_adds(location,term,qgridin_p(1))
!$acc atomic update
             qgridin_2d(1,i-istart+1,j-jstart+1,k-kstart+1,iproc+1)
     &     = qgridin_2d(1,i-istart+1,j-jstart+1,k-kstart+1,iproc+1)
     &     + term
                 exit
               end if
           end do
        end do
      end do
      end
c
c     "grid_pchg_site" places the i-th fractional atomic charge onto
c     the particle mesh Ewald grid
c
      subroutine grid_pchg_sitegpu
      use atmlst
      use charge
      use chunks
      use domdec
      use fft
      use inform,only:deb_Path
      use mpole
      use pme
      use potent
      use sizes
      use utils
      use utilgpu,only:rec_queue,ngangs_rec
      implicit none
      integer istart,iend,jstart,jend,kstart,kend
      integer istat,ied,jstat,jed,kstat,ked
      integer iichg
      integer i,j,k,m,impi,rankloc
      integer mk,mj
      integer iproc,proc,tsize
      integer ii,jj,kk
      integer isite,iatm
      integer offsetx,offsety
      integer offsetz
      integer location,twonlpts_1,twonlpts_12
      integer nlptsit
      real(t_p) q
      real(t_p) v0,u0,t0
      real(t_p) term,term0
      logical,save:: f_in=.true.
c
      if (use_pmecore) then
        rankloc  = rank_bis
      else
        rankloc  = rank
      end if
      if (deb_Path) write(*,'(3x,a)') "grid_pchg_sitegpu"

      kstat = kstart1(rankloc+1)
      ked   = kend1  (rankloc+1)
      jstat = jstart1(rankloc+1)
      jed   = jend1  (rankloc+1)
      istat = istart1(rankloc+1)
      ied   = iend1  (rankloc+1)
      twonlpts_1  = 2*nlpts+1
      twonlpts_12 = twonlpts_1**2
      nlptsit     = (2*nlpts+1)**3 - 1
c
c     put the permanent multipole moments onto the grid
c
!$acc parallel loop gang worker vector_length(32)
!$acc&         present(chgrecglob,igrid,iion,thetai1_p,thetai2_p)
!$acc&         present(thetai3_p,qgridin_2d)
!$acc&         async(rec_queue)
      do impi = 1,nionrecloc
        iichg   = chgrecglob(impi)
        isite   = iion(iichg)
        q       = pchg(iichg)
        offsetx = 1 - (igrid(1,isite) + grdoff - nlpts)
        offsety = 1 - (igrid(2,isite) + grdoff - nlpts)
        offsetz = 1 - (igrid(3,isite) + grdoff - nlpts)
c
c       put the induced dipole moments onto the grid
c
!$acc loop vector
        do kk = 0, nlptsit
           k  = igrid(3,isite) + grdoff + kk/twonlpts_12 - nlpts
           j  = igrid(2,isite) + grdoff 
     &        + mod(kk/twonlpts_1,twonlpts_1) - nlpts
           i  = igrid(1,isite) + grdoff 
     &        + mod(kk,twonlpts_1)-nlpts
           mk = k + offsetz
           if (k .lt. 1) k = k + nfft3
           mj = j + offsety
           if (j .lt. 1) j = j + nfft2
           m  = i + offsetx
           if (i .lt. 1) i = i + nfft1
c
           v0    = thetai3_p(1,mk,impi)
           u0    = thetai2_p(1,mj,impi) 
           t0    = thetai1_p(1,m ,impi) 
           term  = q*u0*v0*t0

           if (((k.ge.kstat).and.(k.le.ked)).and.
     &         ((j.ge.jstat).and.(j.le.jed)).and.
     &         ((i.ge.istat).and.(i.le.ied))) then
c            location = 1+(i-istat)*2+(j-jstat)*2*n1mpimax
c    &                   +(k-kstat)*2*n1mpimax*n2mpimax
c            call atomic_adds(location,term,qgridin_p(1))
!$acc atomic update
             qgridin_2d(1,i-istat+1,j-jstat+1,k-kstat+1,1) =
     &       qgridin_2d(1,i-istat+1,j-jstat+1,k-kstat+1,1)
     &       + term
             cycle
           end if

!$acc loop seq
           do iproc = 1, nrec_send
              proc   = prec_send(iproc)
              kstart = kstart1(proc+1)
              kend   = kend1  (proc+1)
              jstart = jstart1(proc+1)
              jend   = jend1  (proc+1)
              istart = istart1(proc+1)
              iend   = iend1  (proc+1)
              if (((k.ge.kstart).and.(k.le.kend)).and.
     &            ((j.ge.jstart).and.(j.le.jend)).and.
     &            ((i.ge.istart).and.(i.le.iend))) then
c             location = 1+(i-istart)*2+(j-jstart)*2*n1mpimax
c    &                    +(k-kstart)*2*n1mpimax*n2mpimax
c    &                    +iproc*2*n1mpimax*n2mpimax*n3mpimax
c             call atomic_adds(location,term,qgridin_p(1))
!$acc atomic update
             qgridin_2d(1,i-istart+1,j-jstart+1,k-kstart+1,iproc+1)
     &     = qgridin_2d(1,i-istart+1,j-jstart+1,k-kstart+1,iproc+1)
     &     + term
                 exit
               end if
           end do
        end do
      end do
      end
c
c
c       "grid_uind_site" places the fractional induced dipoles onto the
c       particle mesh Ewald grid
c
c
      subroutine grid_uind_sitegpu(fuindvec,fuinpvec,qgrid2loc)
      use atmlst,only:polerecglob
      use chunks
      use domdec
      use fft
      use mpole
      use pme
      use potent
      use utils
      use utilgpu,only:rec_queue,ngangs_rec
      implicit none
      integer istart,iend,jstart,jend,kstart,kend
      integer istat,ied,jstat,jed,kstat,ked
      integer i,j,k,impi
      integer mk,mj,mi
      integer iproc,proc,rankloc
      integer ii,jj,kk
      integer iipole,iglob,isite,iatm
      integer offsetx,offsety
      integer offsetz
      integer twonlpts_1,twonlpts_12,nlptsit
      integer location
      real(t_p) v0,u0,t0
      real(t_p) v1,u1,t1
      real(t_p) vut(6)
      real(t_p) term
      real(t_p) term01,term11
      real(t_p) term02,term12
      integer igrid1,igrid2,igrid3
      integer nearpt(3),abound(6)
      real(t_p)  fuinp(3),fuind(3)
      real(t_p),intent(in),dimension(3,npolerecloc)::fuindvec,fuinpvec
      real(t_p) qgrid2loc(:,:,:,:,:)
      logical,save::first_in=.true.
!$acc routine(adjust) seq
      if (use_pmecore) then
         rankloc = rank_bis
      else
         rankloc = rank
      end if
      kstat = kstart1(rankloc+1)
      ked   = kend1  (rankloc+1)
      jstat = jstart1(rankloc+1)
      jed   = jend1  (rankloc+1)
      istat = istart1(rankloc+1)
      ied   = iend1  (rankloc+1)
      twonlpts_1  = 2*nlpts+1
      twonlpts_12 = twonlpts_1**2
      nlptsit     = (2*nlpts+1)**3-1

c
c       put the induced dipole moments onto the grid
c
!$acc parallel loop gang worker vector_length(32)
!$acc&         present(polerecglob,igrid,ipole,thetai1,thetai2)
!$acc&         present(fuindvec,fuinpvec,qgrid2loc)
!$acc&         async(rec_queue)
      do impi = 1,npolerecloc
         isite     = ipole(polerecglob(impi))
         offsetx   = 1 - (igrid(1,isite) + grdoff - nlpts)
         offsety   = 1 - (igrid(2,isite) + grdoff - nlpts)
         offsetz   = 1 - (igrid(3,isite) + grdoff - nlpts)
c
c     Three dimensional loop on the grid collapse by hand
c
!$acc loop vector private(vut)
         do kk = 0, nlptsit
            k  = igrid(3,isite) + grdoff + kk/twonlpts_12 - nlpts
            j  = igrid(2,isite) + grdoff 
     &         + mod(kk/twonlpts_1,twonlpts_1) - nlpts
            i  = igrid(1,isite) + grdoff 
     &         + mod(kk,twonlpts_1)-nlpts
            mk     = k + offsetz
            if (k .lt. 1) k = k + nfft3
            mj     = j + offsety
            if (j .lt. 1) j = j + nfft2
            mi     = i + offsetx
            if (i .lt. 1) i = i + nfft1

            vut(1) = thetai3(1,mk,impi)
            vut(2) = thetai3(2,mk,impi)
            vut(3) = thetai2(1,mj,impi)
            vut(4) = thetai2(2,mj,impi)
            vut(5) = thetai1(1,mi,impi)
            vut(6) = thetai1(2,mi,impi)
            term01 = fuindvec(2,impi)*vut(4)*vut(1)
     &             + fuindvec(3,impi)*vut(3)*vut(2)
            term11 = fuindvec(1,impi)*vut(3)*vut(1)
            term02 = fuinpvec(2,impi)*vut(4)*vut(1)
     &             + fuinpvec(3,impi)*vut(3)*vut(2)
            term12 = fuinpvec(1,impi)*vut(3)*vut(1)
            term01 = term01*vut(5) + term11*vut(6)
            term02 = term02*vut(5) + term12*vut(6) 
c
            if (((k.ge.kstat).and.(k.le.ked)).and.
     &          ((j.ge.jstat).and.(j.le.jed)).and.
     &          ((i.ge.istat).and.(i.le.ied)))    then
c              location = 1+(i-istat)*2+(j-jstat)*2*n1mpimax
c    &                     +(k-kstat)*2*n1mpimax*n2mpimax
c              call atomic_Adds2(location,term01,term02,qgrid2in_p(1))
!$acc atomic update
               qgrid2loc(1,i-istat+1,j-jstat+1,k-kstat+1,1)
     &       = qgrid2loc(1,i-istat+1,j-jstat+1,k-kstat+1,1)
     &       + term01
!$acc atomic update
               qgrid2loc(2,i-istat+1,j-jstat+1,k-kstat+1,1)
     &       = qgrid2loc(2,i-istat+1,j-jstat+1,k-kstat+1,1)
     &       + term02
               cycle
            end if

!$acc loop seq
            do iproc = 1, nrec_send
               proc   = prec_send(iproc)
               kstart = kstart1(proc+1)
               kend   = kend1  (proc+1)
               jstart = jstart1(proc+1)
               jend   = jend1  (proc+1)
               istart = istart1(proc+1)
               iend   = iend1  (proc+1)
               if (((k.ge.kstart).and.(k.le.kend)).and.
     &             ((j.ge.jstart).and.(j.le.jend)).and.
     &             ((i.ge.istart).and.(i.le.iend))) then
c              location = 1+(i-istart)*2+(j-jstart)*2*n1mpimax
c    &                     +(k-kstart)*2*n1mpimax*n2mpimax
c              call atomic_Adds2(location,term01,term02,qgrid2in_p(1))
!$acc atomic update
                 qgrid2loc(1,i-istart+1,j-jstart+1,k-kstart+1,iproc+1)
     &         = qgrid2loc(1,i-istart+1,j-jstart+1,k-kstart+1,iproc+1)
     &         + term01
!$acc atomic update
                 qgrid2loc(2,i-istart+1,j-jstart+1,k-kstart+1,iproc+1)
     &         = qgrid2loc(2,i-istart+1,j-jstart+1,k-kstart+1,iproc+1)
     &         + term02
                  exit
               end if
            end do
         end do
      end do
      first_in=.false.
      end
c
c
c       "fphi_mpole_site" extracts the permanent multipole potential on the i-th site from
c       the particle mesh Ewald grid
c
c
      subroutine fphi_mpole_sitegpu
      use atmlst
      use domdec
      use fft
      use inform ,only:deb_Path
      use mpole
      use pme
      use potent
      use utilgpu
      implicit none
      integer istart,iend,jstart,jend,kstart,kend
      integer istarti,iendi,jstarti,jendi,kstarti,kendi
      integer i,j,k,impi,rankloc
      integer iproc,proc
      integer isite,iatm,iipole
      integer i0,j0,k0
      integer it1,it2,it3
      integer igrd0,jgrd0,kgrd0
      real(t_p) v0,v1,v2,v3
      real(t_p) u0,u1,u2,u3
      real(t_p) t0,t1,t2,t3,tq
      real(t_p) tu00,tu10,tu01,tu20,tu11
      real(t_p) tu02,tu21,tu12,tu30,tu03
      real(t_p) tuv000,tuv100,tuv010,tuv001
      real(t_p) tuv200,tuv020,tuv002,tuv110
      real(t_p) tuv101,tuv011,tuv300,tuv030
      real(t_p) tuv003,tuv210,tuv201,tuv120
      real(t_p) tuv021,tuv102,tuv012,tuv111
c
      if (deb_Path) write(*,'(4X,A)') "fphi_mpole_sitegpu"
      if (use_pmecore) then
        rankloc  = rank_bis
      else
        rankloc  = rank
      end if
      kstarti = kstart1(rankloc+1)
      kendi   = kend1  (rankloc+1)
      jstarti = jstart1(rankloc+1)
      jendi   = jend1  (rankloc+1)
      istarti = istart1(rankloc+1)
      iendi   = iend1  (rankloc+1)
c
!$acc parallel loop num_gangs(ngangs_rec)
!$acc&         present(polerecglob,ipole,qgridin_2d,igrid,thetai1,
!$acc&  thetai2,thetai3,kstart1,kend1,jstart1,jend1,kstart1,kend1,
!$acc&  fphirec)
!$acc&         async(rec_queue)
      do impi = 1,npolerecloc
         iipole = polerecglob(impi)
         isite  = ipole(iipole)
         igrd0  = igrid(1,isite)
         jgrd0  = igrid(2,isite)
         kgrd0  = igrid(3,isite)
         tuv000 = 0.0_ti_p
         tuv001 = 0.0_ti_p
         tuv010 = 0.0_ti_p
         tuv100 = 0.0_ti_p
         tuv200 = 0.0_ti_p
         tuv020 = 0.0_ti_p
         tuv002 = 0.0_ti_p
         tuv110 = 0.0_ti_p
         tuv101 = 0.0_ti_p
         tuv011 = 0.0_ti_p
         tuv300 = 0.0_ti_p
         tuv030 = 0.0_ti_p
         tuv003 = 0.0_ti_p
         tuv210 = 0.0_ti_p
         tuv201 = 0.0_ti_p
         tuv120 = 0.0_ti_p
         tuv021 = 0.0_ti_p
         tuv102 = 0.0_ti_p
         tuv012 = 0.0_ti_p
         tuv111 = 0.0_ti_p
         k0     = kgrd0
!$acc loop seq
         do it3 = 1, bsorder
            k0   = k0 + 1
c           k   = k0 + 1 + (nfft3-isign(nfft3,k0))/2
            k   = k0 + 1 + ishft(nfft3-isign(nfft3,k0),-1)
            v0   = thetai3(1,it3,impi)
            v1   = thetai3(2,it3,impi)
            v2   = thetai3(3,it3,impi)
            v3   = thetai3(4,it3,impi)
            tu00 = 0.0_ti_p
            tu10 = 0.0_ti_p
            tu01 = 0.0_ti_p
            tu20 = 0.0_ti_p
            tu11 = 0.0_ti_p
            tu02 = 0.0_ti_p
            tu30 = 0.0_ti_p
            tu21 = 0.0_ti_p
            tu12 = 0.0_ti_p
            tu03 = 0.0_ti_p
            j0   = jgrd0
!$acc loop seq
            do it2 = 1, bsorder
               j0 = j0 + 1
c              j  = j0 + 1 + (nfft2-isign(nfft2,j0))/2
               j  = j0 + 1 + ishft(nfft2-isign(nfft2,j0),-1)
               u0 = thetai2(1,it2,impi)
               u1 = thetai2(2,it2,impi)
               u2 = thetai2(3,it2,impi)
               u3 = thetai2(4,it2,impi)
               t0 = 0.0_ti_p
               t1 = 0.0_ti_p
               t2 = 0.0_ti_p
               t3 = 0.0_ti_p
               i0 = igrd0
!$acc loop seq
               do it1 = 1, bsorder
                  i0 = i0 + 1
c                 i = i0 + 1 + (nfft1-isign(nfft1,i0))/2
                  i = i0 + 1 + ishft(nfft1-isign(nfft1,i0),-1)
c
                  tq     = 0.0_ti_p
                  if (((k.ge.kstarti).and.(k.le.kendi)).and.
     $                ((j.ge.jstarti).and.(j.le.jendi)).and.
     $                ((i.ge.istarti).and.(i.le.iendi))) then
                  tq = qgridin_2d(1,i-istarti+1,j-jstarti+1,
     $                              k-kstarti+1,1)
                    goto 10
                  end if
!$acc loop seq
                  do iproc = 1, nrec_send
                    proc   = prec_send(iproc)
                    kstart = kstart1(proc+1)
                    kend   = kend1  (proc+1)
                    jstart = jstart1(proc+1)
                    jend   = jend1  (proc+1)
                    istart = istart1(proc+1)
                    iend   = iend1  (proc+1)
                    if (((k.ge.kstart).and.(k.le.kend)).and.
     $                  ((j.ge.jstart).and.(j.le.jend)).and.
     $                  ((i.ge.istart).and.(i.le.iend))) then
                    tq  = qgridin_2d(1,i-istart+1,j-jstart+1,k-kstart+1,
     $                 iproc+1)
                      goto 10
                    end if
                  end do
                  cycle
c
 10               continue
                  t0 = t0 + tq*thetai1(1,it1,impi)
                  t1 = t1 + tq*thetai1(2,it1,impi)
                  t2 = t2 + tq*thetai1(3,it1,impi)
                  t3 = t3 + tq*thetai1(4,it1,impi)
               end do
               tu00 = tu00 + t0*u0
               tu10 = tu10 + t1*u0
               tu01 = tu01 + t0*u1
               tu20 = tu20 + t2*u0
               tu11 = tu11 + t1*u1
               tu02 = tu02 + t0*u2
               tu30 = tu30 + t3*u0
               tu21 = tu21 + t2*u1
               tu12 = tu12 + t1*u2
               tu03 = tu03 + t0*u3
            end do
            tuv000 = tuv000 + tu00*v0
            tuv100 = tuv100 + tu10*v0
            tuv010 = tuv010 + tu01*v0
            tuv001 = tuv001 + tu00*v1
            tuv200 = tuv200 + tu20*v0
            tuv020 = tuv020 + tu02*v0
            tuv002 = tuv002 + tu00*v2
            tuv110 = tuv110 + tu11*v0
            tuv101 = tuv101 + tu10*v1
            tuv011 = tuv011 + tu01*v1
            tuv300 = tuv300 + tu30*v0
            tuv030 = tuv030 + tu03*v0
            tuv003 = tuv003 + tu00*v3
            tuv210 = tuv210 + tu21*v0
            tuv201 = tuv201 + tu20*v1
            tuv120 = tuv120 + tu12*v0
            tuv021 = tuv021 + tu02*v1
            tuv102 = tuv102 + tu10*v2
            tuv012 = tuv012 + tu01*v2
            tuv111 = tuv111 + tu11*v1
         end do
         fphirec(1,impi) = tuv000
         fphirec(2,impi) = tuv100
         fphirec(3,impi) = tuv010
         fphirec(4,impi) = tuv001
         fphirec(5,impi) = tuv200
         fphirec(6,impi) = tuv020
         fphirec(7,impi) = tuv002
         fphirec(8,impi) = tuv110
         fphirec(9,impi) = tuv101
         fphirec(10,impi) = tuv011
         fphirec(11,impi) = tuv300
         fphirec(12,impi) = tuv030
         fphirec(13,impi) = tuv003
         fphirec(14,impi) = tuv210
         fphirec(15,impi) = tuv201
         fphirec(16,impi) = tuv120
         fphirec(17,impi) = tuv021
         fphirec(18,impi) = tuv102
         fphirec(19,impi) = tuv012
         fphirec(20,impi) = tuv111
      end do
      end

c
c     "grid_pchg_force"
c     get first derivatives of the reciprocal space energy
c
      subroutine grid_pchg_force
      use atmlst ,only: chgrecglob
      use boxes  ,only: recip
      use charge
      use chgpot
      use domdec ,only: rank,rank_bis,nproc,nrec,ndir,comm_rec
     &           ,COMM_TINKER,nrec_send,locrec1,prec_send,nbloc
      use deriv  ,only: decrec
      use fft
      use inform ,only: deb_Path
      use pme    ,only: igrid,nfft1,nfft2,nfft3,bsorder
     &           ,thetai1_p,thetai2_p,thetai3_p,qgridin_2d
      use potent ,only: use_pmecore
      use utilgpu,only: ngangs_rec,rec_queue
      implicit none
      integer iichg,iglob,iloc,iatm,isite,rankloc
     &       ,igrd0,jgrd0,kgrd0,i,j,k,i0,j0,k0,it1,it2,it3
     &       ,istart,iend,jstart,jend,kstart,kend
     &       ,istartl,iendl,jstartl,jendl,kstartl,kendl
     &       ,iproc,proc,commloc,nprocloc

      real(t_p) f,fi,dn1,dn2,dn3,de1,de2,de3,t1,t2,t3
     &       ,dt1,dt2,dt3,term

      if (use_pmecore) then
        nprocloc = nrec
        commloc  = comm_rec
        rankloc = rank_bis
      else
        nprocloc = nproc
        commloc  = COMM_TINKER
        rankloc = rank
      end if
      if (deb_Path) write(*,'(3x,a)') "grid_pchg_force"

      istartl = istart1(rankloc+1)
      jstartl = jstart1(rankloc+1)
      kstartl = kstart1(rankloc+1)
      iendl   =   iend1(rankloc+1)
      jendl   =   jend1(rankloc+1)
      kendl   =   kend1(rankloc+1)
      f       = electric / dielec
      dn1     = real(nfft1,t_p)
      dn2     = real(nfft2,t_p)
      dn3     = real(nfft3,t_p)

!$acc parallel loop num_gangs(1) vector async(rec_queue) 
!$acc&         present(chgrecglob,iion,locrec1,igrid,
!$acc&     pchg,thetai1_p,thetai2_p,thetai3_p,qgridin_2d,
!$acc&     decrec)
      do isite = 1, nionrecloc
        iichg = chgrecglob(isite)
        iglob = iion(iichg)
        iloc  = locrec1(iglob)
        iatm  = iglob
        igrd0 = igrid(1,iatm)
        jgrd0 = igrid(2,iatm)
        kgrd0 = igrid(3,iatm)
        fi    = f * pchg(iichg)
        de1   = 0.0_ti_p
        de2   = 0.0_ti_p
        de3   = 0.0_ti_p
        k0    = kgrd0
!$acc loop seq
        do it3 = 1, bsorder
           k0  = k0 + 1
           k   = k0 + 1 + (nfft3-sign(nfft3,k0))/2
           t3  = thetai3_p(1,it3,isite)
           dt3 = dn3 * thetai3_p(2,it3,isite)
           j0  = jgrd0
!$acc loop seq
           do it2 = 1, bsorder
              j0  = j0 + 1
              j   = j0 + 1 + (nfft2-sign(nfft2,j0))/2
              t2  = thetai2_p(1,it2,isite)
              dt2 = dn2 * thetai2_p(2,it2,isite)
              i0  = igrd0
!$acc loop seq
              do it1 = 1, bsorder
                 i0  = i0 + 1
                 i   = i0 + 1 + (nfft1-sign(nfft1,i0))/2
                 t1  = thetai1_p(1,it1,isite)
                 dt1 = dn1 * thetai1_p(2,it1,isite)
c
                 if (((k.ge.kstartl).and.(k.le.kendl)).and.
     $               ((j.ge.jstartl).and.(j.le.jendl)).and.
     $               ((i.ge.istartl).and.(i.le.iendl)))  then
                    term = qgridin_2d(1,i-istartl+1,j-jstartl+1,
     $                                  k-kstartl+1,1)
                    goto 100
                 end if
c
!$acc loop seq
                 do iproc = 1, nrec_send
                    proc   = prec_send(iproc)
                    kstart = kstart1(proc+1)
                    kend   = kend1(proc+1)
                    jstart = jstart1(proc+1)
                    jend   = jend1(proc+1)
                    istart = istart1(proc+1)
                    iend   = iend1(proc+1)
                    if (((k.ge.kstart).and.(k.le.kend)).and.
     $                  ((j.ge.jstart).and.(j.le.jend)).and.
     $                  ((i.ge.istart).and.(i.le.iend))) then
                       term=qgridin_2d(1,i-istart+1,j-jstart+1
     $                                  ,k-kstart+1,iproc+1)
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
     &                                    +recip(1,3)*de3)
        decrec(2,iloc) =decrec(2,iloc)+fi*(recip(2,1)*de1+recip(2,2)*de2
     &                                    +recip(2,3)*de3)
        decrec(3,iloc) =decrec(3,iloc)+fi*(recip(3,1)*de1+recip(3,2)*de2
     &                                    +recip(3,3)*de3)
      end do

      end subroutine
c
c     "fphi_uind_site" extracts the induced dipole potential at the i-th site from
c     the particle mesh Ewald grid
c
c
      subroutine fphi_uind_sitegpu1(fdip_phi1,fdip_phi2,fdip_sum_phi)
      use atmlst
      use domdec
      use fft
      use mpole, only : npolerecloc, ipole
      use pme
      use potent
      use utilgpu,only:rec_queue,ngangs_rec
      implicit none
      integer istart,iend,jstart,jend,kstart,kend
      integer istat,ied,jstat,jed,kstat,ked
      integer i,j,k,impi
      integer iproc,proc
      integer isite,iatm,iipole
      integer i0,j0,k0
      integer it1,it2,it3
      integer igrd0,jgrd0,kgrd0
      real(t_p) v0,v1,v2,v3
      real(t_p) u0,u1,u2,u3
      real(t_p) t0,t1,t2,t3
      real(t_p) t0_1,t0_2,t1_1,t1_2
      real(t_p) t2_1,t2_2,tq_1,tq_2
      real(t_p) tu00,tu10,tu01,tu20,tu11
      real(t_p) tu02,tu30,tu21,tu12,tu03
      real(t_p) tu00_1,tu01_1,tu10_1
      real(t_p) tu00_2,tu01_2,tu10_2
      real(t_p) tu20_1,tu11_1,tu02_1
      real(t_p) tu20_2,tu11_2,tu02_2
      real(t_p) tuv100_1,tuv010_1,tuv001_1
      real(t_p) tuv100_2,tuv010_2,tuv001_2
      real(t_p) tuv200_1,tuv020_1,tuv002_1
      real(t_p) tuv110_1,tuv101_1,tuv011_1
      real(t_p) tuv200_2,tuv020_2,tuv002_2
      real(t_p) tuv110_2,tuv101_2,tuv011_2
      real(t_p) tuv000,tuv100,tuv010,tuv001
      real(t_p) tuv200,tuv020,tuv002,tuv110
      real(t_p) tuv101,tuv011,tuv300,tuv030
      real(t_p) tuv003,tuv210,tuv201,tuv120
      real(t_p) tuv021,tuv102,tuv012,tuv111
      real(t_p),dimension(:,:)::fdip_phi1,fdip_phi2
      real(t_p) fdip_sum_phi(:,:)

      if (use_pmecore) then
        kstat = kstart1(rank_bis+1)
        ked   = kend1  (rank_bis+1)
        jstat = jstart1(rank_bis+1)
        jed   = jend1  (rank_bis+1)
        istat = istart1(rank_bis+1)
        ied   = iend1  (rank_bis+1)
      else
        kstat = kstart1(rank+1)
        ked   = kend1  (rank+1)
        jstat = jstart1(rank+1)
        jed   = jend1  (rank+1)
        istat = istart1(rank+1)
        ied   = iend1  (rank+1)
      end if
c
!$acc parallel loop num_gangs(ngangs_rec)
!$acc&         present(polerecglob,ipole,igrid,thetai1,thetai2,
!$acc&   kstart1,kend1,jstart1,jend1,istart1,iend1,fdip_phi1,
!$acc&   fdip_phi2,fdip_sum_phi)
!$acc&         async(rec_queue)
      do impi = 1,npolerecloc
         iipole   = polerecglob(impi)
         iatm     = ipole(iipole)
         igrd0    = igrid(1,iatm)
         jgrd0    = igrid(2,iatm)
         kgrd0    = igrid(3,iatm)
         tuv100_1 = 0.0_ti_p
         tuv010_1 = 0.0_ti_p
         tuv001_1 = 0.0_ti_p
         tuv200_1 = 0.0_ti_p
         tuv020_1 = 0.0_ti_p
         tuv002_1 = 0.0_ti_p
         tuv110_1 = 0.0_ti_p
         tuv101_1 = 0.0_ti_p
         tuv011_1 = 0.0_ti_p
         tuv100_2 = 0.0_ti_p
         tuv010_2 = 0.0_ti_p
         tuv001_2 = 0.0_ti_p
         tuv200_2 = 0.0_ti_p
         tuv020_2 = 0.0_ti_p
         tuv002_2 = 0.0_ti_p
         tuv110_2 = 0.0_ti_p
         tuv101_2 = 0.0_ti_p
         tuv011_2 = 0.0_ti_p
         tuv000   = 0.0_ti_p
         tuv001   = 0.0_ti_p
         tuv010   = 0.0_ti_p
         tuv100   = 0.0_ti_p
         tuv200   = 0.0_ti_p
         tuv020   = 0.0_ti_p
         tuv002   = 0.0_ti_p
         tuv110   = 0.0_ti_p
         tuv101   = 0.0_ti_p
         tuv011   = 0.0_ti_p
         tuv300   = 0.0_ti_p
         tuv030   = 0.0_ti_p
         tuv003   = 0.0_ti_p
         tuv210   = 0.0_ti_p
         tuv201   = 0.0_ti_p
         tuv120   = 0.0_ti_p
         tuv021   = 0.0_ti_p
         tuv102   = 0.0_ti_p
         tuv012   = 0.0_ti_p
         tuv111   = 0.0_ti_p
         k0       = kgrd0
!$acc    loop seq
         do it3 = 1, bsorder
            k0 = k0 + 1
c           k      = k0 + 1 + (nfft3-isign(nfft3,k0))/2
            k      = k0 + 1 + ishft(nfft3-isign(nfft3,k0),-1)
            v0     = thetai3(1,it3,impi)
            v1     = thetai3(2,it3,impi)
            v2     = thetai3(3,it3,impi)
            v3     = thetai3(4,it3,impi)
            tu00_1 = 0.0_ti_p
            tu01_1 = 0.0_ti_p
            tu10_1 = 0.0_ti_p
            tu20_1 = 0.0_ti_p
            tu11_1 = 0.0_ti_p
            tu02_1 = 0.0_ti_p
            tu00_2 = 0.0_ti_p
            tu01_2 = 0.0_ti_p
            tu10_2 = 0.0_ti_p
            tu20_2 = 0.0_ti_p
            tu11_2 = 0.0_ti_p
            tu02_2 = 0.0_ti_p
            tu00   = 0.0_ti_p
            tu10   = 0.0_ti_p
            tu01   = 0.0_ti_p
            tu20   = 0.0_ti_p
            tu11   = 0.0_ti_p
            tu02   = 0.0_ti_p
            tu30   = 0.0_ti_p
            tu21   = 0.0_ti_p
            tu12   = 0.0_ti_p
            tu03   = 0.0_ti_p
            j0 = jgrd0
!$acc    loop seq
            do it2 = 1, bsorder
               j0 = j0 + 1
c              j    = j0 + 1 + (nfft2-isign(nfft2,j0))/2
               j    = j0 + 1 + ishft(nfft2-isign(nfft2,j0),-1)
               u0   = thetai2(1,it2,impi)
               u1   = thetai2(2,it2,impi)
               u2   = thetai2(3,it2,impi)
               u3   = thetai2(4,it2,impi)
               t0_1 = 0.0_ti_p
               t1_1 = 0.0_ti_p
               t2_1 = 0.0_ti_p
               t0_2 = 0.0_ti_p
               t1_2 = 0.0_ti_p
               t2_2 = 0.0_ti_p
               t3   = 0.0_ti_p
               i0   = igrd0
!$acc    loop seq
               do it1 = 1, bsorder
                  i0 = i0 + 1
c                 i    = i0 + 1 + (nfft1-isign(nfft1,i0))/2
                  i    = i0 + 1 + ishft(nfft1-isign(nfft1,i0),-1)
c      
                  tq_1 = 0_ti_p
                  tq_2 = 0_ti_p
                  if (((k.ge.kstat).and.(k.le.ked)).and.
     &                ((j.ge.jstat).and.(j.le.jed)).and.
     &                ((i.ge.istat).and.(i.le.ied))) then
                tq_1 = qgrid2in_2d(1,i-istat+1,j-jstat+1,k-kstat+1,1)
                tq_2 = qgrid2in_2d(2,i-istat+1,j-jstat+1,k-kstat+1,1)
                    goto 10
                  end if
!$acc    loop seq
                  do iproc = 1, nrec_send
                    proc   = prec_send(iproc)
                    kstart = kstart1(proc+1)
                    kend   = kend1  (proc+1)
                    jstart = jstart1(proc+1)
                    jend   = jend1  (proc+1)
                    istart = istart1(proc+1)
                    iend   = iend1  (proc+1)
                    if (((k.ge.kstart).and.(k.le.kend)).and.
     &                  ((j.ge.jstart).and.(j.le.jend)).and.
     &                  ((i.ge.istart).and.(i.le.iend))) then
                tq_1 = qgrid2in_2d(1,i-istart+1,j-jstart+1,k-kstart+1,
     &                     iproc+1)
                tq_2 = qgrid2in_2d(2,i-istart+1,j-jstart+1,k-kstart+1,
     &                     iproc+1)
                      goto 10
                    end if
                  end do
                  cycle
 10               continue
c
                  t0_1 = t0_1 +  tq_1      *thetai1(1,it1,impi)
                  t1_1 = t1_1 +  tq_1      *thetai1(2,it1,impi)
                  t2_1 = t2_1 +  tq_1      *thetai1(3,it1,impi)
                  t0_2 = t0_2 +        tq_2*thetai1(1,it1,impi)
                  t1_2 = t1_2 +        tq_2*thetai1(2,it1,impi)
                  t2_2 = t2_2 +        tq_2*thetai1(3,it1,impi)
                  t3   = t3   + (tq_1+tq_2)*thetai1(4,it1,impi)
               end do
               tu00_1 = tu00_1 + t0_1*u0
               tu10_1 = tu10_1 + t1_1*u0
               tu01_1 = tu01_1 + t0_1*u1
               tu20_1 = tu20_1 + t2_1*u0
               tu11_1 = tu11_1 + t1_1*u1
               tu02_1 = tu02_1 + t0_1*u2
               tu00_2 = tu00_2 + t0_2*u0
               tu10_2 = tu10_2 + t1_2*u0
               tu01_2 = tu01_2 + t0_2*u1
               tu20_2 = tu20_2 + t2_2*u0
               tu11_2 = tu11_2 + t1_2*u1
               tu02_2 = tu02_2 + t0_2*u2
               t0     = t0_1 + t0_2
               t1     = t1_1 + t1_2
               t2     = t2_1 + t2_2
               tu00   = tu00 + t0*u0
               tu10   = tu10 + t1*u0
               tu01   = tu01 + t0*u1
               tu20   = tu20 + t2*u0
               tu11   = tu11 + t1*u1
               tu02   = tu02 + t0*u2
               tu30   = tu30 + t3*u0
               tu21   = tu21 + t2*u1
               tu12   = tu12 + t1*u2
               tu03   = tu03 + t0*u3
            end do
            tuv100_1 = tuv100_1 + tu10_1*v0
            tuv010_1 = tuv010_1 + tu01_1*v0
            tuv001_1 = tuv001_1 + tu00_1*v1
            tuv200_1 = tuv200_1 + tu20_1*v0
            tuv020_1 = tuv020_1 + tu02_1*v0
            tuv002_1 = tuv002_1 + tu00_1*v2
            tuv110_1 = tuv110_1 + tu11_1*v0
            tuv101_1 = tuv101_1 + tu10_1*v1
            tuv011_1 = tuv011_1 + tu01_1*v1
            tuv100_2 = tuv100_2 + tu10_2*v0
            tuv010_2 = tuv010_2 + tu01_2*v0
            tuv001_2 = tuv001_2 + tu00_2*v1
            tuv200_2 = tuv200_2 + tu20_2*v0
            tuv020_2 = tuv020_2 + tu02_2*v0
            tuv002_2 = tuv002_2 + tu00_2*v2
            tuv110_2 = tuv110_2 + tu11_2*v0
            tuv101_2 = tuv101_2 + tu10_2*v1
            tuv011_2 = tuv011_2 + tu01_2*v1
            tuv000   = tuv000 + tu00*v0
            tuv100   = tuv100 + tu10*v0
            tuv010   = tuv010 + tu01*v0
            tuv001   = tuv001 + tu00*v1
            tuv200   = tuv200 + tu20*v0
            tuv020   = tuv020 + tu02*v0
            tuv002   = tuv002 + tu00*v2
            tuv110   = tuv110 + tu11*v0
            tuv101   = tuv101 + tu10*v1
            tuv011   = tuv011 + tu01*v1
            tuv300   = tuv300 + tu30*v0
            tuv030   = tuv030 + tu03*v0
            tuv003   = tuv003 + tu00*v3
            tuv210   = tuv210 + tu21*v0
            tuv201   = tuv201 + tu20*v1
            tuv120   = tuv120 + tu12*v0
            tuv021   = tuv021 + tu02*v1
            tuv102   = tuv102 + tu10*v2
            tuv012   = tuv012 + tu01*v2
            tuv111   = tuv111 + tu11*v1
         end do
         fdip_phi1   ( 2,impi) = tuv100_1
         fdip_phi1   ( 3,impi) = tuv010_1
         fdip_phi1   ( 4,impi) = tuv001_1
         fdip_phi1   ( 5,impi) = tuv200_1
         fdip_phi1   ( 6,impi) = tuv020_1
         fdip_phi1   ( 7,impi) = tuv002_1
         fdip_phi1   ( 8,impi) = tuv110_1
         fdip_phi1   ( 9,impi) = tuv101_1
         fdip_phi1   (10,impi) = tuv011_1
         fdip_phi2   ( 2,impi) = tuv100_2
         fdip_phi2   ( 3,impi) = tuv010_2
         fdip_phi2   ( 4,impi) = tuv001_2
         fdip_phi2   ( 5,impi) = tuv200_2
         fdip_phi2   ( 6,impi) = tuv020_2
         fdip_phi2   ( 7,impi) = tuv002_2
         fdip_phi2   ( 8,impi) = tuv110_2
         fdip_phi2   ( 9,impi) = tuv101_2
         fdip_phi2   (10,impi) = tuv011_2
         fdip_sum_phi( 1,impi) = tuv000
         fdip_sum_phi( 2,impi) = tuv100
         fdip_sum_phi( 3,impi) = tuv010
         fdip_sum_phi( 4,impi) = tuv001
         fdip_sum_phi( 5,impi) = tuv200
         fdip_sum_phi( 6,impi) = tuv020
         fdip_sum_phi( 7,impi) = tuv002
         fdip_sum_phi( 8,impi) = tuv110
         fdip_sum_phi( 9,impi) = tuv101
         fdip_sum_phi(10,impi) = tuv011
         fdip_sum_phi(11,impi) = tuv300
         fdip_sum_phi(12,impi) = tuv030
         fdip_sum_phi(13,impi) = tuv003
         fdip_sum_phi(14,impi) = tuv210
         fdip_sum_phi(15,impi) = tuv201
         fdip_sum_phi(16,impi) = tuv120
         fdip_sum_phi(17,impi) = tuv021
         fdip_sum_phi(18,impi) = tuv102
         fdip_sum_phi(19,impi) = tuv012
         fdip_sum_phi(20,impi) = tuv111
      end do
c
      return
      end
c
c     "fphi_uind_sitegpu2" extracts the induced dipole potential at the i-th site from
c     the particle mesh Ewald grid second version
c
c
      subroutine fphi_uind_sitegpu2(fdip_phi1,fdip_phi2)
      use atoms
      use atmlst
      use boxes
      use domdec
      use fft
      use mpole, only : npolerecloc, ipole
      use pme
      use potent
      use utils
      use utilgpu,only:rec_queue,ngangs_rec
      implicit none
      real(t_p),dimension(:,:)::fdip_phi1,fdip_phi2

      integer istart,iend,jstart,jend,kstart,kend
      integer istat,ied,jstat,jed,kstat,ked
      integer i,j,k,impi,ifr
      integer iproc,proc
      integer isite,iatm,iipole
      integer i0,j0,k0
      integer it1,it2,it3
      integer igrd0,jgrd0,kgrd0
      integer ind
      real(t_p) xi,yi,zi
      real(t_p) w,fr,eps
      real(t_p) v0,v1,v2,v3
      real(t_p) u0,u1,u2,u3
      real(t_p) t0,t1,t2,t3
      real(t_p) t0_1,t0_2,t1_1,t1_2
      real(t_p) t2_1,t2_2,tq_1,tq_2
      real(t_p) tu00,tu10,tu01,tu20,tu11
      real(t_p) tu02,tu30,tu21,tu12,tu03
      real(t_p) tu00_1,tu01_1,tu10_1
      real(t_p) tu00_2,tu01_2,tu10_2
      real(t_p) tu20_1,tu11_1,tu02_1
      real(t_p) tu20_2,tu11_2,tu02_2
      real(t_p) tuv100_1,tuv010_1,tuv001_1
      real(t_p) tuv100_2,tuv010_2,tuv001_2
      real(t_p) tuv200_1,tuv020_1,tuv002_1
      real(t_p) tuv110_1,tuv101_1,tuv011_1
      real(t_p) tuv200_2,tuv020_2,tuv002_2
      real(t_p) tuv110_2,tuv101_2,tuv011_2
      real(t_p) tuv000,tuv100,tuv010,tuv001
      real(t_p) tuv200,tuv020,tuv002,tuv110
      real(t_p) tuv101,tuv011,tuv300,tuv030
      real(t_p) tuv003,tuv210,tuv201,tuv120
      real(t_p) tuv021,tuv102,tuv012,tuv111
      real(t_p) tuv1(9),tuv2(9)

      if (use_pmecore) then
        kstat = kstart1(rank_bis+1)
        ked   = kend1  (rank_bis+1)
        jstat = jstart1(rank_bis+1)
        jed   = jend1  (rank_bis+1)
        istat = istart1(rank_bis+1)
        ied   = iend1  (rank_bis+1)
      else
        kstat = kstart1(rank+1)
        ked   = kend1  (rank+1)
        jstat = jstart1(rank+1)
        jed   = jend1  (rank+1)
        istat = istart1(rank+1)
        ied   = iend1  (rank+1)
      end if
c
!$acc parallel loop num_gangs(ngangs_rec)
!$acc&         present(polerecglob,ipole,igrid,thetai1,thetai2,
!$acc&   kstart1,kend1,jstart1,jend1,istart1,iend1)
!$acc&         present(fdip_phi1,fdip_phi2,qgrid2in_p)
!$acc&         private(tuv1,tuv2)
!$acc&         async(rec_queue)
      do impi = 1,npolerecloc
         iipole = polerecglob(impi)
         iatm   = ipole(iipole)
         igrd0  = igrid(1,iatm)
         jgrd0  = igrid(2,iatm)
         kgrd0  = igrid(3,iatm)

!$acc loop seq
         do it1 = 1,9
            tuv1(it1)  = 0.0_ti_p
            tuv2(it1)  = 0.0_ti_p
         end do
c        k0       = kgrd0
!$acc    loop seq
         do it3 = 1, bsorder
c           k0     = k0 + 1
c           k      = k0 + 1 + (nfft3-isign(nfft3,kgrd0+it3))/2
            k      = kgrd0 +it3 +1 + 
     &               ishft(nfft3-isign(nfft3,kgrd0+it3),-1)
            v0     = thetai3(1,it3,impi)
            v1     = thetai3(2,it3,impi)
            v2     = thetai3(3,it3,impi)
            tu00_1 = 0.0_ti_p
            tu01_1 = 0.0_ti_p
            tu10_1 = 0.0_ti_p
            tu20_1 = 0.0_ti_p
            tu11_1 = 0.0_ti_p
            tu02_1 = 0.0_ti_p
            tu00_2 = 0.0_ti_p
            tu01_2 = 0.0_ti_p
            tu10_2 = 0.0_ti_p
            tu20_2 = 0.0_ti_p
            tu11_2 = 0.0_ti_p
            tu02_2 = 0.0_ti_p
c           j0 = jgrd0
!$acc    loop seq
            do it2 = 1, bsorder
c              j0 = j0 + 1
c              j    = jgrd0 +it2 + 1 + (nfft2-isign(nfft2,jgrd0+it2))/2
               j    = jgrd0 +it2 + 1 +
     &                      ishft(nfft2-isign(nfft2,jgrd0+it2),-1)
               u0   = thetai2(1,it2,impi)
               u1   = thetai2(2,it2,impi)
               u2   = thetai2(3,it2,impi)
               t0_1 = 0.0_ti_p
               t1_1 = 0.0_ti_p
               t2_1 = 0.0_ti_p
               t0_2 = 0.0_ti_p
               t1_2 = 0.0_ti_p
               t2_2 = 0.0_ti_p
               i0   = igrd0
!$acc loop seq
               do it1 = 1, bsorder
c                 i0  = i0 + 1
c                 i    = igrd0 +it1 +1 + (nfft1-isign(nfft1,igrd0 +it1))/2
                  i    = igrd0 +it1 +1 +
     &                         ishft(nfft1-isign(nfft1,igrd0+it1),-1)
c
                  tq_1 = 0_ti_p
                  tq_2 = 0_ti_p
                  if (((k.ge.kstat).and.(k.le.ked)).and.
     &                ((j.ge.jstat).and.(j.le.jed)).and.
     &                ((i.ge.istat).and.(i.le.ied))) then
c               ind = 1+(i-istat)*2+(j-jstat)*2*n1mpimax+(k-kstat)*
c    &                2*n1mpimax*n2mpimax
c               tq_1 = qgrid2in_p(ind)
c               tq_2 = qgrid2in_p(ind+1)
                tq_1 = qgrid2in_2d(1,i-istat+1,j-jstat+1,k-kstat+1,1)
                tq_2 = qgrid2in_2d(2,i-istat+1,j-jstat+1,k-kstat+1,1)
                    goto 10
                  end if
!$acc loop seq
                  do iproc = 1, nrec_send
                    proc   = prec_send(iproc)
                    kstart = kstart1(proc+1)
                    kend   = kend1  (proc+1)
                    jstart = jstart1(proc+1)
                    jend   = jend1  (proc+1)
                    istart = istart1(proc+1)
                    iend   = iend1  (proc+1)
                    if (((k.ge.kstart).and.(k.le.kend)).and.
     &                  ((j.ge.jstart).and.(j.le.jend)).and.
     &                  ((i.ge.istart).and.(i.le.iend))) then
                tq_1 = qgrid2in_2d(1,i-istart+1,j-jstart+1,k-kstart+1,
     &                     iproc+1)
                tq_2 = qgrid2in_2d(2,i-istart+1,j-jstart+1,k-kstart+1,
     &                     iproc+1)
                      goto 10
                    end if
                  end do
                  cycle
 10               continue
c
                  t0   = thetai1(1,it1,impi)
                  t1   = thetai1(2,it1,impi)
                  t2   = thetai1(3,it1,impi)
                  t0_1 = t0_1 + tq_1*t0
                  t1_1 = t1_1 + tq_1*t1
                  t2_1 = t2_1 + tq_1*t2
                  t0_2 = t0_2 + tq_2*t0
                  t1_2 = t1_2 + tq_2*t1
                  t2_2 = t2_2 + tq_2*t2
               end do
               tu00_1 = tu00_1 + t0_1*u0
               tu10_1 = tu10_1 + t1_1*u0
               tu01_1 = tu01_1 + t0_1*u1
               tu20_1 = tu20_1 + t2_1*u0
               tu11_1 = tu11_1 + t1_1*u1
               tu02_1 = tu02_1 + t0_1*u2
               tu00_2 = tu00_2 + t0_2*u0
               tu10_2 = tu10_2 + t1_2*u0
               tu01_2 = tu01_2 + t0_2*u1
               tu20_2 = tu20_2 + t2_2*u0
               tu11_2 = tu11_2 + t1_2*u1
               tu02_2 = tu02_2 + t0_2*u2
            end do
            tuv1(1) = tuv1(1) + tu10_1*v0
            tuv1(2) = tuv1(2) + tu01_1*v0
            tuv1(3) = tuv1(3) + tu00_1*v1
            tuv1(4) = tuv1(4) + tu20_1*v0
            tuv1(5) = tuv1(5) + tu02_1*v0
            tuv1(6) = tuv1(6) + tu00_1*v2
            tuv1(7) = tuv1(7) + tu11_1*v0
            tuv1(8) = tuv1(8) + tu10_1*v1
            tuv1(9) = tuv1(9) + tu01_1*v1
            tuv2(1) = tuv2(1) + tu10_2*v0
            tuv2(2) = tuv2(2) + tu01_2*v0
            tuv2(3) = tuv2(3) + tu00_2*v1
            tuv2(4) = tuv2(4) + tu20_2*v0
            tuv2(5) = tuv2(5) + tu02_2*v0
            tuv2(6) = tuv2(6) + tu00_2*v2
            tuv2(7) = tuv2(7) + tu11_2*v0
            tuv2(8) = tuv2(8) + tu10_2*v1
            tuv2(9) = tuv2(9) + tu01_2*v1
         end do
c        call set_dip(impi,tuv1,fdip_phi1(1,1))
c        call set_dip(impi,tuv2,fdip_phi2(1,1))
!$acc loop seq
         do i = 2,10
            fdip_phi1(i,impi) = tuv1(i-1)
            fdip_phi2(i,impi) = tuv2(i-1)
         end do
      end do
      end

#ifdef _CUDA
      subroutine grid_mpole_sitecu(fmpvec)
      use atmlst    ,only:polerecglob
      use chunks
      use domdec
      use fft
      use mpole
      use neigh  ,only:celle_glob
      use pme
      use potent
      use utils
      use pmestuffcu,only:grid_mpole_sitecu_core_1p
      use utilcu    ,only:PME_BLOCK_DIM1,check_launch_kernel
      use utilgpu   ,only:rec_stream,rec_queue,nSMP
      implicit none
      integer istat,ied,jstat,jed,kstat,ked
      integer i,j,k,impi
      integer twonlpts_1,twonlpts_12,nlptsit,rankloc
      real(t_p) fmpvec(10,max(npolerecloc,1))
      integer,save::gS
      logical,save::first_in=.true.
!$acc routine(adjust) seq
      if (use_pmecore) then
         rankloc = rank_bis
      else
         rankloc = rank
      end if
      kstat = kstart1(rankloc+1)
      ked   = kend1  (rankloc+1)
      jstat = jstart1(rankloc+1)
      jed   = jend1  (rankloc+1)
      istat = istart1(rankloc+1)
      ied   = iend1  (rankloc+1)
      twonlpts_1  = 2*nlpts+1
      twonlpts_12 = twonlpts_1**2
      nlptsit     = (2*nlpts+1)**3-1

      if (first_in) then
         call cudaMaxgridsize("grid_mpole_core",gS)
      end if
      call set_pme_texture
c
c     put the permanent multipole moments onto the grid
c     Work only with 1 MPI process
c
!$acc host_data use_device(polerecglob,celle_glob,ipole,igrid,
!$acc&        thetai1,thetai2,thetai3,fmpvec,qgridin_2d)

      call grid_mpole_sitecu_core_1p<<<gS,PME_BLOCK_DIM1,0,rec_stream>>>
     &     (celle_glob,ipole,igrid,thetai1,thetai2,thetai3
     &     ,kstat,ked,jstat,jed,istat,ied
     &     ,nfft1,nfft2,nfft3,npolerecloc
     &     ,nlpts,twonlpts_1,twonlpts_12,nlptsit,grdoff
     &     ,bsorder,n1mpimax,n2mpimax,n3mpimax,nrec_send
     &     ,fmpvec,qgridin_2d,first_in)
      call check_launch_kernel(" grid_mpole_sitecu_core_1p")

!$acc end host_data
      if (first_in) first_in=.false.
      end subroutine
c
c     "grid_pchg_sitecu" places the i-th fractional atomic charge onto
c     the particle mesh Ewald grid (CUDA Routine)
c
      subroutine grid_pchg_sitecu
      use atmlst  ,only:chgrecglob
      use charge
      use chunks
      use domdec
      use fft
      use inform  ,only:deb_Path
      use mpole
      use neigh   ,chg_s=>celle_chg,glob_s=>celle_glob
      use pme
      use pmestuffcu,only:grid_pchg_sitecu_core_1p
      use potent
      use sizes
      use utils
      use utilcu  ,only:PME_BLOCK_DIM1,check_launch_kernel
      use utilgpu ,only:rec_queue,rec_stream,ngangs_rec
      implicit none
      integer istart,iend,jstart,jend,kstart,kend
      integer istat,ied,jstat,jed,kstat,ked
      integer iichg
      integer i,j,k,m,impi,rankloc
      integer mk,mj
      integer iproc,proc,tsize
      integer ii,jj,kk
      integer isite,iatm
      integer offsetx,offsety
      integer offsetz
      integer location,twonlpts_1,twonlpts_12
      integer nlptsit
      real(t_p) q
      real(t_p) v0,u0,t0
      real(t_p) term,term0
      logical,save::f_in=.true.
      integer,save::gS
c
      if (use_pmecore) then
        rankloc  = rank_bis
      else
        rankloc  = rank
      end if
      if (deb_Path) write(*,'(3x,a)') "grid_pchg_sitecu"

      kstat = kstart1(rankloc+1)
      ked   = kend1  (rankloc+1)
      jstat = jstart1(rankloc+1)
      jed   = jend1  (rankloc+1)
      istat = istart1(rankloc+1)
      ied   = iend1  (rankloc+1)
      twonlpts_1  = 2*nlpts+1
      twonlpts_12 = twonlpts_1**2
      nlptsit     = (2*nlpts+1)**3 - 1

      if (f_in) then
         call cudaMaxgridsize("grid_pchg_core",gS)
      end if
      call set_pme_texture
c
c     put the permanent multipole moments onto the grid
c
      if (nproc.eq.1.or.nrec.eq.1) then
!$acc host_data use_device(chg_s,glob_s,igrid,pchg,thetai1,thetai2,
!$acc&    thetai3,qgridin_2d)
      call grid_pchg_sitecu_core_1p<<<gS,PME_BLOCK_DIM1,0,rec_stream>>>
     &     (chg_s,glob_s,igrid,pchg,thetai1,thetai2,thetai3
     &     ,kstat,ked,jstat,jed,istat,ied
     &     ,nfft1,nfft2,nfft3,nionrecloc
     &     ,nlpts,twonlpts_1,twonlpts_12,nlptsit,grdoff
     &     ,bsorder,n1mpimax,n2mpimax,n3mpimax,nrec_send
     &     ,qgridin_2d,f_in)
      call check_launch_kernel(" grid_uind_sitecu_core_1p")
!$acc end host_data
      else
         print*, 'FATAL ERROR, grid_pchg_sitecu is sequential specific'
         call fatal
      end if
      end subroutine

      subroutine grid_uind_sitecu(fuindvec,fuinpvec,qgrid2loc)
      use atmlst    ,only:polerecglob
      use chunks
      use domdec
      use fft
      use mpole
      use neigh  ,only:celle_glob
      use pme
      use potent
      use utils
      use pmestuffcu,only:grid_uind_sitecu_core_1p
      use utilcu    ,only:PME_BLOCK_DIM1,check_launch_kernel
      use utilgpu   ,only:rec_stream,rec_queue,nSMP
      implicit none
      integer istat,ied,jstat,jed,kstat,ked
      integer i,j,k,impi
      integer twonlpts_1,twonlpts_12,nlptsit,rankloc
      real(t_p),intent(in),dimension(3,npolerecloc)::fuindvec,fuinpvec
      real(t_p) qgrid2loc(:,:,:,:,:)
      integer,save::gS
      logical,save::first_in=.true.
!$acc routine(adjust) seq
      if (use_pmecore) then
         rankloc = rank_bis
      else
         rankloc = rank
      end if
      kstat = kstart1(rankloc+1)
      ked   = kend1  (rankloc+1)
      jstat = jstart1(rankloc+1)
      jed   = jend1  (rankloc+1)
      istat = istart1(rankloc+1)
      ied   = iend1  (rankloc+1)
      twonlpts_1  = 2*nlpts+1
      twonlpts_12 = twonlpts_1**2
      nlptsit     = (2*nlpts+1)**3-1

      if (first_in) then
         call cudaMaxgridsize("grid_uind_core",gS)
      end if
      call set_pme_texture
c
c     put the induced dipole moments onto the grid
c     Work only with 1 MPI process
c
!$acc host_data use_device(celle_glob,ipole,igrid,thetai1,thetai2
!$acc&        ,thetai3,fuindvec,fuinpvec,qgrid2loc)

      call grid_uind_sitecu_core_1p<<<gS,PME_BLOCK_DIM1,0,rec_stream>>>
     &     (celle_glob,ipole,igrid,thetai1,thetai2,thetai3
     &     ,kstat,ked,jstat,jed,istat,ied
     &     ,nfft1,nfft2,nfft3,npolerecloc
     &     ,nlpts,twonlpts_1,twonlpts_12,nlptsit,grdoff
     &     ,bsorder,n1mpimax,n2mpimax,n3mpimax,nrec_send
     &     ,fuindvec,fuinpvec,qgrid2loc,first_in)
      call check_launch_kernel(" grid_uind_sitecu_core_1p")

!$acc end host_data
      if (first_in) first_in=.false.
      end subroutine

      subroutine fphi_mpole_sitecu
      use atoms
      use atmlst
      use boxes
      use domdec
      use fft
      use mpole
      use pme
      use potent
      use pmestuffcu,only:fphi_mpole_core
      use utilcu    ,only:PME_BLOCK_DIM,check_launch_kernel
      use utilgpu   ,only:rec_stream,rec_queue,nSMP
      implicit none
      integer istat,ied,jstat,jed,kstat,ked
      integer i,rankloc
      integer,save::gS
      logical,save::first_in=.true.
c
      if (use_pmecore) then
        rankloc  = rank_bis
      else
        rankloc  = rank
      end if
      kstat = kstart1(rankloc+1)
      ked   = kend1  (rankloc+1)
      jstat = jstart1(rankloc+1)
      jed   = jend1  (rankloc+1)
      istat = istart1(rankloc+1)
      ied   = iend1  (rankloc+1)

      if (first_in) then
         call cudaMaxgridsize("fphi_mpole_core",gS)
         first_in=.false.
      end if
      call set_PME_texture

!$acc host_data use_device(polerecglob,ipole,igrid,fphirec)
      call fphi_mpole_core<<<gS,PME_BLOCK_DIM,0,rec_stream>>>
     &     (kstat,ked,jstat,jed,istat,ied,nfft1,nfft2,nfft3
     &     ,npolerecloc,nlocrec,bsorder,nrec_send,nproc
     &     ,polerecglob,ipole,igrid,fphirec)
      call check_launch_kernel("fphi_mpole_core")
!$acc end host_data

      end subroutine

      subroutine grid_pchg_forcecu
      use atoms  ,only: n
      use atmlst ,only: chgrecglob
      use boxes  ,only: recip
      use charge
      use chgpot
      use domdec ,only: rank,rank_bis,nproc,nrec,ndir,comm_rec
     &           ,COMM_TINKER,nrec_send,locrec1,prec_send,nbloc
      use deriv  ,only: decrec
      use fft
      use inform ,only: deb_Path
      use pme    ,only: igrid,nfft1,nfft2,nfft3,bsorder
     &           ,thetai1_p,thetai2_p,thetai3_p,qgridin_2d
      use pmestuffcu,only: grid_pchg_force_core,grid_pchg_force_core1
      use potent ,only: use_pmecore
      use utilcu ,only: PME_BLOCK_DIM,check_launch_kernel
      use utilgpu,only: ngangs_rec,rec_queue,rec_stream
      implicit none
      integer iichg,iglob,iloc,iatm,isite,rankloc
     &       ,igrd0,jgrd0,kgrd0,i,j,k,i0,j0,k0,it1,it2,it3
     &       ,istat,ied,jstat,jed,kstat,ked
     &       ,iproc,proc,commloc,nprocloc
      integer,save:: gS,gS1
      real(t_p) f,fi,dn1,dn2,dn3,de1,de2,de3,t1,t2,t3
     &       ,dt1,dt2,dt3,term
      logical,save:: f_in=.true.

      if (use_pmecore) then
        nprocloc = nrec
        commloc  = comm_rec
        rankloc = rank_bis
      else
        nprocloc = nproc
        commloc  = COMM_TINKER
        rankloc  = rank
      end if
      if (deb_Path) write(*,'(3x,a)') 'grid_pchg_forcecu'

      istat = istart1(rankloc+1)
      jstat = jstart1(rankloc+1)
      kstat = kstart1(rankloc+1)
      ied   =   iend1(rankloc+1)
      jed   =   jend1(rankloc+1)
      ked   =   kend1(rankloc+1)
      f     = electric / dielec
      dn1   = real(nfft1,t_p)
      dn2   = real(nfft2,t_p)
      dn3   = real(nfft3,t_p)

      if (f_in) then
         call cudaMaxgridsize("grid_pchg_force_core",gS)
         call cudaMaxgridsize("grid_pchg_force_core1",gS1)
         !gS = gS - 2*nSMP
         f_in=.false.
      end if
      call set_pme_texture

!$acc host_data use_device(chgrecglob,iion,locrec1,igrid,pchg
!$acc&         ,thetai1_p,thetai2_p,thetai3_p,decrec)
      if (nproc.eq.1.or.nrec.eq.1) then
         call grid_pchg_force_core1<<<gS1,PME_BLOCK_DIM,0,rec_stream>>>
     &        (chgrecglob,iion,locrec1,igrid,pchg
     &        ,thetai1_p,thetai2_p,thetai3_p
     &        ,decrec
     &        ,kstat,ked,jstat,jed,istat,ied
     &        ,nrec_send,nionrecloc,n,nfft1,nfft2,nfft3
     &        ,f,dn1,dn2,dn3)
         call check_launch_kernel(" grid_pchg_force_core1")
      else
         call grid_pchg_force_core<<<gS,PME_BLOCK_DIM,0,rec_stream>>>
     &        (chgrecglob,iion,locrec1,igrid,pchg
     &        ,thetai1_p,thetai2_p,thetai3_p
     &        ,decrec
     &        ,kstat,ked,jstat,jed,istat,ied
     &        ,nrec_send,nionrecloc,n,nfft1,nfft2,nfft3
     &        ,f,dn1,dn2,dn3)
         call check_launch_kernel(" grid_pchg_force_core")
      end if
!$acc end host_data

      end subroutine

      subroutine fphi_uind_sitecu2(fdip_phi1,fdip_phi2)
      use atoms
      use atmlst
      use boxes
      use domdec
      use fft
      use mpole     ,only:npolerecloc,ipole
      use pme
      use pmestuffcu,only:fphi_uind_sitecu2_core
      use potent
      use utilcu    ,only:PME_BLOCK_DIM,check_launch_kernel
      use utilgpu   ,only:rec_stream,rec_queue,nSMP
      implicit none
      integer istat,ied,jstat,jed,kstat,ked
      integer i
      integer iproc,proc
      integer ind,psize,ierr
      real(t_p),dimension(:,:)::fdip_phi1,fdip_phi2
      integer,pointer :: ipole_p(:)
      integer,save::gS
      logical,save::first_in=.true.

      if (use_pmecore) then
        kstat = kstart1(rank_bis+1)
        ked   = kend1  (rank_bis+1)
        jstat = jstart1(rank_bis+1)
        jed   = jend1  (rank_bis+1)
        istat = istart1(rank_bis+1)
        ied   = iend1  (rank_bis+1)
      else
        kstat = kstart1(rank+1)
        ked   = kend1  (rank+1)
        jstat = jstart1(rank+1)
        jed   = jend1  (rank+1)
        istat = istart1(rank+1)
        ied   = iend1  (rank+1)
      end if

      call set_pme_texture
      psize =  size(kstart1)
      if (first_in) then
      if (nproc.eq.1) then
         call cudamaxgridSize("fphi_uind_sitecu2_core_1p",gS)
         gS = gS - 2*nSMP
      else
         call cudamaxgridSize("fphi_uind_sitecu2",gS)
         gS = gS - 2*nSMP
      end if
         first_in=.false.
      end if

!$acc host_data use_device(qgrid2in_2d,polerecglob,ipole,igrid
!$acc&    ,x,y,z,recip,thetai1,thetai2,thetai3,fdip_phi1,fdip_phi2
!$acc&    ,kstart1,kend1,jstart1,jend1,istart1,iend1,prec_send)

      if (nproc.eq.1.or.nrec.eq.1) then
      call fphi_uind_sitecu2_core_1p<<<gS,PME_BLOCK_DIM,0,rec_stream>>>
     &     (kstat,ked,jstat,jed,istat,ied
     &     ,npolerecloc,nlocrec,bsorder,n1mpimax,n2mpimax,n3mpimax
     &     ,nfft1,nfft2,nfft3
     &     ,qgrid2in_2d,polerecglob,ipole,igrid,x,y,z,recip
     &     ,thetai1,thetai2,thetai3
     &     ,fdip_phi1,fdip_phi2,first_in)
      call check_launch_kernel("fphi_uind_sitecu2_core_1p")
      else
      call fphi_uind_sitecu2_core<<<gS,PME_BLOCK_DIM,0,rec_stream>>>
     &     (kstat,ked,jstat,jed,istat,ied
     &     ,npolerecloc,nlocrec,bsorder,n1mpimax,n2mpimax,n3mpimax
     &     ,nrec_send,nproc,prec_send,psize,nfft1,nfft2,nfft3
     &     ,kstart1,kend1,jstart1,jend1,istart1,iend1
     &     ,qgrid2in_2d,polerecglob,ipole,igrid,x,y,z,recip
     &     ,thetai1,thetai2,thetai3
     &     ,fdip_phi1,fdip_phi2,first_in)
      call check_launch_kernel("fphi_uind_sitecu2_core")
      end if

!$acc end host_data

      end subroutine

      subroutine fphi_uind_sitecu1(fdip_phi1,fdip_phi2,fdip_sum_phi)
      use atoms
      use atmlst
      use boxes
      use domdec
      use fft
      use mpole     ,only:npolerecloc,ipole
      use pme
      use pmestuffcu,only:fphi_uind_sitecu1_core
      use potent
      use utilcu    ,only:PME_BLOCK_DIM,check_launch_kernel
      use utilgpu   ,only:rec_stream,rec_queue,nSMP
      implicit none
      integer istat,ied,jstat,jed,kstat,ked
      integer i
      integer iproc,proc
      integer isite,iatm,iipole
      real(t_p),dimension(:,:)::fdip_phi1,fdip_phi2
      real(t_p),dimension(:,:)::fdip_sum_phi
      integer,save::gS
      logical,save::first_in=.true.

      gS = 4*nSMP
      if (use_pmecore) then
        kstat = kstart1(rank_bis+1)
        ked   = kend1  (rank_bis+1)
        jstat = jstart1(rank_bis+1)
        jed   = jend1  (rank_bis+1)
        istat = istart1(rank_bis+1)
        ied   = iend1  (rank_bis+1)
      else
        kstat = kstart1(rank+1)
        ked   = kend1  (rank+1)
        jstat = jstart1(rank+1)
        jed   = jend1  (rank+1)
        istat = istart1(rank+1)
        ied   = iend1  (rank+1)
      end if

      if (first_in) then
         call cudamaxgridSize("fphi_uind_sitecu1",gS)
         first_in=.false.
      end if
      call set_pme_texture

!$acc host_data use_device(qgrid2in_2d,polerecglob,ipole,igrid
!$acc&   ,fdip_phi1,fdip_phi2,fdip_sum_phi
!$acc&   ,kstart1,kend1,jstart1,jend1,istart1,iend1,prec_send)

      call fphi_uind_sitecu1_core<<<gS,PME_BLOCK_DIM,0,rec_stream>>>
     &     (kstat,ked,jstat,jed,istat,ied,nfft1,nfft2,nfft3
     &     ,npolerecloc,nlocrec,bsorder,nrec_send,nproc
     &     ,polerecglob,ipole,igrid,fdip_phi1,fdip_phi2,fdip_sum_phi)
      call check_launch_kernel("fphi_uind_sitecu1_core")

!$acc end host_data

      end subroutine

      subroutine cudaMaxGridSize(kernelname,gS)
      use cudafor
      use echargecu  ,only: ecreal1d_core_cu
      use empole1cu  ,only: emreal1c_core_cu
     &               ,emreal3_cu,emrealshortlong3_cu
     &               ,emrealshortlong1c_core_cu
      use epolar1cu  ,only: epreal1c_core_cu
     &               ,epreal3_cu,mpreal1c_core_cu
      use pmestuffcu
      use tmatxb_pmecu
      use utilcu
      use utilgpu
      implicit none
      character(*),intent(in)::kernelname
      integer,intent(out)::gS
      integer sbytes,ierr

 200  format("error",I5," returned when calling Cuda0ccupancy"//
     &       "MaxActiveBlocksPerMulti in cudamaxGridSize subroutine",
     &       /,A100)

 201  format("cudaMaxGridSize argument ( ", A,
     &       " ) does not match any known Cuda Kernel", /,
     &       " Returning 0 as gridSize ")

      if      (kernelname.eq."fphi_uind_sitecu2") then
         ierr = CUDAOCCUPANCYMAXACTIVEBLOCKSPERMULTIPROCESSOR
     &          (gS,fphi_uind_sitecu2_core,PME_BLOCK_DIM,0)
         gS   = gS*nSMP
         if (ierr.ne.cudaSuccess)
     &      print 200, ierr,cudageterrorstring(ierr)
      else if (kernelname.eq."fphi_uind_sitecu2_core_1p") then
         ierr = CUDAOCCUPANCYMAXACTIVEBLOCKSPERMULTIPROCESSOR
     &          (gS,fphi_uind_sitecu2_core_1p,PME_BLOCK_DIM,0)
         if (ierr.ne.cudaSuccess)
     &      print 200, ierr,cudageterrorstring(ierr)
         gS   = gS*nSMP
      else if (kernelname.eq."fphi_uind_sitecu1") then
         ierr = CUDAOCCUPANCYMAXACTIVEBLOCKSPERMULTIPROCESSOR
     &          (gS,fphi_uind_sitecu1_core,PME_BLOCK_DIM,0)
         if (ierr.ne.cudaSuccess)
     &      print 200, ierr,cudageterrorstring(ierr)
         gS   = gS*nSMP
      else if (kernelname.eq."fphi_mpole_core") then
         ierr = CUDAOCCUPANCYMAXACTIVEBLOCKSPERMULTIPROCESSOR
     &          (gS,fphi_mpole_core,PME_BLOCK_DIM,0)
         if (ierr.ne.cudaSuccess)
     &      print 200, ierr,cudageterrorstring(ierr)
         gS   = gS*nSMP
      else if (kernelname.eq."grid_mpole_core") then
         ierr = CUDAOCCUPANCYMAXACTIVEBLOCKSPERMULTIPROCESSOR
     &          (gS,grid_mpole_sitecu_core_1p,PME_BLOCK_DIM1,0)
         if (ierr.ne.cudaSuccess)
     &      print 200, ierr,cudageterrorstring(ierr)
         gS   = gS*nSMP
      else if (kernelname.eq."grid_pchg_core") then
         ierr = CUDAOCCUPANCYMAXACTIVEBLOCKSPERMULTIPROCESSOR
     &          (gS,grid_pchg_sitecu_core_1p,PME_BLOCK_DIM1,0)
         if (ierr.ne.cudaSuccess)
     &      print 200, ierr,cudageterrorstring(ierr)
         gS   = gS*nSMP
      else if (kernelname.eq."grid_pchg_force_core") then
         ierr = CUDAOCCUPANCYMAXACTIVEBLOCKSPERMULTIPROCESSOR
     &          (gS,grid_pchg_force_core,PME_BLOCK_DIM1,0)
         if (ierr.ne.cudaSuccess)
     &      print 200, ierr,cudageterrorstring(ierr)
         gS   = gS*nSMP
      else if (kernelname.eq."grid_pchg_force_core1") then
         ierr = CUDAOCCUPANCYMAXACTIVEBLOCKSPERMULTIPROCESSOR
     &          (gS,grid_pchg_force_core1,PME_BLOCK_DIM1,0)
         if (ierr.ne.cudaSuccess)
     &      print 200, ierr,cudageterrorstring(ierr)
         gS   = gS*nSMP
      else if (kernelname.eq."grid_uind_core") then
         ierr = CUDAOCCUPANCYMAXACTIVEBLOCKSPERMULTIPROCESSOR
     &          (gS,grid_uind_sitecu_core_1p,PME_BLOCK_DIM1,0)
         if (ierr.ne.cudaSuccess)
     &      print 200, ierr,cudageterrorstring(ierr)
         gS   = gS*nSMP
      else if (kernelname.eq."ecreal1d_core_cu") then
         ierr = CUDAOCCUPANCYMAXACTIVEBLOCKSPERMULTIPROCESSOR
     &          (gS,ecreal1d_core_cu,BLOCK_DIM,0)
         if (ierr.ne.cudaSuccess)
     &      print 200, ierr,cudageterrorstring(ierr)
         gS   = gS*nSMP
      else if (kernelname.eq."ecreal3d_core_cu") then
         ierr = CUDAOCCUPANCYMAXACTIVEBLOCKSPERMULTIPROCESSOR
     &          (gS,ecreal1d_core_cu,BLOCK_DIM,0)
         if (ierr.ne.cudaSuccess)
     &      print 200, ierr,cudageterrorstring(ierr)
         gS   = gS*nSMP
      else if (kernelname.eq."emreal1c_core_cu") then
         ierr = CUDAOCCUPANCYMAXACTIVEBLOCKSPERMULTIPROCESSOR
     &          (gS,emreal1c_core_cu,BLOCK_DIM,0)
         if (ierr.ne.cudaSuccess)
     &      print 200, ierr,cudageterrorstring(ierr)
         gS   = gS*nSMP
      else if (kernelname.eq."emrealshortlong1c_core_cu") then
         ierr = CUDAOCCUPANCYMAXACTIVEBLOCKSPERMULTIPROCESSOR
     &          (gS,emrealshortlong1c_core_cu,BLOCK_DIM,0)
         if (ierr.ne.cudaSuccess)
     &      print 200, ierr,cudageterrorstring(ierr)
         gS   = gS*nSMP
      else if (kernelname.eq."emreal3_cu") then
         ierr = CUDAOCCUPANCYMAXACTIVEBLOCKSPERMULTIPROCESSOR
     &          (gS,emreal3_cu,BLOCK_DIM,0)
         if (ierr.ne.cudaSuccess)
     &      print 200, ierr,cudageterrorstring(ierr)
         gS   = gS*nSMP
      else if (kernelname.eq."emrealshortlong3_cu") then
         ierr = CUDAOCCUPANCYMAXACTIVEBLOCKSPERMULTIPROCESSOR
     &          (gS,emrealshortlong3_cu,BLOCK_DIM,0)
         if (ierr.ne.cudaSuccess)
     &      print 200, ierr,cudageterrorstring(ierr)
         gS   = gS*nSMP
      else if (kernelname.eq."epreal1c_core_cu") then
         ierr = CUDAOCCUPANCYMAXACTIVEBLOCKSPERMULTIPROCESSOR
     &          (gS,epreal1c_core_cu,BLOCK_DIM,0)
         if (ierr.ne.cudaSuccess)
     &      print 200, ierr,cudageterrorstring(ierr)
         gS   = gS*nSMP
      else if (kernelname.eq."mpreal1c_core_cu") then
         ierr = CUDAOCCUPANCYMAXACTIVEBLOCKSPERMULTIPROCESSOR
     &          (gS,mpreal1c_core_cu,BLOCK_DIM,0)
         if (ierr.ne.cudaSuccess)
     &      print 200, ierr,cudageterrorstring(ierr)
         gS   = gS*nSMP
      else if (kernelname.eq."epreal3_cu") then
         ierr = CUDAOCCUPANCYMAXACTIVEBLOCKSPERMULTIPROCESSOR
     &          (gS,epreal3_cu,BLOCK_DIM,0)
         if (ierr.ne.cudaSuccess)
     &      print 200, ierr,cudageterrorstring(ierr)
         gS   = gS*nSMP
      else if (kernelname.eq."tmatxb_pme_core_cu") then
         ierr = CUDAOCCUPANCYMAXACTIVEBLOCKSPERMULTIPROCESSOR
     &          (gS,tmatxb_pme_core_cu,BLOCK_DIM,0)
         if (ierr.ne.cudaSuccess)
     &      print 200, ierr,cudageterrorstring(ierr)
         gS   = gS*nSMP
      else if (kernelname.eq."otf_dc_tmatxb_pme_core_cu") then
         ierr = CUDAOCCUPANCYMAXACTIVEBLOCKSPERMULTIPROCESSOR
     &          (gS,otfdc_tmatxb_pme_core_cu,BLOCK_DIM,0)
         if (ierr.ne.cudaSuccess)
     &      print 200, ierr,cudageterrorstring(ierr)
         gS   = gS*nSMP
      else
         print 201, kernelname
         gS = 0
      end if

      end

      subroutine attach_pmecu_pointer(config)
      use atoms    ,only: x,y,z
      use atmlst   ,only: polerecglob
      use boxes    ,only: recip
      use domdec
      use fft      ,only: kstart1,kend1,istart1,iend1,jstart1,jend1
      use pme      ,only: igrid,qgrid2in_2d,qgridin_2d
     &             ,thetai1,thetai2,thetai3
      use pmestuffcu
      implicit none
      integer,intent(in)::config
      enum,bind(C)
      enumerator pmeD,thetaD
      end enum

      select case (config)
         case (pmeD)
!$acc host_data use_device(thetai1,thetai2,thetai3,igrid
!$acc&  ,x,y,z,qgrid2in_2d,qgridin_2d,kstart1,kend1,jstart1,jend1
!$acc&  ,istart1,iend1,prec_send)
         x_t         => x
         y_t         => y
         z_t         => z
         qgridin_t   => qgridin_2d
         qgrid2in_t  => qgrid2in_2d
         kstart1_t   => kstart1
         kend1_t     => kend1
         istart1_t   => istart1
         iend1_t     => iend1
         jstart1_t   => jstart1
         jend1_t     => jend1
         prec_send_t => prec_send
!$acc end host_data
         case (thetaD)
!$acc host_data use_device(thetai1,thetai2,thetai3)
         thetai1_t => thetai1
         thetai2_t => thetai2
         thetai3_t => thetai3
!$acc end host_data
         case default
         print*, 'Unknown configuration for thetai'
      end select
      end subroutine

      subroutine set_PME_texture
      use atoms    ,only: x,y,z
      use atmlst   ,only: polerecglob
      use boxes    ,only: recip
      use domdec
      use fft      ,only: kstart1,kend1,istart1,iend1,jstart1,jend1
      use pme      ,only: igrid,qgrid2in_2d,qgridin_2d
     &             ,thetai1,thetai2,thetai3
      use pmestuffcu
      implicit none
      logical,save :: first_in=.true.


      if (first_in) then
!$acc host_data use_device(thetai1,thetai2,thetai3,igrid
!$acc&  ,x,y,z,qgrid2in_2d,qgridin_2d,kstart1,kend1,jstart1,jend1
!$acc&  ,istart1,iend1,prec_send)
         x_t         => x
         y_t         => y
         z_t         => z
         qgridin_t   => qgridin_2d
         qgrid2in_t  => qgrid2in_2d
         kstart1_t   => kstart1
         kend1_t     => kend1
         istart1_t   => istart1
         iend1_t     => iend1
         jstart1_t   => jstart1
         jend1_t     => jend1
         prec_send_t => prec_send
!$acc end host_data
      end if

      if (nproc>1.or.first_in) then
!$acc host_data use_device(thetai1,thetai2,thetai3)
         thetai1_t => thetai1
         thetai2_t => thetai2
         thetai3_t => thetai3
!$acc end host_data
         first_in = .false.
      end if

      end subroutine
#endif

      subroutine set_pme_eps
#if (defined(MIXED)||defined(SINGLE))
      use pme
#ifdef _CUDA
      use pmestuffcu,only:setcu_pme_eps
#endif
      use tinheader ,only:prec1_eps
      implicit none
      pme_eps = 2*max(max(nfft1,nfft2),nfft3)*prec1_eps
# ifdef _CUDA
      call setcu_pme_eps(pme_eps)
# endif
#endif
      end
c
c     "cmp_to_fmp_site" transforms the ith atomic multipole from Cartesian
c     to fractional coordinates
c
c
      subroutine cmp_to_fmp_sitegpu(cmp,fmp)
      use mpole
      use sizes
      use utilgpu   ,only:rec_queue,ctf
      use tinheader ,only:ti_p
      implicit none
      integer i,j,k
      real(t_p),intent(in) ::cmp(10,max(npolerecloc,1))
      real(t_p),intent(out)::fmp(10,max(npolerecloc,1))
      logical,save:: init=.true.
c     real(t_p) time0,time1
c
c     find the matrix to convert Cartesian to fractional
c
      if (init) then
         call cart_to_fracgpu
         init=.false.
      end if

!$acc parallel loop collapse(2) async(rec_queue)
!$acc&         present(cmp,fmp,ctf)
      do i= 1,npolerecloc
         do j = 1, 10
c
c        apply the transformation to get the fractional multipoles
c
            if (j.eq.1) then
               fmp(1,i) = ctf(1,1) * cmp(1,i)
            else if (j.lt.5) then
               fmp(j,i) = 0.0_ti_p
!$acc loop seq
               do k = 2, 4
                  fmp(j,i) = fmp(j,i) + ctf(j,k)*cmp(k,i)
               end do
            else
               fmp(j,i) = 0.0_ti_p
!$acc loop seq
               do k = 5, 10
                  fmp(j,i) = fmp(j,i) + ctf(j,k)*cmp(k,i)
               end do
            end if
         end do
      end do
      end
c
c
      subroutine fphi_to_cphi_sitegpu(fphi,cphi)
      use mpole
      use sizes
      use utilgpu   ,only:rec_queue,ftc
      use tinheader ,only:ti_p
      implicit none
      integer i,j,k
      real(t_p),intent(in )::fphi(20,max(npolerecloc,1))
      real(t_p),intent(out)::cphi(10,max(npolerecloc,1))
      logical,save::init=.true.
c
c     find the matrix to convert fractional to Cartesian
c
      if (init) then
         call frac_to_cartgpu
         init=.false.
      end if

!$acc parallel loop collapse(2) async(rec_queue)
!$acc&         present(fphi,cphi,ftc)
      do i= 1,npolerecloc
         do j = 1, 10
c
c        apply the transformation to get the Cartesian potential
c
            if (j.eq.1) then
               cphi(1,i) = ftc(1,1) * fphi(1,i)
            else if (j.lt.5) then
               cphi(j,i) = 0.0_ti_p
!$acc loop seq
               do k = 2, 4
                  cphi(j,i) = cphi(j,i) + ftc(j,k)*fphi(k,i)
               end do
            else
               cphi(j,i) = 0.0_ti_p
!$acc loop seq
               do k = 5, 10
                  cphi(j,i) = cphi(j,i) + ftc(j,k)*fphi(k,i)
               end do
            end if
         end do
      end do
      end
