c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ## subroutine elj3gpu -- Lennard-Jones vdw energy & analysis  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "elj3" calculates the Lennard-Jones 6-12 van der Waals energy
c     and also partitions the energy among the atoms
c
c
#include "tinker_precision.h"
      module elj3gpu_inl
        contains
#include "image.f.inc"
#include "switch_respa.f.inc"
#include "pair_elj.f.inc"
      end module

      subroutine elj3gpu
      use analyz
      use atoms
      use domdec
      use energi
      use inform
      use interfaces,only:elj3c_p
      use iounit
      use tinheader ,only:ti_p,re_p
      use timestat
      use vdwpot
      use mpi
      implicit none
      integer i
      real(t_p) elrc,aelrc
      real(t_p) time0,time1
c
      call timer_enter( timer_elj3 )
      call elj3c_p
      call timer_exit( timer_elj3,quiet_timers )
c
c     apply long range van der Waals correction if desired
c
      if (use_vcorr) then
         call evcorr (elrc)
         ev = ev + elrc
         aelrc = elrc / real(n,t_p)
         do i = 1, nbloc
            aev(i) = aev(i) + aelrc
         end do
         if (verbose .and. elrc.ne.0.0_ti_p) then
            write (iout,10)  elrc
   10       format (/,' Long Range vdw Correction :',9x,f12.4)
         end if
      end if
      end

c
c     #############################################################
c     ##                                                         ##
c     ## subroutine elj3cgpu -- Lennard-Jones analysis via list  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "elj3cgpu" calculates the Lennard-Jones van der Waals energy
c     and also partitions the energy among the atoms using a
c     pairwise neighbor list
c
c
      subroutine elj3cgpu
      use action
      use analyz
      use atmlst
      use atmtyp
      use atoms
      use bound
      use couple
      use domdec
      use energi
      use elj3gpu_inl
      use group
      use inter
      use inform
      use iounit
      use interfaces ,only: elj3_scaling
      use molcul
      use neigh
      use shunt
      use usage
      use utilgpu
      use vdw
      use vdwpot
      use vdw_locArray

      implicit none
      integer i,j,k,iglob,kbis,iivdw,kk,kv,kvloc
      integer ii,iv,it,ivloc,kglob
      integer nn12,nn13,nn14,ntot
      integer nevt
      integer nnvlst,nnvlst1,nvloop8,nvloop16
      integer kt,countsel

      real(t_p) e,etemp,p6
      real(t_p) xr,yr,zr,xi,yi,zi
      real(t_p) half,three
      real(t_p) time0,time1
      real(t_p) taper
      real(t_p) rik,rik2_1,rik3,rik4,rik5,rik2
      real(t_p) rdn,rdn1,redk,redi,rv,eps
      logical   usei
      parameter( half=0.5_ti_p,three=3.0_ti_p )

      character*10 mode
c
c     zero out the van der Waals energy and partitioning terms
c
      if(deb_Path) print*, 'elj3cgpu'

      call prmem_request(xred,nbloc,queue=dir_queue)
      call prmem_request(yred,nbloc,queue=dir_queue)
      call prmem_request(zred,nbloc,queue=dir_queue)
c
c     set the coefficients for the switching function
c     update the number of gangs required for gpu
c
      mode   = 'VDW'
      call switch (mode)
      call update_gang(nvdwbloc)
c     nev_   = 0
c     ev     = 0.0_re_p

!$acc data present(vdwglob,vdwlocnl,vdwglobnl,ired,kred,loc,
!$acc&   ivdw,jvdw,vlst,nvlst,epsilon,epsilon4,radmin,radmin4)
!$acc&     present(xred,yred,zred)
!$acc&     present(ev,nev_)

c
c     apply any reduction factor to the atomic coordinates
c
!$acc parallel loop async(dir_queue)
      do k = 1 , nvdwbloc
         iglob   = ivdw (vdwglob (k))
         i       = loc  (iglob)
         iv      = ired (iglob)
         rdn     = kred (iglob)
         xred(i) =  rdn * x(iglob)  + (1.0_ti_p-rdn) * x(iv)
         yred(i) =  rdn * y(iglob)  + (1.0_ti_p-rdn) * y(iv)
         zred(i) =  rdn * z(iglob)  + (1.0_ti_p-rdn) * z(iv)
      enddo
c
c     find the van der Waals energy via neighbor list search
c
!$acc parallel loop gang vector_length(32) async(dir_queue)
      do ii = 1, nvdwlocnl
         iivdw = vdwglobnl(ii)
         iglob = ivdw(iivdw)
         i     = loc(iglob)
         iv    = ired(iglob)
         ivloc = loc(iv)
         redi  = merge (kred(iglob),1.0_ti_p,(i.ne.ivloc))
         it    = jvdw(iglob)
         xi    = xred(i)
         yi    = yred(i)
         zi    = zred(i)
c        usei  = (use(iglob) .or. use(iv))
c
c     decide whether to compute the current interaction
c
!$acc loop vector
         do kk = 1, nvlst(ii)
            kglob = vlst(kk,ii)
            kbis  = loc(kglob)
            kv    = ired(kglob)
            kvloc = loc(kv)
            kt    = jvdw(kglob)
            xr    = xi - xred(kbis)
            yr    = yi - yred(kbis)
            zr    = zi - zred(kbis)
            if (use_bounds) call image_inl (xr,yr,zr)
            rik2  = xr*xr + yr*yr + zr*zr
c
c     check for an interaction distance less than the cutoff
c
            if (rik2.le.off2) then
               redk = merge (kred(kglob),1.0_ti_p,(kbis.ne.kvloc))
               rv   = radmin (kt,it)
               eps  = epsilon(kt,it)

               !compute the energy contribution for this interaction
               call elj3_couple(rik2,xr,yr,zr,rv,eps,cut2
     &                         ,cut,off,e)
c
c     increment the total van der Waals energy and derivatives
c
               ev   = ev   + e
               nev_ = nev_ + 1.0_ti_p
c!$acc atomic update
c               aev(i) = aev(i) + 0.5_ti_p * e

            end if
         end do
      end do

      call elj3_scaling(xred,yred,zred)

!$acc serial async present(nev,nev_)
      nev = int(nev_)
!$acc end serial

!$acc end data

      end

#ifdef _CUDA
      subroutine elj3c_cu
      use action    ,only: nev,nev_
      use atmlst    ,only: vdwglobnl,vdwglob
      use atoms     ,only: x,y,z,n
      use domdec    ,only: loc,rank,nbloc,nproc
     &              ,xbegproc,xendproc,ybegproc,yendproc,zbegproc
     &              ,zendproc,glob
      use ehal1cu   ,only: set_vdw_texture
      use eljcu
      use energi    ,only: ev
      use inform    ,only: deb_Path
      use interfaces,only: elj3_scaling
      use neigh     ,only: cellv_glob,cellv_loc,cellv_jvdw
     &              ,vblst,ivblst
      use tinheader ,only: ti_p
      use shunt     ,only: c0,c1,c2,c3,c4,c5,off2,off,cut2,cut
      use utilcu    ,only: check_launch_kernel
     &              ,BLOCK_DIM
      use utilgpu   ,only: def_queue,dir_queue,rec_queue,dir_stream
     &              ,rec_stream,rec_event,stream_wait_async
     &              ,warp_size,def_stream,inf
     &              ,ered_buff,nred_buff,reduce_energy_action
     &              ,zero_en_red_buffer,prmem_request
      use vdw       ,only: ired,kred,jvdw,ivdw,radmin_c
     &              ,epsilon_c,nvdwbloc,nvdwlocnl
     &              ,nvdwlocnlb,nvdwclass
     &              ,nvdwlocnlb_pair,nvdwlocnlb2_pair
      use vdwpot    ,only: dhal,ghal
      use vdw_locArray
      implicit none
      integer i,k
      integer iglob,iivdw,iv,grid
      integer ierrSync,lst_start
#ifdef TINKER_DEBUG
#endif
      real(t_p)  xbeg,xend,ybeg,yend,zbeg,zend
      real(t_p)  rdn,rdn1

      logical,save:: first_in=.true.
      logical,parameter:: dyn_gS=.true.
      integer,save:: gS
      character*10 mode

      call prmem_request(xred    ,nvdwlocnlb,queue=def_queue)
      call prmem_request(yred    ,nvdwlocnlb,queue=def_queue)
      call prmem_request(zred    ,nvdwlocnlb,queue=def_queue)
      call prmem_request(xredc   ,nbloc     ,queue=def_queue)
      call prmem_request(yredc   ,nbloc     ,queue=def_queue)
      call prmem_request(zredc   ,nbloc     ,queue=def_queue)
      call prmem_request(loc_ired,nvdwlocnlb,queue=def_queue)
      call prmem_request(loc_kred,nvdwlocnlb,queue=def_queue)

c
      if(deb_Path) write (*,*) 'elj3c_cu'
      def_queue = dir_queue
      def_stream = dir_stream
      xbeg = xbegproc(rank+1)
      xend = xendproc(rank+1)
      ybeg = ybegproc(rank+1)
      yend = yendproc(rank+1)
      zbeg = zbegproc(rank+1)
      zend = zendproc(rank+1)
      lst_start = 2*nvdwlocnlb_pair+1

      if (first_in) then
         first_in = .false.
      end if

      if(dyn_gS) gS = nvdwlocnlb2_pair/8

#ifdef _OPENACC
      if (dir_queue.ne.rec_queue)
     &   call stream_wait_async(rec_stream,dir_stream,rec_event)
#endif

#ifdef TINKER_DEBUG
#endif

c
c     apply any reduction factor to the atomic coordinates
c
!$acc parallel loop default(present) async(def_queue)
      do k = 1,nvdwlocnlb
         if (k.le.nvdwlocnl) then
            iglob    = cellv_glob(k)
            iv       = ired (iglob)
            rdn      = kred (iglob)
            rdn1     = 1.0_ti_p - rdn
            cellv_loc(k) = loc(iglob)
            loc_ired(k)  = loc(iv)
            if (iglob.eq.iv) then
               loc_kred(k) = rdn
            else
               loc_kred(k) = 1.0_ti_p
            end if
            xred(k)  = rdn * x(iglob) + rdn1 * x(iv)
            yred(k)  = rdn * y(iglob) + rdn1 * y(iv)
            zred(k)  = rdn * z(iglob) + rdn1 * z(iv)
         else
            ! Exclusion buffer to prevent interaction compute
            cellv_loc(k) = nbloc
            loc_ired(k)  = nbloc
            xred(k) = inf
            yred(k) = inf
            zred(k) = inf
         end if
      end do

!$acc parallel loop default(present) async(def_queue)
      do k = 1,nvdwbloc
         iglob    = ivdw(vdwglob(k))
         i        = loc  (iglob)
         iv       = ired (iglob)
         rdn      = kred (iglob)
         rdn1     = 1.0_ti_p - rdn
         xredc(i)  = rdn * x(iglob) + rdn1 * x(iv)
         yredc(i)  = rdn * y(iglob) + rdn1 * y(iv)
         zredc(i)  = rdn * z(iglob) + rdn1 * z(iv)
      end do

      call zero_en_red_buffer(def_queue)
c
c     set the coefficients for the switching function
c
      !print*, nvdwlocnlb_pair
      mode = 'VDW'
      call switch (mode)

c
c     Call Vdw kernel in CUDA using C2 nblist
c
!$acc host_data use_device(xred,yred,zred,cellv_glob,cellv_loc
!$acc&    ,loc_ired,ivblst,vblst,cellv_jvdw,epsilon_c
!$acc&    ,radmin_c,ired,kred,ered_buff,nred_buff
#ifdef TINKER_DEBUG
#endif
!$acc&    )

      call set_vdw_texture
     &     (kred,radmin_c,epsilon_c,xred,yred,zred
     &     ,vblst,ivblst,loc_ired,cellv_jvdw,cellv_glob,cellv_loc
     &     ,n,nvdwlocnl,nvdwlocnlb,nvdwclass,nvdwlocnlb_pair)
      call elj3_cu<<<gS,BLOCK_DIM,0,def_stream>>>
     &             (xred,yred,zred,cellv_glob,cellv_loc,loc_ired
     &             ,ivblst,vblst(lst_start),cellv_jvdw
     &             ,epsilon_c,radmin_c,ired,kred
     &             ,ered_buff,nred_buff
     &             ,nvdwlocnlb2_pair,n,nbloc,nvdwlocnl,nvdwlocnlb
     &             ,nvdwclass
     &             ,c0,c1,c2,c3,c4,c5,cut2,cut,off2,off,ghal,dhal
     &             ,xbeg,xend,ybeg,yend,zbeg,zend
#ifdef TINKER_DEBUG
#endif
     &             )
      call check_launch_kernel(" elj3_cu ")

!$acc end host_data

      call reduce_energy_action(ev,nev,def_queue)

#ifdef TINKER_DEBUG
#endif

      call elj3_scaling(xredc,yredc,zredc)

!$acc serial async present(nev,nev_)
      nev = int(nev_) + nev
!$acc end serial

      end subroutine
#endif
c
c     Scaling interaction correction subroutines for Lennard-Jones
c
      subroutine elj3_scaling(xred,yred,zred)

      use action    ,only: nev_
      use atmlst    ,only: vdwglobnl
      use domdec    ,only: loc,rank
      use elj3gpu_inl
      use energi    ,only: ev
      use inform    ,only: deb_Path
      use tinheader ,only: ti_p
      use tintypes  ,only: real3
      use shunt     ,only: c0,c1,c2,c3,c4,c5,off2,off,cut2,cut
      use vdw       ,only: ired,kred,jvdw,ivdw,radmin,radmin4,
     &                     epsilon,epsilon4
      use vdwpot    ,only: vcorrect_ik,vcorrect_scale,n_vscale,dhal,ghal
      use utilgpu   ,only: dir_queue
      use virial
      implicit none
      integer i,j,k,kk,ksave
      integer kt,kglob,kbis,kvloc,kv,ki
      integer iglob,iivdw
      integer ii,iv,it,ivloc
      integer nnvlst,nnvlst2
      integer nn12,nn13,nn14,ntot
      integer interac
      real(t_p)  xi,yi,zi,redi,e,de
      real(t_p)  half,one
      real(t_p)  rdn,rdn1,redk
      real(t_p)  rik2
      type(real3) ded
      real(r_p)  devx,devy,devz
      real(t_p)  invrho,rv7orho
      real(t_p)  dtau,gtau,tau,tau7,rv7
      real(t_p)  rv2,eps2
      real(t_p)  xpos,ypos,zpos
      real(t_p)  vscale,vscale4
      logical    do_scale4
      character*10 mode

      real(t_p),intent(in):: xred(:)
      real(t_p),intent(in):: yred(:)
      real(t_p),intent(in):: zred(:)
      parameter(half=0.5_ti_p,
     &           one=1.0_ti_p)

      if (deb_Path) write(*,'(2x,a)') "elj3_scaling"

      ! Scaling factor correction loop
!$acc parallel loop async(dir_queue) gang vector
!$acc&         present(xred,yred,zred)
!$acc&         present(loc,ired,kred,ivdw,loc,jvdw,vir,radmin,
!$acc&  radmin4,epsilon,epsilon4,vcorrect_ik,vcorrect_scale)
!$acc&         present(ev,nev_)
      do ii = 1,n_vscale
         iglob  = vcorrect_ik(ii,1)
         kglob  = vcorrect_ik(ii,2)
         vscale = vcorrect_scale(ii)
         i      = loc(iglob)
         kbis   = loc(kglob)

         ivloc  = loc (ired(iglob))
         kvloc  = loc (ired(kglob))
         it     = jvdw(iglob)
         kt     = jvdw(kglob)

         redi   = merge (kred(iglob),1.0_ti_p,(i.ne.ivloc))
         redk   = merge (kred(kglob),1.0_ti_p,(kbis.ne.kvloc))

         do_scale4 = .false.
         vscale4   = 0

         if (vscale.lt.0) then 
            vscale4 = -vscale
            vscale = 1
         end if
c
c     compute the energy contribution for this interaction
c
         xpos   = xred(i) - xred(kbis)
         ypos   = yred(i) - yred(kbis)
         zpos   = zred(i) - zred(kbis)
         call image_inl(xpos,ypos,zpos)
c
c     decide whether to compute the current interaction
c     and check for an interaction distance less than the cutoff
c
         rik2   = xpos**2 + ypos**2 + zpos**2
         if (rik2>off2) cycle
c
c     replace 1-4 interactions
c
 20      continue
         if (do_scale4) then
            rv2  = radmin4 (kt,it)
            eps2 = epsilon4(kt,it)
         else
            rv2  =  radmin (kt,it)
            eps2 = epsilon (kt,it)
         end if

         call elj3_couple(rik2,xpos,ypos,zpos,rv2,eps2*vscale
     &                    ,cut2,cut,off,e)

         if (.not.do_scale4) then
         e    = -e
         end if

         ev   =   ev + e
         !if(rank.eq.0.and.mod(ii,1).eq.0) print*,iglob,kglob,vscale,e
         if (vscale.eq.1.0_ti_p) nev_=nev_-1
         if (vscale4.lt.0)       nev_=nev_+1

         ! deal with 1-4 Interactions
         if (vscale4.gt.0) then
            vscale    =  vscale4
            do_scale4 = .true.
            vscale4   = 0
            goto 20
         end if
      end do

      end subroutine
