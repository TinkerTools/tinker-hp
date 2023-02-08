c
c     sorbonne university
c     washington university in saint louis
c     university of texas at austin
c
c     ###################################################################
c     ##                                                               ##
c     ##  subroutine ehal1gpu  --  buffered 14-7 energy & derivatives  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     "ehal1" calculates the buffered 14-7 van der waals energy and
c     its first derivatives with respect to cartesian coordinates
c
c
#include "tinker_macro.h"
      module ehal1gpu_inl
        real(r_p) elrc,vlrc
!$acc declare create(elrc,vlrc)
#include "atomicOp.h.f"
        contains
#include "convert.f.inc"
#include "image.f.inc"
#include "atomicOp.inc.f"
#include "groups.inc.f"
#include "pair_ehal.inc.f"
      end module
#define __tver__ (__use_ene__+__use_grd__+__use_vir__)
#define __tfea__ (__use_mpi__+__use_softcore__+__use_groups__)

      subroutine ehal1gpu

      use domdec     ,only: ndir,rank
      use ehal1gpu_inl,only: elrc,vlrc
      use energi     ,only: ev
      use interfaces ,only: ehal1c_p
      use utilgpu    ,only: def_queue,dir_queue
      use potent
      use virial     ,only: vir
      use vdwpot     ,only: use_vcorr
      implicit none
      character*11 mode
      logical      fullrange
c
      ! PME-Core case
      if (use_pmecore.and.rank.ge.ndir) return
      def_queue = dir_queue
      fullrange = .not.(use_vdwshort.or.use_vdwlong)
c
c     choose the method for summing over pairwise interactions
c
      call ehal1c_p
c
c     apply long range van der Waals correction if desired
c
      if (use_vcorr.and.(fullrange.or.use_vdwlong)) then
         mode = "VDW"
         call evcorr1gpu (mode,elrc,vlrc)
!$acc serial present(ev,vir,elrc,vlrc) async(def_queue)
         ev = ev + elrc
         vir(1,1) = vir(1,1) + vlrc
         vir(2,2) = vir(2,2) + vlrc
         vir(3,3) = vir(3,3) + vlrc
!$acc end serial
      end if

      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine ehal1c  --  buffered 14-7 vdw derivs via list  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ehal1c" calculates the buffered 14-7 van der Waals energy and
c     its first derivatives with respect to Cartesian coordinates
c     using a pairwise neighbor list
c
c
      subroutine ehal1c_ac
      use atmlst    ,only: vdwglobnl,vdwglob
      use atoms     ,only: x,y,z,n
      use couple    ,only: i12,n12
      use cutoff    ,only: vdwshortcut,shortheal
      use deriv     ,only: dev=>de_ws2,delambdav
      use domdec    ,only: loc,rank,nbloc
      use ehal1gpu_inl
      use energi    ,only: ev=>ev_r
      use group     ,only: use_group,ngrp,grplist,wgrp
      use inform    ,only: deb_Path
      use mutant    ,only: scalpha,scexp,vlambda,vcouple,mut=>mutInt
     &              ,nmut
      use neigh     ,only: vlst,nvlst,shortvlst,nshortvlst
      use potent    ,only: use_vdwshort,use_vdwlong,use_lambdadyn
      use tinheader ,only: ti_p,re_p,one1,two1
      use tinMemory ,only: prmem_request
      use sizes     ,only: maxvalue,tinkerdebug
      use shunt     ,only: c0,c1,c2,c3,c4,c5,off2,off,cut2,cut
      use vdw       ,only: ired,kred,jvdw,ivdw,radmin,radmin4,epsilon
     &              ,epsilon4,nvdwbloc,nvdwlocnl,nvdwclass,skipvdw12
      use vdwpot    ,only: vcorrect_ik,vcorrect_scale,n_vscale,dhal,ghal
      use vdw_locArray
      use utilgpu   ,only: def_queue,dir_queue,rec_queue
#ifdef _OPENACC
     &                    ,dir_stream
     &                    ,rec_stream,rec_event,stream_wait_async
#endif
      use virial
      implicit none
      integer i,j,k,kk,ksave,ver,ver1,fea
      integer kt,kglob,kbis,kvloc,kv,ki
      integer iglob,iivdw
      integer ii,iv,it,ivloc
      integer nnvlst,nnvlst2
      integer nn12,nn13,nn14,ntot
      integer in12,ai12(maxvalue)
      integer,save:: ncall=0
      integer,pointer:: lst(:,:),nlst(:)
      real(t_p)  xi,yi,zi,xpos,ypos,zpos,redi,e,de,rinv
      real(t_p)  rdn,rdn1,redk,dedx,dedy,dedz
      real(t_p)  invrik,rik,rik2,rik3,rik4,rik5,rik6,rik7
      real(t_p)  invrho,rv7orho,dtau,gtau,tau,tau7,rv7
      real(t_p)  rv2,eps2,dtaper,taper,fgrp
      real(t_p)  vscale,vscale4
      real(t_p)  loff2,dt1lam,dt2lam,delambdav_
      logical    do_scale4,ik12
      integer(1) muti,mutik
      character*11 mode
      parameter( ver=__use_grd__+__use_ene__+__use_vir__
     &         , ver1=ver+__use_sca__)

c
      if(deb_Path) write (*,*) 'ehal1c_ac',use_vdwshort,use_vdwlong
      ncall = ncall + 1

      fea    = __use_mpi__
      if (use_group) fea = fea + __use_groups__
      if (nmut.ne.0) fea = fea + __use_softcore__
      if (use_lambdadyn) fea = fea + __use_lambdadyn__

#ifdef _OPENACC
      if (dir_queue.ne.rec_queue) then
!!$acc wait(rec_queue) async(dir_queue)
         call stream_wait_async(rec_stream,dir_stream,rec_event)
      end if
#endif
      call prmem_request(xred,nbloc,queue=def_queue)
      call prmem_request(yred,nbloc,queue=def_queue)
      call prmem_request(zred,nbloc,queue=def_queue)
c
c     apply any reduction factor to the atomic coordinates
c
!$acc parallel loop default(present)
!$acc&         async(def_queue)
      do k = 1,nvdwbloc
         iglob   = ivdw (vdwglob(k))
         i       = loc  (iglob)
         iv      = ired (iglob)
         rdn     = kred (iglob)
         rdn1    = 1.0_ti_p - rdn
         xred(i) = rdn * x(iglob) + rdn1 * x(iv)
         yred(i) = rdn * y(iglob) + rdn1 * y(iv)
         zred(i) = rdn * z(iglob) + rdn1 * z(iv)
      enddo
c
c     set the configuration for the main loop
c
      if    (use_vdwshort) then
         mode = 'SHORTVDW'
         loff2= 0.0
         fea  = fea + __use_shortRange__
          lst =>  shortvlst
         nlst => nshortvlst
      else if (use_vdwlong) then
         mode = 'VDW'
         loff2= (vdwshortcut-shortheal)**2
         fea  = fea + __use_longRange__
          lst =>  vlst
         nlst => nvlst
      else
         mode = 'VDW'
         loff2= 0.0
          lst =>  vlst
         nlst => nvlst
      end if
      call switch (mode)
      rinv    = 1.0/(cut-off)

      if (use_lambdadyn) then
         dt1lam = (1d0+dhal)**7 *2d0*scalpha*(1d0-vlambda)
         dt2lam = (1d0+ghal)    *2d0*scalpha*(1d0-vlambda)
      end if
c
c     Compute vdw pairwise scaling interactions
c
!$acc parallel loop gang vector async(def_queue)
!$acc&     present(xred,yred,zred,loc,ivdw,loc,jvdw,vir,dev
!$acc&  ,radmin,radmin4,epsilon,epsilon4,vcorrect_ik,vcorrect_scale
!$acc&  ,grplist,wgrp,mut)
!$acc&     present(ev,delambdav,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
!$acc&     reduction(+:ev,delambdav,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
      do ii = 1,n_vscale
         iglob  = vcorrect_ik(ii,1)
         kglob  = vcorrect_ik(ii,2)
         vscale = vcorrect_scale(ii)
         i      = loc(iglob)
         kbis   = loc(kglob)

         it     = jvdw(iglob)
         kt     = jvdw(kglob)
         mutik  = mut(iglob) + mut(kglob)

         do_scale4 = .false.
         vscale4   = 0

         if ( vscale.lt.0 ) then 
            vscale4 = -vscale
            vscale  = 1
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
         if (rik2<loff2.or.rik2>off2) cycle

         ! Annihilate
         if (vcouple.eq.one1.and.mutik.eq.two1) mutik=one1
 
         !replace 1-4 interactions
 20      continue
         if (do_scale4) then
            rv2  = radmin4 (kt,it)
            eps2 = epsilon4(kt,it)
         else
            rv2  =  radmin (kt,it)
            eps2 = epsilon (kt,it)
         end if

         if (use_group) 
     &      call groups2_inl(fgrp,iglob,kglob,ngrp,grplist,wgrp)

         ! Compute pairwaise
         call duo_hal
     &       (xpos,ypos,zpos,rik2,rv2,eps2,vscale,fgrp,rinv,cut2
     &       ,vdwshortcut,off,shortheal,ghal,dhal,scexp,vlambda
     &       ,scalpha,dt1lam,dt2lam,mutik,use_lambdadyn,use_group
     &       ,delambdav_,e,dedx,dedy,dedz
     &       ,ver1,fea)

         if (.not.do_scale4) then
         e    = -e
         dedx = -dedx; dedy = -dedy; dedz = -dedz;
         delambdav_ = -delambdav_
         end if

         ev  =   ev + tp2enr(e)

         if (use_lambdadyn) delambdav = delambdav + delambdav_

         call atomic_sub( dev(1,kbis),dedx )
         call atomic_sub( dev(2,kbis),dedy )
         call atomic_sub( dev(3,kbis),dedz )

         call atomic_add( dev(1,i),dedx )
         call atomic_add( dev(2,i),dedy )
         call atomic_add( dev(3,i),dedz )
c
c     increment the total van der Waals energy 
c
         g_vxx  = g_vxx + real(xpos * dedx,r_p)
         g_vxy  = g_vxy + real(ypos * dedx,r_p)
         g_vxz  = g_vxz + real(zpos * dedx,r_p)
         g_vyy  = g_vyy + real(ypos * dedy,r_p)
         g_vyz  = g_vyz + real(zpos * dedy,r_p)
         g_vzz  = g_vzz + real(zpos * dedz,r_p)

         ! deal with 1-4 Interactions
         if (vscale4.gt.0) then
            vscale    =  vscale4
            do_scale4 = .true.
            vscale4   = 0
            goto 20
         end if
      end do
c
c     find van der Waals energy and derivatives via neighbor list
c
!$acc parallel loop gang vector_length(32) default(present)
!$acc&         present(ev,delambdav,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
!$acc&         private(ai12)
!$acc&         async(def_queue)
      MAINLOOP:
     &do ii = 1, nvdwlocnl
         iivdw = vdwglobnl(ii)
         iglob = ivdw(iivdw)
         i     = loc(iglob)
         it    = jvdw(iglob)
         xi    = xred(i)
         yi    = yred(i)
         zi    = zred(i)
         muti  = mut(iglob)

         if (skipvdw12) then
            in12 = n12(iglob)
!$acc loop vector
            do j = 1,in12
               ai12(j) = i12(j,iglob)
            end do
         end if

!$acc loop vector 
         do k = 1, nlst(ii)
            kglob  = lst(k,ii)
            kbis   = loc (kglob)
            kt     = jvdw(kglob)
            mutik  = muti + mut (kglob)
            if (skipvdw12) then
               ik12 = .false.
!$acc loop seq
               do j = 1, in12
                  if (ai12(j).eq.kglob) ik12=.true.
               end do
               if (ik12) cycle
            end if
            !vscale  = 1.0_ti_p
            xpos   = xi - xred(kbis)
            ypos   = yi - yred(kbis)
            zpos   = zi - zred(kbis)
            call image_inl(xpos,ypos,zpos)
c
c     decide whether to compute the current interaction
c     and check for an interaction distance less than the cutoff
c
            rik2   = xpos**2 + ypos**2 + zpos**2
            if (rik2<loff2.or.rik2>off2) cycle

            ! Annihilate
            if (vcouple.eq.1.and.mutik.eq.two1) mutik=one1

            rv2  =  radmin (kt,it)
            eps2 = epsilon (kt,it)

            if (use_group) 
     &         call groups2_inl(fgrp,iglob,kglob,ngrp,grplist,wgrp)
            
            ! Compute pairwaise
            call duo_hal
     &          (xpos,ypos,zpos,rik2,rv2,eps2,1.0,fgrp,rinv,cut2
     &          ,vdwshortcut,off,shortheal,ghal,dhal,scexp,vlambda
     &          ,scalpha,dt1lam,dt2lam,mutik,use_lambdadyn,use_group
     &          ,delambdav_,e,dedx,dedy,dedz
     &          ,ver,fea)

            !Increment interaction energy
            ev    =   ev + tp2enr(e)

            delambdav = delambdav + delambdav_

            !Increment the van der Waals derivatives
            call atomic_sub( dev(1,kbis),dedx )
            call atomic_sub( dev(2,kbis),dedy )
            call atomic_sub( dev(3,kbis),dedz )

            call atomic_add( dev(1,i),dedx )
            call atomic_add( dev(2,i),dedy )
            call atomic_add( dev(3,i),dedz )

            !Increment the total van der Waals virial
            g_vxx = g_vxx + xpos * dedx
            g_vxy = g_vxy + ypos * dedx
            g_vxz = g_vxz + zpos * dedx
            g_vyy = g_vyy + ypos * dedy
            g_vyz = g_vyz + zpos * dedy
            g_vzz = g_vzz + zpos * dedz

         end do
      end do MAINLOOP

      call vdw_gradient_reduce

      end subroutine

c      subroutine apply_vdw_reduction_factor_
c      use atmlst    ,only: vdwglob
c      use atoms     ,only: x,y,z
c      use domdec    ,only: loc
c      use tinMemory ,only: prmem_request
c      use utilgpu   ,only: inf,def_queue
c      use vdw       ,only: ired,ivdw,kred
c     &              ,nvdwlocnl,nvdwlocnlb,nvdwbloc
c      use vdw_locArray
c      implicit none
c      integer   k,i,iglob,iv
c      real(t_p) rdn,rdn1
c
c      call prmem_request(xred    ,nvdwbloc,queue=def_queue)
c      call prmem_request(yred    ,nvdwbloc,queue=def_queue)
c      call prmem_request(zred    ,nvdwbloc,queue=def_queue)
cc
cc     apply any reduction factor to the atomic coordinates
cc
c!$acc parallel loop default(present) async(def_queue)
c      do k = 1,nvdwbloc
c         iglob    = ivdw(vdwglob(k))
c         i        = loc  (iglob)
c         iv       = ired (iglob)
c         rdn      = kred (iglob)
c         rdn1     = 1.0_ti_p - rdn
c         xredc(i)  = rdn * x(iglob) + rdn1 * x(iv)
c         yredc(i)  = rdn * y(iglob) + rdn1 * y(iv)
c         zredc(i)  = rdn * z(iglob) + rdn1 * z(iv)
c      end do
c      end subroutine

      subroutine apply_vdw_reduction_factor
      use atmlst    ,only: vdwglob
      use atoms     ,only: x,y,z
      use domdec    ,only: loc,nbloc
      use ehal1gpu_inl
      use neigh     ,only: sgl_id,slc_id,cellv_jvdw
      use tinMemory ,only: prmem_request
      use utilgpu   ,only: inf,def_queue
      use vdw       ,only: ired,ivdw,kred
     &              ,nvdwlocnl,nvdwlocnlb,nvdwbloc
      use vdw_locArray
      implicit none
      integer   k,i,iglob,iv,nvdwa
      real(t_p) rdn,rdn1,xr,yr,zr

      nvdwa = max(nvdwbloc,nvdwlocnlb)

      ! Reserve memory space
      call prmem_request(xred    ,nvdwlocnlb,queue=def_queue)
      call prmem_request(yred    ,nvdwlocnlb,queue=def_queue)
      call prmem_request(zred    ,nvdwlocnlb,queue=def_queue)
      call prmem_request(xredc   ,nbloc     ,queue=def_queue)
      call prmem_request(yredc   ,nbloc     ,queue=def_queue)
      call prmem_request(zredc   ,nbloc     ,queue=def_queue)
      call prmem_request(loc_ired,nvdwlocnlb,queue=def_queue)
      call prmem_request(loc_kred,nvdwlocnlb,queue=def_queue)
c
c     apply any reduction factor to the atomic coordinates
c
!$acc parallel loop default(present) async(def_queue)
      do k = 1,nvdwa
         if (k.le.nvdwlocnl) then
            iglob    = sgl_id(k)
            iv       = ired (iglob)
            rdn      = kred (iglob)
            rdn1     = 1.0_ti_p - rdn
            slc_id(k)   = loc(iglob)
            loc_ired(k) = loc(iv)
            loc_kred(k) = merge(rdn,1.0_ti_p,iglob.eq.iv)
            xr       = rdn * x(iglob) + rdn1 * x(iv)
            yr       = rdn * y(iglob) + rdn1 * y(iv)
            zr       = rdn * z(iglob) + rdn1 * z(iv)
            call image_inl(xr,yr,zr)
            xred(k)  = xr
            yred(k)  = yr
            zred(k)  = zr
         else if (k.le.nvdwlocnlb) then
            ! Exclusion buffer to prevent interaction compute
            slc_id(k)   = nbloc
            loc_ired(k) = nbloc
            xred(k) = inf
            yred(k) = inf
            zred(k) = inf
         end if

         if (k.le.nvdwbloc) then
            iglob    = ivdw(vdwglob(k))
            i        = loc  (iglob)
            iv       = ired (iglob)
            rdn      = kred (iglob)
            rdn1     = 1.0_ti_p - rdn
            xredc(i)  = rdn * x(iglob) + rdn1 * x(iv)
            yredc(i)  = rdn * y(iglob) + rdn1 * y(iv)
            zredc(i)  = rdn * z(iglob) + rdn1 * z(iv)
         end if
      end do
      end subroutine

      subroutine vdw_gradient_reduce
      use atoms  , only: n
      use atmlst , only: vdwglob
      use ehal1gpu_inl ,only: atomic_add
#ifdef USE_DETERMINISTIC_REDUCTION
     &                 , mdr2md
#endif
      use vdw    , only: ivdw,ired,kred,nvdwbloc
      use domdec , only: loc,rank,nbloc
      use deriv  , dev_=>de_ws2
      use inform , only: deb_Path
      use utilgpu, only: def_queue
      use group
      implicit none
      integer ii,j,i,il,iglob
      real(md_p) ded,de_red
      real(t_p) redi

12    format(2x,A)
      if (deb_Path) write(*,12) 'ehal1_gradient_reduce'

!$acc parallel loop collapse(2) default(present) async(def_queue)
      do ii = 1, nvdwbloc
         do j = 1, 3
            iglob  = ivdw(vdwglob(ii))
            i      = loc (iglob)
            ded    = mdr2md(dev_(j,i))
            ! Eject if not computed
            ! Necessary when (nbloc<n)
            if (ded.eq.0.0_re_p) cycle
            redi   = kred(iglob)
            il     = loc (ired(iglob))
            if ( i.ne.il ) then
               de_red    = redi*ded
               call atomic_add( dev(j,i ),de_red )
               de_red    = (1.0_re_p-redi)*ded
               call atomic_add( dev(j,il),de_red )
            else
               call atomic_add( dev(j,i),ded )
            end if
            dev_(j,i) = 0
         end do
      end do

      end subroutine

#ifdef _CUDA
      subroutine ehal1c_cu
      use atmlst    ,only: vdwglobnl,vdwglob
      use atoms     ,only: x,y,z,n
      use cutoff    ,only: vdwshortcut,shortheal
      use deriv     ,only: dev=>de_ws2,d_x,d_y,d_z,dr_stride,delambdav
      use domdec    ,only: loc,rank,nbloc,nproc,ndir
     &              ,xbegproc,xendproc,ybegproc,yendproc,zbegproc
     &              ,zendproc,glob
      use ehal1cu
      use energi    ,only: calc_e,ev=>ev_r
      use group     ,only: use_group,grplist,wgrp
      use inform    ,only: deb_Path,minmaxone
      use kvdws     ,only: radv,epsv
      use mutant    ,only: scalpha,scexp,vlambda,mut=>mutInt,vcouple
     &              ,nmut
      use neigh     ,only: sgl_id,slc_id,cellv_jvdw
     &              ,na,nab,nb2p,nb2p_0,nb2p_1,nbap,nbap_1
     &              ,vdw_nbl,b2pl,b2pl_1,bapl,bapl_1,abpl,abpl_1
     &              ,b_stat,b_rmid,load_nbl2mod
      use potent    ,only: use_vdwshort,use_vdwlong,use_lambdadyn
      use tinheader ,only: ti_p
      use shunt     ,only: c0,c1,c2,c3,c4,c5,off2,off,cut2,cut
      use utilcu    ,only: check_launch_kernel
      use utilgpu   ,only: def_queue,dir_queue,rec_queue,dir_stream
     &              ,rec_stream,rec_event,stream_wait_async
     &              ,warp_size,def_stream,inf,RED_BUFF_SIZE
     &              ,ered_buff=>ered_buf1,vred_buff,nred_buff,lam_buff
     &              ,reduce_energy_virial,zero_evir_red_buffer
     &              ,reduce_buffer,prmem_request,BLOCK_SIZE,nSMP
     &              ,CUDAOCCUPANCYMAXACTIVEBLOCKSPERMULTIPROCESSOR
      use vdw       ,only: ired,kred,jvdw,ivdw,radmin,radmin4,radmin_c
     &              ,epsilon,epsilon4,epsilon_c
     &              ,nvdwbloc,nvdwlocnl,nvdwlocnlb,nvdwclass
     &              ,nvdwlocnlb_pair,nvdwlocnlb_pair1,nvdwlocnlb2_pair
      use vdwpot    ,only: vcorrect_ik,vcorrect_scale,n_vscale,dhal,ghal
     &              ,radrule_i,epsrule_i
      use vdw_locArray
      use virial    ,only:use_virial,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz
      implicit none
      integer i,k
      integer iglob,iivdw,iv,hal_Gs,sizvc
      integer ierrSync
      real(t_p)  xbeg,xend,ybeg,yend,zbeg,zend
      real(t_p)  rinv,vdwshortcut2,dt1lamb,dt2lamb
      logical v0
      character*10 mode
c
      if(deb_Path) write (*,*) 'ehal1c_cu'
      def_stream = dir_stream
      xbeg  = xbegproc(rank+1)
      xend  = xendproc(rank+1)
      ybeg  = ybegproc(rank+1)
      yend  = yendproc(rank+1)
      zbeg  = zbegproc(rank+1)
      zend  = zendproc(rank+1)
      sizvc = size(vcorrect_scale)

      ! Load vdw_nbl to module global var
      call load_nbl2mod(vdw_nbl)
 
      !set the coefficients for the switching function
      mode = merge('SHORTVDW  ','VDW       ',use_vdwshort)
      call switch (mode)
      vdwshortcut2 = (vdwshortcut-shortheal)**2
      rinv    = 1.0_ti_p/(cut-off)

      !hal_Gs  = nvdwlocnlb_pair/4
      if      (use_vdwshort) then
      i       = CUDAOCCUPANCYMAXACTIVEBLOCKSPERMULTIPROCESSOR
     &          (hal_Gs,hal1short_kcu,VDW_BLOCK_DIM,0)
      else if (use_vdwlong) then
      i       = CUDAOCCUPANCYMAXACTIVEBLOCKSPERMULTIPROCESSOR
     &          (hal_Gs,hal1long_kcu,VDW_BLOCK_DIM,0)
      else
      i       = CUDAOCCUPANCYMAXACTIVEBLOCKSPERMULTIPROCESSOR
     &          (hal_Gs,hal1_v0_kcu,VDW_BLOCK_DIM,0)
      end if
      hal_Gs  = hal_Gs*nSMP*4
      if (i.ne.0) print*,'Issue detected with ehal1:cudaoccupancy',i

      if (use_lambdadyn) then
         dt1lamb = (1d0+dhal)**7 *2d0*scalpha*(1d0-vlambda)
         dt2lamb = (1d0+ghal)    *2d0*scalpha*(1d0-vlambda)
      end if

      if (dir_queue.ne.rec_queue)
     &   call stream_wait_async(rec_stream,dir_stream,rec_event)

      !if (use_virial.or.calc_e) call zero_evir_red_buffer(def_queue)
      call apply_vdw_reduction_factor
c
c     Call Vdw kernel in CUDA using C2 nblist
c
!$acc host_data use_device(xred,yred,zred,sgl_id,slc_id
!$acc&    ,loc_ired,bapl,bapl_1,abpl,abpl_1,b2pl,b2pl_1,b_stat,b_rmid
!$acc&    ,cellv_jvdw,epsilon_c,epsv,mut,radmin_c,radv,ired,kred
!$acc&    ,d_x,d_y,d_z,ered_buff,vred_buff,nred_buff,lam_buff
!$acc&    ,grplist,wgrp
!$acc&    ,vcorrect_ik,vcorrect_scale,loc,jvdw,xredc,yredc,zredc
!$acc&    ,radmin,radmin4,epsilon,epsilon4,dev
!$acc&    )
      if (use_vdwshort) then
      call hal1short_kcu<<<hal_Gs,VDW_BLOCK_DIM,0,def_stream>>>
     &             (xred,yred,zred,sgl_id,slc_id,loc_ired,b2pl_1
     &             ,bapl_1,abpl_1,b_stat,b_rmid,grplist,cellv_jvdw
     &             ,epsilon_c,radmin_c,radv,epsv,wgrp,ired,kred
     &             ,d_x,d_y,d_z,ered_buff,vred_buff,nred_buff,lam_buff
     &             ,nb2p_1,nbap_1,n,nbloc,na,nab
     &             ,nvdwclass,radrule_i,epsrule_i
     &             ,c0,c1,c2,c3,c4,c5,rinv,shortheal,ghal,dhal
     &             ,cut2,cut,vdwshortcut2,vdwshortcut,off2,off
     &             ,scexp,vlambda,scalpha,dt1lamb,dt2lamb,mut
     &             ,use_lambdadyn,use_group
     &             ,xbeg,xend,ybeg,yend,zbeg,zend,rank
     &             ! Scaling Factor Params
     &             ,n_vscale,sizvc,vcorrect_ik,vcorrect_scale,loc
     &             ,jvdw,xredc,yredc,zredc
     &             ,radmin4,epsilon4,radmin,epsilon,dev
     &             )
      call check_launch_kernel(" hal1short_kcu ")
      else if (use_vdwlong) then
      call hal1long_kcu<<<hal_Gs,VDW_BLOCK_DIM,0,def_stream>>>
     &             (xred,yred,zred,sgl_id,slc_id,loc_ired
     &             ,b2pl,bapl,abpl,b_stat,b_rmid,grplist,cellv_jvdw
     &             ,epsilon_c,radmin_c,radv,epsv,wgrp,ired,kred
     &             ,d_x,d_y,d_z,ered_buff,vred_buff,nred_buff,lam_buff
     &             ,nb2p,nbap,n,nbloc,na,nab
     &             ,nvdwclass,radrule_i,epsrule_i
     &             ,c0,c1,c2,c3,c4,c5,rinv,shortheal,ghal,dhal
     &             ,cut2,cut,vdwshortcut2,vdwshortcut,off2,off
     &             ,scexp,vlambda,scalpha,dt1lamb,dt2lamb,mut
     &             ,use_lambdadyn,use_group
     &             ,xbeg,xend,ybeg,yend,zbeg,zend,rank
     &             ! Scaling Factor Params
     &             ,n_vscale,sizvc,vcorrect_ik,vcorrect_scale,loc
     &             ,jvdw,xredc,yredc,zredc
     &             ,radmin4,epsilon4,radmin,epsilon,dev
     &             )
      call check_launch_kernel(" hal1long_kcu ")
      else
         v0 = nproc.eq.1.and.vcouple.eq.0.and.nmut.eq.0
     &      .and..not.calc_e   .and..not.use_virial
     &      .and..not.use_group.and..not.use_lambdadyn

         if (v0) then
      call hal1_v0_kcu<<<hal_Gs,VDW_BLOCK_DIM,0,def_stream>>>
     &             (xred,yred,zred,sgl_id,slc_id,loc_ired
     &             ,b2pl,bapl,abpl,b_stat,b_rmid,grplist,cellv_jvdw
     &             ,epsilon_c,radmin_c,radv,epsv,wgrp,ired,kred
     &             ,d_x,d_y,d_z,ered_buff,vred_buff,nred_buff,lam_buff
     &             ,nb2p,nbap,n,nbloc,na,nab
     &             ,nvdwclass,radrule_i,epsrule_i
     &             ,c0,c1,c2,c3,c4,c5,rinv,shortheal,ghal,dhal
     &             ,cut2,cut,vdwshortcut2,vdwshortcut,off2,off
     &             ,scexp,vlambda,scalpha,dt1lamb,dt2lamb,mut
     &             ,use_lambdadyn,use_group
     &             ,xbeg,xend,ybeg,yend,zbeg,zend,rank
     &             ! Scaling Factor Params
     &             ,n_vscale,sizvc,vcorrect_ik,vcorrect_scale,loc
     &             ,jvdw,xredc,yredc,zredc
     &             ,radmin4,epsilon4,radmin,epsilon,dev
     &             )
      call check_launch_kernel(" hal1_v0_kcu ")
         else
      call hal1_kcu<<<hal_Gs,VDW_BLOCK_DIM,0,def_stream>>>
     &             (xred,yred,zred,sgl_id,slc_id,loc_ired
     &             ,b2pl,bapl,abpl,b_stat,b_rmid,grplist,cellv_jvdw
     &             ,epsilon_c,radmin_c,radv,epsv,wgrp,ired,kred
     &             ,d_x,d_y,d_z,ered_buff,vred_buff,nred_buff,lam_buff
     &             ,nb2p,nbap,n,nbloc,na,nab
     &             ,nvdwclass,radrule_i,epsrule_i
     &             ,c0,c1,c2,c3,c4,c5,rinv,shortheal,ghal,dhal
     &             ,cut2,cut,vdwshortcut2,vdwshortcut,off2,off
     &             ,scexp,vlambda,scalpha,dt1lamb,dt2lamb,mut
     &             ,use_lambdadyn,use_group
     &             ,xbeg,xend,ybeg,yend,zbeg,zend,rank
     &             ! Scaling Factor Params
     &             ,n_vscale,sizvc,vcorrect_ik,vcorrect_scale,loc
     &             ,jvdw,xredc,yredc,zredc
     &             ,radmin4,epsilon4,radmin,epsilon,dev
     &             )
      call check_launch_kernel(" hal1_kcu ")
         end if
      end if
!$acc end host_data
      if (use_virial.or.calc_e) then
         call reduce_energy_virial(ev,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz
     &                            ,g_vzz,ered_buff,def_queue)
      end if

!$acc parallel loop collapse(2) default(present) async(def_queue)
      do i = 1,nvdwlocnl; do k = 1,3;
         dev(k,slc_id(i)) = dev(k,slc_id(i))
     &                       + d_x(i+(k-1)*dr_stride)
         d_x(i+(k-1)*dr_stride) = 0
      end do; end do;
      call vdw_gradient_reduce

      if (use_lambdadyn) then
         call reduce_buffer(lam_buff,RED_BUFF_SIZE,delambdav,def_queue)
!$acc update host(delambdav) async(def_queue)
      end if

      end subroutine
#endif

      subroutine dist(idx1,idx2)
      use atoms
      use domdec
      use ehal1gpu_inl
      integer idx1,idx2
      integer i,j
      real(t_p) xr,yr,zr

!$acc serial async default(present)
      if (loc(idx1).lt.nbloc.and.loc(idx2).lt.nbloc) then
         xr = x(idx1) - x(idx2)
         yr = y(idx1) - y(idx2)
         zr = z(idx1) - z(idx2)
         call image_inl(xr,yr,zr)
         d = sqrt(xr**2+yr**2+zr**2)
         print*, 'dist',idx1,idx2,d,rank
      end if
!$acc end serial

      end subroutine

      subroutine searchpair(nlst,lst,maxlst,int1,int2)
      use atmlst    ,only: vdwglobnl,vdwglob
      use atoms
      use cutoff    ,only: vdwshortcut,shortheal
      use cell
      use deriv     ,only: dev=>de_ws2
      use domdec
      use ehal1gpu_inl
#if defined(_OPENACC) && 1
      use ehal1cu   ,only: VDW_BLOCK_DIM,ehal1short_cu_deb
      use utilcu    ,only: check_launch_kernel
#endif
      use energi    ,only: ev=>ev_r
      use group     ,only: ngrp,wgrp,grplist,use_group
      use inform    ,only: deb_Path,dibuff
      use mutant    ,only: scalpha,scexp,vlambda,vcouple,mut=>mutInt
      use neigh 
      use potent
      use tinheader ,only: ti_p,re_p,one1,two1
      use tinMemory ,only: prmem_request
      use shunt     ,only: c0,c1,c2,c3,c4,c5,off2,off,cut2,cut
      use vdw       ,only: ired,kred,jvdw,ivdw,radmin,radmin_c
     &              ,epsilon_c,epsilon,nvdwbloc,nvdwlocnl
     &              ,nvdwclass,nvdwlocnlb,nvdwlocnlb_pair
     &              ,nvdwlocnlb2_pair,nshortvdwlocnlb2_pair
      use vdwpot    ,only: vcorrect_ik,vcorrect_scale,n_vscale,dhal,ghal
      use vdw_locArray
      use utilgpu   ,only: def_queue,dir_queue,rec_queue
     &              ,ered_buff=>ered_buf1,vred_buff,reduce_energy_virial
     &              ,zero_evir_red_buffer,prmem_request,inf
#ifdef _OPENACC
     &                    ,dir_stream,def_stream
     &                    ,rec_stream,rec_event,stream_wait_async
#endif
      use virial
      implicit none
      integer,intent(in)::maxlst,int1,int2
      integer nlst(nvdwlocnl)
      integer lst(maxlst,nvdwlocnl)

      integer i,j,k,kk,ksave
      integer kt,kglob,kbis,kvloc,kv,ki
      integer iglob,iivdw,hal_gs,lst_start
      integer ii,iv,it,ivloc
      integer,save:: ncall=0
      integer nnvlst,nnvlst2
      integer nn12,nn13,nn14,ntot
      integer(1) muti,mutik
      integer xw,yw,zw,proc
      integer icell_len,icell,imageCell,distImage2Cell
      integer,parameter:: ncell2buff=2
      integer ncell2buffb
      real(t_p) xi,yi,zi,redi,e,de
      real(t_p) rdn,rdn1,redk
      real(t_p) rik2,rinv
      real(t_p) dedx,dedy,dedz
      real(t_p) invrho,rv7orho
      real(t_p) dtau,gtau,tau,tau7,rv7
      real(t_p) rv2,eps2
      real(t_p) xpos,ypos,zpos
      real(t_p) dtaper,taper
      real(t_p) vscale,vscale4
      real(t_p) xr,yr,zr,mbuf,vbuf,bigbuf
      real(t_p) lenx,leny,lenz
      real(t_p) lenx_cell,leny_cell,lenz_cell
      real(t_p) xmin,xmax,ymin,ymax,zmin,zmax
      real(t_p)  vdwshortcut2
      real(t_p)  xbeg,xend,ybeg,yend,zbeg,zend
      logical   do_scale4
      character*10 mode
!$acc routine(distprocpart1)
c
      write (*,*) 'searchpair',int1,int2
      ncall = ncall + 1

#ifdef _OPENACC
      if (dir_queue.ne.rec_queue) then
!!$acc wait(rec_queue) async(dir_queue)
         call stream_wait_async(rec_stream,dir_stream,rec_event)
      end if
#endif

      call prmem_request(xred,nbloc,queue=def_queue)
      call prmem_request(yred,nbloc,queue=def_queue)
      call prmem_request(zred,nbloc,queue=def_queue)
      if (ncall.eq.1) then
!$acc wait
      end if

!$acc data present(xred,yred,zred)
!$acc&     present(loc,ired,kred,x,y,z,vdwglobnl,ivdw,loc,jvdw,
!$acc&  vir,dev,vdwglob,lst,nlst,radmin,epsilon,mut)
!$acc&     present(g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)

c
c     set the coefficients for the switching function
c
      mode = 'SHORTVDW'
      call switch (mode)
      rinv = 1.0/(cut-off)
      vdwshortcut2 = (vdwshortcut-shortheal)**2
c
c     apply any reduction factor to the atomic coordinates
c
!$acc parallel loop default(present)
!$acc&         async(def_queue)
      do k = 1,nvdwbloc
         iglob   = ivdw (vdwglob (k))
         i       = loc  (iglob)
         iv      = ired (iglob)
         rdn     = kred (iglob)
         rdn1    = 1.0_ti_p - rdn
         xred(i) = rdn * x(iglob) + rdn1 * x(iv)
         yred(i) = rdn * y(iglob) + rdn1 * y(iv)
         zred(i) = rdn * z(iglob) + rdn1 * z(iv)
      end do

!$acc parallel loop default(present)
        do i = 1, nbloc
           iglob = glob(i)
           call distprocpart1(iglob,rank,rv2,.true.,x,y,z)
       if (iglob.eq.int1) print*, "found glob",int1,ii,real(rv2,4),rank
       if (iglob.eq.int2) print*, "found glob",int2,ii,real(rv2,4),rank
        end do

!$acc parallel loop
      do ii = 1,nvdwlocnl
         iivdw = vdwglobnl(ii)
         iglob = ivdw(iivdw)
         if (iglob.eq.int1) print*, "found vdwglobnl",int1,ii,rank
         if (iglob.eq.int2) print*, "found vdwglobnl",int2,ii,rank
      end do

      ncell2buffb = ncell2buff
      xmin  = xbegproc(rank+1)
      xmax  = xendproc(rank+1)
      ymin  = ybegproc(rank+1)
      ymax  = yendproc(rank+1)
      zmin  = zbegproc(rank+1)
      zmax  = zendproc(rank+1)

      do i = 1, nbig_recep
        proc = pbig_recep(i)
        if (xbegproc(proc+1).le.xmin) xmin = xbegproc(proc+1)
        if (xendproc(proc+1).ge.xmax) xmax = xendproc(proc+1)
        if (ybegproc(proc+1).le.ymin) ymin = ybegproc(proc+1)
        if (yendproc(proc+1).ge.ymax) ymax = yendproc(proc+1)
        if (zbegproc(proc+1).le.zmin) zmin = zbegproc(proc+1)
        if (zendproc(proc+1).ge.zmax) zmax = zendproc(proc+1)
      end do

      if ((mbuf2).gt.vbuf2) then
        bigbuf = sqrt(mbuf2)/ncell2buffb
      else
        bigbuf = sqrt(vbuf2)/ncell2buffb
      end if
      lenx       = abs(xmax-xmin)
      nx_cell    = max(2*ncell2buffb+1,int(lenx/(bigbuf)))
      lenx_cell  = lenx/nx_cell
      leny       = abs(ymax-ymin)
      ny_cell    = max(2*ncell2buffb+1,int(leny/(bigbuf)))
      leny_cell  = leny/ny_cell
      lenz       = abs(zmax-zmin)
      nz_cell    = max(2*ncell2buffb+1,int(lenz/(bigbuf)))
      lenz_cell  = lenz/nz_cell

!$acc parallel loop async
      do i = 1,ncell_tot
         cell_len(i) = 0
      end do
c!$acc parallel loop async
c      do i = 1,n
c         indcelltemp(i) = 0
c      end do
       if (rank.eq.0) then
c         print*,'nx',nx_cell,ny_cell,nz_cell
c         print*,'lx',lenx,leny,lenz,bigbuf
          print*,'cell h',xcell2,ycell2,zcell2
!$acc serial async present(xcell2,ycell2,zcell2)
          print*,'cell d',xcell2,ycell2,zcell2
!$acc end serial
       end if

!$acc parallel loop async
      do i = 1,nlocnl
         iglob = ineignl(i)
         xr    = x(iglob)
         yr    = y(iglob)
         zr    = z(iglob)
         call image_inl(xr,yr,zr)
         if ((xcell2-xr).lt.eps_cell) xr = xr-0.05*lenx_cell
         if ((ycell2-yr).lt.eps_cell) yr = yr-0.05*leny_cell
         if ((zcell2-zr).lt.eps_cell) zr = zr-0.05*lenz_cell
         xw  = int((xr-xmin)/lenx_cell)
         yw  = int((yr-ymin)/leny_cell)
         zw  = int((zr-zmin)/lenz_cell)
         if (iglob.eq.int1) print*,'box1',rank,xw,yw,zw
         if (iglob.eq.int1) print*,'box1',rank,x(iglob),xr
     &                     ,y(iglob),yr,z(iglob),zr
         if (iglob.eq.int2) print*,'box2',rank,xw,yw,zw
         if (iglob.eq.int2) print*,'box2',rank,x(iglob),xr
     &                     ,y(iglob),yr,z(iglob),zr
         icell = (xw + nx_cell*yw + nx_cell*ny_cell*zw) + 1
         repartcell (iglob) = icell
!$acc atomic
         cell_len (icell) = cell_len(icell)+1
c         icell_len        = cell_len(icell)
c!$acc end atomic
      end do
c
c     find van der Waals energy and derivatives via neighbor list
c
#ifdef _OPENACC
#if 0
!$acc parallel loop num_gangs(2) vector_length(32)
!$acc&         async(def_queue)
      MAINLOOP:
     &do ii = 1, nvdwlocnl
         iivdw = vdwglobnl(ii)
         iglob = ivdw(iivdw)
         i     = loc(iglob)
         it    = jvdw(iglob)
         xi    = xred(i)
         yi    = yred(i)
         zi    = zred(i)
         muti  = mut(iglob)

!$acc loop vector 
         do k = 1, nlst(ii)
            kglob  = lst(k,ii)
            kbis   = loc (kglob)
            kt     = jvdw (kglob)
            mutik  = muti + mut(kglob)
            !vscale = 1.0_ti_p
c
c     compute the energy contribution for this interaction
c
            xpos   = xi - xred(kbis)
            ypos   = yi - yred(kbis)
            zpos   = zi - zred(kbis)
            call image_inl(xpos,ypos,zpos)
c
c     decide whether to compute the current interaction
c     and check for an interaction distance less than the cutoff
c
            rik2   = xpos**2 + ypos**2 + zpos**2
         if (iglob.eq.int1.and.kglob.lt.int1+10) then
             print*,"comp0",kglob,rank,nlst(ii),rik2
             print*,"comp0",kglob,rank,xpos,ypos,zpos
         end if
            if (rik2>off2) cycle
            ! Annihilate
            if (vcouple.eq.1.and.mutik.eq.two1) mutik=one1

            rv2  =  radmin (kt,it)
            eps2 = epsilon (kt,it)

            call ehal1_couple(xpos,ypos,zpos,rik2,rv2,eps2,1.0_ti_p
     &                       ,cut2,rinv,off,ghal,dhal
     &                       ,scexp,vlambda,scalpha,mutik
     &                       ,e,dedx,dedy,dedz
     &                       ,__tver__,__tfea__)

         if (iglob.eq.int1.and.kglob.lt.int1+10) then
             print*,"comp",kglob,rank,dedx,dedy,dedz
         end if

         end do
      end do MAINLOOP
#else
      xmin = xbegproc(rank+1)
      xmax = xendproc(rank+1)
      ymin = ybegproc(rank+1)
      ymax = yendproc(rank+1)
      zmin = zbegproc(rank+1)
      zmax = zendproc(rank+1)
      lst_start = 2*nvdwlocnlb_pair+1
c
c     apply any reduction factor to the atomic coordinates
c
!$acc parallel loop default(present) async(def_queue)
      do k = 1,nvdwlocnlb
         if (k.le.nvdwlocnl) then
            iglob    = sgl_id(k)
            iv       = ired (iglob)
            rdn      = kred (iglob)
            rdn1     = 1.0_ti_p - rdn
            slc_id(k)   = loc(iglob)
            loc_ired(k) = loc(iv)
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
            slc_id(k)   = nbloc
            loc_ired(k) = nbloc
            xred(k) = inf
            yred(k) = inf
            zred(k) = inf
         end if
      end do

      call zero_evir_red_buffer(def_queue)
c
c     Call Vdw kernel in CUDA using C2 nblist
c
!$acc host_data use_device(xred,yred,zred,sgl_id,slc_id
!$acc&    ,loc_ired,ivblst,vblst,cellv_jvdw,grplist,wgrp
!$acc&    ,epsilon_c,radmin_c,mut,ired,kred,dev,ered_buff,vred_buff
#ifdef TINKER_DEBUG
!$acc&    ,inter
#endif
!$acc&    )

      if (.true.) then

      hal_Gs = nshortvdwlocnlb2_pair/4
      hal_Gs = 1
      call ehal1short_cu_deb <<<hal_Gs,VDW_BLOCK_DIM,0,def_stream>>>
     &           (xred,yred,zred,sgl_id,slc_id,loc_ired
     &           ,ivblst,vblst(lst_start),cellv_jvdw
     &           ,epsilon_c,radmin_c
     &           ,ired,kred,grplist,wgrp,dev,ered_buff,vred_buff
     &           ,nvdwlocnlb2_pair,n,nbloc,nvdwlocnl,nvdwlocnlb
     &           ,nvdwclass,rinv
     &           ,c0,c1,c2,c3,c4,c5,cut2,off2,off
     &           ,scexp,vlambda,scalpha,mut
     &           ,shortheal,ghal,dhal,use_vdwshort,use_group
     &           ,xmin,xmax,ymin,ymax,zmin,zmax
     &           ,int1,int2,rank
#ifdef TINKER_DEBUG
     &           ,inter
#endif
     &           )
      call check_launch_kernel(" ehal1short_cu_deb ")

      else if (use_vdwlong) then

      hal_Gs = nvdwlocnlb2_pair/4
      call ehal1long_cu <<<hal_Gs,VDW_BLOCK_DIM,0,def_stream>>>
     &           (xred,yred,zred,sgl_id,slc_id,loc_ired
     &           ,ivblst,vblst(lst_start),cellv_jvdw,epsilon_c,radmin_c
     &           ,ired,kred,dev,ered_buff,vred_buff
     &           ,nvdwlocnlb2_pair,n,nbloc,nvdwlocnl,nvdwlocnlb
     &           ,nvdwclass
     &           ,c0,c1,c2,c3,c4,c5,cut2,cut,off2,off,vdwshortcut2
     &           ,scexp,vlambda,scalpha,mut
     &           ,vdwshortcut,shortheal,ghal,dhal,use_vdwshort
     &           ,xmin,xmax,ymin,ymax,zmin,zmax
#ifdef TINKER_DEBUG
     &           ,inter,rank
#endif
     &           )
      call check_launch_kernel(" ehal1long_cu ")
      end if

!$acc end host_data

      call reduce_energy_virial(ev,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz
     &                         ,ered_buff,def_queue)
#endif
#endif

#if 0
      ! Save positions before removing comments on this kernel
c!$acc parallel loop async
c!$acc&         present(dibuff,xold_nl,yold_nl,zold_nl)
c      do i = 1,nbloc
c         iglob    = dibuff(i)
c         x(iglob) = xold_nl(iglob)
c         y(iglob) = yold_nl(iglob)
c         z(iglob) = zold_nl(iglob)
c         loc(iglob) = i
c      end do
      call reassignrespa(2,2)
      call kvdw(.false.,-1)
      call nblist(0)

      mode = 'SHORTVDW'
      call switch (mode)
c
c     apply any reduction factor to the atomic coordinates
c
!$acc parallel loop default(present)
!$acc&         async(def_queue)
      do k = 1,nvdwbloc
         iglob   = ivdw (vdwglob (k))
         i       = loc  (iglob)
         iv      = ired (iglob)
         rdn     = kred (iglob)
         rdn1    = 1.0_ti_p - rdn
         xred(i) = rdn * x(iglob) + rdn1 * x(iv)
         yred(i) = rdn * y(iglob) + rdn1 * y(iv)
         zred(i) = rdn * z(iglob) + rdn1 * z(iv)
      end do
c
c     find van der Waals energy and derivatives via neighbor list
c
!$acc parallel loop num_gangs(2) vector_length(32)
!$acc&         async(def_queue)
      do ii = 1, nvdwlocnl
         iivdw = vdwglobnl(ii)
         iglob = ivdw(iivdw)
         i     = loc(iglob)
         it    = jvdw(iglob)
         xi    = xred(i)
         yi    = yred(i)
         zi    = zred(i)
         muti  = mut(iglob)

!$acc loop vector 
         do k = 1, nlst(ii)
            kglob  = lst(k,ii)
            kbis   = loc (kglob)
            kt     = jvdw (kglob)
            mutik  = muti + mut(kglob)
            !vscale = 1.0_ti_p
c
c     compute the energy contribution for this interaction
c
            xpos   = xi - xred(kbis)
            ypos   = yi - yred(kbis)
            zpos   = zi - zred(kbis)
            call image_inl(xpos,ypos,zpos)
c
c     decide whether to compute the current interaction
c     and check for an interaction distance less than the cutoff
c
            rik2   = xpos**2 + ypos**2 + zpos**2
         if (iglob.eq.int1.and.kglob.le.int1+10) then
             print*,"comb0",kglob,rank,nlst(ii),rik2
             print*,"comp0",kglob,rank,xpos,ypos,zpos
         end if
            if (rik2>off2) cycle
            ! Annihilate
            if (vcouple.eq.1.and.mutik.eq.two1) mutik=one1

            rv2  =  radmin (kt,it)
            eps2 = epsilon (kt,it)

            call ehal1_couple(xpos,ypos,zpos,rik2,rv2,eps2,1.0_ti_p
     &                       ,cut2,cut,off,ghal,dhal
     &                       ,scexp,vlambda,scalpha,mutik
     &                       ,e,dedx,dedy,dedz
     &                       ,__tver__,__tfea__)

         if (iglob.eq.int1.and.kglob.le.int1+10) then
             print*,"comb",kglob,rank,dedx,dedy,dedz
         end if
         end do
      end do
#endif

!$acc parallel loop default(present)
!$acc&         async(def_queue)
      do k = 1,nvdwbloc
         iglob   = ivdw (vdwglob (k))
         i       = loc  (iglob)
         iv      = ired (iglob)
         rdn     = kred (iglob)
         rdn1    = 1.0_ti_p - rdn
         xred(i) = rdn * x(iglob) + rdn1 * x(iv)
         yred(i) = rdn * y(iglob) + rdn1 * y(iv)
         zred(i) = rdn * z(iglob) + rdn1 * z(iv)
      end do

      call searchVdwScaled(xred,yred,zred,
     &     g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz,int1,int2)

c      call vdw_gradient_reduce
!$acc end data
!$acc wait
      end subroutine

      subroutine searchVdwScaled(xred,yred,zred,
     &           vxx,vxy,vxz,vyy,vyz,vzz,int1,int2)

      use atmlst    ,only: vdwglobnl
      use deriv     ,only: dev=>de_ws2
      use domdec    ,only: loc,rank,nbloc
      use cutoff    ,only: shortheal
      use ehal1gpu_inl
      use energi    ,only: ev
      use inform    ,only: deb_Path
      use mutant    ,only: scexp,scalpha,vlambda,vcouple,mut=>mutInt
      use tinheader ,only: ti_p,one1,two1
      use shunt     ,only: c0,c1,c2,c3,c4,c5,off2,off,cut2,cut
      use vdw       ,only: ired,kred,jvdw,ivdw,radmin,radmin4,
     &                     epsilon,epsilon4,nvdwbloc
      use vdwpot    ,only: vcorrect_ik,vcorrect_scale,n_vscale,dhal,ghal
      use utilgpu   ,only: def_queue
      use virial
      implicit none
      integer i,j,k,kk,ksave
      integer kt,kglob,kbis,kvloc,kv,ki
      integer iglob,iivdw
      integer ii,iv,it,ivloc
      integer nnvlst,nnvlst2
      integer nn12,nn13,nn14,ntot
      integer interac
      integer(1) muti,mutik
      real(t_p)  xi,yi,zi,redi,e,de
      real(t_p)  rdn,rdn1,redk
      real(t_p)  rik2,rinv
      real(t_p)  dedx,dedy,dedz
      real(t_p)  invrho,rv7orho
      real(t_p)  dtau,gtau,tau,tau7,rv7
      real(t_p)  rv2,eps2
      real(t_p)  xpos,ypos,zpos
      real(t_p)  dtaper,taper
      real(t_p)  vscale,vscale4
      logical    do_scale4
      character*10 mode

      integer int1,int2
      real(t_p),intent(in):: xred(nbloc)
      real(t_p),intent(in):: yred(nbloc)
      real(t_p),intent(in):: zred(nbloc)
      real(r_p)  vxx,vxy,vxz
      real(r_p)  vyy,vyz,vzz

      ! Scaling factor correction loop
      if(deb_Path) write(*,*) "searchVdwScaled",n_vscale
      rinv = 1.0/(cut-off)

!$acc parallel loop async(def_queue)
!$acc&     gang vector
!$acc&     present(xred,yred,zred,vxx,vxy,vxz,vyy,vyz,vzz)
!$acc&     present(loc,ivdw,loc,jvdw,vir,dev,radmin,mut,
!$acc&  radmin4,epsilon,epsilon4,vcorrect_ik,vcorrect_scale)
      do ii = 1,n_vscale
         iglob  = vcorrect_ik(ii,1)
         kglob  = vcorrect_ik(ii,2)
         vscale = vcorrect_scale(ii)
         i      = loc(iglob)
         kbis   = loc(kglob)

         it     = jvdw(iglob)
         kt     = jvdw(kglob)
         muti   = mut(iglob)
         mutik  = muti + mut(kglob)

         do_scale4 = .false.
         vscale4   = 0

         if (vscale.lt.0) then 
            vscale4 = -vscale
            vscale  = 1
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

         ! Annihilate
         if (vcouple.eq.1.and.mutik.eq.two1) mutik=one1
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

         call ehal1_couple(xpos,ypos,zpos,rik2,rv2,eps2,vscale
     &                    ,cut2,rinv,off,ghal,dhal
     &                    ,scexp,vlambda,scalpha,mutik
     &                    ,e,dedx,dedy,dedz
     &                    ,__use_ene__,__use_softcore__)
c        call ehal1_couple_short(xpos,ypos,zpos,rik2,rv2,eps2,vscale
c    &                    ,cut2,off
c    &                    ,scexp,vlambda,scalpha,mutik
c    &                    ,shortheal,ghal,dhal,e,dedx,dedy,dedz)

         if (.not.do_scale4) then
         e    = -e
         dedx = -dedx; dedy = -dedy; dedz = -dedz;
         end if

         if (iglob.eq.int1.or.kglob.eq.int1) then
             print*,iglob,kglob,rik2,dedx,dedy,dedz,"cor1",rank
         end if
         if (kglob.eq.int2.or.iglob.eq.int2) then
             print*,iglob,kglob,rik2,dedx,dedy,dedz,"cor2",rank
         end if

         ! deal with 1-4 Interactions
         if (vscale4.gt.0) then
            vscale    =  vscale4
            do_scale4 = .true.
            vscale4   = 0
            goto 20
         end if
      end do
      end
