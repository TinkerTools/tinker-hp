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
#include "tinker_macro.h"
      module ehal3gpu_inl
        implicit none
        contains
#include "convert.f.inc"
#include "image.f.inc"
#include "groups.inc.f"
#include "pair_ehal.inc.f"
      end module

      subroutine ehal3gpu
      use analyz
      use atoms
      use domdec
      use energi
      use inform
      use iounit
      use interfaces ,only: ehal3c_p
      use potent
      use tinheader  ,only:ti_p,re_p
      use timestat
      use vdwpot
      use mpi

      implicit none
      integer i
      real(r_p) elrc,aelrc
      character*11 mode
      logical fullrange
c
      fullrange = .not.(use_vdwshort.or.use_vdwlong)
c
c     choose the method for summing over pairwise interactions
c
      call timer_enter( timer_ehal3 )
      call ehal3c_p
c
c     apply long range van der Waals correction if desired
c
      if (use_vcorr.and.(fullrange.or.use_vdwlong)) then
         mode = "VDW"
!$acc data copyout(elrc) present(ev)
         call evcorr (mode,elrc)
!$acc serial async
         ev = ev + elrc
!$acc end serial
!$acc wait
c        aelrc = elrc / real(n,r_p)
c        do i = 1, nbloc
c           aev(i) = aev(i) + aelrc
c        end do
!$acc end data
         if (rank.eq.0.and.verbose) then
            if (elrc.ne.0.0_ti_p.and.app_id.eq.analyze_a) then
               write (iout,10)  elrc
   10          format (/,' Long Range vdw Correction :',9x,f12.4)
            end if
         end if
      end if
      call timer_exit( timer_ehal3,quiet_timers )
      end
c
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine ehal3cgpu  --  buffered 14-7 analysis via list  ##
c     ##                                                             ##
c     #################################################################
c
c     "ehal3cvec" calculates the buffered 14-7 van der Waals energy
c     and also partitions the energy among the atoms using a
c     pairwise neighbor list
c
c
      subroutine ehal3c_ac
      use action    ,only: nev,nev_
      use atmlst    ,only: vdwglobnl,vdwglob
      use atoms     ,only: x,y,z,n
      use couple    ,only: i12,n12
      use cutoff    ,only: vdwshortcut,shortheal
      use domdec    ,only: loc,rank,nbloc
      use ehal3gpu_inl
      use energi    ,only: ev=>ev_r
      use group     ,only: use_group,grplist,wgrp,ngrp
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
      parameter( ver =__use_ene__+__use_act__
     &         , ver1=ver+__use_sca__ )
c     &         , fea =__use_mpi__+__use_groups__+__use_softcore__)

c
      if(deb_Path) write (*,*) 'ehal3c_ac',use_vdwshort,use_vdwlong
      ncall = ncall + 1

      fea    = __use_mpi__
      if (use_group) fea = fea + __use_groups__
      if (nmut.ne.0) fea = fea + __use_softcore__
      !if (use_lambdadyn) fea = fea + __use_lambdadyn__

#ifdef _OPENACC
      if (dir_queue.ne.rec_queue) then
!!$acc wait(rec_queue) async(dir_queue)
         call stream_wait_async(rec_stream,dir_stream,rec_event)
      end if
#endif
      def_queue = dir_queue
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

      !if (use_lambdadyn) then
      !   dt1lam = (1d0+dhal)**7 *2d0*scalpha*(1d0-vlambda)
      !   dt2lam = (1d0+ghal)    *2d0*scalpha*(1d0-vlambda)
      !end if
c
c     Compute vdw pairwise scaling interactions
c
!$acc parallel loop gang vector async(def_queue)
!$acc&     present(xred,yred,zred,loc,ivdw,loc,jvdw
!$acc&  ,radmin,radmin4,epsilon,epsilon4,vcorrect_ik,vcorrect_scale
!$acc&  ,grplist,wgrp,mut)
!$acc&     present(ev,nev_)
!$acc&     reduction(+:ev,nev_)
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
         if (vscale.eq.1.0) nev_=nev_-1
         end if

         ev   =  ev + tp2enr(e)

         ! deal with 1-4 Interactions
         if (vscale4.gt.0) then
            vscale    = vscale4
            nev_      = nev_+1
            do_scale4 = .true.
            vscale4   = 0
            goto 20
         end if
      end do
c
c     find van der Waals energy and derivatives via neighbor list
c
!$acc parallel loop gang vector_length(32) default(present)
!$acc&         present(radmin,epsilon,wgrp)
!$acc&         present(ev,nev_) reduction(+:ev,nev_)
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

!$acc loop vector reduction(+:ev,nev_)
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
            ev    =   ev  + tp2enr(e)
            nev_  =  nev_ + 1.0

         end do
      end do MAINLOOP

!$acc serial async(def_queue) present(nev,nev_)
      nev = nev + int(nev_)
!$acc end serial

      end subroutine

#ifdef _CUDA
      subroutine ehal3c_cu
      use action    ,only: nev
      use atmlst    ,only: vdwglobnl,vdwglob
      use atoms     ,only: x,y,z,n
      use cutoff    ,only: vdwshortcut,shortheal
      use deriv     ,only: dev=>de_ws2,d_x,d_y,d_z,dr_stride
      use domdec    ,only: loc,rank,nbloc,nproc
     &              ,xbegproc,xendproc,ybegproc,yendproc,zbegproc
     &              ,zendproc,glob
      use ehal1cu
      use energi    ,only: calc_e,ev=>ev_r
      use group     ,only: use_group,grplist,wgrp
      use inform    ,only: deb_Path,minmaxone
      use kvdws     ,only: radv,epsv
      use mutant    ,only: scalpha,scexp,vlambda,mut=>mutInt
      use neigh     ,only: sgl_id,slc_id,cellv_jvdw
     &              ,na,nab,nb2p,nb2p_1,nbap,nbap_1
     &              ,vdw_nbl,b2pl,b2pl_1,bapl,bapl_1,abpl,abpl_1
     &              ,b_stat,b_rmid,load_nbl2mod
      use potent    ,only: use_lambdadyn,use_vdwshort,use_vdwlong
      use tinheader ,only: ti_p
      use shunt     ,only: c0,c1,c2,c3,c4,c5,off2,off,cut2,cut
      use utilcu    ,only: check_launch_kernel
      use utilgpu   ,only: def_queue,dir_queue,rec_queue,dir_stream
     &              ,rec_stream,rec_event,stream_wait_async
     &              ,WARP_SIZE,def_stream,inf
     &              ,ered_buff=>ered_buf1,vred_buff,nred_buff,lam_buff
     &              ,reduce_energy_action
     &              ,prmem_request,BLOCK_SIZE,nSMP
     &              ,CUDAOCCUPANCYMAXACTIVEBLOCKSPERMULTIPROCESSOR
      use vdw       ,only: ired,kred,jvdw,ivdw,radmin,radmin4,radmin_c
     &              ,epsilon,epsilon4,epsilon_c
     &              ,nvdwbloc,nvdwlocnl,nvdwlocnlb,nvdwclass
     &              ,nvdwlocnlb_pair,nvdwlocnlb_pair1,nvdwlocnlb2_pair
      use vdwpot    ,only: vcorrect_ik,vcorrect_scale,n_vscale,dhal,ghal
     &              ,radrule_i,epsrule_i
      use vdw_locArray
      implicit none
      integer i,k
      integer iglob,iivdw,iv,hal_Gs,sizvc
      integer ierrSync
      real(t_p)  xbeg,xend,ybeg,yend,zbeg,zend
      real(t_p)  rinv,vdwshortcut2,dt1lamb,dt2lamb
      character*10 mode
c
      if(deb_Path) write (*,*) 'ehal3c_cu',use_vdwshort,use_vdwlong
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
      mode = merge('SHORTVDW','VDW',use_vdwshort)
      call switch (mode)
      vdwshortcut2 = merge(0.0_ti_p,(vdwshortcut-shortheal)**2
     &                    ,use_vdwshort)
      rinv    = 1.0_ti_p/(cut-off)
      !hal_Gs  = nvdwlocnlb_pair/4
      i       = CUDAOCCUPANCYMAXACTIVEBLOCKSPERMULTIPROCESSOR
     &          (hal_Gs,hal3_kcu,VDW_BLOCK_DIM,0)
      hal_Gs  = hal_Gs*nSMP*4
      if (i.ne.0) print*,'Issue detected with ehal1:cudaoccupancy',i

      if (dir_queue.ne.rec_queue)
     &   call stream_wait_async(rec_stream,dir_stream,rec_event)

      call apply_vdw_reduction_factor
c
c     Call Vdw kernel in CUDA using C2 nblist
c
!$acc host_data use_device(xred,yred,zred,sgl_id,slc_id
!$acc&    ,loc_ired,bapl,abpl,b2pl,b_stat,b_rmid
!$acc&    ,cellv_jvdw,epsilon_c,epsv,mut,radmin_c,radv,ired,kred
!$acc&    ,d_x,d_y,d_z,ered_buff,vred_buff,nred_buff,lam_buff
!$acc&    ,grplist,wgrp
!$acc&    ,vcorrect_ik,vcorrect_scale,loc,jvdw,xredc,yredc,zredc
!$acc&    ,radmin,radmin4,epsilon,epsilon4,dev
!$acc&    )
      if (use_vdwshort) then
!$acc host_data use_device(b2pl_1,bapl_1,abpl_1)
        call hal3s_kcu<<<hal_Gs,VDW_BLOCK_DIM,0,def_stream>>>
     &             (xred,yred,zred,sgl_id,slc_id,loc_ired,b2pl_1,bapl_1
     &             ,abpl_1,b_stat,b_rmid,grplist,cellv_jvdw
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
!$acc end host_data
        call check_launch_kernel(" hal3s_kcu ")
      else if (use_vdwlong) then
              call hal3l_kcu<<<hal_Gs,VDW_BLOCK_DIM,0,def_stream>>>
     &             (xred,yred,zred,sgl_id,slc_id,loc_ired,b2pl,bapl
     &             ,abpl,b_stat,b_rmid,grplist,cellv_jvdw
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
        call check_launch_kernel(" hal3l_kcu ")
      else
        call hal3_kcu<<<hal_Gs,VDW_BLOCK_DIM,0,def_stream>>>
     &             (xred,yred,zred,sgl_id,slc_id,loc_ired,b2pl,bapl
     &             ,abpl,b_stat,b_rmid,grplist,cellv_jvdw
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
      call check_launch_kernel(" hal3_kcu ")
      endif
!$acc end host_data

      call reduce_energy_action(ev,nev,ered_buff,def_queue)

      end subroutine
#endif
