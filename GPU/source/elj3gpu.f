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
#include "tinker_macro.h"
      module elj3gpu_inl
        integer(1) one1,two1
        parameter( one1=1, two1=2 )
        contains
#include "convert.f.inc"
#include "image.f.inc"
#include "switch_respa.f.inc"
#include "groups.inc.f"
#include "pair_lj.inc.f"
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
      character*11 mode
c
      call timer_enter( timer_elj3 )
      call elj3c_p
      call timer_exit( timer_elj3,quiet_timers )
c
c     apply long range van der Waals correction if desired
c
      if (use_vcorr) then
         mode = 'VDW'
         call evcorr (mode,elrc)
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
c     ## subroutine elj3c_ac -- Lennard-Jones analysis via list  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "elj3c_ac" calculates the Lennard-Jones van der Waals energy
c     and also partitions the energy among the atoms using a
c     pairwise neighbor list
c
c
      subroutine elj3c_ac
      use action
      use atmlst
      use atoms
      use bound
      use couple
      use cutoff    ,only: vdwshortcut,shortheal
      use deriv     ,only: dev,delambdav
      use domdec
      use elj3gpu_inl
      use energi    ,only: ev=>ev_r
      use group
      use inform
      use inter
      use iounit
      use molcul
      use mutant
      use neigh
      use potent    ,only: use_vdwshort,use_vdwlong,use_lambdadyn
      use shunt
      use tinMemory
      use tinTypes
      use usage
      use utilgpu
      use vdw
      use vdw_locArray
      use vdwpot
      use virial
      implicit none
      logical     usei,ik12,do_scale4
      integer(1)  muti,mutk,mutik
      integer     i,j,iglob,kglob,kbis,iivdw
     &           ,ii,iv,it,ivloc,kvloc,kk,kv,kt
     &           ,in12,ai12(maxvalue),ver,ver1,fea
      real(t_p)   e,de,p6,p12,eps
     &           ,rinv,rv,rdn,fgrp,loff2
     &           ,xi,yi,zi,xr,yr,zr
     &           ,rik,rik2,rik3,rik4,rik5
     &           ,taper,dtaper
     &           ,galpha,glamb,lambdavt,lambdavt1
     &           ,delambdav_
     &           ,rv2,eps2,vscale,vscale4
      type(real3) ded
      character*11 mode
      integer    ,pointer:: lst(:,:),nlst(:)

      parameter( ver=__use_ene__+__use_act__
     &         ,ver1=        ver+__use_sca__)
c
 1000 format(' Warning, system moved too much since last neighbor list'
     $   ' update, try lowering nlupdate VDW')
      if (deb_Path) write(*,*) 'elj3c_ac'

      ! set elj3c_ac Configuration
      fea = __use_mpi__
      if (nmut.ne.0)     fea = fea + __use_softcore__
      if (use_group)     fea = fea + __use_groups__
 
      if      (use_vdwshort) then
         fea   = fea + __use_shortRange__
         lst   => shortvlst
         nlst  => nshortvlst
         loff2 = 0.0
         !set the coefficients for the switching function
         mode  = 'SHORTVDW'
         call  switch (mode)
      else if (use_vdwlong ) then
         fea   = fea + __use_longRange__
         lst   => vlst
         nlst  => nvlst
         !set the coefficients for the switching function
         mode  = 'VDW'
         call  switch (mode)
         loff2 = (vdwshortcut-shortheal)**2
      else
         lst   => vlst
         nlst  => nvlst
         loff2 = 0.0
         !set the coefficients for the switching function
         mode  = 'VDW'
         call  switch (mode)
      end if
      rinv   = 1.0/(cut-off)
      if (nmut.ne.0) then
         galpha    = scalpha/(2d0**(sck/6d0))
         glamb     = 1.0-vlambda
         lambdavt  = vlambda ** sct
         lambdavt1 = sct*(vlambda**(sct-1.0))
      end if
c
c     Scaling interaction correction for Lennard-Jones
c
!$acc parallel loop async(dir_queue) gang vector
!$acc&         present(loc,ired,x,y,z,loc,jvdw,vir,radmin,mutInt,grplist
!$acc&     ,wgrp,radmin4,epsilon,epsilon4,vcorrect_ik,vcorrect_scale
!$acc&     ,nev_,ev)
!$acc&         private(ded,e)
!$acc&         reduction(+:nev_,ev)
      do ii = 1,n_vscale
         iglob  = vcorrect_ik(ii,1)
         kglob  = vcorrect_ik(ii,2)
         vscale = vcorrect_scale(ii)
         i      = loc(iglob)
         kbis   = loc(kglob)

         it     = jvdw(iglob)
         kt     = jvdw(kglob)

         do_scale4 = .false.
         vscale4   = 0

         if (vscale.lt.0) then 
            vscale4 = -vscale
            vscale = 1
         end if
c
c     compute the energy contribution for this interaction
c
         xr   = x(iglob) - x(kglob)
         yr   = y(iglob) - y(kglob)
         zr   = z(iglob) - z(kglob)
         call image_inl(xr,yr,zr)
c
c     decide whether to compute the current interaction
c     and check for an interaction distance less than the cutoff
c
         rik2   = xr**2 + yr**2 + zr**2
         if (rik2.lt.loff2.or.rik2.gt.off2) cycle

         mutik  = mutInt(iglob) + mutInt(kglob)
         if (vcouple.and.mutik.eq.two1) mutik=one1 ! Annihilation

         if (use_group)
     &      call groups2_inl(fgrp,iglob,kglob,ngrp,grplist,wgrp)
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

         call duo_lj(rik2,xr,yr,zr,rv2,eps2*vscale,cut2
     &              ,rinv,off,shortheal,use_group,fgrp,mutik
     &              ,sck,sct,scs,scalpha,galpha,use_lambdadyn
     &              ,lambda,vlambda,glamb,lambdavt,lambdavt1
     &              ,delambdav_,e,ded,ver1,fea)

         if (.not.do_scale4) then
         e     = -e
         end if

         ev = ev + tp2enr(e)

         if (vscale.eq.1.0) nev_ = nev_ - 1

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
!$acc parallel loop gang vector_length(32) async(dir_queue)
!$acc&         present(vdwglobnl,ivdw,grplist,wgrp
!$acc&   ,loc,ired,kred,x,y,z,jvdw,lst,nlst,mutInt,radmin,epsilon
!$acc&   ,ev,nev_)
!$acc&         private(ai12) reduction(+:nev_,ev)
      do ii = 1, nvdwlocnl
         iivdw = vdwglobnl(ii)
         iglob = ivdw(iivdw)
         i     = loc(iglob)
         !iv    = ired(iglob)
         it    = jvdw(iglob)
         xi    = x(i)
         yi    = y(i)
         zi    = z(i)
         muti  = mutInt(iglob)
c        usei  = (use(iglob) .or. use(iv))
         if (skipvdw12) then
            in12 = n12(iglob)
!$acc loop vector
            do j = 1,in12
               ai12(j) = i12(j,iglob)
            end do
         end if
c
c     decide whether to compute the current interaction
c
!$acc loop vector private(ded) reduction(+:nev_,ev)
         do kk = 1, nlst(ii)
            kglob = lst(kk,ii)
            kbis  = loc(kglob)
            !kv    = ired(kglob)
            kt    = jvdw(kglob)

            if (skipvdw12) then
               ik12 = .false.
!$acc loop seq
               do j = 1, in12
                  if (ai12(j).eq.kglob) ik12=.true.
               end do
               if (ik12) cycle
            end if
            xr    = xi - x(kbis)
            yr    = yi - y(kbis)
            zr    = zi - z(kbis)
            if (use_bounds) call image_inl (xr,yr,zr)
            rik2  = xr*xr + yr*yr + zr*zr
c
c     check for an interaction distance less than the cutoff
c
            if (loff2.le.rik2.and.rik2.le.off2) then
               rv   = radmin (kt,it)
               eps  = epsilon(kt,it)
               mutik= muti + mutInt( kglob )
               if (vcouple.and.mutik.eq.two1) mutik=one1 ! Annihilation

               if (use_group)
     &            call groups2_inl(fgrp,iglob,kglob,ngrp,grplist,wgrp)

               !compute the energy contribution for this interaction
               call duo_lj(rik2,xr,yr,zr,rv,eps,cut2
     &                    ,rinv,off,shortheal,use_group,fgrp,mutik
     &                    ,sck,sct,scs,scalpha,galpha,use_lambdadyn
     &                    ,lambda,vlambda,glamb,lambdavt,lambdavt1
     &                    ,delambdav_,e,ded,ver,fea)
c
c     increment the total van der Waals energy and derivatives
c
               ev  = ev  + tp2enr(e)
               nev_= nev_ + 1.0
            end if
         end do
      end do
!$acc serial async present(nev,nev_)
      nev = int(nev_)
!$acc end serial
      end

#ifdef _CUDA
      subroutine elj3c_cu
      use action    ,only: nev,nev_
      use atmlst    ,only: vdwglobnl,vdwglob
      use atoms     ,only: x,y,z,n
      use couple    ,only: i12
      use cudafor   ,only: dim3
      use cutoff    ,only: shortheal,vdwshortcut
      use deriv     ,only: dev,devx,de_ws2,d_x,d_y,d_z,dr_stride
      use domdec    ,only: loc,rank,nbloc,nproc
     &              ,xbegproc,xendproc,ybegproc,yendproc,zbegproc
     &              ,zendproc,glob
      use eljcu     ,only: lj3_kcu
#ifdef USE_DETERMINISTIC_REDUCTION
      use elj1gpu_inl,only: enr2en,mdr2md
#endif
      use energi    ,only: ev=>ev_r,calc_e
      use group     ,only: use_group,grplist,wgrp
      use inform    ,only: deb_Path,minmaxone
      use kvdws     ,only: radv,epsv
      use kvdwpr    ,only: vdwpr_l
      use mutant    ,only: mutInt,lambda,vlambda,bvlambda
     &              ,scexp,scalpha,sck,sct,scs,nmut
      use neigh
      use potent    ,only: use_lambdadyn,use_vdwshort,use_vdwlong
      use tinheader ,only: ti_p
      use tinMemory ,only: prmem_request
      use timestat  ,only: timer_enter,timer_exit,timer_elj3
      use shunt     ,only: c0,c1,c2,c3,c4,c5,off2,off,cut2,cut
      use utilcu    ,only: check_launch_kernel,VDW_BLOCK_DIM
     &              ,TRP_BLOCK_DIM,transpose_z3fl
      use utilgpu   ,only: def_queue,dir_queue,rec_queue,dir_stream
     &              ,rec_stream,def_stream,WARP_SIZE,RED_BUFF_SIZE,inf
     &              ,ered_buff=>ered_buf1,vred_buff,nred_buff,lam_buff
     &              ,reduce_buffer,reduce_energy_action
      use vdw       ,only: ired,kred,jvdw,ivdw,radmin_c
     &              ,epsilon_c,radmin,radmin4,epsilon,epsilon4
     &              ,nvdwbloc,nvdwlocnl
     &              ,nvdwlocnlb,nvdwclass
     &              ,nvdwlocnlb_pair,nvdwlocnlb2_pair
      use vdwpot    ,only: vcorrect_ik,vcorrect_scale,n_vscale,dhal,ghal
     &              ,radrule_i,epsrule_i
      use vdw_locArray
      use virial    ,only: vir,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz
     &              ,use_virial
      implicit none
      integer i,k,sz1
      integer iglob,iivdw,iv,it,hal_Gs
      real(t_p)  galpha,glamb,lambdavt,lambdavt1
      real(t_p)  xbeg,xend,ybeg,yend,zbeg,zend
      real(t_p)  rinv,loff2
      character*11 mode
c
      if(deb_Path) write (*,*) 'elj3c_cu'

      call load_nbl2mod(vdw_nbl)

      call spatialOrder_pos

      xbeg = xbegproc(rank+1)
      xend = xendproc(rank+1)
      ybeg = ybegproc(rank+1)
      yend = yendproc(rank+1)
      zbeg = zbegproc(rank+1)
      zend = zendproc(rank+1)

      def_queue  = dir_queue
      def_stream = dir_stream

      !set the coefficients for the switching function
      mode   = 'VDW'
      hal_Gs = min(max((nb2p+nbap)/16,1),2**14)
      sz1    = size(vcorrect_ik,1)
      call   switch (mode)
      rinv   = 1.0/(cut-off)
      loff2  = merge((vdwshortcut-shortheal)**2,0.0_ti_p,use_vdwlong)

      if (nmut.ne.0) then
         galpha    = scalpha/(2d0**(sck/6d0))
         glamb     = 1.0-vlambda
         lambdavt  = vlambda ** sct
         lambdavt1 = sct*(vlambda**(sct-1.0))
      end if
c
c     Call Lennard-Jones kernel in CUDA using C2 nblist
c
!$acc host_data use_device(so_x,so_y,so_z,sgl_id,slc_id
!$acc&    ,i12,mutInt,b2pl,bapl,abpl,cellv_jvdw,epsilon_c
!$acc&    ,radmin_c,radv,epsv,grplist,wgrp
!$acc&    ,d_x,d_y,d_z,ered_buff,vred_buff,nred_buff,lam_buff
!$acc&    ,x,y,z,vcorrect_ik,vcorrect_scale,loc,jvdw
!$acc&    ,radmin,epsilon,radmin4,epsilon4,dev
!$acc&    )

      call lj3_kcu<<<hal_Gs,VDW_BLOCK_DIM,0,def_stream>>>
     &            (so_x,so_y,so_z,sgl_id,slc_id
     &            ,b2pl,bapl,abpl,cellv_jvdw,i12,mutInt
     &            ,epsilon_c,radmin_c,radv,epsv
     &            ,d_x,d_y,d_z,ered_buff,vred_buff,nred_buff
     &            ,nb2p,nbap,n,nbloc,nvdwlocnl,nab,nvdwclass
     &            ,radrule_i,epsrule_i
     &            ,cut2,cut,off2,off,loff2,shortheal,rinv
     &            ,use_group,grplist,wgrp
     &            ! lambdaDyn
     &            ,sck,sct,scs,scalpha,galpha,use_lambdadyn,lambda
     &            ,vlambda,glamb,lambdavt,lambdavt1,lam_buff
     &            ! Scaling factor
     &            ,x,y,z,vcorrect_ik,vcorrect_scale,sz1,n_vscale
     &            ,loc,jvdw
     &            ,radmin,epsilon,radmin4,epsilon4,dev
     &            ! Box data
     &            ,xbeg,xend,ybeg,yend,zbeg,zend
     &            )
      call check_launch_kernel(" lj3_kcu ")

!$acc end host_data

      call reduce_energy_action(ev,nev,ered_buff,def_queue)

      end subroutine
#endif
