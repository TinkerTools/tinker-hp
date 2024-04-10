c
c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine elj1  --  Lennard-Jones energy & derivatives  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "elj1" calculates the Lennard-Jones 6-12 van der Waals energy
c     and its first derivatives with respect to Cartesian coordinates
c
c
#include "tinker_macro.h"
      module elj1gpu_inl
        integer(1) one1,two1
#include "atomicOp.h.f"
        parameter( one1=1, two1=2 )
        contains
#include "convert.f.inc"
#include "image.f.inc"
#include "atomicOp.inc.f"
#include "groups.inc.f"
#include "pair_lj.inc.f"
      end module

      subroutine elj1gpu
      use energi
      use deriv     ,only: delambdav,delambdavsave
      use interfaces,only: elj1c_p
      use potent
      use virial
      use vdwpot
      use utilgpu
      implicit none
      real(r_p) elrc,vlrc
      character*11 mode
c
c     choose the method for summing over pairwise interactions
c
      call elj1c_p
c
c     apply long range van der Waals correction if desired
c
      if (use_vcorr) then
         def_queue = dir_queue
         mode = 'VDW'
!$acc data  create(elrc,vlrc) async(def_queue)
!$acc&      present(ev,vir)
         call evcorr1gpu (mode,elrc,vlrc)
!$acc serial async(def_queue)
         ev = ev + elrc
         vir(1,1) = vir(1,1) + vlrc
         vir(2,2) = vir(2,2) + vlrc
         vir(3,3) = vir(3,3) + vlrc
!$acc end serial
!$acc wait
!$acc end data
      end if

      if (use_lambdadyn) then
c
c     save delambdav for short range computation
c
         if (use_vdwshort) then
!$acc serial async(def_queue) present(delambdav,delambdavsave)
            delambdavsave = delambdav
!$acc end serial
         end if
!$acc update host(delambdav,delambdavsave) async
      end if
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine elj1c_ac  --  Lennard-Jones vdw derivs via list  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "elj1c" calculates the Lennard-Jones 12-6 van der Waals energy
c     and its first derivatives using a pairwise neighbor list
c
c
      subroutine elj1c_ac
      use atmlst
      use atoms
      use bound
      use couple
      use cutoff    ,only: vdwshortcut,shortheal
      use deriv     ,only: dev,delambdav
      use domdec
      use elj1gpu_inl
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

      parameter( ver=__use_grd__+__use_ene__+__use_vir__
     &         ,ver1=        ver+__use_sca__)
c
 1000 format(' Warning, system moved too much since last neighbor list'
     $   ' update, try lowering nlupdate VDW')
      if (deb_Path) write(*,*) 'elj1c_ac'

      ! set elj1cgpu Configuration
      fea = __use_mpi__
      if (nmut.ne.0)     fea = fea + __use_softcore__
      if (use_lambdadyn) fea = fea + __use_lambdadyn__
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
         loff2 = (vdwshortcut-shortheal)**2
         !set the coefficients for the switching function
         mode  = 'VDW'
         call  switch (mode)
      else
         lst   => vlst
         nlst  => nvlst
         loff2 = 0.0
         !set the coefficients for the switching function
         mode  = 'VDW'
         call  switch (mode)
      end if
      rinv   = 1.0/(cut-off)

      ! Set constant value for lambda dynamic
      if (nmut.ne.0) then
         galpha    = scalpha/(2d0**(sck/6d0))
         glamb     = 1.0-vlambda
         lambdavt  = vlambda ** sct
         lambdavt1 = sct*(vlambda**(sct-1.0))
      end if
c
c     Scaling interaction correction for Lennard-Jones
c
!$acc parallel loop async(dir_queue) gang vector private(ded)
!$acc&         present(dev,delambdav,ev,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz
!$acc&                ,g_vzz)
!$acc&         present(loc,ired,x,y,z,loc,jvdw,vir,radmin,mutInt
!$acc&                ,grplist,wgrp,radmin4,epsilon,epsilon4
!$acc&                ,vcorrect_ik,vcorrect_scale)
!$acc&         reduction(+:delambdav,ev,g_vxx,g_vxy,g_vxz
!$acc&                ,g_vyy,g_vyz,g_vzz)
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
         if (vcouple.eq.1.and.mutik.eq.two1) mutik=one1 ! Annihilation

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
         ded%x = -ded%x; ded%y = -ded%y; ded%z = -ded%z;
         delambdav_ = -delambdav_
         end if

         ev = ev + tp2enr(e)
         !if(rank.eq.0.and.mod(ii,1).eq.0) print*,iglob,kglob,vscale,e

         if (use_lambdadyn) delambdav = delambdav + delambdav_

         call atomic_sub( dev(1,kbis),ded%x )
         call atomic_sub( dev(2,kbis),ded%y )
         call atomic_sub( dev(3,kbis),ded%z )

         call atomic_add( dev(1,i),ded%x )
         call atomic_add( dev(2,i),ded%y )
         call atomic_add( dev(3,i),ded%z )
c
c     increment the total van der Waals energy 
c
         g_vxx = g_vxx + xr * ded%x
         g_vxy = g_vxy + yr * ded%x
         g_vxz = g_vxz + zr * ded%x
         g_vyy = g_vyy + yr * ded%y
         g_vyz = g_vyz + zr * ded%y
         g_vzz = g_vzz + zr * ded%z

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
!$acc&   ,ev,dev,delambdav,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
!$acc&         private(ai12) reduction(+:delambdav,ev
!$acc&   ,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
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
!$acc loop vector private(ded) reduction(+:delambdav,ev
!$acc&           ,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
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
               if (vcouple.eq.1.and.mutik.eq.two1) mutik=one1 ! Annihilation

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

               if (use_lambdadyn) delambdav = delambdav + delambdav_

               call atomic_sub( dev(1,kbis),ded%x )
               call atomic_sub( dev(2,kbis),ded%y )
               call atomic_sub( dev(3,kbis),ded%z )

               call atomic_add( dev(1,i),ded%x )
               call atomic_add( dev(2,i),ded%y )
               call atomic_add( dev(3,i),ded%z )
c
c     increment the internal virial tensor components
c
               g_vxx = g_vxx + xr * ded%x
               g_vxy = g_vxy + yr * ded%x
               g_vxz = g_vxz + zr * ded%x
               g_vyy = g_vyy + yr * ded%y
               g_vyz = g_vyz + zr * ded%y
               g_vzz = g_vzz + zr * ded%z
            end if
         end do
      end do
      end

#ifdef _CUDA
c=============================================================
c            CUDA Routine for Lennard-Jones 
c=============================================================
      subroutine elj1c_cu
      use atmlst    ,only: vdwglobnl,vdwglob
      use atoms     ,only: x,y,z,n
      use couple    ,only: i12
      use cudafor   ,only: dim3
      use cutoff    ,only: shortheal,vdwshortcut
      use deriv     ,only: dev,devx,de_ws2,d_x,d_y,d_z,delambdav
     &              ,dr_stride
      use domdec    ,only: loc,rank,nbloc,nproc
     &              ,xbegproc,xendproc,ybegproc,yendproc,zbegproc
     &              ,zendproc,glob
      use eljcu     ,only: lj1_kcu
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
     &              ,reduce_energy_virial,reduce_buffer,transpose_az3
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
      if(deb_Path) write (*,*) 'elj1c_cu'
      call timer_enter(timer_elj3)

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
      hal_Gs =merge(nb2p+nbap,min((nb2p+nbap)/8,2**14),nb2p+nbap.lt.257)
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

      call lj1_kcu<<<hal_Gs,VDW_BLOCK_DIM,0,def_stream>>>
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
      call check_launch_kernel(" lj1_kcu ")

!$acc end host_data

      if (calc_e.or.use_virial) then
      call reduce_energy_virial(ev,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz
     &                         ,ered_buff,def_queue)
      end if
      if (use_lambdadyn)
     &   call reduce_buffer(lam_buff,RED_BUFF_SIZE,delambdav,def_queue)

      call transpose_az3(d_x,devx,slc_id,na,dr_stride,def_queue)

      call timer_exit(timer_elj3)
      end subroutine
#endif
