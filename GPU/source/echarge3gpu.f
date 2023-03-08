c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine echarge3  --  charge-charge energy & analysis  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "echarge3" calculates the charge-charge interaction energy
c     and partitions the energy among the atoms
c
c
#include "tinker_macro.h"
      module echarge3gpu_inl
        use tinTypes , only: real3,real3_red,rpole_elt
        implicit none
#include "atomicOp.h.f"
        include "erfcore_data.f.inc"
        contains
#include "convert.f.inc"
#include "image.f.inc"
#include "atomicOp.inc.f"
#if defined(SINGLE) | defined(MIXED)
        include "erfcscore.f.inc"
#else
        include "erfcdcore.f.inc"
#endif
#include "groups.inc.f"
#include "switch_respa.f.inc"
#include "pair_charge.f.inc"
      end module

      subroutine echarge3gpu
      implicit none  
c
c     choose the method for summing over pairwise interactions
c
      call echarge3cgpu
c
      return
      end
c
c
c     ###################################################################
c     ##                                                               ##
c     ##  subroutine echarge3cgpu  --  Ewald charge analysis via list  ##
c     ##                                                               ##
c     ###################################################################
c
c
c
c     "echarge3cgpu" calculates the charge-charge interaction energy
c     and partitions the energy among the atoms using a particle
c     mesh Ewald summation
c
c
      subroutine echarge3cgpu
      use action
      use analyz
      use atmlst
      use atoms
      use bound
      use boxes
      use charge
      use chgpot
      use couple
      use domdec
      use echarge3gpu_inl
      use energi
      use erf_mod
      use ewald
      use group
      use inform
      use inter
      use interfaces,only: ecreal3d_p
      use iounit
      use math
      use molcul
      use mpi
      use neigh
      use potent
      use shunt
      use timestat
      use usage
      use utilgpu
      use utilvec

      implicit none
      integer i,k,ii,kk,ksave
      integer iglob,iichg,inl
      integer nnchg,nnchg1,nnchg2
      integer nn12,nn13,nn14,ntot
      real(t_p) f,fi
      real(en_p) e,efull
      real(t_p) fs
      real(t_p) xi,yi,zi
      real(t_p) xd,yd,zd
      real(t_p) rew,erf1,fik
      real(t_p) r,r2
      real(t_p) term,term1,term2
      real(t_p) one,half
      real(t_p) time0,time1
      character*10 mode
      parameter(half=0.5_ti_p)
      parameter(one=1.0_ti_p)

      if(deb_Path) write (*,*) 'echarge3cgpu'
      call timer_enter(timer_echarge)
c
c     zero out the Ewald summation energy and partitioning
c
      nec = 0
      ec  = 0.0_re_p
c     aec = 0.0_ti_p

      if (nion .eq. 0)  return
c
c     set conversion factor, cutoff and switching coefficients
c     update the number of gangs required for gpu
c
      f    = electric / dielec
      mode = 'EWALD'
      call switch (mode)
      call update_gang(nionlocnl)
c
c     compute the Ewald self-energy term over all the atoms
c
      fs   = -f * aewald / sqrtpi

!$acc data present(ec,ec_r,ecrec,nec,nec_) async
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.lt.ndir))
     $   then

        if (use_cself) then

          call timer_enter(timer_other)
          f  = electric / dielec
          fs = -f * aewald / sqrtpi
!$acc   parallel loop async(dir_queue)
!$acc&           present(chgglob,pchg)
          do ii = 1, nionloc
             iichg = chgglob(ii)
             ec    = ec + fs * pchg(iichg)**2
          end do
c
c       compute the cell dipole boundary correction term
c
          if (boundary .eq. 'VACUUM') then
             term = (2.0_ti_p/3.0_ti_p) * f * (pi/volbox)
             xd   = 0.0_ti_p
             yd   = 0.0_ti_p
             zd   = 0.0_ti_p
!$acc   parallel loop async
!$acc&           present(chgglob,iion,pchg,x,y,z)
             do k = 1,nionloc
                iichg = chgglob(ii)
                iglob = iion(iichg)
                xd    = xd + pchg(iichg)*x(iglob)
                yd    = yd + pchg(iichg)*y(iglob)
                zd    = zd + pchg(iichg)*z(iglob)
             enddo
!$acc serial async
             ec   = ec + term * ( xd**2 + yd**2 + zd**2 )
             nec_ = nec_ + 1
!$acc end serial
          end if
          call timer_exit( timer_other,quiet_timers )

        end if
c
c     compute the real space portion of the Ewald summation
c
        call timer_enter(timer_real)
        call ecreal3d_p
        call timer_exit( timer_real )

      end if
c
c     compute the reciprocal space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.gt.ndir-1))
     $   then
         call timer_enter( timer_rec )
         call ecrecipgpu
         call timer_exit ( timer_rec,quiet_timers )
      end if

!$acc serial async
      !print*, 'echarge3gpu',ec, ec_r, ecrec
       ec = ec + enr2en( ec_r ) + ecrec
      nec = int(nec_) + nionloc + nec
!$acc end serial

!$acc end data

      call timer_exit(timer_echarge)
      end

      subroutine ecreal3dgpu
      use action
      use atmlst
      use atoms
      use bound
      use boxes
      use charge
      use chgpot
      use couple
      use cutoff
      use deriv
      use domdec
      use echarge3gpu_inl
      use energi ,only: ec=>ec_r
      use ewald
      use group
      use iounit
      use inform
      use inter
      use math
      use molcul
      use mutant
      use neigh
      use potent
      use shunt
      use timestat
      use usage
      use tinheader ,only: ti_p
      use tinTypes
      use utilgpu
      use virial
      use mpi
      implicit none
      integer    i,j,k,iichg,iglob,ii,kkk,kglob,kkchg,ver,ver1,fea
      integer   ,pointer:: lst(:,:),nlst(:)
      integer(1) mutik,muti
      real(t_p)  e,delambdae_
      real(t_p)  f,fi,fi_,fik,fik_,r,r2,rew,rb,rb2
      real(t_p)  xi,yi,zi,xr,yr,zr
      real(t_p)  loff2,scut,scale_f,scale,fgrp
      type(real3)   ded
      character*11  mode
      parameter( scale_f= 1.0
     &         , ver    = __use_ene__+__use_act__
     &         , ver1   =       ver  +__use_sca__
     &         )

      if (deb_Path) write(*,'(2x,a)') 'ecreal3dgpu'

      fea = __use_mpi__
      if (use_lambdadyn) fea = fea + __use_lambdadyn__

      !set conversion factor, cutoff and switching coefficients
      f    = electric / dielec
      if (use_cshortreal) then
         mode  = 'SHORTEWALD'
         call switch (mode)
         loff2 = 0
         scut  = off
         lst   => shortelst
         nlst  => nshortelst
         fea   = fea + __use_shortRange__
      else if (use_clong) then
         mode  = 'EWALD'
         call switch (mode)
         loff2 = (chgshortcut-shortheal)**2
         scut  = chgshortcut
         lst   => elst
         nlst  => nelst
         fea   = fea + __use_longRange__
      else
         mode  = 'EWALD'
         call switch (mode)
         loff2 = 0
         scut  = chgshortcut
         lst   => elst
         nlst  => nelst
      end if
c
c     Scaling interaction correction subroutines
c
!$acc parallel loop gang vector async(dir_queue)
!$acc&         present(ccorrect_ik,ccorrect_scale,loc,x,y,z,mutInt,
!$acc&   grplist,wgrp,pchg,pchg_orig,ec,nec_) reduction(+:ec,nec_)
!$acc&         private(ded)
      do ii = 1, n_cscale
         iglob = ccorrect_ik(ii,1)
         kglob = ccorrect_ik(ii,2)
         scale = ccorrect_scale(2*ii+1)
         if (use_lambdadyn) then
         fi    = pchg(chglist(iglob))
         fik   = f*fi*pchg(chglist(kglob))
         fi    = pchg_orig(chglist(iglob))
         fik_  = f*fi*pchg_orig(chglist(kglob))
         mutik = mutInt(iglob)+mutInt(kglob)
         else
         fik   = f*ccorrect_scale(2*ii+2)
         end if
         xi    = x(iglob)
         yi    = y(iglob)
         zi    = z(iglob)
c
c     compute the energy contribution for this interaction
c
         xr    = xi - x(kglob)
         yr    = yi - y(kglob)
         zr    = zi - z(kglob)
         call image_inl (xr,yr,zr)
         r2 = xr*xr + yr*yr + zr*zr

         if (use_group) then
            call groups2_inl(fgrp,iglob,kglob,ngrp,grplist,wgrp)
            scale = scale*fgrp
         end if
c
c     find energy for interactions within real space cutoff
c
         if (r2.gt.loff2 .and. r2.le.off2) then
            call charge_couple(r2,xr,yr,zr,ebuffer,fik_,fik,aewald
     &                        ,scale,mutik,use_lambdadyn,shortheal
     &                        ,scut,elambda,delambdae_,e,ded,ver1,fea)
 
           !increment the overall energy and derivative expressions
            ec = ec + tp2enr(e)
            if (scale.eq.-1.0) nec_ = nec_-1
         end if
      end do
c
c     compute the real space Ewald energy and first derivatives
c
!$acc parallel loop gang vector_length(32) async(dir_queue)
!$acc&         present(chgglobnl,iion,loc,x,y,z,pchg,nlst,lst,mutInt,
!$acc&    ec,nec_) reduction(+:ec,nec_)
      do ii = 1, nionlocnl
         iichg = chgglobnl(ii)
         iglob = iion(iichg)
         if (use_lambdadyn) muti = mutInt(iglob)
         xi    = x(iglob)
         yi    = y(iglob)
         zi    = z(iglob)
         fi    = f * pchg(iichg)
!$acc loop vector private(ded) reduction(+:ec,nec_)
         do kkk = 1, nlst(ii)
            kkchg = lst(kkk,ii)
            if (kkchg.eq.0) cycle
            kglob = iion(kkchg)
            if (use_lambdadyn) mutik = muti+mutInt(kglob)
c
c     compute the energy contribution for this interaction
c
            xr = xi - x(kglob)
            yr = yi - y(kglob)
            zr = zi - z(kglob)
            call image_inl (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
c
c     find energy for interactions within real space cutoff
c
            if (r2.gt.loff2 .and. r2.le.off2) then
               if (use_group) then
                  call groups2_inl(fgrp,iglob,kglob,ngrp,grplist,wgrp)
                  fgrp = fgrp-1.0
               else
                  fgrp = 0.0
               end if
               fik  = fi*pchg(kkchg)
               call charge_couple(r2,xr,yr,zr,ebuffer,fik_,fik,aewald
     &                    ,fgrp-1.0,mutik,use_lambdadyn,shortheal,scut
     &                    ,elambda,delambdae_,e,ded,ver,fea)
 
              !increment the overall energy and derivative expressions
               ec   = ec + tp2enr(e)
               nec_ = nec_ + 1
            end if
         end do
      end do
      end

#ifdef _CUDA
      subroutine ecreal3d_cu
      use action
      use atmlst
      use atoms
      use bound
      use boxes
      use charge
      use chgpot
      use couple
      use cutoff
      use deriv
      use domdec
      use echargecu
      use energi , only: ec=>ec_r,calc_e
      use ewald
      use group
      use iounit
      use inform
      use inter
      use math
      use molcul
      use mpi
      use mutant
      use neigh  , iion_s=>celle_glob,chg_s=>celle_chg
     &           , loc_s=>celle_loc, ieblst_s=>ieblst, eblst_s=>eblst
     &           , x_s=>celle_x, y_s=>celle_y, z_s=>celle_z
      use potent
      use shunt
      use timestat
      use usage
      use utilcu  ,only: BLOCK_DIM,check_launch_kernel 
      use utilgpu ,only: def_queue,dir_queue,dir_stream,def_stream
     &            , vred_buff,lam_buff,ered_buff=>ered_buf1,nred_buff
     &            , reduce_energy_virial, reduce_energy_action
      use virial
      implicit none
      integer i,szcik
      real(t_p) f
      real(t_p) loff2,scut,scale_f
      real(t_p) p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
      logical,save:: first_in=.true.
      integer,save:: gS
      character*10 mode
      logical,parameter::dyn_gS=.true.

      if (deb_Path) write(*,'(2X,A)') 'ecreal3d_cu'

      if (first_in) then
         call cudaMaxGridSize("ecreal3_kcu",gS)
         first_in=.false.
      end if
      if (dyn_gS) gs = max(nionlocnlb2_pair/8,1)

      p_xbeg = xbegproc(rank+1)
      p_xend = xendproc(rank+1)
      p_ybeg = ybegproc(rank+1)
      p_yend = yendproc(rank+1)
      p_zbeg = zbegproc(rank+1)
      p_zend = zendproc(rank+1)

      !set conversion factor, cutoff and switching coefficients
      f    = electric / dielec
      if (use_cshortreal) then
         mode  = 'SHORTEWALD'
         call switch (mode)
         loff2 = 0
         scut  = off
      else if (use_clong) then
         mode  = 'EWALD'
         call switch (mode)
         loff2 = (chgshortcut-shortheal)**2
         scut  = chgshortcut
      else
         mode  = 'EWALD'
         call switch (mode)
         loff2 = 0
         scut  = chgshortcut
      end if

      def_stream = dir_stream
      def_queue  = dir_queue
      szcik      = size(ccorrect_ik,1)

!$acc host_data use_device(iion_s,chg_s,loc_s,ieblst_s,eblst_s
!$acc&    ,x_s,y_s,z_s,pchg,pchg_orig,mutInt,grplist,wgrp
!$acc&    ,dec,ered_buff,vred_buff,lam_buff,nred_buff
!$acc&    ,ccorrect_ik,ccorrect_scale,chglist,loc,x,y,z
!$acc&    )

      call ecreal3_kcu<<<gS,BLOCK_DIM,0,def_stream>>>
     &     ( iion_s,chg_s,loc_s,ieblst_s
     &     , eblst_s(2*nionlocnlb_pair+1)
     &     , nionlocnlb,nionlocnlb2_pair,nionbloc,n
     &     , x_s,y_s,z_s,pchg,pchg_orig,mutInt
     &     , off2,loff2,scut,shortheal,f,aewald,ebuffer
     &     , elambda,use_lambdadyn
     &     , dec,ered_buff,vred_buff,lam_buff,nred_buff
     &     , use_group, grplist, wgrp
     &     , ccorrect_ik,ccorrect_scale,chglist,x,y,z,loc,n_cscale,szcik
     &     , p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
     &     )
      call check_launch_kernel(" ecreal1d_core_cu")

!$acc end host_data

      call reduce_energy_action(ec,nec,ered_buff,def_queue)

      end subroutine
#endif
c
c     GPU version of ecrecip
c
      subroutine ecrecipgpu
      use atoms
      use atmlst
      use bound
      use boxes
      use charge
      use chgpot
      use domdec
      use energi
      use ewald
      use fft
      use inform    ,only: deb_Path
      use interfaces,only: grid_pchg_site_p
      use math
      use mpi
      use pme
      use pme1
      use potent
      use timestat
      use tinheader ,only: ti_p
      use utils     ,only: set_to_zero1
      use utilgpu   ,only: rec_queue
      implicit none
      integer i,j,k
      integer iichg,iglob
      integer k1,k2,k3
      integer m1,m2,m3
      integer nf1,nf2,nf3
      integer nff,npoint
      integer status(MPI_STATUS_SIZE),tag,ierr,proc
      integer nprocloc,rankloc,commloc
      real(t_p) e,f,denom
      real(t_p) term,expterm
      real(t_p) pterm,volterm
      real(t_p) hsq,struc2
      real(t_p) h1,h2,h3
      real(t_p) r1,r2,r3
      integer,dimension(nproc*nproc)::req,reqbcast
      logical,save:: f_in=.true.
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
c
      if (deb_Path)
     &   write(*,'(2x,a)') 'ecrecipgpu'

      if (use_pmecore) then
        nprocloc = nrec
        rankloc = rank_bis
        commloc = comm_rec
      else
        nprocloc = nproc
        rankloc = rank
        commloc = MPI_COMM_WORLD
      end if
c
c     dynamic allocation of local arrays
c
      call mallocMpiGrid

!$acc data present(qgridin_2d,qgridout_2d,qfac_2d,
!$acc&  kstart2,kend2,jstart2,jend2,istart2,iend2,
!$acc&  bsmod1,bsmod2,bsmod3,use_bounds,
!$acc&  octahedron,ecrec)
c
c     MPI : Begin reception
c
      do i = 1, nrec_recep
         tag = nprocloc*rankloc + prec_recep(i) + 1
!$acc host_data use_device(qgridmpi)
         call MPI_IRECV(qgridmpi(1,1,1,1,i),2*n1mpimax*n2mpimax*
     $        n3mpimax,MPI_TPREC,prec_recep(i),tag,commloc,req(tag),
     $        ierr)
!$acc end host_data
      end do
c
      call bspline_fill_sitegpu(1)
c
      call timer_enter( timer_grid1 )
      call set_to_zero1(qgridin_2d(1,1,1,1,1),
     &     2*n1mpimax*n2mpimax*n3mpimax*(nrec_send+1),rec_queue)
c
      call grid_pchg_site_p
      call timer_exit ( timer_grid1,quiet_timers )

      ! FFt Grid Communication
      call commGridFront( qgridin_2d,r_comm )
      call commGridFront( qgridin_2d,r_wait )

#ifdef _OPENACC
      ! perform the 3-D FFT forward transformation
      call cufft2d_frontmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
     $     n3mpimax)
#else
      call fft2d_frontmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
     $     n3mpimax)
#endif
c
c     use scalar sum to get the reciprocal space energy
c
      call timer_enter( timer_scalar )
      f   = 0.5_ti_p * electric / dielec
      pterm = (pi/aewald)**2
      volterm = pi * volbox
      nff = nfft1 * nfft2
      nf1 = (nfft1+1) / 2
      nf2 = (nfft2+1) / 2
      nf3 = (nfft3+1) / 2
!$acc serial async(rec_queue)
      if ((istart2(rankloc+1).eq.1).and.(jstart2(rankloc+1).eq.1).and.
     $   (kstart2(rankloc+1).eq.1)) then
            qfac_2d(1,1,1) = 0.0_ti_p
      end if
!$acc end serial

!$acc parallel loop collapse(3) async(rec_queue)
      do k3 = kstart2(rankloc+1),kend2(rankloc+1)
         do k2 = jstart2(rankloc+1),jend2(rankloc+1)
            do k1 = istart2(rankloc+1),iend2(rankloc+1)
              m1 = k1 - 1
              m2 = k2 - 1
              m3 = k3 - 1
              if (k1 .gt. nf1)  m1 = m1 - nfft1
              if (k2 .gt. nf2)  m2 = m2 - nfft2
              if (k3 .gt. nf3)  m3 = m3 - nfft3
              if ((m1.eq.0).and.(m2.eq.0).and.(m3.eq.0)) goto 10
              r1 = real(m1,t_p)
              r2 = real(m2,t_p)
              r3 = real(m3,t_p)
              h1 = recip(1,1)*r1 + recip(1,2)*r2 + recip(1,3)*r3
              h2 = recip(2,1)*r1 + recip(2,2)*r2 + recip(2,3)*r3
              h3 = recip(3,1)*r1 + recip(3,2)*r2 + recip(3,3)*r3
              hsq = h1*h1 + h2*h2 + h3*h3
              term = -pterm * hsq
              expterm = 0.0_ti_p
              if (term .gt. -50.0_ti_p) then
                 denom = volterm*hsq*bsmod1(k1)*bsmod2(k2)*bsmod3(k3)
                 expterm = exp(term) / denom
                 if (.not. use_bounds) then
                    expterm = expterm*(1.0_ti_p-cos(pi*xbox*sqrt(hsq)))
                 else if (octahedron) then
                    if (mod(m1+m2+m3,2) .ne. 0)  expterm = 0.0_ti_p
                 end if
                 struc2 = qgridout_2d(1,k1-istart2(rankloc+1)+1,
     $  k2-jstart2(rankloc+1)+1,k3-kstart2(rankloc+1)+1)**2 + 
     $  qgridout_2d(2,k1-istart2(rankloc+1)+1,k2-jstart2(rankloc+1)+1,
     $  k3-kstart2(rankloc+1)+1)**2
                 e = f * expterm * struc2
                 ecrec = ecrec + e
              end if
 10           continue
            end do
         end do
      end do
c
c     account for zeroth grid point for nonperiodic system
c
!$acc serial async(rec_queue)
      if ((istart2(rankloc+1).eq.1).and.(jstart2(rankloc+1).eq.1).and.
     $    (kstart2(rankloc+1).eq.1)) then
         if (.not. use_bounds) then
            expterm = 0.5_ti_p * pi / xbox
            struc2 = qgridout_2d(1,1,1,1)**2 + qgridout_2d(2,1,1,1)**2
            e = f * expterm * struc2
            ecrec = ecrec + e
         end if
      end if
!$acc end serial
      call timer_exit ( timer_scalar,quiet_timers )

!$acc end data
      if(f_in) f_in = .false.

      end
