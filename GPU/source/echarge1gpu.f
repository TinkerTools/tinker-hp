c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine echarge1  --  charge-charge energy & derivs  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "echarge1" calculates the charge-charge interaction energy
c     and first derivatives with respect to Cartesian coordinates
c
c
#include "tinker_precision.h"
#include "tinker_types.h"
      module echarge1gpu_inl
        use utilgpu , only: real3,real3_red,rpole_elt
        implicit none
        include "erfcore_data.f.inc"
        contains
#include "convert.f.inc"
#include "image.f.inc"
#if defined(SINGLE) | defined(MIXED)
        include "erfcscore.f.inc"
#else
        include "erfcdcore.f.inc"
#endif
#include "switch_respa.f.inc"
#include "pair_charge.f.inc"
      end module

      subroutine echarge1gpu
      implicit none
c
c     choose the method for summing over pairwise interactions
c
      call echarge1cgpu
c
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine echarge1c  --  Ewald charge derivs via list  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "echarge1c" calculates the charge-charge interaction energy
c     and first derivatives with respect to Cartesian coordinates
c     using a particle mesh Ewald summation and a neighbor list
c
c
      subroutine echarge1cgpu
      use atmlst
      use atoms
      use bound
      use boxes
      use charge
      use chgpot
      use deriv
      use echarge1gpu_inl
      use energi
      use ewald
      use domdec
      use iounit
      use inform
      use interfaces,only: ecreal1d_p,ecrealshortlong1d_p
      use inter
      use math
      use potent
      use timestat
      use tinheader
      use usage
      use utilgpu
      use virial
      use mpi
      implicit none
      integer ii,i,iglob,iichg
      real(r_p) e
      real(t_p) de,term
      real(t_p) f,fs
      real(t_p) xd,yd,zd
      real(t_p) dedx,dedy,dedz
      real(8) time0,time1

      if (nion.eq.0) return

      if (deb_Path)  write(*,'(1x,a)') 'echarge1cgpu'
      call timer_enter(timer_echarge)
c
c     zero out the Ewald summation energy and derivatives
c
!$acc data create(e) async(rec_queue)
!$acc serial async(dir_queue)
      e     = 0.0
!$acc end serial
c
c     compute the Ewald self-energy term over all the atoms
c
      if (use_cself) then
        call timer_enter(timer_other)
        f  = electric / dielec
        fs = -f * aewald / sqrtpi
!$acc parallel loop async(dir_queue)
!$acc&         present(chgglob,pchg)
        do ii = 1, nionloc
           iichg = chgglob(ii)
           e = e + fs * pchg(iichg)**2
        end do
c
c     compute the cell dipole boundary correction term
c
        if (boundary .eq. 'VACUUM') then
           xd = 0.0
           yd = 0.0
           zd = 0.0
!$acc parallel loop async(dir_queue) default(present)
           do ii = 1, nionloc
             iichg = chgglob(ii)
             iglob = iion(iichg)
             i  = loc(iglob)
             xd = xd + pchg(iichg)*x(iglob)
             yd = yd + pchg(iichg)*y(iglob)
             zd = zd + pchg(iichg)*z(iglob)
           end do
           term = (2.0_ti_p/3.0_ti_p) * f * (pi/volbox)
!$acc serial async(dir_queue)
           e = e + term * (xd*xd+yd*yd+zd*zd)
!$acc end serial
!$acc parallel loop async(dir_queue) default(present)
           do ii = 1, nionloc
              iichg = chgglob(ii)
              iglob = iion(iichg)
              i     = loc(iglob)
              de    = 2.0 * term * pchg(iichg)
              dedx  = de * xd
              dedy  = de * yd
              dedz  = de * zd
              dec(1,i) = dec(1,i) + dedx
              dec(2,i) = dec(2,i) + dedy
              dec(3,i) = dec(3,i) + dedz
           end do
        end if
        call timer_exit( timer_other,quiet_timers )
      end if
c
c     compute the real space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.le.ndir-1))
     $   then
         if (use_creal) then
            call timer_enter(timer_real)
            if (use_cshortreal.or.use_clong) then
               call ecrealshortlong1d_p
            else
               call ecreal1d_p
            end if
            call timer_exit( timer_real )
         end if
      end if
c
c     compute the reciprocal space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.gt.ndir-1))
     $   then
         if (use_crec) then
            call timer_enter(timer_rec)
            call ecrecip1gpu
            call timer_exit(timer_rec )
         end if
      end if
c
#ifdef _OPENACC
      if (dir_queue.ne.rec_queue)
     &   call stream_wait_async(dir_stream,rec_stream)
#endif
!$acc serial present(ecrec,ec,ec_r) async(rec_queue)
      ec = ec + e + ecrec + enr2en( ec_r )
!$acc end serial

!$acc end data
      call timer_exit(timer_echarge)
      end
c
c     "ecreal1dgpu" evaluates the real space portion of the Ewald sum
c     energy and forces due to atomic charge interactions, using a neighbor list
c
      subroutine ecreal1dgpu
      use atmlst
      use atoms
      use bound
      use boxes
      use charge
      use chgpot
      use couple
      use deriv
      use domdec
      use echarge1gpu_inl
      use energi ,only: ec=>ec_r
      use ewald
      use iounit
      use inform
      use inter
      use math
      use molcul
      use neigh
      use potent
      use shunt
      use timestat
      use usage
      use tinTypes
      use utilgpu
      use virial
      use mpi
      implicit none
      integer i,j,k,iichg,iglob
      integer ii,kkk,kglob,kkchg
#ifdef TINKER_DEBUG
      integer ninte(n)
#endif
      real(t_p) e
      real(t_p) f,fi,fik,pchgk
      real(t_p) r,r2,rew
      real(t_p) rb,rb2

      real(t_p) xi,yi,zi
      real(t_p) xr,yr,zr

      real(t_p),parameter:: scale_f=1.0
      type(real3) ded
      type(mdyn3_r) dedc

      character*10 mode

      if (deb_Path) write(*,'(2x,a)') 'ecreal1dgpu'
#ifdef TINKER_DEBUG
      ninte = 0
!$acc enter data copyin(ninte)
#endif
c
c     set conversion factor, cutoff and switching coefficients
c
      f    = electric / dielec
      mode = 'EWALD'
      call switch (mode)
c
c     compute the real space Ewald energy and first derivatives
c
!$acc parallel loop vector_length(32) async(dir_queue)
!$acc&         present(chgglobnl,iion,loc,x,y,z,pchg,nelst,elst,
!$acc&    dec,ec,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
!$acc&         private(ded,dedc)
#ifdef TINKER_DEBUG
!$acc&    present(ninte)
#endif
      do ii = 1, nionlocnl
         iichg = chgglobnl(ii)
         iglob = iion(iichg)
         i     = loc(iglob)
         xi    = x(iglob)
         yi    = y(iglob)
         zi    = z(iglob)
         fi    = f * pchg(iichg)
!$acc loop vector
         do kkk = 1, nelst(ii)
            kkchg = elst(kkk,ii)
            if (kkchg.eq.0) cycle
            kglob = iion(kkchg)
c
c     compute the energy contribution for this interaction
c
            xr = xi - x(kglob)
            yr = yi - y(kglob)
            zr = zi - z(kglob)
c
c     find energy for interactions within real space cutoff
c
            call image_inl (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. off2) then
#ifdef TINKER_DEBUG
!$acc atomic
               ninte(iglob) = ninte(iglob) + 1
               if (iglob.eq.1) print*,kglob
#endif
               fik   = fi*pchg(kkchg)
               k     = loc(kglob)
               call charge_couple(r2,xr,yr,zr,ebuffer
     &                           ,fik,aewald,scale_f
     &                           ,e,ded,dedc,0)
c
c     increment the overall energy and derivative expressions
c
               ec       = ec + tp2enr(e)
!$acc atomic
               dec(1,i) = dec(1,i) + dedc%x
!$acc atomic
               dec(2,i) = dec(2,i) + dedc%y
!$acc atomic
               dec(3,i) = dec(3,i) + dedc%z
!$acc atomic
               dec(1,k) = dec(1,k) - dedc%x
!$acc atomic
               dec(2,k) = dec(2,k) - dedc%y
!$acc atomic
               dec(3,k) = dec(3,k) - dedc%z
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

      call ecreal_scaling
#ifdef TINKER_DEBUG
!$acc wait
!$acc exit data copyout(ninte)
c     do i = 1,100
c        print*,i,ninte(i)
c     end do
      print*, 'charge interactions', sum(ninte)
#endif
      end
c
c     "ecrealshortlong1dgpu" evaluates either short or long range real space portion of the Ewald sum
c     energy and forces due to atomic charge interactions, using a neighbor list
c
      subroutine ecrealshortlong1dgpu
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
      use echarge1gpu_inl
      use energi
      use ewald
      use iounit
      use inform
      use inter
      use interfaces
      use math
      use molcul
      use neigh
      use potent
      use shunt
      use timestat
      use tinTypes
      use usage
      use utilgpu
      use virial
      use mpi
      implicit none
      integer i,j,k,iichg,iglob,kglob,kkchg
      integer ii,kkk,range_cfg
      integer ninte
      integer,pointer,save::lst(:,:),nlst(:)

      real(t_p) e
      real(t_p) f,fi,fik
      real(t_p) r,r2,rew

      real(t_p) xi,yi,zi
      real(t_p) xr,yr,zr

      real(t_p),parameter:: scale_f=1.0
      type(real3) ded
      type(mdyn3_r) dedc

      real(t_p) cshortcut2,coff
      character*10 mode
c
      if (deb_Path) write(*,'(2x,A)') 'ecrealshortlong1dgpu'
c
c     set conversion factor, cutoff and switching coefficients
c
#ifdef TINKER_DEBUG
      ninte = 0
#endif
      f = electric / dielec
      if (use_cshortreal) then
         mode       = 'SHORTEWALD'
         call switch (mode)
         range_cfg  = m_short
         cshortcut2 = 0
         coff       = off
         lst        => shortelst
         nlst       => nshortelst
      else if (use_clong) then
         mode       = 'EWALD'
         call switch (mode)
         range_cfg  = m_long
         cshortcut2 = (chgshortcut-shortheal)**2
         coff       = chgshortcut
         lst        => elst
         nlst       => nelst
      else
 12      format( "Unknown config for ecrealshortlong1dgpu" )
         print 12
         call fatal
      end if
c
c     compute the real space Ewald energy and first derivatives
c
!$acc parallel loop vector_length(32) async(dir_queue)
!$acc&         present(chgglobnl,iion,loc,x,y,z,pchg,nlst,lst,
!$acc&    dec,ec,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
!$acc&         private(ded,dedc)
#ifdef TINKER_DEBUG
!$acc&    reduction(+:ninte)
#endif
      do ii = 1, nionlocnl
         iichg = chgglobnl(ii)
         iglob = iion(iichg)
         i     = loc(iglob)
         xi    = x(iglob)
         yi    = y(iglob)
         zi    = z(iglob)
         fi    = f * pchg(iichg)
!$acc loop vector
         do kkk = 1, nlst(ii)
            kkchg = lst(kkk,ii)
            if (kkchg.eq.0) cycle
            kglob = iion(kkchg)
c
c     compute the energy contribution for this interaction
c
            xr = xi - x(kglob)
            yr = yi - y(kglob)
            zr = zi - z(kglob)
c
c     find energy for interactions within real space cutoff
c
            call image_inl (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2.ge.cshortcut2 .and. r2.le.off2) then
#ifdef TINKER_DEBUG
               ninte = ninte + 1
#endif
               fik   = fi*pchg(kkchg)
               k     = loc(kglob)
               call charge_couple_shortlong
     &                     (r2,xr,yr,zr,ebuffer
     &                     ,fik,aewald,scale_f,coff,shortheal
     &                     ,e,ded,dedc,0,range_cfg)
c
c     increment the overall energy and derivative expressions
c
               ec       = ec + tp2enr(e)
!$acc atomic
               dec(1,i) = dec(1,i) + dedc%x
!$acc atomic
               dec(2,i) = dec(2,i) + dedc%y
!$acc atomic
               dec(3,i) = dec(3,i) + dedc%z
!$acc atomic
               dec(1,k) = dec(1,k) - dedc%x
!$acc atomic
               dec(2,k) = dec(2,k) - dedc%y
!$acc atomic
               dec(3,k) = dec(3,k) - dedc%z
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

      call ecrealshortlong_scaling
#ifdef TINKER_DEBUG
!$acc wait
      print*, 'charge interactions', ninte
#endif
      end

#ifdef _CUDA
      subroutine ecreal1d_cu
      use atmlst
      use atoms
      use bound
      use boxes
      use charge
      use chgpot
      use couple
      use deriv
      use domdec
      use echargecu
      use energi , only: ec=>ec_r
      use ewald
      use iounit
      use inform
      use inter
      use math
      use molcul
      use mpi
      use neigh  , iion_s=>celle_glob,chg_s=>celle_chg
     &           , loc_s=>celle_loc, ieblst_s=>ieblst, eblst_s=>eblst
     &           , x_s=>celle_x, y_s=>celle_y, z_s=>celle_z
      use potent
      use shunt
      use timestat
      use usage
      use utilcu
      use utilgpu ,only: def_queue,dir_queue,vred_buff
     &            , reduce_energy_virial,zero_evir_red_buffer
     &            , ered_buff=>ered_buf1,dir_stream,def_stream
      use virial
      implicit none
      integer i,gs1
#ifdef TINKER_DEBUG
      integer inter(n)
#endif
      real(t_p) f
      real(t_p) p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
      logical,save:: first_in=.true.
      integer,save:: gS
      character*10 mode
      logical,parameter::dyn_gS=.true.

      if (deb_Path) write(*,'(2X,A)') 'ecreal1d_cu'
#ifdef TINKER_DEBUG
!$acc enter data create(inter) async(dir_queue)
!$acc parallel loop async(dir_queue) present(inter)
      do i = 1,n
         inter(i) = 0
      end do
#endif
c
c     set conversion factor, cutoff and switching coefficients
c
      f    = electric / dielec
      mode = 'EWALD'
      call switch (mode)

      if (first_in) then
         call cudaMaxGridSize("ecreal1d_core_cu",gS)
         first_in=.false.
      end if
      if(dyn_gS) gs = nionlocnlb2_pair/8

      p_xbeg = xbegproc(rank+1)
      p_xend = xendproc(rank+1)
      p_ybeg = ybegproc(rank+1)
      p_yend = yendproc(rank+1)
      p_zbeg = zbegproc(rank+1)
      p_zend = zendproc(rank+1)

      def_stream = dir_stream
      def_queue  = dir_queue
      call zero_evir_red_buffer(def_queue)
      call set_ChgData_CellOrder(.false.)

!$acc host_data use_device(iion_s,chg_s,loc_s,ieblst_s,eblst_s,
!$acc&   x_s,y_s,z_s,pchg,dec,ered_buff,vred_buff,
!$acc&   ccorrect_ik,ccorrect_scale,loc,x,y,z
#ifdef TINKER_DEBUG
!$acc&    ,inter
#endif
!$acc&    )

      call ecreal1d_core_cu<<<gS,BLOCK_DIM,0,def_stream>>>
     &     ( iion_s,chg_s,loc_s,ieblst_s
     &     , eblst_s(2*nionlocnlb_pair+1)
     &     , nionlocnlb,nionlocnlb2_pair,nionbloc,n
     &     , x_s,y_s,z_s,pchg
     &     , off2,f,aewald, ebuffer
     &     , dec, ered_buff, vred_buff
     &     , p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
#ifdef TINKER_DEBUG
     &     , inter
#endif
     &     )
      call check_launch_kernel(" ecreal1d_core_cu")

      gs1 = n_cscale/(BLOCK_DIM)
      call ecreal_scaling_cu<<<gs1,BLOCK_DIM,0,def_stream>>>
     &            ( ccorrect_ik,ccorrect_scale,loc,x,y,z
     &            , dec,ered_buff,vred_buff,n,nbloc,n_cscale
     &            , f,aewald,ebuffer,off2 )
      call check_launch_kernel(" ecreal_scaling_cu")

!$acc end host_data

      call reduce_energy_virial(ec,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz
     &                         ,ered_buff,def_queue)

#ifdef TINKER_DEBUG
 34   format(2I10,3F12.4)
 35   format(A30,2I16)
!$acc wait
!$acc exit data copyout(inter)
c     do i =1,100
c        print 34,i,inter(i)
c     end do
      print 35,'total charge interactions ', sum(inter)
#endif

      end subroutine
#endif
c
c     Scaling interaction correction subroutines
c
      subroutine ecreal_scaling
      use atmlst
      use atoms
      use bound
      use boxes
      use charge
      use chgpot
      use couple
      use deriv
      use domdec
      use echarge1gpu_inl
      use energi  ,only: ec=>ec_r
      use ewald
      use iounit
      use inform
      use inter
      use math
      use molcul
      use neigh
      use potent
      use shunt
      use timestat
      use tinTypes
      use usage
      use utilgpu
      use virial
      use mpi
      implicit none
      integer i,j,k,iichg,iglob
      integer ii,kkk,kglob,kkchg

      real(t_p) e
      real(t_p) f,fi,fik
      real(t_p) r,r2,rew
      real(t_p) rb,rb2

      real(t_p) xi,yi,zi
      real(t_p) xr,yr,zr

      real(t_p) scale_f
      type(real3) ded
      type(mdyn3_r) dedc

      character*10 mode

      if (deb_Path) write(*,'(3x,A)') 'ecreal_scaling'
c
c     set conversion factor, cutoff and switching coefficients
c
      f    = electric / dielec
      mode = 'EWALD'
      call switch (mode)
c
c     compute the real space Ewald energy and first derivatives
c
!$acc parallel loop vector_length(32) async(dir_queue)
!$acc&         present(ccorrect_ik,ccorrect_scale,loc,x,y,z,
!$acc&    dec,ec,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
!$acc&         private(ded,dedc)
      do ii = 1, n_cscale
         iglob = ccorrect_ik(ii,1)
         kglob = ccorrect_ik(ii,2)
         scale_f =   ccorrect_scale(2*ii+1)
         fik     = f*ccorrect_scale(2*ii+2)
         xi    = x(iglob)
         yi    = y(iglob)
         zi    = z(iglob)
c
c     compute the energy contribution for this interaction
c
         xr    = xi - x(kglob)
         yr    = yi - y(kglob)
         zr    = zi - z(kglob)
c
c     find energy for interactions within real space cutoff
c
         call image_inl (xr,yr,zr)
         r2 = xr*xr + yr*yr + zr*zr
         if (r2 .le. off2) then
            i     = loc(iglob)
            k     = loc(kglob)
            call charge_couple(r2,xr,yr,zr,ebuffer
     &                        ,fik,aewald,scale_f
     &                        ,e,ded,dedc,1)
c
c     increment the overall energy and derivative expressions
c
            ec       = ec + tp2enr(e)
!$acc atomic
            dec(1,i) = dec(1,i) + dedc%x
!$acc atomic
            dec(2,i) = dec(2,i) + dedc%y
!$acc atomic
            dec(3,i) = dec(3,i) + dedc%z
!$acc atomic
            dec(1,k) = dec(1,k) - dedc%x
!$acc atomic
            dec(2,k) = dec(2,k) - dedc%y
!$acc atomic
            dec(3,k) = dec(3,k) - dedc%z
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

      end

      subroutine ecrealshortlong_scaling
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
      use echarge1gpu_inl
      use energi
      use ewald
      use iounit
      use inform
      use inter
      use interfaces
      use math
      use molcul
      use neigh
      use potent
      use shunt
      use timestat
      use tinTypes
      use usage
      use utilgpu
      use virial
      use mpi
      implicit none
      integer i,j,k,iichg,iglob
      integer ii,kkk,kglob,kkchg
      integer range_cfg

      real(t_p) e
      real(t_p) f,fi,fik
      real(t_p) r,r2,rew
      real(t_p) rb,rb2

      real(t_p) xi,yi,zi
      real(t_p) xr,yr,zr

      real(t_p) scale_f
      type(real3) ded
      type(mdyn3_r) dedc

      real(t_p) cshortcut2,coff
      character*10 mode

      if (deb_Path) write(*,'(3x,A)') 'ecrealshortlond_scaling'
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      if (use_cshortreal) then
         range_cfg  = m_short
         cshortcut2 = 0
         coff       = off
      else if (use_clong) then
         range_cfg  = m_long
         cshortcut2 = (chgshortcut-shortheal)**2
         coff       = chgshortcut
      else
 12      format( "Unknown config for ecrealshortlong1dgpu" )
         print 12
         call fatal
      end if
c
c     compute the real space Ewald energy and first derivatives
c
!$acc parallel loop vector_length(32) async(dir_queue)
!$acc&         present(ccorrect_ik,ccorrect_scale,loc,x,y,z,
!$acc&    dec,ec,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
!$acc&         private(ded,dedc)
      do ii = 1, n_cscale
         iglob = ccorrect_ik(ii,1)
         kglob = ccorrect_ik(ii,2)
         scale_f =   ccorrect_scale(2*ii+1)
         fik     = f*ccorrect_scale(2*ii+2)
         xi    = x(iglob)
         yi    = y(iglob)
         zi    = z(iglob)
c
c     compute the energy contribution for this interaction
c
         xr    = xi - x(kglob)
         yr    = yi - y(kglob)
         zr    = zi - z(kglob)
c
c     find energy for interactions within real space cutoff
c
         call image_inl (xr,yr,zr)
         r2 = xr*xr + yr*yr + zr*zr
         if (r2.ge.cshortcut2 .and. r2.le.off2) then
            i     = loc(iglob)
            k     = loc(kglob)
            call charge_couple_shortlong
     &                 (r2,xr,yr,zr,ebuffer
     &                 ,fik,aewald,scale_f,coff,shortheal
     &                 ,e,ded,dedc,1,range_cfg)
c
c     increment the overall energy and derivative expressions
c
            ec       = ec + e
!$acc atomic
            dec(1,i) = dec(1,i) + dedc%x
!$acc atomic
            dec(2,i) = dec(2,i) + dedc%y
!$acc atomic
            dec(3,i) = dec(3,i) + dedc%z
!$acc atomic
            dec(1,k) = dec(1,k) - dedc%x
!$acc atomic
            dec(2,k) = dec(2,k) - dedc%y
!$acc atomic
            dec(3,k) = dec(3,k) - dedc%z
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
      
      end
c
c
c     ####################################################################
c     ##                                                                ##
c     ##  subroutine ecrecip1gpu  --  PME recip charge energy & derivs  ##
c     ##                                                                ##
c     ####################################################################
c
c
c     "ecrecip1" evaluates the reciprocal space portion of the particle
c     mesh Ewald summation energy and gradient due to partial charges
c
c     literature reference:
c
c     U. Essmann, L. Perera, M. L Berkowitz, T. Darden, H. Lee and
c     L. G. Pedersen, "A Smooth Particle Mesh Ewald Method", Journal
c     of Chemical Physics, 103, 8577-8593 (1995)
c
c     modifications for nonperiodic systems suggested by Tom Darden
c     during May 2007
c
c
      subroutine ecrecip1gpu
      use atmlst
      use atoms
      use bound
      use boxes
      use charge
      use chgpot
      use deriv
      use domdec
      use energi
      use ewald
      use fft
      use inform
      use interfaces ,only: grid_pchg_site_p,grid_pchg_force_p
      use math
      use pme
      use pme1
      use potent
      use timestat
      use utils
      use utilgpu
      use virial
      use mpi
      implicit none
      integer i,j,k
      integer iichg,iglob
      integer isite,iatm
      integer i0,j0,k0
      integer it1,it2,it3
      integer k1,k2,k3
      integer m1,m2,m3
      integer nf1,nf2,nf3
      integer nff
      integer status(MPI_STATUS_SIZE),tag,ierr
      integer nprocloc,rankloc,commloc,proc
      integer ist2,jst2,kst2,ien2,jen2,ken2
      integer istart,iend,jstart,jend,kstart,kend
      integer iloc,iproc
      integer igrd0,jgrd0,kgrd0
      real(t_p) e,term,expterm
      real(t_p) vterm,pterm
      real(t_p) volterm
      real(t_p) f,fi,denom
      real(t_p) hsq,struc2
      real(t_p) de1,de2,de3
      real(t_p) dn1,dn2,dn3
      real(t_p) t1,t2,t3
      real(t_p) dt1,dt2,dt3
      real(t_p) h1,h2,h3
      real(t_p) r1,r2,r3
      real(t_p),save:: vxx,vxy,vxz,vyy,vyz,vzz
      integer  , allocatable :: req(:),reqbcast(:)
c     real(t_p), allocatable :: qgridmpi(:,:,:,:,:)
      real(8) time0,time1
      logical,save::f_in=.true.
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6) return

      if (deb_Path) write(*,'(2x,A)') 'ecrecip1gpu'
c
      if (use_pmecore) then
        nprocloc = nrec
        commloc  = comm_rec
        rankloc  = rank_bis
      else
        nprocloc = nproc
        commloc  = COMM_TINKER
        rankloc  = rank
      end if

      if (f_in) then
!$acc enter data create(vxx,vxy,vxz,vyy,vyz,vzz)
         f_in = .false.
      end if
c
c     dynamic allocation of local arrays
c
      call mallocMpiGrid
      allocate (req(nprocloc*nprocloc))
      allocate (reqbcast(nprocloc*nprocloc))

!$acc data present(vxx,vxy,vxz,vyy,vyz,vzz)
!$acc&     present(qgridin_2d,qgridout_2d,recip,ecrec,
!$acc&        g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz,vir,
!$acc&        bsmod1,bsmod2,bsmod3)
c
!$acc serial async(rec_queue)
      vxx = 0
      vxy = 0
      vxz = 0
      vyy = 0
      vyz = 0
      vzz = 0
!$acc end serial

c
      call timer_enter( timer_grid1 )
      call bspline_fill_sitegpu(1)
      call set_to_zero1(qgridin_2d(1,1,1,1,1),
     &     2*n1mpimax*n2mpimax*n3mpimax*(nrec_send+1),rec_queue)

      call grid_pchg_site_p
      call timer_exit ( timer_grid1,quiet_timers )
c
c     MPI : Begin reception
c
      if (nrec_send.gt.0) then
!$acc wait(rec_queue)
      end if
      call timer_enter( timer_recreccomm )
!$acc host_data use_device(qgridmpi,qgridin_2d)
      do i = 1, nrec_recep
        tag = nprocloc*rankloc + prec_recep(i) + 1
        call MPI_IRECV(qgridmpi(1,1,1,1,i),
     $       2*n1mpimax*n2mpimax*n3mpimax,MPI_TPREC,
     $       prec_recep(i),tag,commloc,req(tag),ierr)
      end do
c
c     MPI : begin sending
c
      do i = 1, nrec_send
        proc = prec_send(i)
        tag  = nprocloc*prec_send(i) + rankloc + 1
        call MPI_ISEND(qgridin_2d(1,1,1,1,i+1),
     $       2*n1mpimax*n2mpimax*n3mpimax,MPI_TPREC,
     $       proc,tag,commloc,req(tag),ierr)
      end do
!$acc end host_data

      do i = 1, nrec_recep
        tag = nprocloc*rankloc + prec_recep(i) + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
      do i = 1, nrec_send
        tag = nprocloc*prec_send(i) + rankloc + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
c
c     do the reduction 'by hand'
c
      do i = 1, nrec_recep
         call aaddgpuAsync(2*n1mpimax*n2mpimax*n3mpimax,
     $        qgridin_2d(1,1,1,1,1),qgridmpi(1,1,1,1,i),
     $        qgridin_2d(1,1,1,1,1))
      end do
      call timer_exit ( timer_recreccomm,quiet_timers )
c
c     perform the 3-D FFT forward transformation
c
#ifdef _OPENACC
      call cufft2d_frontmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
     $                      n3mpimax)
#else
      call   fft2d_frontmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
     $                      n3mpimax)
#endif
#ifdef TINKER_DEBUG
      if (rankloc.eq.0) then
!$acc update host(qgridin_2d,qgridout_2d)
      print*,'gridi norm2',comput_norm(qgridin_2d,size(qgridin_2d),2)
      print*,'grido norm2',comput_norm(qgridout_2d,size(qgridout_2d),2)
      end if
#endif
c
c     use scalar sum to get reciprocal space energy and virial
c
      call timer_enter( timer_scalar )
      ist2 = istart2(rankloc+1)
      jst2 = jstart2(rankloc+1)
      kst2 = kstart2(rankloc+1)
      ien2 =   iend2(rankloc+1)
      jen2 =   jend2(rankloc+1)
      ken2 =   kend2(rankloc+1)
c     if ((istart2(rankloc+1).eq.1).and.(jstart2(rankloc+1).eq.1).and.
c    $    (kstart2(rankloc+1).eq.1)) then
c          qfac_2d(1,1,1) = 0.0d0
c     end if
      f       = 0.5d0 * electric / dielec
      pterm   = (pi/aewald)**2
      volterm = pi * volbox
      nff     = nfft1 * nfft2
      nf1     = (nfft1+1) / 2
      nf2     = (nfft2+1) / 2
      nf3     = (nfft3+1) / 2
!$acc parallel loop collapse(3) async(rec_queue)
      do k3 = kst2,ken2
        do k2 = jst2,jen2
          do k1 = ist2,ien2
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
                  expterm = expterm * (1.0_ti_p-cos(pi*xbox*sqrt(hsq)))
               else if (octahedron) then
                  if (mod(m1+m2+m3,2).ne.0)  expterm = 0.0_ti_p
               end if
               struc2 = qgridout_2d(1,k1-ist2+1,k2-jst2+1,k3-kst2+1)**2
     $                + qgridout_2d(2,k1-ist2+1,k2-jst2+1,k3-kst2+1)**2
                e = f * expterm * struc2
               ecrec = ecrec + e
               vterm = (2.0_ti_p/hsq) * (1.0_ti_p-term) * e
               vxx   = vxx + h1*h1*vterm - e
               vxy   = vxy + h1*h2*vterm
               vxz   = vxz + h1*h3*vterm
               vyy   = vyy + h2*h2*vterm - e
               vyz   = vyz + h3*h2*vterm
               vzz   = vzz + h3*h3*vterm - e
            end if
            qgridout_2d(1,k1-ist2+1,k2-jst2+1,k3-kst2+1) = expterm *
     $      qgridout_2d(1,k1-ist2+1,k2-jst2+1,k3-kst2+1)

            qgridout_2d(2,k1-ist2+1,k2-jst2+1,k3-kst2+1) = expterm *
     $      qgridout_2d(2,k1-ist2+1,k2-jst2+1,k3-kst2+1)
 10         continue
          end do
        end do
      end do
c
c     account for zeroth grid point for nonperiodic system
c
      if ((istart2(rankloc+1).eq.1).and.(jstart2(rankloc+1).eq.1).and.
     $   (kstart2(rankloc+1).eq.1)) then
        if (.not. use_bounds) then
           expterm = 0.5_ti_p * pi / xbox
!$acc serial async(rec_queue)
           struc2 = qgridout_2d(1,1,1,1)**2 + qgridout_2d(2,1,1,1)**2
               e = f * expterm * struc2
           ecrec = ecrec + e
           qgridout_2d(1,1,1,1) = expterm * qgridout_2d(1,1,1,1)
           qgridout_2d(2,1,1,1) = expterm * qgridout_2d(2,1,1,1)
!$acc end serial
        end if
      end if
!$acc serial async(rec_queue)
      vir(1,1) = vir(1,1) + vxx
      vir(2,1) = vir(2,1) + vxy
      vir(3,1) = vir(3,1) + vxz
      vir(1,2) = vir(1,2) + vxy
      vir(2,2) = vir(2,2) + vyy
      vir(3,2) = vir(3,2) + vxz
      vir(1,3) = vir(1,3) + vxz
      vir(2,3) = vir(2,3) + vyz
      vir(3,3) = vir(3,3) + vzz
!$acc end serial
      call timer_exit ( timer_scalar,quiet_timers )
c
c     perform the 3-D FFT backward transformation
c
#ifdef _OPENACC
      call cufft2d_backmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
     $                     n3mpimax)
#else
      call   fft2d_backmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
     $                     n3mpimax)
#endif
c
c     MPI : Begin reception
c
      if (nrec_recep.gt.0) then
!$acc wait(rec_queue)
      end if
      call timer_enter( timer_recreccomm )
!$acc host_data use_device(qgridin_2d)
      do i = 1, nrec_send
         proc = prec_send(i)
         tag  = nprocloc*rankloc + prec_send(i) + 1
         call MPI_IRECV(qgridin_2d(1,1,1,1,i+1),
     $        2*n1mpimax*n2mpimax*n3mpimax,MPI_TPREC,
     $        prec_send(i),tag,commloc,reqbcast(tag),ierr)
      end do
c
c     MPI : begin sending
c
      do i = 1, nrec_recep
         tag = nprocloc*prec_recep(i) + rankloc + 1
         call MPI_ISEND(qgridin_2d,2*n1mpimax*n2mpimax*n3mpimax,
     $        MPI_TPREC,prec_recep(i),tag,commloc,reqbcast(tag),
     $        ierr)
      end do
!$acc end host_data
c
      do i = 1, nrec_send
         tag = nprocloc*rankloc + prec_send(i) + 1
         call MPI_WAIT(reqbcast(tag),status,ierr)
      end do
c
      do i = 1, nrec_recep
         tag = nprocloc*prec_recep(i) + rankloc + 1
         call MPI_WAIT(reqbcast(tag),status,ierr)
      end do
      call timer_exit ( timer_recreccomm,quiet_timers )

#ifdef TINKER_DEBUG
      if (rankloc.eq.0) then
!$acc update host(qgridin_2d,qgridout_2d)
      print*,'grido norm2',comput_norm(qgridout_2d,size(qgridout_2d),2)
      print*,'gridi norm2',comput_norm(qgridin_2d,size(qgridin_2d),2)
      end if
#endif
c
c     get first derivatives of the reciprocal space energy
c
      call timer_enter( timer_grid2 )
      call grid_pchg_force_p
      call timer_exit ( timer_grid2,quiet_timers )

!$acc end data
c     deallocate (qgridmpi)
      deallocate (req)
      deallocate (reqbcast)
      end
