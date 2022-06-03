c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine epolar3  --  induced dipole energy & analysis  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "epolar3" calculates the induced dipole polarization energy,
c     and partitions the energy among atoms
c
c
#include "tinker_precision.h"
#include "tinker_types.h"
      module epolar3gpu_inl
        include "erfcore_data.f.inc"
        contains
#include "convert.f.inc"
#include "image.f.inc"
#if defined(SINGLE) | defined(MIXED)
        include "erfcscore.f.inc"
#else
        include "erfcdcore.f.inc"
#endif
#include "pair_polar.f.inc"
      end module

      subroutine epolar3gpu
      implicit none
c
c
c     choose the method for summing over polarization interactions
c
      call epolar3cgpu
      end

c     #####################################################################
c     ##                                                                 ##
c     ##  subroutine epolar3cgpu  --  Ewald polarization analysis; list  ##
c     ##                                                                 ##
c     #####################################################################
c
c
c     "epolar3c" calculates the polarization energy and analysis with
c     respect to Cartesian coordinates using particle mesh Ewald and
c     a neighbor list
c
c
      subroutine epolar3cgpu
      use action
      use analyz
      use atmlst
      use atoms
      use boxes
      use chgpot
      use domdec
      use energi     ,only: ep,ep_r,eprec
      use ewald
      use epolar3gpu_inl
      use inform     ,only: deb_Path
      use interfaces ,only: epreal3d_p
      use math
      use mpole
      use polar
      use polpot
      use potent
      use mpi
      use utilgpu
      use sizes
      use timestat

      implicit none
      integer i,k,ii,ierr
      integer iipole,iglob
      real(t_p) e,f
      real(t_p) term,fterm
      real(t_p) dix,diy,diz
      real(t_p) uix,uiy,uiz,uii
      real(t_p) xd,yd,zd
      real(t_p) xu,yu,zu
c
      if (npole .eq. 0)  return
      if(deb_Path) write(*,*) 'epolar3cgpu'
c
c     zero out the dipole polarization energy and components
c
!$acc serial async(rec_queue) present(nep,nep_,ep,eprec)
      nep   = 0
      nep_  = 0.0
      ep    = 0
      eprec = 0
!$acc end serial
c     aep = 0.0_ti_p
c
c     set the energy unit conversion factor
c
      f = electric / dielec
c
c     check the sign of multipole components at chiral sites
c
      if (.not. use_mpole)  call chkpolegpu(.false.)
c
c     rotate the multipole components into the global frame
c
      if (.not. use_mpole)  call rotpolegpu
c
c     compute the induced dipoles at each polarizable atom
c
      if (use_polarshortreal) then
         call newinduce_shortrealgpu
      else if (use_pmecore) then
         if (polalg.eq.5) then
            call dcinduce_pme
         else
            call newinduce_pmegpu
         end if
      else
         if (polalg.eq.5) then
            call dcinduce_pme2gpu
         else
            call newinduce_pme2gpu
         end if
      end if
c
c     compute the real space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.le.ndir-1))
     &   then
         call timer_enter( timer_real )
         def_queue = dir_queue
         if (use_preal) then
           if (use_polarshortreal) then
             call epreal3d_p
           else
             call epreal3d_p
           end if
         end if

         if (use_pself) then
!$acc data present(poleglob,ipole,loc,rpole,uind,ep,ep_r,nep_)

c
c     compute the Ewald self-energy term over all the atoms
c
          term = 2.0_ti_p * aewald * aewald
          fterm = -f * aewald / sqrtpi
!$acc parallel loop async(def_queue)
          do ii = 1, npoleloc
             iipole = poleglob(ii)
             iglob  = ipole(iipole)
             dix    = rpole(2,iipole)
             diy    = rpole(3,iipole)
             diz    = rpole(4,iipole)
             uix    = uind(1,iipole)
             uiy    = uind(2,iipole)
             uiz    = uind(3,iipole)
             uii    = dix*uix + diy*uiy + diz*uiz
             e      = fterm * term * uii / 3.0_ti_p
             ep_r   = ep_r + tp2enr(e)
             nep_   = nep_ + 1
          end do
c
c         compute the cell dipole boundary correction term
c
          if (boundary .eq. 'VACUUM') then
             xd = 0.0_ti_p
             yd = 0.0_ti_p
             zd = 0.0_ti_p
             xu = 0.0_ti_p
             yu = 0.0_ti_p
             zu = 0.0_ti_p
!$acc parallel loop async(def_queue)
             do ii = 1, npoleloc
                iipole = poleglob(ii)
                iglob = ipole(iipole)
                i = loc(iglob)
                xd = xd + rpole(2,iipole) + rpole(1,iipole)*x(iglob)
                yd = yd + rpole(3,iipole) + rpole(1,iipole)*y(iglob)
                zd = zd + rpole(4,iipole) + rpole(1,iipole)*z(iglob)
                xu = xu + uind(1,iipole)
                yu = yu + uind(2,iipole)
                zu = zu + uind(3,iipole)
             end do
             !implicit wait added here by the compiler
             call MPI_ALLREDUCE(MPI_IN_PLACE,xd,1,MPI_TPREC,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,yd,1,MPI_TPREC,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,zd,1,MPI_TPREC,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,xu,1,MPI_TPREC,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,yu,1,MPI_TPREC,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,zu,1,MPI_TPREC,MPI_SUM,
     $          COMM_TINKER,ierr)
             if (rank.eq.0) then
               term = (2.0_ti_p/3.0_ti_p) * f * (pi/volbox)
!$acc serial async
               ep   = ep + term*(xd*xu+yd*yu+zd*zu)
               nep_ = nep_+ 1.0
!$acc end serial
             end if
          end if

!$acc end data
         end if
         call timer_exit( timer_real,quiet_timers )
      end if
c
c     compute the reciprocal space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.gt.ndir-1))
     &   then
         if (use_prec) then
            call timer_enter( timer_rec )
            call eprecipgpu
            call timer_exit( timer_rec,quiet_timers )
         end if
      end if
!$acc serial async(rec_queue) present(nep_,nep,ep,eprec,ep_r)
         ep  =  ep + eprec + enr2en(ep_r)
         nep = nep + int(nep_)
!$acc end serial
      end
c
c
c     ######################################################################
c     ##                                                                  ##
c     ##  subroutine epreal3dgpu  --  real space polar analysis via list  ##
c     ##                                                                  ##
c     ######################################################################
c
c
c     "epreal3d" calculates the induced dipole polarization energy
c     and analysis using particle mesh Ewald and a neighbor list
c
c
      subroutine epreal3dgpu
      use action  ,only: nep_
      use analyz
      use atoms   ,only: x,y,z
      use chgpot  ,only: dielec,electric
      use atmlst  ,only: poleglobnl
      use atmtyp
      use bound
      use couple
      use cutoff  ,only: shortheal
      use domdec
      use energi  ,only: ep=>ep_r
      use ewald
      use epolar3gpu_inl
      use erf_mod
      use inform
      use inter
      use iounit
      use math
      use molcul
      use mpole
      use neigh
      use polar
      use polgrp
      use polpot
      use potent
      use shunt
      use mpi
      use utilgpu
      use timestat

      implicit none

      integer i,iglob,j,k,iploc,kploc
      integer nnelst
      integer ii,iipole
      integer kpole,kglob,kbis
      integer,pointer :: lst(:,:),nlst(:)
      real(t_p) alsq2,alsq2n
      real(t_p) r2,pgamma,damp
      real(t_p) f,e
      real(t_p) pdi,pti
      real(t_p) one
      real(t_p) pscale,dscale,uscale

      type(rpole_elt):: ip,kp
      type(real6) :: dpui,dpuk
      type(real3) :: pos,posi
      parameter(one=1.0_ti_p)
      character*10:: mode


      if (deb_Path) 
     &   write(*,'(2x,a)') 'epreal3dgpu'

      f = 0.5_ti_p * electric / dielec
      if (use_polarshortreal) then
         mode = 'SHORTEWALD'
         call switch (mode)
          lst =>  shortelst
         nlst => nshortelst
      else
         mode = 'MPOLE     '
         call switch (mode)
          lst =>  elst
         nlst => nelstc
      end if

      alsq2  = 2.0_ti_p * aewald**2
      alsq2n = 0.0_ti_p
      if (aewald .gt. 0.0_ti_p)  alsq2n = 1.0_ti_p / (sqrtpi*aewald)
      call timer_enter( timer_epreal )

!$acc enter data attach(nlst,lst) async(def_queue)

!$acc parallel loop gang vector_length(32)
!$acc&         present(poleglobnl,ipole,rpole,thole,pdamp,
!$acc&  loc,x,y,z,uind,uinp,lst,nlst,polelocnl,ep,nep_)
!$acc&         private(ip,dpui,posi)
!$acc&         async(def_queue)
      MAINLOOP:
     &do ii = 1, npolelocnl
         iipole  = poleglobnl(ii)
         iglob   = ipole (iipole)
         i       = loc(iglob)
         nnelst  = nlst (ii)
         !No neighbours
         if (nnelst.eq.0) cycle MAINLOOP
         posi%x  = x(iglob)
         posi%y  = y(iglob)
         posi%z  = z(iglob)

         pdi     = pdamp(iipole)
         pti     = thole(iipole)

         ip%c    = rpole( 1, iipole)
         ip%dx   = rpole( 2, iipole)
         ip%dy   = rpole( 3, iipole)
         ip%dz   = rpole( 4, iipole)
         ip%qxx  = rpole( 5, iipole)
         ip%qxz  = rpole( 7, iipole)
         ip%qxy  = rpole( 6, iipole)
         ip%qyy  = rpole( 9, iipole)
         ip%qyz  = rpole(10, iipole)
         ip%qzz  = rpole(13, iipole)

         dpui%x  = uind ( 1, iipole)
         dpui%y  = uind ( 2, iipole)
         dpui%z  = uind ( 3, iipole)
         dpui%xx = uinp ( 1, iipole)
         dpui%yy = uinp ( 2, iipole)
         dpui%zz = uinp ( 3, iipole)
c
c     loop on the neighbors
c
!$acc loop vector private(dpuk,kp,pos)
         do k = 1, nnelst
            kpole    = lst(k,ii)
            kglob    = ipole(kpole)
            kbis     = loc  (kglob)
            kploc    = polelocnl(kpole)
            pos%x    = x(kglob) - posi%x
            pos%y    = y(kglob) - posi%y
            pos%z    = z(kglob) - posi%z

            call image_inl(pos%x,pos%y,pos%z)
            r2       = pos%x**2 + pos%y**2 + pos%z**2
            if (r2>off2) cycle
c
            kp%c     = rpole( 1, kpole)
            kp%dx    = rpole( 2, kpole)
            kp%dy    = rpole( 3, kpole)
            kp%dz    = rpole( 4, kpole)
            kp%qxx   = rpole( 5, kpole)
            kp%qxy   = rpole( 6, kpole)
            kp%qxz   = rpole( 7, kpole)
            kp%qyy   = rpole( 9, kpole)
            kp%qyz   = rpole(10, kpole)
            kp%qzz   = rpole(13, kpole)

            dpuk%x   = uind ( 1, kpole)
            dpuk%y   = uind ( 2, kpole)
            dpuk%z   = uind ( 3, kpole)
            dpuk%xx  = uinp ( 1, kpole)
            dpuk%yy  = uinp ( 2, kpole)
            dpuk%zz  = uinp ( 3, kpole)

            pgamma   = min( pti,thole(kpole) )
            damp     = pdi * pdamp (kpole)
c
c     Compute polar interaction
c
            call epolar3_couple(dpui,ip,dpuk,kp,r2,pos,
     &               aewald,alsq2,alsq2n,pgamma,damp,f,
     &               off,shortheal,1.0_ti_p,
     &               e,use_polarshortreal,.false.)
c
c     increment energy and interaction
c
            ep   =  ep  + tp2enr(e)
            nep_ = nep_ + 1
         enddo

      end do  MAINLOOP
!$acc exit data detach(nlst,lst) async(def_queue)
c
      call epreal3c_correct_scale

      call timer_exit( timer_epreal  )
      end

      subroutine epreal3d_cu
      use action  ,only: nep
      use atmlst  ,only: poleglobnl
      use atoms   ,only: x,y,z,n
      use chgpot  ,only: dielec,electric
      use cutoff  ,only: shortheal
      use domdec  ,only: xbegproc,ybegproc,zbegproc
     &            ,nproc,rank,xendproc,yendproc,zendproc
     &            ,nbloc,loc
      use energi  ,only: ep=>ep_r
#ifdef _CUDA
      use epolar1cu,only:epreal3_cu
#endif
      use ewald   ,only: aewald
      use inform  ,only: deb_Path
      use interfaces ,only: epreal3c_correct_scale
      use math    ,only: sqrtpi
      use mpole   ,only: npolelocnl,ipole,rpole,polelocnl
     &            ,npolelocnlb,npolebloc
     &            ,npolelocnlb_pair,npolelocnlb2_pair
     &            ,nspnlb2=>nshortpolelocnlb2_pair
      use neigh   ,only:ipole_s=>celle_glob,pglob_s=>celle_pole
     &            , loc_s=>celle_loc, plocnl_s=>celle_plocnl
     &            , ieblst_s=>ieblst, eblst_s=>eblst
     &            , iseblst_s=>ishorteblst, seblst_s=>shorteblst
     &            , x_s=>celle_x, y_s=>celle_y, z_s=>celle_z
      use polar   ,only: uind,uinp,thole,pdamp
      use potent  ,only: use_polarshortreal
      use shunt   ,only: off2,off
      use tinheader,only: ti_p
#ifdef _CUDA
      use cudafor
      use utilcu  ,only: BLOCK_DIM,check_launch_kernel
#endif
      use utilgpu ,only: def_queue,dir_queue,rec_queue
     &            ,real3,real6,real3_red,rpole_elt
     &            ,ered_buff=>ered_buf1,nred_buff,reduce_energy_action
     &            ,RED_BUFF_SIZE,zero_en_red_buffer
#ifdef  _OPENACC
     &            ,dir_stream,def_stream,nSMP
#endif

      implicit none

      integer i,iglob,j,k,iploc,kploc
      integer ii,iipole,ierrSync
      integer lst_beg
      integer,save:: gS=0
      real(t_p) alsq2,alsq2n,f
      real(t_p) p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
      logical,save::first_in=.true.
      logical,parameter::dyn_gS=.false.
      character*10:: mode

      if(deb_Path)write(*,'(2x,a)') 'epreal3d_cu'
c
c     set conversion factor, cutoff and switching coefficients
c
      f      = 0.5_ti_p * electric / dielec
      alsq2  = 2.0_ti_p * aewald**2
      alsq2n = 0.0_ti_p
      if (aewald .gt. 0.0_ti_p) alsq2n = 1.0_ti_p / (sqrtpi * aewald)

      p_xbeg = xbegproc(rank+1)
      p_xend = xendproc(rank+1)
      p_ybeg = ybegproc(rank+1)
      p_yend = yendproc(rank+1)
      p_zbeg = zbegproc(rank+1)
      p_zend = zendproc(rank+1)
      lst_beg= 2*npolelocnlb_pair+1

      if (use_polarshortreal) then
         mode = 'SHORTEWALD'
         call switch (mode)
      else
         mode = 'MPOLE     '
         call switch (mode)
      end if

#ifdef _CUDA
      def_stream = dir_stream
      if (first_in) then
         ! Compute though occupancy the right gridSize to launch the kernel with
         call cudaMaxGridSize("epreal3_cu",gS)
         if (deb_Path) print*,'epreal3_cu blockSize',gS
         first_in = .false.
      end if
      if (dyn_gS) gS = max(npolelocnlb2_pair/8,1)

      call zero_en_red_buffer(def_queue)

!$acc host_data use_device(ipole_s,pglob_s,loc_s,plocnl_s,ieblst_s,
!$acc&    iseblst_s,eblst_s,seblst_s,x_s,y_s,z_s,rpole,pdamp,thole,
!$acc&    uind,uinp,ered_buff,nred_buff)
      
      if (use_polarshortreal) then
      call epreal3_cu<<<gS,BLOCK_DIM,0,def_stream>>>
     &     (ipole_s,pglob_s,loc_s,plocnl_s
     &     ,iseblst_s,seblst_s(lst_beg)
     &     ,x_s,y_s,z_s,rpole,pdamp,thole,uind,uinp
     &     ,ered_buff,nred_buff
     &     ,npolelocnlb,nspnlb2,npolebloc,n
     &     ,off2,f,alsq2,alsq2n,aewald,off,shortheal
     &     ,use_polarshortreal
     &     ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend)
      else
      call epreal3_cu<<<gS,BLOCK_DIM,0,def_stream>>>
     &     (ipole_s,pglob_s,loc_s,plocnl_s
     &     ,ieblst_s,eblst_s(lst_beg)
     &     ,x_s,y_s,z_s,rpole,pdamp,thole,uind,uinp
     &     ,ered_buff,nred_buff
     &     ,npolelocnlb,npolelocnlb2_pair,npolebloc,n
     &     ,off2,f,alsq2,alsq2n,aewald,off,shortheal
     &     ,use_polarshortreal
     &     ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend)
      end if
      call check_launch_kernel(" epreal3_cu")

!$acc end host_data

      call reduce_energy_action(ep,nep,ered_buff,def_queue)
c
      call epreal3c_correct_scale
#else
      print 100
 100  format('epreal3_cu is a specific device routine',/,
     &       'you are not supposed to get inside with your compile',
     &       'type.')
      call fatal
#endif

      end


      subroutine epreal3c_correct_scale
      use action  ,only: nep
      use atmlst  ,only: poleglobnl
      use atoms   ,only: x,y,z
      use chgpot  ,only: dielec,electric
      use cutoff  ,only: shortheal
      use domdec  ,only: rank,loc
      use energi  ,only: ep=>ep_r
      use epolar3gpu_inl
      use ewald   ,only: aewald
      use inform  ,only: deb_Path
      use math    ,only: sqrtpi
      use mpole   ,only: npolelocnl,ipole,rpole,polelocnl
      use polar   ,only: uind,uinp,thole,pdamp
      use polpot  ,only: n_dpuscale,dpucorrect_ik,dpucorrect_scale
      use potent  ,only: use_polarshortreal
      use shunt   ,only: off2,off
      use tinheader,only: ti_p
      use utilgpu ,only: def_queue,real3,real6,real3_red,rpole_elt
      use atoms   ,only: x,y,z

      implicit none

      integer i,k,iglob,kglob,iploc,kploc
      integer nnelst
      integer ii,iipole,kpole
      integer j,kbis
      real(t_p) alsq2,alsq2n
      real(t_p) r2,pgamma,damp
      real(t_p) f,e
      real(t_p) pdi,pti
      real(t_p) one
      real(t_p) pscale,dscale,uscale

      type(rpole_elt):: ip,kp
      type(real6) :: dpui,dpuk
      type(real3) :: pos,posi

      parameter(one=1.0_ti_p)

      if(deb_Path)
     &   write(*,'(2x,a)') 'epreal3c_correct_scale'
c
c     set conversion factor, cutoff and switching coefficients
c
      f      = 0.5_ti_p * electric / dielec
      alsq2  = 2.0_ti_p * aewald**2
      alsq2n = 0.0_ti_p
      if (aewald .gt. 0.0_ti_p) alsq2n = 1.0_ti_p / (sqrtpi * aewald)

!$acc parallel loop gang vector_length(32)
!$acc&         present(ipole,rpole,thole,pdamp,loc,
!$acc&     x,y,z,uind,uinp,polelocnl,dpucorrect_ik,
!$acc&     dpucorrect_scale,nep,ep)
!$acc&     private(pos,ip,kp,dpui,dpuk)
!$acc&         reduction(+:nep)
!$acc&         async(def_queue)
      do ii = 1, n_dpuscale
         iipole   = dpucorrect_ik(2*(ii-1)+1)
         kpole    = dpucorrect_ik(2*(ii-1)+2)

         !dscale   = dpucorrect_scale(3*(ii-1)+1)
         pscale   = dpucorrect_scale(3*(ii-1)+2)
         !uscale   = dpucorrect_scale(3*(ii-1)+3)
         if (pscale.eq.0.0) cycle

         iglob    = ipole(iipole)
         kglob    = ipole(kpole)
         i        = loc  (iglob)
         k        = loc  (kglob)
         iploc    = polelocnl(iipole)
         kploc    = polelocnl(kpole)

         pos%x    = x(kglob) - x(iglob)
         pos%y    = y(kglob) - y(iglob)
         pos%z    = z(kglob) - z(iglob)

         call image_inl(pos%x,pos%y,pos%z)
         ! cutoff
         r2       = pos%x**2 + pos%y**2 + pos%z**2
         if (r2>off2) cycle

         pdi      = pdamp(iipole)
         pti      = thole(iipole)

         ip%c     = rpole( 1,iipole)
         ip%dx    = rpole( 2,iipole)
         ip%dy    = rpole( 3,iipole)
         ip%dz    = rpole( 4,iipole)
         ip%qxx   = rpole( 5,iipole)
         ip%qxy   = rpole( 6,iipole)
         ip%qxz   = rpole( 7,iipole)
         ip%qyy   = rpole( 9,iipole)
         ip%qyz   = rpole(10,iipole)
         ip%qzz   = rpole(13,iipole)

         dpui%x   = uind ( 1,iipole)
         dpui%y   = uind ( 2,iipole)
         dpui%z   = uind ( 3,iipole)
         dpui%xx  = uinp ( 1,iipole)
         dpui%yy  = uinp ( 2,iipole)
         dpui%zz  = uinp ( 3,iipole)

         kp%c     = rpole( 1, kpole)
         kp%dx    = rpole( 2, kpole)
         kp%dy    = rpole( 3, kpole)
         kp%dz    = rpole( 4, kpole)
         kp%qxx   = rpole( 5, kpole)
         kp%qxy   = rpole( 6, kpole)
         kp%qxz   = rpole( 7, kpole)
         kp%qyy   = rpole( 9, kpole)
         kp%qyz   = rpole(10, kpole)
         kp%qzz   = rpole(13, kpole)

         dpuk%x   = uind ( 1, kpole)
         dpuk%y   = uind ( 2, kpole)
         dpuk%z   = uind ( 3, kpole)
         dpuk%xx  = uinp ( 1, kpole)
         dpuk%yy  = uinp ( 2, kpole)
         dpuk%zz  = uinp ( 3, kpole)

         pgamma   = min( pti,thole(kpole) )
         damp     = pdi * pdamp (kpole)
c
c     Compute polar interaction
c
         call epolar3_couple(dpui,ip,dpuk,kp,r2,pos,
     &            aewald,alsq2,alsq2n,pgamma,damp,f,
     &            off,shortheal,pscale,e,
     &            use_polarshortreal,.true.)
c
c     increment energy
c
         ep       = ep + tp2enr(e)
         if (pscale.eq.1.0) nep = nep-1
      end do
c
      end

c
c     ######################################################################
c     ##                                                                  ##
c     ##  subroutine eprecipgpu  --  PME recip space polarization energy  ##
c     ##                                                                  ##
c     ######################################################################
c
c
c     "eprecip" evaluates the reciprocal space portion of particle
c     mesh Ewald summation energy due to dipole polarization
c
c     literature reference:
c
c     C. Sagui, L. G. Pedersen and T. A. Darden, "Towards an Accurate
c     Representation of Electrostatics in Classical Force Fields:
c     Efficient Implementation of Multipolar Interactions in
c     Biomolecular Simulations", Journal of Chemical Physics, 120,
c     73-87 (2004)
c
c     modifications for nonperiodic systems suggested by Tom Darden
c     during May 2007
c
c
      subroutine eprecipgpu
      use atmlst
      use atoms
      use bound
      use boxes
      use chgpot
      use domdec
      use energi
      use ewald
      use fft
      use inform ,only: deb_Path
      use math
      use mpole
      use pme
      use polar
      use polpot
      use potent
      use utilgpu,only:rec_queue
      use mpi
      use timestat
      implicit none
      integer ierr,iipole,proc
      integer status(MPI_STATUS_SIZE),tag,commloc
      integer nprocloc,rankloc
      integer i,j,k,iglob
      integer k1,k2,k3
      integer m1,m2,m3
      integer ntot,nff
      integer nf1,nf2,nf3
      real(r_p) e
      real(t_p) f,h1,h2,h3
      real(t_p) volterm,denom
      real(t_p) hsq,expterm
      real(t_p) term,pterm
      real(t_p) struc2
      real(t_p) a(3,3),ftc(10,10)
      real(t_p) fuind
c
      if (aewald .lt. 1.0d-6)  return
c
      if (deb_Path) write(*,'(2x,a)') 'eprecipgpu'
      call timer_enter( timer_eprecip )

      if (use_pmecore) then
        nprocloc = nrec
        rankloc  = rank_bis
        commloc  =  comm_rec
      else
        nprocloc = nproc
        rankloc  = rank
        commloc  = MPI_COMM_WORLD
      end if
c
c     return if the Ewald coefficient is zero
c
      f = electric / dielec
!$acc enter data create(a,e) async(rec_queue)
c
c     convert Cartesian induced dipoles to fractional coordinates
c
!$acc data async(rec_queue)
!$acc&     present(ipole,fphirec,polerecglob,uind,qgrid2in_2d,
!$acc&   istart2,jstart2,kstart2,use_bounds,eprec)

!$acc serial async(rec_queue) present(e,a)
      e = 0.0_re_p
!$acc end serial

!$acc parallel loop async(rec_queue) present(a)
      do i = 1, 3
         a(1,i) = real(nfft1,t_p) * recip(i,1)
         a(2,i) = real(nfft2,t_p) * recip(i,2)
         a(3,i) = real(nfft3,t_p) * recip(i,3)
      end do

!$acc parallel loop collapse(2) async(rec_queue)
!$acc&         present(e,a)
      do i = 1, npolerecloc
         do j = 1, 3
            iipole = polerecglob(i)
            fuind = a(j,1)*uind(1,iipole) 
     &            + a(j,2)*uind(2,iipole)
     &            + a(j,3)*uind(3,iipole)
            e     = e + fuind*fphirec(j+1,i)
         end do
      end do

!$acc serial present(e) async(rec_queue)
      e     = 0.5_re_p * electric * e
      eprec = eprec + e
c
c     account for zeroth grid point for nonperiodic system
c
      if ((istart2(rankloc+1).eq.1).and.(jstart2(rankloc+1).eq.1)
     &   .and.(kstart2(rankloc+1).eq.1)) then
         if (.not. use_bounds) then
            expterm = 0.5_re_p * real(pi,r_p) / xbox
            struc2  = qgrid2in_2d(1,1,1,1,1)**2 +
     &                qgrid2in_2d(2,1,1,1,1)**2
            e       = f * expterm * struc2
            eprec   = eprec + e
         end if
      end if
!$acc end serial

!$acc end data
!$acc exit data delete(a,e) async(rec_queue)
      call timer_exit( timer_eprecip )
      end
