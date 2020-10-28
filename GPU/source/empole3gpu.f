c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine empole3  --  mpole/polar energy & analysis  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "empole3" calculates the electrostatic energy due to
c     atomic multipole and dipole polarizability interactions,
c     and partitions the energy among the atoms
c
c
#include "tinker_precision.h"
      module empole3gpu_inl
        use utilgpu , only: real3,real3_red,rpole_elt
        include "erfcore_data.f.inc"
        contains
#include "image.f.inc"
#if defined(SINGLE) | defined(MIXED)
        include "erfcscore.f.inc"
#else
        include "erfcdcore.f.inc"
#endif
#include "switch_respa.f.inc"
#include "pair_mpole1.f.inc"
      end module

      subroutine empole3gpu
      use energi
      use potent
      use mpi

      implicit none
c
c     choose the method for summing over multipole interactions
c
      call empole3cgpu
c
      return
      end
c
c     ########################################################################
c     ##                                                                    ##
c     ##  subroutine empole3cgpu     --  Ewald multipole analysis via list  ##
c     ##                                                                    ##
c     ########################################################################
c
c
c     "empole3d" calculates the atomic multipole interaction energy
c     using particle mesh Ewald summation and a neighbor list, and
c     partitions the energy among the atoms
c
c
      subroutine empole3cgpu
      use sizes
      use action
      use analyz
      use atmlst
      use atoms
      use boxes
      use chgpot
      use domdec
      use energi
      use ewald
      use inform     ,only: deb_Path
      use interfaces ,only: reorder_nblist,commpoleglob
     &               ,emrealshortlong3d_p,emreal3d_p
      use math
      use mpole
      use mpi
      use neigh
      use potent
      use shunt
      use utilgpu
      use timestat
      implicit none
      integer k,i,ii,ierr
      integer iipoledefault
      integer iipole
      integer iglob
      real(t_p) f,emdir
      real(t_p) term,fterm
      real(t_p) xd,yd,zd
      real(t_p) cii,ci,dix,diy,diz
      real(t_p) qixx,qixy,qixz,qiyy,qiyz,qizz
      real(t_p) dii,qii,e
c
      if (npole .eq. 0)  return
      if(deb_Path) write(*,*) 'empole3cgpu'
c
c     zero out the multipole and polarization energies
c
!$acc serial present(nem,em,emrec,nem_) async
      nem   = 0
      em    = 0.0_ti_p
      emrec = 0.0_ti_p
      nem_  = 0.0
!$acc end serial
cold  aem   = 0_ti_p

c
c     set the energy unit conversion factor
c
      f = electric / dielec
c
c     Communicate poleglob
c
#ifdef _OPENACC
      call orderPole
#endif
      if (mlst2_enable) call set_ElecData_cellOrder(.false.)
c
c     check the sign of multipole components at chiral sites
c
      call chkpolegpu(.false.)
c
c     rotate the multipole components into the global frame
c
      call rotpolegpu
c
c     Reorder neigbhor list
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.le.ndir-1))
     &   then
         if (use_mreal) then
            if (mlst_enable.and.use_mpoleshortreal) then
               call switch('SHORTEWALD')
               call reorder_nblist(shortelst,nshortelst,nshortelstc
     &                    ,npolelocnl,off2,ipole,poleglobnl)
            end if
            if (mlst_enable) then
               call switch('EWALD     ')
               call reorder_nblist(elst,nelst,nelstc
     &                    ,npolelocnl,off2,ipole,poleglobnl)
            end if
         end if
      end if

#ifdef _OPENACC
      ! Start stream overlapping
      if (dir_queue.ne.rec_queue) then
         call stream_wait_async(dir_stream,rec_stream,dir_event)
      end if
#endif
c
c     compute the real space part of the Ewald summation
c
      def_queue = dir_queue
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.le.ndir-1))
     $   then
        call timer_enter (timer_real )
        if (use_mreal) then
           if (use_mpoleshortreal) then
              call emrealshortlong3d_p
           else if (use_mpolelong) then
              call emrealshortlong3d_p
           else
              call emreal3d_p
           end if
        end if
c
c     compute the self-energy part of the Ewald summation
c
        if (use_mself) then
           term  = 2.0_ti_p * aewald * aewald
           fterm = -f    * aewald / sqrtpi

!$acc parallel loop async(def_queue)
!$acc&         default(present) present(em,nem_)
           do ii = 1, npoleloc
              iipole = poleglob(ii)
              iglob  = ipole(iipole)
              i      = loc(iglob)
              ci     = rpole(1,iipole)
              dix    = rpole(2,iipole)
              diy    = rpole(3,iipole)
              diz    = rpole(4,iipole)
              qixx   = rpole(5,iipole)
              qixy   = rpole(6,iipole)
              qixz   = rpole(7,iipole)
              qiyy   = rpole(9,iipole)
              qiyz   = rpole(10,iipole)
              qizz   = rpole(13,iipole)
              cii    = ci*ci
              dii    = dix*dix + diy*diy + diz*diz
              qii    = 2.0_ti_p*(qixy*qixy+qixz*qixz+qiyz*qiyz)
     &                         + qixx*qixx+qiyy*qiyy+qizz*qizz
              e   = fterm * (cii + term*(dii/3.0_ti_p+
     &                     2.0_ti_p*term*qii/5.0_ti_p))
              em  = em + e
              nem_ = nem_ + 1.0
           end do
c
c          compute the cell dipole boundary correction term
c
           if (boundary .eq. 'VACUUM') then
              xd = 0.0_ti_p
              yd = 0.0_ti_p
              zd = 0.0_ti_p
!$acc parallel loop default(present) async(def_queue)
              do ii = 1, npoleloc
                 iipole = poleglob(ii)
                 iglob  = ipole(iipole)
                 dix    = rpole(2,iipole)
                 diy    = rpole(3,iipole)
                 diz    = rpole(4,iipole)
                 xd     = xd + dix + rpole(1,iipole)*x(iglob)
                 yd     = yd + diy + rpole(1,iipole)*y(iglob)
                 zd     = zd + diz + rpole(1,iipole)*z(iglob)
              end do
              call MPI_ALLREDUCE(MPI_IN_PLACE,xd,1,MPI_TPREC,MPI_SUM,
     $             COMM_TINKER,ierr)
              call MPI_ALLREDUCE(MPI_IN_PLACE,yd,1,MPI_TPREC,MPI_SUM,
     $             COMM_TINKER,ierr)
              call MPI_ALLREDUCE(MPI_IN_PLACE,zd,1,MPI_TPREC,MPI_SUM,
     $             COMM_TINKER,ierr)
              if (rank.eq.0) then
                term = (2.0_ti_p/3.0_ti_p) * f * (pi/volbox)
                e  = term * (xd*xd+yd*yd+zd*zd)
!$acc serial async(def_queue) copyin(e) present(em,nem_)
                em = em + e
                nem_ = nem_ + 1
!$acc end serial
              end if
           end if
        end if
        call timer_exit( timer_real,quiet_timers )
      end if
c
c     compute the reciprocal space part of the Ewald summation
c
      def_queue = rec_queue
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.gt.ndir-1))
     &   then
         call timer_enter( timer_rec )
         call emrecipgpu
         call timer_exit( timer_rec,quiet_timers )
      end if
c
c     Finalize async overlapping
c
#ifdef _OPENACC
      if (dir_queue.ne.rec_queue) then
         call stream_wait_async(dir_stream,rec_stream)
      end if
#endif
c
c     Sum both contribution
c
!$acc serial present(em,emrec,nem,nem_) async(rec_queue)
      em  = em + emrec
      nem = nem + int(nem_)
!$acc end serial
c
      end
c
c     ######################################################################
c     ##                                                                  ##
c     ##  subroutine emreal3dgpu  --  real space mpole analysis via list  ##
c     ##                                                                  ##
c     ######################################################################
c
c
c     "emreal3d" evaluates the real space portion of the Ewald sum
c     energy due to atomic multipole interactions, and partitions
c     the energy among the atoms using a pairwise neighbor list
c
c     literature reference:
c
c     W. Smith, "Point Multipoles in the Ewald Summation (Revisited)",
c     CCP5 Newsletter, 46, 18-30, 1998  (see http://www.ccp5.org/)
c
c
      subroutine emreal3dgpu
      use action ,only:nem_
      use atmlst ,only:poleglobnl
      use atoms  ,only:x,y,z,n
      use chgpot ,only:electric,dielec
      use domdec ,only:rank,nbloc,loc
      use empole3gpu_inl ,only: image_inl,mpole3_couple
      use energi ,only:em
      use ewald  ,only:aewald
      use inform ,only:deb_Path
      use interfaces,only: emreal3_correct_interactions
      use math   ,only:sqrtpi
      !use mpi
      use mplpot ,only:n_mscale,mcorrect_ik,mcorrect_scale
      use mpole  ,only:rpole,ipole,polelocnl,npolelocnl,npole
      use neigh  ,only:nelst,elst
      use shunt  ,only:off2
      use utilgpu,only:dir_queue,rec_queue,def_queue,warning,
     &                 real3,real3_red,rpole_elt,maxscaling
      use timestat
      implicit none

      integer i,j,k,iglob,kglob,kbis
      integer ii,kk,iipole,kkpole
      integer nnelst
#ifdef TINKER_DEBUG
#endif
      real(t_p) zero
      real(t_p) r2,f,e
      real(t_p) alsq2,alsq2n
      real(t_p) xi,yi,zi
      real(t_p) xr,yr,zr
      type(rpole_elt) ip,kp
      type(real3_red) frc_r
      type(real3) ttmi,ttmk,frc
      real(t_p) mscale
      parameter(zero=0.0)
      character*10 mode

      if (npole .eq. 0)  return
      if(deb_Path) write(*,'(3x,a)') 'emreal3dgpu'

#ifdef TINKER_DEBUG
#endif
c
c     set conversion factor, cutoff and switching coefficients
c
      call timer_enter( timer_emreal )
      mode   = 'EWALD'
      call switch (mode)

      f      = electric / dielec
      alsq2  = 2 * aewald**2
      alsq2n = zero
      if (aewald .gt. zero)
     &   alsq2n = 1.0 / (sqrtpi*aewald)
c
c     compute the real space portion of the Ewald summation
c
!$acc parallel loop gang vector_length(32)
!$acc&         present(nem_,em,poleglobnl,ipole,loc,
!$acc&  x,y,z,rpole,elst,nelst)
!$acc&         private(ip)
#ifdef TINKER_DEBUG
#endif
!$acc&         async(def_queue)
      do ii = 1, npolelocnl
         iipole = poleglobnl(ii)
         iglob  = ipole(iipole)
         i      = loc(iglob)
         nnelst = nelst(ii)
         if (nnelst.eq.0) cycle
         xi     = x(iglob)
         yi     = y(iglob)
         zi     = z(iglob)
         ip%c   = rpole(01,iipole)
         ip%dx  = rpole(02,iipole)
         ip%dy  = rpole(03,iipole)
         ip%dz  = rpole(04,iipole)
         ip%qxx = rpole(05,iipole)
         ip%qxy = rpole(06,iipole)
         ip%qxz = rpole(07,iipole)
         ip%qyy = rpole(09,iipole)
         ip%qyz = rpole(10,iipole)
         ip%qzz = rpole(13,iipole)
c
c     evaluate all sites within the cutoff distance
c
!$acc loop vector private(kp)
         do kk = 1, nnelst
            kkpole = elst(kk,ii)
            kglob  = ipole(kkpole)
            kbis   = loc(kglob)
            xr     = x(kglob) - xi
            yr     = y(kglob) - yi
            zr     = z(kglob) - zi
            call image_inl(xr,yr,zr)
            r2     = xr*xr + yr*yr + zr*zr
            if (r2 .gt. off2) cycle       !apply cutoff

            kp%c     = rpole( 1,kkpole)
            kp%dx    = rpole( 2,kkpole)
            kp%dy    = rpole( 3,kkpole)
            kp%dz    = rpole( 4,kkpole)
            kp%qxx   = rpole( 5,kkpole)
            kp%qxy   = rpole( 6,kkpole)
            kp%qxz   = rpole( 7,kkpole)
            kp%qyy   = rpole( 9,kkpole)
            kp%qyz   = rpole(10,kkpole)
            kp%qzz   = rpole(13,kkpole)

            ! compute mpole one interaction
            call mpole3_couple(r2,xr,yr,zr,ip,kp,zero,
     &                         aewald,f,alsq2n,alsq2,e,.false.)

            ! update energy
            em    =  em + e
            nem_  = nem_+ 1

         end do
      end do

      call emreal3_correct_interactions
      call timer_exit( timer_emreal )

      end

      subroutine emrealshortlong3d
      use action ,only:nem_
      use atmlst ,only:poleglobnl
      use atoms  ,only:x,y,z,n
      use cutoff ,only:mpoleshortcut,shortheal
      use chgpot ,only:electric,dielec
      use domdec ,only:rank,nbloc,loc
      use empole3gpu_inl ,only: image_inl,mpole3_couple_shortlong
      use energi ,only:em
      use ewald  ,only:aewald
      use inform ,only:deb_Path
      use interfaces,only: emreal3_correct_interactions_shortlong
     &              ,long_mode,short_mode
      use math   ,only:sqrtpi
      !use mpi
      use mplpot ,only:n_mscale,mcorrect_ik,mcorrect_scale
      use mpole  ,only:rpole,ipole,polelocnl,npolelocnl,npole
      use neigh  ,only:nelst,elst,shortelst,nshortelst
      use potent ,only:use_mpoleshortreal,use_mpolelong
      use shunt  ,only:off,off2
      use timestat
      use tinheader,only: ti_p
      use utilgpu,only:dir_queue,rec_queue,def_queue,warning,
     &                 real3,real3_red,rpole_elt,maxscaling
      implicit none

      integer i,j,k,iglob,kglob,kbis
      integer ii,kk,iipole,kkpole
      integer nnelst,mode
      integer,pointer::lst(:,:),nlst(:)
#ifdef TINKER_DEBUG
#endif
      real(t_p) zero,r_cut
      real(t_p) r2,f,e
      real(t_p) alsq2,alsq2n
      real(t_p) xi,yi,zi
      real(t_p) xr,yr,zr
      type(rpole_elt) ip,kp
      real(t_p) mscale
      real(t_p) mpoleshortcut2
      parameter(zero=0.0)

      if (npole .eq. 0)  return
      if(deb_Path)
     &   write(*,'(3x,a)') 'emrealshortlong3d'

#ifdef TINKER_DEBUG
#endif
c
c     set conversion factor, cutoff and switching coefficients
c
      f      = electric / dielec
      alsq2  = 2 * aewald**2
      alsq2n = zero
      if (aewald .gt. zero)
     &   alsq2n = 1.0 / (sqrtpi*aewald)
      call timer_enter( timer_emreal )

      ! Configure data to be use in next loop
      if (use_mpoleshortreal) then
         call switch('SHORTEWALD')
         mode  = short_mode
         r_cut = off
         mpoleshortcut2 = 0.0_ti_p
          lst =>  shortelst
         nlst => nshortelst
      else
         call switch('EWALD     ')
         mode  = long_mode
         r_cut = mpoleshortcut
         mpoleshortcut2 = (mpoleshortcut-shortheal)**2
          lst =>  elst
         nlst => nelst
      end if
!$acc enter data attach(nlst,lst) async(def_queue)
c
c     compute the real space portion of the Ewald summation
c
!$acc parallel loop gang vector_length(32)
!$acc&         present(em,nem_,poleglobnl,ipole,loc,x,y,z,rpole,
!$acc&  lst,nlst)
!$acc&         private(ip)
#ifdef TINKER_DEBUG
#endif
!$acc&         async(def_queue)
      do ii = 1, npolelocnl
         iipole = poleglobnl(ii)
         iglob  = ipole(iipole)
         i      = loc(iglob)
         nnelst = nlst(ii)
         if (nnelst.eq.0) cycle
         xi     = x(iglob)
         yi     = y(iglob)
         zi     = z(iglob)
         ip%c   = rpole(01,iipole)
         ip%dx  = rpole(02,iipole)
         ip%dy  = rpole(03,iipole)
         ip%dz  = rpole(04,iipole)
         ip%qxx = rpole(05,iipole)
         ip%qxy = rpole(06,iipole)
         ip%qxz = rpole(07,iipole)
         ip%qyy = rpole(09,iipole)
         ip%qyz = rpole(10,iipole)
         ip%qzz = rpole(13,iipole)
c
c     evaluate all sites within the cutoff distance
c
!$acc loop vector private(kp)
         do kk = 1, nnelst
            kkpole = lst(kk,ii)
            kglob  = ipole(kkpole)
            kbis   = loc(kglob)
            xr     = x(kglob) - xi
            yr     = y(kglob) - yi
            zr     = z(kglob) - zi
            call image_inl(xr,yr,zr)
            r2     = xr*xr + yr*yr + zr*zr
            if (r2.lt.mpoleshortcut2.or.r2.gt.off2) cycle       !apply cutoff

            kp%c     = rpole( 1,kkpole)
            kp%dx    = rpole( 2,kkpole)
            kp%dy    = rpole( 3,kkpole)
            kp%dz    = rpole( 4,kkpole)
            kp%qxx   = rpole( 5,kkpole)
            kp%qxy   = rpole( 6,kkpole)
            kp%qxz   = rpole( 7,kkpole)
            kp%qyy   = rpole( 9,kkpole)
            kp%qyz   = rpole(10,kkpole)
            kp%qzz   = rpole(13,kkpole)

            ! compute mpole one interaction
            call mpole3_couple_shortlong(r2,xr,yr,zr,ip,kp,zero,
     &                    r_cut,shortheal,aewald,f,alsq2n,alsq2,
     &                    e,.false.,mode)

            ! update energy
            em       = em  + e
            nem_     = nem_+ 1

#ifdef TINKER_DEBUG
#endif
         end do
      end do
!$acc exit data detach(nlst,lst) async(def_queue)

      call emreal3_correct_interactions_shortlong
      call timer_exit( timer_emreal )

      end

      ! CUDA Fortran version of emreal1c using C2 nblist
      subroutine emreal3d_cu
      use action ,only:nem
      use atmlst ,only:poleglobnl
      use atoms  ,only:x,y,z,n
      use chgpot ,only:electric,dielec
      use cell
      use domdec  ,only: xbegproc,ybegproc,zbegproc
     &            ,nproc,rank,xendproc,yendproc,zendproc
     &            ,nbloc,loc
#ifdef _CUDA
      use empole1cu ,only: emreal3_cu
#endif
      use energi ,only:em
      use ewald  ,only:aewald
      use inform ,only:deb_Path
      use interfaces,only: emreal_correct_interactions
      use math   ,only:sqrtpi
      !use mpi
      use mplpot ,only:n_mscale,mcorrect_ik,mcorrect_scale
      use mpole  ,only:rpole,ipole,polelocnl,npolelocnl,npole
     &           ,npolelocnlb,npolelocnlb_pair,npolebloc
     &           ,npolelocnlb2_pair
      use neigh  ,only:ipole_s=>celle_glob,pglob_s=>celle_pole
     &           , loc_s=>celle_loc, ieblst_s=>ieblst, eblst_s=>eblst
     &           , x_s=>celle_x, y_s=>celle_y, z_s=>celle_z
      use shunt  ,only:off2
#ifdef _CUDA
      use utilcu ,only:BLOCK_DIM,check_launch_kernel
#endif
      use utilgpu,only:dir_queue,rec_queue,def_queue,warning
     &           ,real3,real3_red,rpole_elt,RED_BUFF_SIZE
     &           ,ered_buff,nred_buff,reduce_energy_action
     &           ,zero_en_red_buffer
#ifdef  _OPENACC
     &           ,dir_stream,def_stream
#endif
      use timestat
      implicit none

      integer i
#ifdef TINKER_DEBUG
#endif
      real(t_p) r2,f,e
      real(t_p) alsq2,alsq2n
      real(t_p) p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
      logical,save:: first_in=.true.
      integer,save:: gS
      character*10 mode

      if (npole .eq. 0)  return
      if(deb_Path) write(*,'(2x,a)') 'emreal3d_cu'
      call timer_enter( timer_emreal )

#ifdef TINKER_DEBUG
#endif
c
c     set conversion factor, cutoff and switching coefficients
c
      mode   = 'EWALD'
      call switch (mode)

      f      = electric / dielec
      alsq2  = 2 * aewald**2
      alsq2n = 0
      if (aewald .gt. 0.0)
     &   alsq2n = 1.0 / (sqrtpi*aewald)

      p_xbeg = xbegproc(rank+1)
      p_xend = xendproc(rank+1)
      p_ybeg = ybegproc(rank+1)
      p_yend = yendproc(rank+1)
      p_zbeg = zbegproc(rank+1)
      p_zend = zendproc(rank+1)

      if (first_in) then
#ifdef _CUDA
         call cudaMaxGridSize("emreal3_cu",gS)
#endif
         first_in = .false.
      end if

#ifdef _CUDA
      def_stream = dir_stream
      call zero_en_red_buffer(def_queue)

!$acc host_data use_device(ipole_s,pglob_s,loc_s,ieblst_s,eblst_s
!$acc&    ,x_s,y_s,z_s,rpole,ered_buff,nred_buff
#ifdef TINKER_DEBUG
#endif
!$acc&    )

      !call CUDA kernel to compute the real space portion of the Ewald summation
      call emreal3_cu<<<gS,BLOCK_DIM,0,def_stream>>>
     &        ( ipole_s, pglob_s, loc_s, ieblst_s
     &        , eblst_s(2*npolelocnlb_pair+1)
     &        , npolelocnlb, npolelocnlb2_pair, npolebloc, n
     &        , x_s, y_s, z_s, rpole
     &        , off2, f, alsq2, alsq2n, aewald
     &        , ered_buff, nred_buff
     &        , p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
#ifdef TINKER_DEBUG
#endif
     &        )
      call check_launch_kernel(" emreal3_cu")

!$acc end host_data

      call reduce_energy_action(em,nem,def_queue)
#else
 100  format('emreal3c_core2 is a specific device routine !!',/,
     &       'you are not supposed to get inside with your compile ',
     &       'mode.')
      print 100
      call fatal
#endif

#ifdef TINKER_DEBUG
#endif

      call emreal3_correct_interactions
      call timer_exit( timer_emreal )

      end

      subroutine emrealshortlong3d_cu
      use action ,only:nem
      use atmlst ,only:poleglobnl
      use atoms  ,only:x,y,z,n
      use cutoff ,only:mpoleshortcut,shortheal
      use chgpot ,only:electric,dielec
      use cell
      use domdec  ,only: xbegproc,ybegproc,zbegproc
     &            ,nproc,rank,xendproc,yendproc,zendproc
     &            ,nbloc,loc
#ifdef _CUDA
      use empole1cu ,only: emrealshortlong3_cu
#endif
      use energi ,only:em
      use ewald  ,only:aewald
      use inform ,only:deb_Path
      use interfaces,only: emreal_correct_interactions_shortlong
     &              ,long_mode,short_mode
      use math   ,only:sqrtpi
      !use mpi
      use mplpot ,only:n_mscale,mcorrect_ik,mcorrect_scale
      use mpole  ,only:rpole,ipole,polelocnl,npolelocnl,npole
     &           ,npolelocnlb,npolelocnlb_pair,npolebloc
     &           ,npolelocnlb2_pair,nshortpolelocnlb2_pair
      use neigh  ,only:ipole_s=>celle_glob,pglob_s=>celle_pole
     &           , loc_s=>celle_loc, ieblst_s=>ieblst, eblst_s=>eblst
     &           , x_s=>celle_x, y_s=>celle_y, z_s=>celle_z
     &           , seblst_s=>shorteblst,iseblst_s=>ishorteblst
      use potent ,only:use_mpoleshortreal,use_mpolelong
      use shunt  ,only:off,off2
      use timestat
      use tinheader,only: ti_p
#ifdef _CUDA
      use utilcu ,only:BLOCK_DIM,check_launch_kernel
#endif
      use utilgpu,only:dir_queue,rec_queue,def_queue,warning
     &           ,real3,real3_red,rpole_elt,RED_BUFF_SIZE
     &           ,ered_buff,nred_buff,reduce_energy_action
     &           ,zero_en_red_buffer
#ifdef  _OPENACC
     &           ,dir_stream,def_stream
#endif
      implicit none

      integer i,mode
#ifdef TINKER_DEBUG
#endif
      real(t_p) r2,f,e
      real(t_p) alsq2,alsq2n
      real(t_p) p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
      real(t_p) mpoleshortcut2,r_cut
      logical,save:: first_in=.true.
      integer,save:: gS

      if (npole .eq. 0)  return
      if(deb_Path) 
     &   write(*,'(3x,a)') 'emrealshortlong3d_cu'
      call timer_enter( timer_emreal )

#ifdef TINKER_DEBUG
#endif
c
c     set conversion factor, cutoff and switching coefficients
c

      f      = electric / dielec
      alsq2  = 2 * aewald**2
      alsq2n = 0
      if (aewald .gt. 0.0)
     &   alsq2n = 1.0 / (sqrtpi*aewald)

      p_xbeg = xbegproc(rank+1)
      p_xend = xendproc(rank+1)
      p_ybeg = ybegproc(rank+1)
      p_yend = yendproc(rank+1)
      p_zbeg = zbegproc(rank+1)
      p_zend = zendproc(rank+1)

      if (use_mpoleshortreal) then
         call switch('SHORTEWALD')
         mode = short_mode
         mpoleshortcut2 = 0.0_ti_p
         r_cut = off
      else if (use_mpolelong) then
         call switch('EWALD     ')
         mode = long_mode
         mpoleshortcut2 = (mpoleshortcut-shortheal)**2
         r_cut = mpoleshortcut
      end if

      if (first_in) then
#ifdef _CUDA
        call cudaMaxGridSize("emrealshortlong3_cu",gS)
#endif
         first_in = .false.
      end if

#ifdef _CUDA
      def_stream = dir_stream
      call zero_en_red_buffer(def_queue)


      if (use_mpoleshortreal) then
!$acc host_data use_device(ipole_s,pglob_s,loc_s,iseblst_s,seblst_s
!$acc&    ,x_s,y_s,z_s,rpole,ered_buff,nred_buff
#ifdef TINKER_DEBUG
#endif
!$acc&    )
      !call CUDA kernel to compute the real space portion of the Ewald summation
      call emrealshortlong3_cu<<<gS,BLOCK_DIM,0,def_stream>>>
     &        ( ipole_s, pglob_s, loc_s, iseblst_s
     &        , seblst_s(2*npolelocnlb_pair+1)
     &        , npolelocnlb, nshortpolelocnlb2_pair, npolebloc, n
     &        , x_s, y_s, z_s, rpole
     &        , shortheal, r_cut, mpoleshortcut2, off2, f
     &        , alsq2, alsq2n, aewald, mode
     &        , ered_buff, nred_buff
     &        , p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
#ifdef TINKER_DEBUG
#endif
     &        )
!$acc end host_data

      else
!$acc host_data use_device(ipole_s,pglob_s,loc_s,ieblst_s,eblst_s
!$acc&    ,x_s,y_s,z_s,rpole,ered_buff,nred_buff
#ifdef TINKER_DEBUG
#endif
!$acc&    )
      !call CUDA kernel to compute the real space portion of the Ewald summation
      call emrealshortlong3_cu<<<gS,BLOCK_DIM,0,def_stream>>>
     &        ( ipole_s, pglob_s, loc_s, ieblst_s
     &        , eblst_s(2*npolelocnlb_pair+1)
     &        , npolelocnlb, npolelocnlb2_pair, npolebloc, n
     &        , x_s, y_s, z_s, rpole
     &        , shortheal, r_cut, mpoleshortcut2, off2, f
     &        , alsq2, alsq2n, aewald, mode
     &        , ered_buff, nred_buff
     &        , p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
#ifdef TINKER_DEBUG
#endif
     &        )
!$acc end host_data
      end if
      call check_launch_kernel(" emrealshortlong3_cu")

      call reduce_energy_action(em,nem,def_queue)
#else
 100  format('emrealshortlong3d_cu is a specific device routine !!',
     &       /,'you are not supposed to get inside with your compile ',
     &       'mode.')
      print 100
      call fatal
#endif

#ifdef TINKER_DEBUG
#endif

      call emreal3_correct_interactions_shortlong
      call timer_exit( timer_emreal )
      end


      ! Procede to corrections of scaling interactions
      subroutine
     &   emreal3_correct_interactions
      use action ,only:nem
      use atoms  ,only:x,y,z
      use chgpot ,only:electric,dielec
      use domdec ,only:rank,nbloc,loc
      use empole3gpu_inl ,only: image_inl,mpole3_couple
      use energi ,only:em
      use ewald  ,only:aewald
      use inform ,only:deb_Path
      use math   ,only:sqrtpi
      use mplpot ,only:n_mscale,mcorrect_ik,mcorrect_scale
      use mpole  ,only:rpole,ipole,polelocnl,npolelocnl
      use shunt  ,only:off2
      use tinheader,only: ti_p
      use utilgpu,only:real3,real3_red,rpole_elt,def_queue
      implicit none

      integer   i,j,k,iglob,kglob,kbis
      integer   ii,kk,iipole,kkpole
      real(t_p) zero
      real(t_p) r2,f,e
      real(t_p) alsq2,alsq2n
      real(t_p) xr,yr,zr
      type(rpole_elt) ip,kp
      real(t_p) mscale
      logical   corre
      parameter(zero=0.0)

      if (deb_Path)
     &   write(*,'(2x,a)') 'emreal3_correct_interactions'
c
c     set conversion factor, cutoff and switching coefficients
c
      f      = electric / dielec
      alsq2  = 2 * aewald**2
      alsq2n = zero
      if (aewald .gt. zero)
     &   alsq2n = 1.0 / (sqrtpi*aewald)
c
c     Apply correction to scaling factor interactions
c
!$acc parallel loop vector_length(32)
!$acc&         present(em,nem,ipole,loc,x,y,z,rpole
!$acc&     ,mcorrect_ik,mcorrect_scale)
!$acc&         private(ip,kp)
!$acc&         reduction(+:nem)
!$acc&         async(def_queue)
      do ii = 1, n_mscale
         iipole = mcorrect_ik(ii,1)
         kkpole = mcorrect_ik(ii,2)
         mscale = mcorrect_scale(ii)
         iglob  = ipole(iipole)
         kglob  = ipole(kkpole)
         i      = loc(iglob)
         kbis   = loc(kglob)

         xr     = x(kglob) - x(iglob)
         yr     = y(kglob) - y(iglob)
         zr     = z(kglob) - z(iglob)
         call image_inl(xr,yr,zr)
         r2     = xr*xr + yr*yr + zr*zr
         if (r2 .gt. off2) cycle       !apply cutoff

         ip%c   = rpole(01,iipole)
         ip%dx  = rpole(02,iipole)
         ip%dy  = rpole(03,iipole)
         ip%dz  = rpole(04,iipole)
         ip%qxx = rpole(05,iipole)
         ip%qxy = rpole(06,iipole)
         ip%qxz = rpole(07,iipole)
         ip%qyy = rpole(09,iipole)
         ip%qyz = rpole(10,iipole)
         ip%qzz = rpole(13,iipole)
         kp%c   = rpole(01,kkpole)
         kp%dx  = rpole(02,kkpole)
         kp%dy  = rpole(03,kkpole)
         kp%dz  = rpole(04,kkpole)
         kp%qxx = rpole(05,kkpole)
         kp%qxy = rpole(06,kkpole)
         kp%qxz = rpole(07,kkpole)
         kp%qyy = rpole(09,kkpole)
         kp%qyz = rpole(10,kkpole)
         kp%qzz = rpole(13,kkpole)

         call mpole3_couple(r2,xr,yr,zr,ip,kp,mscale,
     &                      aewald,f,alsq2n,alsq2,e,.true.)

         ! update energy
         em     = em  + e
         if (mscale.eq.1.0_ti_p) nem = nem - 1
      end do

      end

      ! Procede to corrections of scaling interactions
      subroutine
     &   emreal3_correct_interactions_shortlong
      use action ,only:nem
      use atoms  ,only:x,y,z
      use cutoff ,only:mpoleshortcut,shortheal
      use chgpot ,only:electric,dielec
      use domdec ,only:rank,nbloc,loc
      use empole3gpu_inl ,only: image_inl,mpole3_couple_shortlong
      use energi ,only:em
      use ewald  ,only:aewald
      use inform ,only: deb_Path
      use interfaces ,only:long_mode,short_mode
      use math   ,only:sqrtpi
      use mplpot ,only:n_mscale,mcorrect_ik,mcorrect_scale
      use mpole  ,only:rpole,ipole,polelocnl,npolelocnl
      use potent ,only:use_mpoleshortreal,use_mpolelong
      use shunt  ,only:off2,off
      use tinheader,only: ti_p
      use utilgpu,only:real3,real3_red,rpole_elt,def_queue
      implicit none

      integer   i,j,k,iglob,kglob,kbis
      integer   ii,kk,iipole,kkpole
      integer   mode
      real(t_p) zero
      real(t_p) r2,f,e
      real(t_p) alsq2,alsq2n
      real(t_p) xr,yr,zr
      type(rpole_elt) ip,kp
      real(t_p) mpoleshortcut2,r_cut
      real(t_p) mscale
      logical   corre
      parameter(zero=0.0)

      if (deb_Path)
     &   write(*,'(2x,a)') 'emreal3_correct_interactions_shortlong'
c
c     set conversion factor, cutoff and switching coefficients
c
      f      = electric / dielec
      alsq2  = 2 * aewald**2
      alsq2n = zero
      if (aewald .gt. zero)
     &   alsq2n = 1.0 / (sqrtpi*aewald)

      if (use_mpoleshortreal) then
         mpoleshortcut2 = 0.0_ti_p
         mode  = short_mode
         r_cut = off
      else if (use_mpolelong) then
         mpoleshortcut2 = (mpoleshortcut-shortheal)**2
         mode  = long_mode
         r_cut = mpoleshortcut
      else
         write(*,*) 'unknown mode for '
     &            ,'emreal3_correct_interactions_shortlong'
         call fatal
      end if

c
c     Apply correction to scaling factor interactions
c
!$acc parallel loop present(em,nem,ipole,loc,x,y,z,rpole,
!$acc&  mcorrect_ik,mcorrect_scale)
!$acc&         private(ip,kp)
!$acc&         reduction(+:nem)
!$acc&         async(def_queue)
      do ii = 1, n_mscale
         iipole = mcorrect_ik(ii,1)
         kkpole = mcorrect_ik(ii,2)
         mscale = mcorrect_scale(ii)
         iglob  = ipole(iipole)
         kglob  = ipole(kkpole)
         i      = loc(iglob)
         kbis   = loc(kglob)

         xr     = x(kglob) - x(iglob)
         yr     = y(kglob) - y(iglob)
         zr     = z(kglob) - z(iglob)
         call image_inl(xr,yr,zr)
         r2     = xr*xr + yr*yr + zr*zr
         if (r2.lt.mpoleshortcut2.or.r2.gt.off2) cycle  !apply cutoff

         ip%c   = rpole(01,iipole)
         ip%dx  = rpole(02,iipole)
         ip%dy  = rpole(03,iipole)
         ip%dz  = rpole(04,iipole)
         ip%qxx = rpole(05,iipole)
         ip%qxy = rpole(06,iipole)
         ip%qxz = rpole(07,iipole)
         ip%qyy = rpole(09,iipole)
         ip%qyz = rpole(10,iipole)
         ip%qzz = rpole(13,iipole)
         kp%c   = rpole(01,kkpole)
         kp%dx  = rpole(02,kkpole)
         kp%dy  = rpole(03,kkpole)
         kp%dz  = rpole(04,kkpole)
         kp%qxx = rpole(05,kkpole)
         kp%qxy = rpole(06,kkpole)
         kp%qxz = rpole(07,kkpole)
         kp%qyy = rpole(09,kkpole)
         kp%qyz = rpole(10,kkpole)
         kp%qzz = rpole(13,kkpole)

         ! compute mpole one interaction
         call mpole3_couple_shortlong(r2,xr,yr,zr,ip,kp,mscale,
     &                   r_cut,shortheal,aewald,f,alsq2n,alsq2,
     &                       e,.true.,mode)

         ! update energy
         em     = em  + e
         if (mscale.eq.1.0) nem = nem - 1
      end do

      end

c
c
c
c     ###################################################################
c     ##                                                               ##
c     ##  subroutine emrecipgpu  --  PME recip space multipole energy  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     "emrecipgpu" evaluates on device the reciprocal space portion of the particle
c     mesh Ewald energy due to atomic multipole interactions
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
      subroutine emrecipgpu
      use atoms
      use atmlst
      use bound
      use boxes
      use chgpot
      use domdec
      use energi
      use ewald
      use fft
      use inform    ,only: deb_Path
      use interfaces,only: torquegpu,fphi_mpole_site_p
     &              ,grid_mpole_site_p,bspline_fill_sitegpu
      use polar_temp,only: cmp=>fphid,fmp=>fphip !Use register pool
      use math
      use mpole
      use pme
      use potent
      use utilgpu ,only:dir_queue,rec_queue
     &            ,prmem_request
      use utils   ,only:set_to_zero1,rpole_scale
     &            ,rpole_ind_extract,comput_norm
     &            ,convolution_product
      use mpi
      use timestat
      implicit none
      integer ierr,iipole,proc,iglob
      integer status(MPI_STATUS_SIZE),tag,commloc
      integer i,j,k
      integer k1,k2,k3
      integer m1,m2,m3
      integer nff
      integer nf1,nf2,nf3
      integer nprocloc,rankloc
      real(t_p) e,r1,r2,r3
      real(t_p) f,h1,h2,h3
      real(t_p) volterm,denom
      real(t_p) hsq,expterm
      real(t_p) term,pterm
      real(t_p) struc2
      integer  reqsend(nproc), reqrec(nproc)
      integer req2send(nproc),req2rec(nproc)
      real(t_p) qgridmpi(2,n1mpimax,n2mpimax,n3mpimax,
     &                   nrec_recep)
c
      if (deb_Path)
     &   write(*,'(2x,a)') 'emrecipgpu'
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
c
      if (use_pmecore) then
        nprocloc = nrec
        rankloc  = rank_bis
        commloc  =  comm_rec
      else
        nprocloc = nproc
        rankloc  = rank
        commloc  = MPI_COMM_WORLD
      end if
      f = electric / dielec
c
c     dynamic allocation of global arrays
c
      j = max(npolerecloc,1)
      call prmem_request(cphirec,10,j,async=.true.)
      call prmem_request(fphirec,20,j,async=.true.)
      call prmem_request(cmp    ,10,j,async=.true.)
      call prmem_request(fmp    ,10,j,async=.true.)
!$acc enter data create(qgridmpi,e) async(rec_queue)

!$acc data present(polerecglob,ipole,rpole,kstart2,kend2,
!$acc&  jstart2,jend2,istart2,iend2,bsmod1,bsmod2,bsmod3,
!$acc&  repart,qgridin_2d,qgridout_2d,use_bounds,qfac_2d,
!$acc&  octahedron,thetai1,thetai2,thetai3,igrid,fmp,cmp,
!$acc&  fphirec,cphirec) async(rec_queue)
c
c     zero out the PME grid
c
      call timer_enter( timer_grid1 )
      call set_to_zero1(qgridin_2d,size(qgridin_2d),rec_queue)

!$acc parallel loop collapse(2) async(rec_queue)
      do i = 1,npolerecloc
         do j = 1,20
            if (j.lt.10)
     &      cphirec(j,i) = 0.0_ti_p
            fphirec(j,i) = 0.0_ti_p
         end do
      end do
c
c     copy the multipole moments into local storage areas
c
!$acc parallel loop collapse(2) async(rec_queue)
      do i = 1, npolerecloc
         do j = 1, 10
            iipole   = polerecglob(i)
            cmp(j,i) = rpole_scale(j)*rpole(rpole_ind_extract(j),iipole)
         end do
      end do
      call cmp_to_fmp_sitegpu(cmp,fmp)
      call bspline_fill_sitegpu
c
c     assign permanent multipoles to PME grid and perform
c     the 3-D FFT forward transformation
c
      call grid_mpole_site_p(fmp)
      call timer_exit( timer_grid1,quiet_timers )
c
c     MPI : Begin reception
c
      call timer_enter( timer_recreccomm )
      do i = 1, nrec_recep
         tag = nprocloc*rankloc + prec_recep(i) + 1
!$acc host_data use_device(qgridmpi)
         call MPI_IRECV(qgridmpi(1,1,1,1,i),2*n1mpimax*n2mpimax*
     $        n3mpimax,MPI_TPREC,prec_recep(i),tag,
     $        commloc,reqrec(i),ierr)
!$acc end host_data
      end do
c
c     MPI : begin sending
c
      do i = 1, nrec_send
!$acc wait(rec_queue)
         proc = prec_send(i)
         tag  = nprocloc*prec_send(i) + rankloc + 1
!$acc host_data use_device(qgridin_2d)
         call MPI_ISEND(qgridin_2d(1,1,1,1,i+1),
     $        2*n1mpimax*n2mpimax*n3mpimax,MPI_TPREC,
     $        proc,tag,commloc,reqsend(i),ierr)
!$acc end host_data
      end do
c
      do i = 1, nrec_recep
         call MPI_WAIT(reqrec(i),status,ierr)
      end do
      do i = 1, nrec_send
         call MPI_WAIT(reqsend(i),status,ierr)
      end do
c
c     do the reduction 'by hand'
c
      do i = 1, nrec_recep
         call aaddgpuAsync(2*n1mpimax*n2mpimax*n3mpimax,
     $        qgridin_2d(1,1,1,1,1),
     $        qgridmpi(1,1,1,1,i),qgridin_2d(1,1,1,1,1))
      end do
      call timer_exit( timer_recreccomm,quiet_timers )
c
c     Perform 3-D FFT forward transform
c
#ifdef _OPENACC
      call cufft2d_frontmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
     $     n3mpimax)
#else
      call fft2d_frontmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
     $     n3mpimax)
#endif
c
c     make the scalar summation over reciprocal lattice
c
      pterm = (pi/aewald)**2
      volterm = pi * volbox
      nff = nfft1 * nfft2
      nf1 = (nfft1+1) / 2
      nf2 = (nfft2+1) / 2
      nf3 = (nfft3+1) / 2
c
      call timer_enter( timer_scalar )
!$acc serial async(rec_queue)
      if ((istart2(rankloc+1).eq.1).and.(jstart2(rankloc+1).eq.1).and.
     $   (kstart2(rankloc+1).eq.1)) then
           qfac_2d(1,1,1) = 0.0_ti_p
      end if
!$acc end serial
c
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
                  expterm = expterm * (1.0_ti_p-cos(pi*xbox*sqrt(hsq)))
               else if (octahedron) then
                  if (mod(m1+m2+m3,2) .ne. 0)  expterm = 0.0_ti_p
               end if
            qfac_2d(k1-istart2(rankloc+1)+1,k2-jstart2(rankloc+1)+1,k3-
     $       kstart2(rankloc+1)+1) = expterm
            end if
 10         continue
          end do
        end do
      end do
c
c     account for zeroth grid point for nonperiodic system
c
!$acc serial async(rec_queue) present(emrec,e)
      if ((istart2(rankloc+1).eq.1).and.(jstart2(rankloc+1).eq.1)
     $   .and.(kstart2(rankloc+1).eq.1)) then
        if (.not. use_bounds) then
           expterm = 0.5_ti_p * pi / xbox
           struc2 = qgrid2in_2d(1,1,1,1,1)**2 +
     $       qgrid2in_2d(2,1,1,1,1)**2
           e  = f * expterm * struc2
           emrec = emrec + e
        end if
      end if
!$acc end serial
c
c     complete the transformation of the charge grid
c
      call convolution_product(qfac_2d,size(qgridout_2d),qgridout_2d
     &                        ,rec_queue)
      call timer_exit( timer_scalar,quiet_timers )
c
c     perform 3-D FFT backward transform and get potential
c
#ifdef _OPENACC
      call cufft2d_backmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
     $ n3mpimax)
#else
      call fft2d_backmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
     $ n3mpimax)
#endif
c
c     MPI : Begin reception
c
      call timer_enter( timer_recreccomm )
      do i = 1, nrec_send
         proc = prec_send(i)
         tag = nprocloc*rankloc + prec_send(i) + 1
!$acc host_data use_device(qgridin_2d)
         call MPI_IRECV(qgridin_2d(1,1,1,1,i+1),
     $        2*n1mpimax*n2mpimax*n3mpimax,MPI_TPREC,
     $        prec_send(i),tag,commloc,req2rec(i),ierr)
!$acc end host_data
      end do
c
c     MPI : begin sending
c
      do i = 1, nrec_recep
!$acc wait(rec_queue)
         tag = nprocloc*prec_recep(i) + rankloc + 1
!$acc host_data use_device(qgridin_2d)
         call MPI_ISEND(qgridin_2d,
     $        2*n1mpimax*n2mpimax*n3mpimax,MPI_TPREC,
     $        prec_recep(i),tag,commloc,req2send(i),ierr)
!$acc end host_data
      end do
c
      do i = 1, nrec_send
         call MPI_WAIT(req2rec(i),status,ierr)
      end do
      do i = 1, nrec_recep
         call MPI_WAIT(req2send(i),status,ierr)
      end do
      call timer_exit( timer_recreccomm,quiet_timers )

      call timer_enter( timer_grid2 )
      call fphi_mpole_site_p
      call timer_exit( timer_grid2,quiet_timers )
c
c     sum over multipoles and increment total multipole energy
c
!$acc serial async(rec_queue) present(e)
      e = 0.0_ti_p
!$acc end serial

!$acc parallel loop collapse(2) present(e) async(rec_queue)
      do i = 1, npolerecloc
         do k = 1, 10
            e = e + fmp(k,i)*fphirec(k,i)
         end do
      end do

!$acc serial async(rec_queue) present(emrec,e)
      e  = 0.5_ti_p * f * e
      emrec = emrec + e
!$acc end serial

!$acc end data
c
!$acc exit data delete(qgridmpi,e) async(rec_queue)
      end
