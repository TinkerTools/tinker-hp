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
#include "tinker_precision.h"
#include "tinker_types.h"
      module echarge3gpu_inl
        use tinTypes , only: real3,real3_red,rpole_elt
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
      use vec
      use vec_elec
      use vec_charge

      implicit none
      integer i,k,ii,kk,ksave
      integer iglob,iichg,inl
      integer nnchg,nnchg1,nnchg2
      integer nn12,nn13,nn14,ntot
      integer nect
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
!DIR$ ATTRIBUTES ALIGN:64::iichgvec,iglobvec
      integer iichgvec(nionlocloop),iglobvec(nionlocloop)
!DIR$ ATTRIBUTES ALIGN:64::ivec
      integer ivec(nionlocloop)
!DIR$ ATTRIBUTES ALIGN:64::pchgvec
      real(t_p) pchgvec(nionlocloop)
      character*10 mode
      parameter(half=0.5_ti_p)
      parameter(one=1.0_ti_p)

      if(rank.eq.0.and.tinkerdebug) write (*,*) 'echarge3cgpu'
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

!$acc data create(e) present(ec,ec_r,nec,nec_) async
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
             e = e + fs * pchg(iichg)**2
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
      !print*,'echarge3gpu',nec_,nec,nionloc
       ec = ec + e + enr2en( ec_r )
      nec = int(nec_) + nionloc + nec
!$acc end serial

!$acc end data

      call timer_exit(timer_echarge)
      end

      subroutine ecreal3dgpu
      use action ,only: nec_
      use atmlst
      use atoms
      use bound
      use boxes
      use charge
      use chgpot
      use couple
      use deriv
      use domdec
      use echarge3gpu_inl
      use energi  ,only:ec=>ec_r
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
      use utilgpu
      use virial
      use mpi
      implicit none
      integer i,j,k,iichg,iglob
      integer ii,kkk,kglob,kkchg
      real(t_p) e
      real(t_p) f,fi,fik,pchgk
      real(t_p) r,r2,rew
      real(t_p) rb,rb2

      real(t_p) xi,yi,zi
      real(t_p) xr,yr,zr

      real(t_p),parameter:: scale_f=1.0

      character*10 mode

      if (deb_path) write(*,'(2x,a)') 'ecreal3dgpu'
c
c     set conversion factor, cutoff and switching coefficients
c
      f    = electric / dielec
      mode = 'ewald'
      call switch (mode)
c
c     compute the real space ewald energy and first derivatives
c
!$acc parallel loop vector_length(32) async(dir_queue)
!$acc&         present(chgglobnl,iion,loc,x,y,z,pchg,nelst,elst,
!$acc&    ec,nec_)
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
#endif
               fik   = fi*pchg(kkchg)
               k     = loc(kglob)
               call charge3_couple(r2,xr,yr,zr,ebuffer
     &                            ,fik,aewald,scale_f,e,0)
c
c     increment the overall energy and derivative expressions
c
               ec       = ec + tp2enr(e)
               nec_     = nec_ + 1
            end if
         end do
      end do

      call ecreal3_scaling
      end

#ifdef _CUDA
      subroutine ecreal3d_cu
      use action ,only: nec_,nec
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
      use energi ,only: ec=>ec_r
      use ewald
      use iounit
      use inform
      use inter
      use math
      use molcul
      use neigh  , iion_s=>celle_glob,chg_s=>celle_chg
     &           , loc_s=>celle_loc, ieblst_s=>ieblst, eblst_s=>eblst
     &           , x_s=>celle_x, y_s=>celle_y, z_s=>celle_z
      use potent
      use shunt
      use timestat
      use usage
      use utilcu
      use utilgpu ,only: def_queue,dir_queue,nred_buff
     &            , zero_en_red_buffer,reduce_energy_action
     &            , ered_buff=>ered_buf1,dir_stream,def_stream
      use virial
      use mpi
      implicit none
      integer i,j,k,iichg,iglob
      integer ii,kkk,kglob,kkchg
      real(t_p) e
      real(t_p) f,fi,fik,pchgk
      real(t_p) r,r2,rew
      real(t_p) rb,rb2

      real(t_p) xi,yi,zi
      real(t_p) xr,yr,zr
      real(t_p) p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend

      logical,save:: first_in=.true.
      logical,parameter:: dyn_gS=.true.
      integer,save:: gS

      real(t_p),parameter:: scale_f=1.0

      character*10 mode

      if (deb_path) write(*,'(2x,a)') 'ecreal3d_cu'
c
c     set conversion factor, cutoff and switching coefficients
c
      f    = electric / dielec
      mode = 'ewald'
      call switch (mode)

      if (first_in) then
         call cudaMaxGridSize("ecreal3d_core_cu",gS)
         first_in = .false.
      end if
      if(dyn_gS) gs = nionlocnlb2_pair/4

      p_xbeg = xbegproc(rank+1)
      p_xend = xendproc(rank+1)
      p_ybeg = ybegproc(rank+1)
      p_yend = yendproc(rank+1)
      p_zbeg = zbegproc(rank+1)
      p_zend = zendproc(rank+1)

      def_stream = dir_stream
      def_queue  = dir_queue

      call zero_en_red_buffer(def_queue)
      call set_ChgData_CellOrder(.false.)

c
c     compute the real space ewald energy and first derivatives
c
!$acc host_data use_device(iion_s,chg_s,loc_s,ieblst_s,eblst_s,
!$acc&   x_s,y_s,z_s,pchg,ered_buff,nred_buff
!$acc&    )
      call ecreal3d_core_cu<<<gS,BLOCK_DIM,0,def_stream>>>
     &     ( iion_s,chg_s,loc_s,ieblst_s
     &     , eblst_s(2*nionlocnlb_pair+1)
     &     , nionlocnlb,nionlocnlb2_pair,nionbloc,n
     &     , x_s,y_s,z_s,pchg
     &     , off2,f,aewald, ebuffer
     &     , ered_buff, nred_buff
     &     , p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
#ifdef TINKER_DEBUG
#endif
     &     )
      call check_launch_kernel(" ecreal1d_core_cu")
!$acc end host_data

      call reduce_energy_action(ec,nec,ered_buff,def_queue)

      call ecreal3_scaling
      end
#endif
c
c     Scaling interaction correction subroutines
c
      subroutine ecreal3_scaling
      use action
      use atmlst
      use atoms
      use bound
      use boxes
      use charge
      use chgpot
      use couple
      use deriv
      use domdec
      use echarge3gpu_inl
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
      character*10 mode

      if (deb_Path) write(*,'(3x,a)') 'ecreal3_scaling'
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
!$acc&    ec,nec_)
      do ii = 1, n_cscale
         iglob = ccorrect_ik(ii,1)
         kglob = ccorrect_ik(ii,2)
         scale_f =   ccorrect_scale(2*ii+1)
         fik     = f*ccorrect_scale(2*ii+2)
         xi    = x(iglob)
         yi    = y(iglob)
         zi    = z(iglob)

         !compute the energy contribution for this interaction
         xr    = xi - x(kglob)
         yr    = yi - y(kglob)
         zr    = zi - z(kglob)

         !find energy for interactions within real space cutoff
         call image_inl (xr,yr,zr)
         r2 = xr*xr + yr*yr + zr*zr

         if (r2 .le. off2) then
            i     = loc(iglob)
            k     = loc(kglob)
            call charge3_couple(r2,xr,yr,zr,ebuffer
     &                         ,fik,aewald,scale_f
     &                         ,e,1)
            ec    = ec  + tp2enr(e)
            if (scale_f.eq.-1.0_ti_p) nec_=nec_-1
         end if
      end do

      end
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
      use interfaces,only: grid_pchg_site_p
      use math
      use mpi
      use pme
      use pme1
      use potent
      use timestat
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
      integer,dimension(nproc*proc)::req(:),reqbcast(:)
      logical,save:: f_in=.true.
c
      if (rank.eq.0.and.tinkerdebug)
     &   write(*,'(2x,a)') 'ecrecipgpu'
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return

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
!$acc&  octahedron,ec)
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
c
c     MPI : begin sending
c
      if (nrec_send.gt.0) then
!$acc wait
      end if
      call timer_enter( timer_recreccomm )
      do i = 1, nrec_send
         proc = prec_send(i)
         tag = nprocloc*prec_send(i) + rankloc + 1
!$acc host_data use_device(qgridin_2d)
         call MPI_ISEND(qgridin_2d(1,1,1,1,i+1),
     $        2*n1mpimax*n2mpimax*n3mpimax,MPI_TPREC,proc,tag,
     $        commloc,req(tag),ierr)
!$acc end host_data
      end do
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
     $        qgridin_2d(1,1,1,1,1),
     $        qgridmpi(1,1,1,1,i),qgridin_2d(1,1,1,1,1))
      end do
      call timer_exit ( timer_recreccomm,quiet_timers )
c
c     perform the 3-D FFT forward transformation
c
#ifdef _OPENACC
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
                 ec = ec + e
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
            ec = ec + e
         end if
      end if
!$acc end serial
      call timer_exit ( timer_scalar,quiet_timers )

!$acc end data
      if(f_in) f_in = .false.

      end
