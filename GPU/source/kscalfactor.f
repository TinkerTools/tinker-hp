c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     workS[-p]  workspace required to build polar scaling with
c                low amount of memory
c     maxgangs   maximum gangs number used to process polar scaling
c                kernel
c     vecl       length of each gang
#include "tinker_precision.h"
      module kscalfactor_inl
        integer maxgangs,vecl,worksize1
        integer,allocatable:: workS(:,:),workS_p(:,:)
        parameter(maxgangs=1024,vecl=128,
     &            worksize1=maxgangs*vecl)
        contains
#include "image.f.inc"
#include "midpointimage.f.inc"

        function more_than_num(list,n,number)
!$acc routine
        implicit none
        integer,intent(in) ::list(*)
        integer,intent(in) ::n,number
        integer            ::more_than_num
        integer i,elt_old

        elt_old = 0
        more_than_num = 0
        do i = 1,n
           if (number.ge.list(i).or.list(i).eq.elt_old) cycle
           more_than_num = more_than_num+1
           elt_old = list(i)
        end do
        end function
        function more_than_num_p(list,n,number)
!$acc routine
        implicit none
        integer,pointer,intent(in) :: list(:,:)
        integer,intent(in)         :: n,number
        integer                    :: more_than_num_p
        integer i,elt_old

        elt_old = 0
        more_than_num_p = 0
        do i = 1,n
           if (number.ge.list(i,number).or.list(i,number).eq.elt_old)
     &        cycle
           more_than_num_p = more_than_num_p+1
           elt_old = list(i,number)
        end do
        end function
      end module

      subroutine scaling_preop(istep)
      use couple
      use chgpot ,only: n_cscale
      use domdec ,only: nlocnl,rank
      use inform ,only: deb_Path
      use interfaces ,only: polar_scaling,polar_scaling1,polar_scaling_p
      use mplpot ,only: n_mscale
      use mpole  ,only: npolelocnl,npolebloc
      use polgrp
      use polpot ,only: n_uscale,n_dpscale,n_dpuscale
      use potent ,only: use_polar
      use tinMemory
      use utilcomm,only: no_commdir
      use vdwpot ,only: n_vscale
      implicit none
      integer,intent(in)::istep
      logical,save :: first_in=.true.
      integer(mipk) szo_workspace
      integer(mipk),parameter :: max_sfworkSpace=2*(int(2,mipk)**30)
      integer nplocnl_e
      real(8) wsSize

      if (first_in) then
!$acc enter data create(max_factorn,
!$acc&   n_vscale,n_mscale,n_uscale,n_dpscale,n_dpuscale,n_cscale)
         first_in    = .false.
      end if

      !Initialize scaling structure data on device
      if (istep.eq.0) then
         call fuse_scaling_factor_n
         if (use_polar) call fuse_scaling_factor_p
      end if

 12   format(' Rank ',I3,'; Scaling factor workspace size ',F10.3,
     &       ' Mio', 2I4)
      if (allocated(pair_factorn)) then
         print*, "ERROR!! call scaling_preop an inconvenient way"
         print*, " Alteration of the source code "
         call fatal
      end if

      nplocnl_e = merge(npolebloc,npolelocnl,(no_commdir))
      szo_workspace= nlocnl*(max_factorn+2)*szoi +
     &               nlocnl*max_factorn*szoTp
      if (use_polar) then
      szo_workspace= szo_workspace +
     &                 nlocnl*(max_factorn+max_factorp)*szoi  +
     &               3*nlocnl*(max_factorn+max_factorp)*szoTp +
     &               4*nplocnl_e*szoi
      end if
      !print*,"scal1",szo_workspace/(2**20)*1.0

      ! Selet polar scaling subroutine
      if (szo_workspace.gt.max_sfworkSpace) then
         s_sfWork = szo_workspace -
     &              ( npolelocnl*(max_factorn+max_factorp)*szoi  +
     &              3*npolelocnl*(max_factorn+max_factorp)*szoTp )
         polar_scaling_p => polar_scaling1
      else
         s_sfWork = szo_workspace
         polar_scaling_p => polar_scaling
      end if

      allocate (pair_factorn(max_factorn,nlocnl))
      allocate (scal_factorn(max_factorn,nlocnl))
      allocate (   n_factorn(nlocnl))
      allocate (scan_factorn(nlocnl))

!$acc enter data create(pair_factorn,scal_factorn,n_factorn,
!$acc&     scan_factorn) async
      if (use_polar) then
         if (associated(polar_scaling_p,polar_scaling)) then
         allocate (pair_factorp(   max_factorp+max_factorn ,npolelocnl))
         allocate (scal_factorp(3*(max_factorp+max_factorn),npolelocnl))
!$acc enter data create(pair_factorp,scal_factorp) async
         end if
         allocate (   n_factorp(nplocnl_e),   n_factordp(nplocnl_e))
         allocate (scan_factorp(nplocnl_e),scan_factordp(nplocnl_e))
!$acc enter data create(n_factorp,scan_factorp,n_factordp,
!$acc&    scan_factordp) async
      end if

      if (deb_Path) then
         wsSize = s_sfWork/(1024*1024)
         print 12, rank,wsSize,max_factorn,max_factorp
      end if

      end subroutine

c
c     Gather in a structure array all the n scaling factor along with
c     their scaling type (1-2, 1-3, 1-4, 1-5)
c
      subroutine fuse_scaling_factor_n
      use atoms  ,only:n
      use couple
      use inform ,only:deb_Path
      use tinMemory ,only:shmem_request,mhostacc
      use kscalfactor_inl
      use utilgpu,only:maxscaling
c#ifdef _OPENACC
c      use thrust ,only:thrust_exclusive_scan
c#endif
      implicit none
      integer i,j,k,kold,ii,iscalbeg,count
      integer nn12,nn13,nn14,nn15,ntot

      ninteract_scaling_n = 0
      max_factorn         = 0

      call shmem_request(numscal_n,winnumscal_n,[n],config=mhostacc)
      call shmem_request(scalbeg_n,winscalbeg_n,[n],config=mhostacc)

      ! compute optimal size for allocation
c!$acc parallel loop present(n12,n13,n14,n15,i12,i13,i14,i15,
c!$acc&   numscal_n)
      do i = 1, n
         ntot =        more_than_num(i12(1,i),n12(i),i)
         ntot = ntot + more_than_num_p(i13,n13(i),i)
         ntot = ntot + more_than_num_p(i14,n14(i),i)
         ntot = ntot + more_than_num_p(i15,n15(i),i)

         numscal_n(i)        = ntot
         ninteract_scaling_n = ninteract_scaling_n + ntot
         max_factorn         = max(max_factorn,ntot)
      end do

 12   format(' fuse_scaling_factor_n  ( max n scaled ',I5,' )')
      if (deb_Path) write(*,12) max_factorn

      if (max_factorn.gt.maxscaling)
     &   print*,'n scaling gather buffer too short for direct space
     & computation'

c#ifdef _OPENACC
c!$acc host_data use_device(numscal_n,scalbeg_n)
c      call thrust_exclusive_scan(numscal_n, n, scalbeg_n)
c!$acc end host_data
c#else
      scalbeg_n(1) = 0
      do i = 2, n
         scalbeg_n(i) = scalbeg_n(i-1) + numscal_n(i-1)
      end do
c#endif

      !print*,max_factorn,ninteract_scaling_n,ninteract_scaling_n/n
      call shmem_request(allscal_n,winallscal_n,[ninteract_scaling_n]
     &     ,config=mhostacc)
      call shmem_request(typscal_n,wintypscal_n,[ninteract_scaling_n]
     &     ,config=mhostacc)


c!$acc parallel loop
c!$acc&         present(numscal_n,scalbeg_n,typscal_n,allscal_n,
c!$acc&   i12,i13,i14,i15,n12,n13,n14,n15)
      do i = 1, n
         nn12  = n12(i)
         nn13  = n13(i) + nn12 
         nn14  = n14(i) + nn13 
         ntot  = n15(i) + nn14 
         iscalbeg = scalbeg_n(i)
         kold  = 0
         count = 0

c!$acc loop seq
         do j=1,ntot
            if      (j.le.nn12) then
               k  = i12(j,i)
               if (i.ge.k.or.k.eq.kold) cycle
               kold  = k
               count = count+1
               allscal_n(iscalbeg+count) = k
               typscal_n(iscalbeg+count) = 2
            else if (j.le.nn13) then
               k  = i13 (j-nn12,i)
               if (i.ge.k.or.k.eq.kold) cycle
               kold  = k
               count = count+1
               allscal_n(iscalbeg+count) = k
               typscal_n(iscalbeg+count) = 3
            else if (j.le.nn14) then
               k  = i14 (j-nn13,i)
               if (i.ge.k.or.k.eq.kold) cycle
               kold  = k
               count = count+1
               allscal_n(iscalbeg+count) = k
               typscal_n(iscalbeg+count) = 4
            else
               k  = i15 (j-nn14,i)
               if (i.ge.k.or.k.eq.kold) cycle
               kold  = k
               count = count+1
               allscal_n(iscalbeg+count) = k
               typscal_n(iscalbeg+count) = 5
            end if
         end do
      end do

      !Update n scaling interactions data on device
!$acc update device(numscal_n,scalbeg_n,typscal_n,allscal_n)
      end subroutine

c
c     Gather in a structure array all the p scaling factor along with
c     their scaling type (1-1, 1-2, 1-3, 1-4)
c
      subroutine fuse_scaling_factor_p
      use atoms  ,only:n
      use inform ,only:deb_Path
      use polgrp
      use tinMemory ,only:shmem_request,mhostacc
      use kscalfactor_inl
      use utilgpu,only:maxscaling1
c#ifdef _OPENACC
c      use thrust
c#endif
      implicit none
      integer i,j,k,kold,ii,iscalbeg,count
      integer nnp11,nnp12,nnp13,nnp14

      ninteract_scaling_p = 0
      max_factorp         = 0

      ! compute optimal size for allocation
      call shmem_request(numscal_p,winnumscal_p,[n],config=mhostacc)
      call shmem_request(scalbeg_p,winscalbeg_p,[n],config=mhostacc)

c!$acc parallel loop present(np11,np12,np13,np14,ip11,ip12,
c!$acc&  ip13,ip14,numscal_p)
      do i = 1, n
         nnp14 =         more_than_num_p(ip11,np11(i),i)
         nnp14 = nnp14 + more_than_num_p(ip12,np12(i),i)
         nnp14 = nnp14 + more_than_num_p(ip13,np13(i),i)
         nnp14 = nnp14 + more_than_num_p(ip14,np14(i),i)

         numscal_p(i) = nnp14
         max_factorp         = max(max_factorp,nnp14)
         ninteract_scaling_p = ninteract_scaling_p + nnp14
      end do

12    format(' fuse_scaling_factor_p  ( max p scaled ',I5,' )')
      if (deb_Path) write(*,12) max_factorp

      if (max_factorp.gt.maxscaling1)
     &   print*,'p scaling gather buffer too short for direct space
     & computation'

c#ifdef _OPENACC
c!$acc host_data use_device(numscal_p,scalbeg_p)
c      call thrust_exclusive_scan(numscal_p, n, scalbeg_p)
c!$acc end host_data
c#else
      scalbeg_p(1) = 0
      do i = 2, n
         scalbeg_p(i) = scalbeg_p(i-1) + numscal_p(i-1)
      end do
c#endif

      !print*,max_factorp,ninteract_scaling_p,ninteract_scaling_p/n
      call shmem_request(allscal_p,winallscal_p,[ninteract_scaling_p]
     &     ,config=mhostacc)
      call shmem_request(typscal_p,wintypscal_p,[ninteract_scaling_p]
     &     ,config=mhostacc)


c!$acc parallel loop
c!$acc&         present(numscal_p,scalbeg_p,allscal_p,typscal_p,
c!$acc&   ip11,ip12,ip13,ip14,np11,np12,np13,np14)
      do i = 1, n
         nnp11    = np11(i)
         nnp12    = np12(i) + nnp11
         nnp13    = np13(i) + nnp12
         nnp14    = np14(i) + nnp13
         iscalbeg = scalbeg_p(i)
         kold  = 0
         count = 0

c!$acc loop seq
         do j=1,nnp14
            if      (j.le.nnp11) then
               k  = ip11(j,i)
               if (i.ge.k.or.k.eq.kold) cycle
               kold  = k
               count = count+1
               allscal_p(iscalbeg+count) = k
               typscal_p(iscalbeg+count) = 1
            else if (j.le.nnp12) then
               k  = ip12(j-nnp11,i)
               if (i.ge.k.or.k.eq.kold) cycle
               kold  = k
               count = count+1
               allscal_p(iscalbeg+count) = k
               typscal_p(iscalbeg+count) = 2
            else if (j.le.nnp13) then
               k  = ip13(j-nnp12,i)
               if (i.ge.k.or.k.eq.kold) cycle
               kold  = k
               count = count+1
               allscal_p(iscalbeg+count) = k
               typscal_p(iscalbeg+count) = 3
            else
               k  = ip14(j-nnp13,i)
               if (i.ge.k.or.k.eq.kold) cycle
               kold  = k
               count = count+1
               allscal_p(iscalbeg+count) = k
               typscal_p(iscalbeg+count) = 4
            end if
         end do
      end do

      !Update p scaling interactions data on device
!$acc update device(numscal_p,scalbeg_p,typscal_p,allscal_p)
      end subroutine


      subroutine vdw_scaling
      use atmlst ,only: vdwglobnl,vdwglob
      use couple
      use domdec ,only: rank,nlocnl
      use tinheader
      use tinMemory ,only: prmem_request
      use vdw    ,only: nvdwlocnl,ivdw,nvdwloc
      use vdwpot
#ifdef _OPENACC
      use thrust
#endif
      use utilgpu,only: maxscaling
#ifdef _OPENACC
     &                 ,rec_stream
#endif
      implicit none
      integer i,ii,j,k,cap,iscalbeg
      integer*1:: iscal
      integer iivdw,iglob,kglob,iscan
      integer nn12,nn13,nn14,nn15,ntot
      real(t_p) vscalevec(5),vscale

!$acc data present(vdwglob,ivdw,pair_factorn,scal_factorn,
!$acc& scan_factorn,n_factorn,n_vscale,
!$acc& allscal_n,typscal_n,scalbeg_n,numscal_n)
!$acc&     create(vscalevec)
!$acc&     async

!$acc serial async
      n_vscale  = 0
      vscalevec = (/0.0_ti_p,v2scale,v3scale,v4scale,v5scale/)
!$acc end serial

!$acc parallel loop async
      do ii = 1,nlocnl
         n_factorn(ii)=0
      end do

c
c     Filter scaling factor for vscale
c
!$acc parallel loop gang vector async
      do ii = 1,nvdwloc
         iglob = ivdw(vdwglob(ii))
         k     = 0
c        nn12  = n12(iglob)
c        nn13  = n13(iglob)
c        nn14  = n14(iglob)
c        nn15  = n15(iglob)
c        ntot  = nn12+nn13+nn14+nn15
         ntot     = numscal_n(iglob)
         iscalbeg = scalbeg_n(iglob)

         do j = 1, ntot
            kglob = allscal_n(iscalbeg+j)
            iscal = typscal_n(iscalbeg+j)
            if (iglob.ge.kglob) cycle
            if (iscal.eq.4) then
               k = k+1
               pair_factorn(k,ii) = kglob
               scal_factorn(k,ii) = -v4scale
            else if (vscalevec(iscal).ne.1) then
               k = k+1
               pair_factorn(k,ii) = kglob
               scal_factorn(k,ii) = 1-vscalevec(iscal)
            end if
         end do

c        if (v2scale.ne.1) then
c           do j = 1,nn12
c              kglob = i12(j,iglob)
c              if (iglob.lt.kglob) then
c                 k = k + 1
c                 pair_factorn(k,ii) = kglob
c                 scal_factorn(k,ii) = 1-v2scale
c              end if
c           end do
c        end if
c        if (v3scale.ne.1) then
c           do j = 1,nn13
c              kglob = i13(j,iglob)
c              if (iglob.lt.kglob) then
c                 k = k + 1
c                 pair_factorn(k,ii) = kglob
c                 scal_factorn(k,ii) = 1-v3scale
c              end if
c           end do
c        end if
c        ! 1-4 interaction need absolutely to be take into account
c           do j = 1,nn14
c              kglob = i14(j,iglob)
c              if (iglob.lt.kglob) then
c                 k = k + 1
c                 pair_factorn(k,ii) = kglob
c                 scal_factorn(k,ii) = -v4scale
c              end if
c           end do
c        if (v5scale.ne.1) then
c           do j = 1,nn15
c              kglob = i15(j,iglob)
c              if (iglob.lt.kglob) then
c                 k = k + 1
c                 pair_factorn(k,ii) = kglob
c                 scal_factorn(k,ii) = 1-v5scale
c              end if
c           end do
c        end if
         n_factorn(ii) = k
         n_vscale      = n_vscale + k
      end do
!$acc update host(n_vscale) async

c
c     Procede to scan n_factorn
c
#ifdef _OPENACC
!$acc wait
!$acc host_data use_device(n_factorn,scan_factorn)
      call thrust_exclusive_scan(n_factorn, nvdwlocnl, scan_factorn
     &                          ,rec_stream)
!$acc end host_data
#else
      scan_factorn(1) = 0
      do ii = 2,nvdwlocnl
         scan_factorn(ii) = n_factorn(ii-1) + scan_factorn(ii-1)
      end do
#endif

      ! Allocate memory to store v_scaling
      call prmem_request(vcorrect_ik,n_vscale,2,cdim=.false.)
      call prmem_request(vcorrect_scale,n_vscale)
c
c     Fill scaling the scaling_factor intercations along with their
c     values
c
!$acc parallel loop present(vcorrect_ik,vcorrect_scale)
!$acc&         gang vector
!$acc&         async
      do ii = 1, nvdwloc
         iglob = ivdw(vdwglob(ii))
         iscan = scan_factorn(ii)
         do j = 1, n_factorn(ii)
            vcorrect_ik(iscan+j,1)  = iglob
            vcorrect_ik(iscan+j,2)  = pair_factorn(j,ii)
            vcorrect_scale(iscan+j) = scal_factorn(j,ii)
         end do
      end do

!$acc end data
      end

c
c     Construct vdw scaling factor interation to postprocess
c
      subroutine vdw_scalingnl
      use atoms  ,only: x,y,z
      use atmlst ,only: vdwglobnl
      use couple
      use domdec ,only: rank,nlocnl,
     &                  xbegproc,ybegproc,zbegproc,
     &                  xendproc,yendproc,zendproc
#ifdef _OPENACC
      use utilcu ,only: cu_update_skipvdw12
#endif
      use inform ,only: deb_Path,tindPath
      use kscalfactor_inl
      use tinheader
      use tinMemory ,only: prmem_request
      use vdw       ,only: nvdwlocnl,ivdw,skipvdw12
      use vdwpot
#ifdef _OPENACC
      use thrust
      use utilgpu,only: rec_stream
#endif
      implicit none
      integer i,ii,j,k,cap
      integer*1:: iscal
      integer iscalbeg,iglob,kglob,iscan
      integer nn12,nn13,nn14,nn15,ntot
      real(t_p) pos1,pos2,pos3
      real(t_p) xi,yi,zi,xk,yk,zk
      real(t_p) zbeg,zend,ybeg,yend,xbeg,xend
      real(t_p) vscalevec(5),vscale

      if (deb_Path)
     &   write(*,'(2x,a)') 'vdw_scalingnl'

!$acc data create(vscalevec) 
!$acc&     present(vdwglobnl,ivdw,pair_factorn,scal_factorn,
!$acc& scan_factorn,n_factorn,n_vscale,
!$acc& allscal_n,typscal_n,scalbeg_n,numscal_n)
!$acc&     async

      xbeg = xbegproc(rank+1)
      xend = xendproc(rank+1)
      ybeg = ybegproc(rank+1)
      yend = yendproc(rank+1)
      zbeg = zbegproc(rank+1)
      zend = zendproc(rank+1)

      skipvdw12 = .false.
      if (v2scale.eq.0.0_ti_p) then
         skipvdw12 = .true.
      end if
#ifdef _OPENACC
      call cu_update_skipvdw12(skipvdw12)
#endif

!$acc serial async
      n_vscale  = 0
      vscalevec = (/0.0_ti_p,1.0_ti_p-v2scale,1.0_ti_p-v3scale,
     &                               -v4scale,1.0_ti_p-v5scale/)
!$acc end serial

!$acc parallel loop async
      do ii = 1,nlocnl
         n_factorn(ii)=0
      end do
c
c     Filter scaling factor for vscale
c
!$acc parallel loop gang vector async
      do ii = 1,nvdwlocnl
         iglob = ivdw(vdwglobnl(ii))
         k     = 0
         ntot     = numscal_n(iglob)
         iscalbeg = scalbeg_n(iglob)
         xi    = x(iglob)
         yi    = y(iglob)
         zi    = z(iglob)
!$acc loop seq
         do j = 1, ntot
            kglob = allscal_n(iscalbeg+j)
            iscal = typscal_n(iscalbeg+j)
            vscale= vscalevec(iscal)
            if (skipvdw12.and.iscal.eq.2.and.vscale.eq.1.0_ti_p) cycle
            if (iscal.ne.4.and.vscale.eq.0) cycle

            xk   =  x(kglob)
            yk   =  y(kglob)
            zk   =  z(kglob)
            pos1 =  xi - xk
            pos2 =  yi - yk
            pos3 =  zi - zk
            call midpointimage_inl(xk,yk,zk,pos1,pos2,pos3)
            if ((zk.lt.zbeg).or.(zk.ge.zend)
     &      .or.(yk.lt.ybeg).or.(yk.ge.yend)
     &      .or.(xk.lt.xbeg).or.(xk.ge.xend)) cycle
            k  = k+1
            pair_factorn(k,ii) = kglob
            scal_factorn(k,ii) = vscale
         end do

         n_factorn(ii) = k
         n_vscale      = n_vscale + k
      end do
!$acc update host(n_vscale) async
!$acc wait

c
c     Procede to scan
c
#ifdef _OPENACC
!$acc host_data use_device(n_factorn,scan_factorn)
      call thrust_exclusive_scan(n_factorn, nvdwlocnl, scan_factorn
     &                          ,rec_stream)
!$acc end host_data
#else
      scan_factorn(1) = 0
      do ii = 2,nvdwlocnl
         scan_factorn(ii) = n_factorn(ii-1) + scan_factorn(ii-1)
      end do
#endif
      if (btest(tinkerdebug,tindPath)) print*, 'n_vscale',n_vscale,rank

      !print*,n_vscale,nvdwlocnl
      ! Allocate memory to store v_scaling
      call prmem_request(vcorrect_ik,n_vscale,2,cdim=.false.)
      call prmem_request(vcorrect_scale,n_vscale)
c
c     Fill scaling the scaling_factor intercations along with their
c     types and value
c
!$acc parallel loop present(vcorrect_ik,vcorrect_scale)
!$acc&         gang vector
!$acc&         async
      do ii = 1, nvdwlocnl
         iglob = ivdw(vdwglobnl(ii))
         iscan = scan_factorn(ii)
         do j = 1, n_factorn(ii)
            vcorrect_ik(iscan+j,1)  = iglob
            vcorrect_ik(iscan+j,2)  = pair_factorn(j,ii)
            vcorrect_scale(iscan+j) = scal_factorn(j,ii)
         end do
      end do

!$acc end data

      end

c
c     Construct mpole scaling factor to postprocess
c
      subroutine mpole_scaling
      use atoms  ,only: x,y,z
      use atmlst ,only: poleglobnl
      use couple
      use domdec ,only: rank,nlocnl,
     &                  xbegproc,ybegproc,zbegproc,
     &                  xendproc,yendproc,zendproc
      use inform ,only: deb_Path
      use kscalfactor_inl
      use tinheader
      use tinMemory ,only: prmem_request
      use mpole  ,only: npolelocnl, ipole
      use mplpot
#ifdef _OPENACC
      use thrust
      use utilgpu,only: rec_stream
#endif
      implicit none
      integer i,ii,j,k,cap
      integer*1:: iscal
      integer iscalbeg,iglob,kglob,iscan
      integer nn12,nn13,nn14,nn15,ntot
      real(t_p) pos1,pos2,pos3
      real(t_p) xi,yi,zi,xk,yk,zk
      real(t_p) zbeg,zend,ybeg,yend,xbeg,xend
      real(t_p) mscalevec(5),mscale

!$acc data create(mscalevec) 
!$acc&     present(poleglobnl,ipole,pair_factorn,scal_factorn,
!$acc& scan_factorn,n_factorn,n_mscale,
!$acc& allscal_n,typscal_n,scalbeg_n,numscal_n)
!$acc&     async

      if (deb_Path)
     &   write(*,'(2x,a)') 'mpole_scaling'

      xbeg = xbegproc(rank+1)
      xend = xendproc(rank+1)
      ybeg = ybegproc(rank+1)
      yend = yendproc(rank+1)
      zbeg = zbegproc(rank+1)
      zend = zendproc(rank+1)

!$acc serial async
      n_mscale  = 0
      mscalevec = (/0.0_ti_p,1.0_ti_p-m2scale,1.0_ti_p-m3scale,
     &                       1.0_ti_p-m4scale,1.0_ti_p-m5scale/)
!$acc end serial

!$acc parallel loop async
      do ii = 1,nlocnl
         n_factorn(ii)=0
      end do
c
c     Filter scaling factor for vscale
c
!$acc parallel loop gang vector async
      do ii = 1,npolelocnl
         iglob = ipole(poleglobnl(ii))
         k     = 0
         ntot     = numscal_n(iglob)
         iscalbeg = scalbeg_n(iglob)
         xi    = x(iglob)
         yi    = y(iglob)
         zi    = z(iglob)
!$acc loop seq
         do j = 1, ntot
            kglob  = allscal_n(iscalbeg+j)
            iscal  = typscal_n(iscalbeg+j)
            mscale = mscalevec(iscal)
            if (mscale.eq.0.0_ti_p) cycle

            xk   =  x(kglob)
            yk   =  y(kglob)
            zk   =  z(kglob)
            pos1 =  xi - xk
            pos2 =  yi - yk
            pos3 =  zi - zk
            call midpointimage_inl(xk,yk,zk,pos1,pos2,pos3)
            if ((zk.lt.zbeg).or.(zk.ge.zend)
     &      .or.(yk.lt.ybeg).or.(yk.ge.yend)
     &      .or.(xk.lt.xbeg).or.(xk.ge.xend)) cycle
            k  = k+1
            pair_factorn(k,ii) = kglob
            scal_factorn(k,ii) = mscale
         end do

         n_factorn(ii) = k
         n_mscale      = n_mscale   +k
      end do

!$acc update host(n_mscale) async

c
c     Procede to scan
c
#ifdef _OPENACC
!$acc wait
!$acc host_data use_device(n_factorn,scan_factorn)
      call thrust_exclusive_scan(n_factorn, npolelocnl, scan_factorn
     &                          ,rec_stream)
!$acc end host_data
#else
      scan_factorn(1) = 0
      do ii = 2,npolelocnl
         scan_factorn(ii) = n_factorn(ii-1) + scan_factorn(ii-1)
      end do
#endif
c     print*,n_mscale,npolelocnl

      ! Allocate memory to store m_scaling
      call prmem_request(mcorrect_ik,n_mscale,2,cdim=.false.)
      call prmem_request(mcorrect_scale,n_mscale)
c
c     Fill scaling the scaling_factor intercations along with their
c     types and value
c
!$acc parallel loop present(mcorrect_ik,mcorrect_scale)
!$acc&         gang vector
!$acc&         async
      do ii = 1, npolelocnl
         iglob = ipole(poleglobnl(ii))
         iscan = scan_factorn(ii)
         do j = 1, n_factorn(ii)
            mcorrect_ik(iscan+j,1)  = iglob
            mcorrect_ik(iscan+j,2)  = pair_factorn(j,ii)
            mcorrect_scale(iscan+j) = scal_factorn(j,ii)
         end do
      end do

!$acc end data

      end

      subroutine polar_scaling
      use atoms  ,only: x,y,z
      use atmlst ,only: poleglobnl
      use couple
      use domdec ,only: rank,nlocnl,
     &                  xbegproc,ybegproc,zbegproc,
     &                  xendproc,yendproc,zendproc
      use inform ,only: deb_Path,tindPath
      use kscalfactor_inl
      use tinheader ,only: ti_p
      use tinMemory ,only: prmem_request
      use mpole  ,only: npolelocnl, ipole
      use polgrp
      use polpot
#ifdef _OPENACC
      use thrust
#endif
      use utilcomm,only: no_commdir
      use utils  ,only: set_to_zero1
      use utilgpu,only: rec_queue
#ifdef _OPENACC
     &           , rec_stream
#endif
      !use timestat
      implicit none
      integer i,ii,j,k,kdp,kdpu,cpt
      integer*1:: iscal
      integer iscalbeg,iscalbegp,iglob,kglob,kglob_s
      integer iscan,iscandp,iscandpu
      integer ntot,nnp14
      !integer n_dscale,n_pscale,kd,kp
      integer  :: buff(max_factorp)
      real(t_p) pos1,pos2,pos3
      real(t_p) xi,yi,zi,xk,yk,zk
      real(t_p) zbeg,zend,ybeg,yend,xbeg,xend
      real(t_p) uscalevec(5),uscale
      real(t_p) pscalevec(5),pscale,pscale41
      real(t_p) dscalevec(5),dscale

      if (deb_Path) print '(2X,A)','polar_scaling'

!$acc data create(uscalevec,pscalevec,dscalevec)
!$acc&     present(poleglobnl,ipole,
!$acc&  scan_factorn,n_factorn,pair_factorp,scal_factorp,
!$acc&  scan_factorp,n_factorp,n_uscale,n_dpscale,n_dpuscale,
!$acc&  allscal_n,typscal_n,scalbeg_n,numscal_n,allscal_p,typscal_p,
!$acc&  scalbeg_p,numscal_p,n_factordp,scan_factordp)
!$acc&     async

      xbeg = xbegproc(rank+1)
      xend = xendproc(rank+1)
      ybeg = ybegproc(rank+1)
      yend = yendproc(rank+1)
      zbeg = zbegproc(rank+1)
      zend = zendproc(rank+1)
      pscale41 = 1.0_ti_p-p4scale*p41scale
      kglob_s=0
c
c     Construct induced dipoles scalar product scaling interactions
c
!$acc serial async
      n_uscale   = 0
      !n_dscale   = 0
      !n_pscale   = 0
      n_dpscale  = 0
      n_dpuscale = 0
      uscalevec  = (/1.0_ti_p-u1scale,1.0_ti_p-u2scale,
     &               1.0_ti_p-u3scale,1.0_ti_p-u4scale,0.0_ti_p/)
      dscalevec  = (/1.0_ti_p-d1scale,1.0_ti_p-d2scale,
     &               1.0_ti_p-d3scale,1.0_ti_p-d4scale,0.0_ti_p/)
      pscalevec  = (/0.0_ti_p,1.0_ti_p-p2scale,1.0_ti_p-p3scale,
     &                        1.0_ti_p-p4scale,1.0_ti_p-p5scale/)
!$acc end serial

      call set_to_zero1(scal_factorp(1,1),
     &            3*(max_factorn+max_factorp)*npolelocnl,rec_queue)
c
c     Filter scaling factor for vscale
c
!$acc parallel loop gang vector private(buff) async
      do ii = 1,npolelocnl
         iglob = ipole(poleglobnl(ii))
         k     = 0
         !kd    = 0
         !kp    = 0
         kdp   = 0
         kdpu  = 0
         cpt   = 1
         ntot     = numscal_n(iglob)
         nnp14    = numscal_p(iglob)
         iscalbeg = scalbeg_p(iglob)
         iscalbegp= scalbeg_n(iglob)
         xi    = x(iglob)
         yi    = y(iglob)
         zi    = z(iglob)
!$acc loop seq
         do j = 1, nnp14
            kglob  = allscal_p(iscalbeg+j)
            iscal  = typscal_p(iscalbeg+j)
            uscale = uscalevec(iscal)
            dscale = dscalevec(iscal)
            if (uscale.eq.0.0_ti_p.and.dscale.eq.0.0_ti_p) cycle

            xk   =  x(kglob)
            yk   =  y(kglob)
            zk   =  z(kglob)
            pos1 =  xi - xk
            pos2 =  yi - yk
            pos3 =  zi - zk
            call midpointimage_inl(xk,yk,zk,pos1,pos2,pos3)
            if ((zk.lt.zbeg).or.(zk.ge.zend)
     &      .or.(yk.lt.ybeg).or.(yk.ge.yend)
     &      .or.(xk.lt.xbeg).or.(xk.ge.xend)) cycle
            if (uscale.ne.0.0_ti_p) then
               k  = k+1
               kdpu = kdpu+1
               kglob_s = kglob
               pair_factorp(kdpu,ii) = kglob
               scal_factorp(3*(kdpu-1)+3,ii) = uscale
            end if
            if (dscale.ne.0.0_ti_p) then
               !kd  = kd+1
               kdp = kdp+1
               if (kglob.ne.kglob_s) then !look if not found previously
                  kdpu = kdpu+1
                  pair_factorp(kdpu,ii) = kglob
                  scal_factorp(3*(kdpu-1)+1,ii) = dscale
               else
                  scal_factorp(3*(kdpu-1)+1,ii) = dscale
               end if
            end if
         end do

         ! gather (4n - 1p) interactions
         do while (cpt.le.nnp14.and.typscal_p(iscalbeg+cpt).eq.1)
            buff(cpt) = allscal_p(iscalbeg+cpt)
            cpt       = cpt + 1
         end do
!$acc loop seq
         do j = 1, ntot
            kglob  = allscal_n(iscalbegp+j)
            iscal  = typscal_n(iscalbegp+j)
            pscale = pscalevec(iscal)
            if (iscal.ne.4.and.pscale.eq.0.0_ti_p) cycle

            xk   =  x(kglob)
            yk   =  y(kglob)
            zk   =  z(kglob)
            pos1 =  xi - xk
            pos2 =  yi - yk
            pos3 =  zi - zk
            call midpointimage_inl(xk,yk,zk,pos1,pos2,pos3)
            if ((zk.lt.zbeg).or.(zk.ge.zend)
     &      .or.(yk.lt.ybeg).or.(yk.ge.yend)
     &      .or.(xk.lt.xbeg).or.(xk.ge.xend)) cycle
            if (iscal.eq.4) then ! deal with 4-1 interaction if 4n
               do i = 1,cpt-1
                  if (buff(i).eq.kglob) then
                     pscale = pscale41
                     exit
                  end if
               end do
               if (pscale.eq.0.0_ti_p) cycle
            end if
            !kp = kp+1
            do i = 1,kdpu  !look among p interactions previously gathered
               if (pair_factorp(i,ii).eq.kglob) exit
            end do
            if (i.le.kdpu) then
               scal_factorp(3*(i-1)+2,ii) = pscale
               ! check if dscale==0.0 to increment couple (d,p) counter
               if (scal_factorp(3*(i-1)+1,ii).eq.0.0_ti_p) then
                  kdp=kdp+1
               end if
            else !add pscale to the triplet
               kdp = kdp+1
               kdpu= kdpu+1
               pair_factorp(kdpu,ii) = kglob
               scal_factorp(3*(kdpu-1)+2,ii) = pscale
            end if
         end do

         n_factorp(ii) = k
         n_factordp(ii)= kdp
         n_factorn(ii) = kdpu
         n_uscale      = n_uscale   + k
         !n_dscale      = n_dscale   + kd
         n_dpscale     = n_dpscale  + kdp
         n_dpuscale    = n_dpuscale + kdpu
         !n_pscale      = n_pscale   + kp
      end do

!$acc update host(n_uscale,n_dpscale,n_dpuscale)
!$acc&       async
!$acc wait
c
c     Procede to scan
c
#ifdef _OPENACC
!$acc host_data use_device(n_factorp,n_factorn,n_factordp,
!$acc&     scan_factorp,scan_factorn,scan_factordp)
      call thrust_exclusive_scan(n_factorn, npolelocnl, scan_factorn
     &                          ,rec_stream)
      if (.not.no_commdir) then
      call thrust_exclusive_scan(n_factorp, npolelocnl, scan_factorp
     &                          ,rec_stream)
      call thrust_exclusive_scan(n_factordp, npolelocnl, scan_factordp
     &                          ,rec_stream)
      end if
!$acc end host_data
#else
      scan_factorp (1) = 0
      scan_factorn (1) = 0
      scan_factordp(1) = 0
      do ii = 2,npolelocnl
         scan_factorp(ii)  = n_factorp (ii-1) + scan_factorp (ii-1)
         scan_factorn(ii)  = n_factorn (ii-1) + scan_factorn (ii-1)
         scan_factordp(ii) = n_factordp(ii-1) + scan_factordp(ii-1)
      end do
#endif
      if (no_commdir) then  ! Skip u and dp if commdirdir is disabled
         n_uscale  = 0
         n_dpscale = 0
!$acc update host(n_uscale,n_dpscale) async
      end if

      if (btest(tinkerdebug,tindPath)) write(*,'(A,1X,5I10)')
     &         'polar_scaling',n_uscale,n_dpscale,n_dpuscale,rank

      ! Allocate memory to store u_scaling,dp_scaling & dpu_scaling
      call prmem_request(  ucorrect_ik,2*n_uscale  )
      call prmem_request( dpcorrect_ik,2*n_dpscale )
      call prmem_request(dpucorrect_ik,2*n_dpuscale)
      call prmem_request(  ucorrect_scale, n_uscale)
      call prmem_request( dpcorrect_scale,2*n_dpscale)
      call prmem_request(dpucorrect_scale,3*n_dpuscale)
c
c     Fill scaling the scaling_factor intercations along with their
c     types and value
c
!$acc parallel loop present(ucorrect_ik,ucorrect_scale,dpcorrect_ik,
!$acc&  dpcorrect_scale,dpucorrect_ik,dpucorrect_scale)
!$acc&         gang vector
!$acc&         async
      do ii = 1, npolelocnl
         iglob    = ipole(poleglobnl(ii))
         iscan    = scan_factorp (ii)
         iscandp  = scan_factordp(ii)
         iscandpu = scan_factorn (ii)
         k     = 0
         kdp   = 0
         kdpu  = 0
         ! gather (d,p,u) scaling for polarisation compute
         do j = 1, n_factorn(ii)
            kglob  = pair_factorp(j,ii)
            dscale = scal_factorp(3*(j-1)+1,ii)
            pscale = scal_factorp(3*(j-1)+2,ii)
            uscale = scal_factorp(3*(j-1)+3,ii)

            dpucorrect_ik   (2*(iscandpu+j-1)+1) = iglob
            dpucorrect_ik   (2*(iscandpu+j-1)+2) = kglob
            dpucorrect_scale(3*(iscandpu+j-1)+1) = dscale
            dpucorrect_scale(3*(iscandpu+j-1)+2) = pscale
            dpucorrect_scale(3*(iscandpu+j-1)+3) = uscale
            if (no_commdir) cycle
            ! filter (d,p) scaling for efld0_direct
            if (dscale.ne.0.0_ti_p.or.pscale.ne.0.0_ti_p) then
               kdp  = kdp+1
               dpcorrect_ik   (2*(iscandp+kdp-1)+1) = iglob
               dpcorrect_ik   (2*(iscandp+kdp-1)+2) = kglob
               dpcorrect_scale(2*(iscandp+kdp-1)+1) = dscale
               dpcorrect_scale(2*(iscandp+kdp-1)+2) = pscale
            end if
            ! filter u scaling for tmatxb
            if (uscale.ne.0.0_ti_p) then
               k   = k+1
               ucorrect_ik(2*(iscan+k-1)+1) = iglob
               ucorrect_ik(2*(iscan+k-1)+2) = kglob
               ucorrect_scale(iscan+k)      = uscale
            end if
         end do
      end do

!$acc end data

      end

      ! This routine performs the same build as polar_scaling
      ! But this one does not require any temporary buffer
      ! Double processing is performed instead
      subroutine polar_scaling1
      use atoms  ,only: x,y,z
      use atmlst ,only: poleglobnl
      use couple
      use domdec ,only: rank,nlocnl,
     &                  xbegproc,ybegproc,zbegproc,
     &                  xendproc,yendproc,zendproc
      use inform ,only: deb_Path,tindPath
      use kscalfactor_inl
      use tinheader ,only: ti_p
      use tinMemory ,only: prmem_request
      use mpole  ,only: npolelocnl, ipole
      use polgrp
      use polpot
      use sizes  ,only: tinkerdebug
#ifdef _OPENACC
      use thrust
#endif
      use utils  ,only: set_to_zero1
      use utilgpu,only: rec_queue
#ifdef _OPENACC
     &           , rec_stream
#endif
      !use timestat
      implicit none
      integer i,i1,ii,j,k,kdp,kdpu,kdpu_s,cpt
      integer*1:: iscal
      integer iscalbeg,iscalbegp,iglob,kglob,kglob_s
      integer iscan,iscandp,iscandpu
      integer ntot,nnp14
      integer ib,ic
      !integer n_dscale,n_pscale,kd,kp
      integer  :: buff(max_factorp)
      integer  :: buffp(max_factorp)
      integer  :: buffdp(0:5)
      logical fdpu,fdp
      real(t_p) pos1,pos2,pos3
      real(t_p) xi,yi,zi,xk,yk,zk
      real(t_p) zbeg,zend,ybeg,yend,xbeg,xend
      real(t_p) uscalevec(5),uscale
      real(t_p) pscalevec(5),pscale,pscale41
      real(t_p) dscalevec(5),dscale

      if (deb_Path) print '(2X,A)','polar_scaling1'

!$acc data create(uscalevec,pscalevec,dscalevec)
!$acc&     present(poleglobnl,ipole,x,y,z,
!$acc&  n_uscale,n_dpscale,n_dpuscale,
!$acc&  allscal_n,typscal_n,scalbeg_n,numscal_n,
!$acc&  allscal_p,typscal_p,scalbeg_p,numscal_p)
!$acc&     async

      xbeg = xbegproc(rank+1)
      xend = xendproc(rank+1)
      ybeg = ybegproc(rank+1)
      yend = yendproc(rank+1)
      zbeg = zbegproc(rank+1)
      zend = zendproc(rank+1)
      pscale41 = 1.0_ti_p-p4scale*p41scale
      kglob_s=0
c
c     Construct induced dipoles scalar product scaling interactions
c
!$acc serial async
      n_uscale   = 0
      !n_dscale   = 0
      !n_pscale   = 0
      n_dpscale  = 0
      n_dpuscale = 0
      uscalevec  = (/1.0_ti_p-u1scale,1.0_ti_p-u2scale,
     &               1.0_ti_p-u3scale,1.0_ti_p-u4scale,0.0_ti_p/)
      dscalevec  = (/1.0_ti_p-d1scale,1.0_ti_p-d2scale,
     &               1.0_ti_p-d3scale,1.0_ti_p-d4scale,0.0_ti_p/)
      pscalevec  = (/0.0_ti_p,1.0_ti_p-p2scale,1.0_ti_p-p3scale,
     &                        1.0_ti_p-p4scale,1.0_ti_p-p5scale/)
!$acc end serial

c
c     Analyze polar scaling interactions
c
!$acc parallel loop gang vector
!$acc&         private(buff,buffp,buffdp,cpt,i)
!$acc&         present(n_factorn,n_factorp,n_factordp)
!$acc&         async
      do ii = 1,npolelocnl
         !Fetch and init data
         iglob = ipole(poleglobnl(ii))
         k     = 0
         kdp   = 0
         kdpu  = 0
         cpt   = 1
         ntot     = numscal_n(iglob)
         nnp14    = numscal_p(iglob)
         iscalbeg = scalbeg_p(iglob)
         iscalbegp= scalbeg_n(iglob)
         xi    = x(iglob)
         yi    = y(iglob)
         zi    = z(iglob)
!$acc loop seq
         do j = 0, 5
            buffdp(j)=0
         end do
!$acc loop seq
         do j = 1, nnp14  ! Loop on p interactions
            kglob  = allscal_p(iscalbeg+j)
            iscal  = typscal_p(iscalbeg+j)
            uscale = uscalevec(iscal)
            dscale = dscalevec(iscal)
            if (uscale.eq.0.0_ti_p.and.dscale.eq.0.0_ti_p) cycle

            xk   =  x(kglob)
            yk   =  y(kglob)
            zk   =  z(kglob)
            pos1 =  xi - xk
            pos2 =  yi - yk
            pos3 =  zi - zk
            call midpointimage_inl(xk,yk,zk,pos1,pos2,pos3)
            if ((zk.lt.zbeg).or.(zk.ge.zend)
     &      .or.(yk.lt.ybeg).or.(yk.ge.yend)
     &      .or.(xk.lt.xbeg).or.(xk.ge.xend)) cycle
            if (uscale.ne.0.0_ti_p) then
               k  = k+1       !Find uscale insteraction
               kdpu = kdpu+1  !Find dpu triplet
               kglob_s = kglob
               buffp(kdpu) = kglob
            end if
            if (dscale.ne.0.0_ti_p) then
               kdp = kdp+1    !Find new dp scale interaction
               if (kglob.ne.kglob_s) then !look if not found previously
                  kdpu = kdpu+1
                  buffp(kdpu) = kglob
               end if
               buffdp((kdpu-1)/32)= ior(buffdp((kdpu-1)/32),
     &                  ishft(1,iand(kdpu-1,31)))  !Mark flag for dp scale
            end if
         end do

         kdpu_s = kdpu
         ! gather (4n - 1p) interactions
         do while (cpt.le.nnp14.and.typscal_p(iscalbeg+cpt).eq.1)
            buff(cpt) = allscal_p(iscalbeg+cpt)
            cpt       = cpt + 1
         end do

!$acc loop seq
         do j = 1, ntot      ! Loop on n interactions
            kglob  = allscal_n(iscalbegp+j)
            iscal  = typscal_n(iscalbegp+j)
            pscale = pscalevec(iscal)
            if (iscal.ne.4.and.pscale.eq.0.0_ti_p) cycle

            xk   =  x(kglob)
            yk   =  y(kglob)
            zk   =  z(kglob)
            pos1 =  xi - xk
            pos2 =  yi - yk
            pos3 =  zi - zk
            call midpointimage_inl(xk,yk,zk,pos1,pos2,pos3)
            if ((zk.lt.zbeg).or.(zk.ge.zend)
     &      .or.(yk.lt.ybeg).or.(yk.ge.yend)
     &      .or.(xk.lt.xbeg).or.(xk.ge.xend)) cycle
            if (iscal.eq.4) then ! deal with 4-1 interaction if 4n
               do i = 1,cpt-1
                  if (buff(i).eq.kglob) then
                     pscale = pscale41
                     exit
                  end if
               end do
               if (pscale.eq.0.0_ti_p) cycle
            end if
            fdp = .false.
            fdpu = .false.
            do i = 1,kdpu_s  !look among p interactions previously gathered in dpu triplet
               if (buffp(i).eq.kglob) then
                  fdpu = .true.
                  ! Check if Flag has been marked for dp couple
                  if (btest(buffdp((i-1)/32),iand(i-1,31))) fdp=.true.
                  exit
               end if
            end do
            if (.not.fdp) then
               kdp = kdp+1
            end if
            if (.not.fdpu) kdpu= kdpu+1
         end do

         n_factorp(ii) = k
         n_factordp(ii)= kdp
         n_factorn(ii) = kdpu
         n_uscale      = n_uscale   + k
         n_dpscale     = n_dpscale  + kdp
         n_dpuscale    = n_dpuscale + kdpu
      end do

!$acc update host(n_uscale,n_dpscale,n_dpuscale)
!$acc&       async

c
c     Procede to scan
c
#ifdef _OPENACC
!$acc wait
!$acc host_data use_device(n_factorp,n_factorn,n_factordp,
!$acc&     scan_factorp,scan_factorn,scan_factordp)
      call thrust_exclusive_scan(n_factorp, npolelocnl, scan_factorp
     &                          ,rec_stream)
      call thrust_exclusive_scan(n_factorn, npolelocnl, scan_factorn
     &                          ,rec_stream)
      call thrust_exclusive_scan(n_factordp, npolelocnl, scan_factordp
     &                          ,rec_stream)
!$acc end host_data
#else
      scan_factorp (1) = 0
      scan_factorn (1) = 0
      scan_factordp(1) = 0
      do ii = 2,npolelocnl
         scan_factorp(ii)  = n_factorp (ii-1) + scan_factorp (ii-1)
         scan_factorn(ii)  = n_factorn (ii-1) + scan_factorn (ii-1)
         scan_factordp(ii) = n_factordp(ii-1) + scan_factordp(ii-1)
      end do
#endif
      if (btest(tinkerdebug,tindPath)) print '(A,4I10)',
     &   'polar_scaling1',n_uscale,n_dpscale,n_dpuscale,rank
c
      ! Allocate memory to store u_scaling,dp_scaling & dpu_scaling
      call prmem_request(  ucorrect_ik,2*n_uscale  )
      call prmem_request( dpcorrect_ik,2*n_dpscale )
      call prmem_request(dpucorrect_ik,2*n_dpuscale)
      call prmem_request(  ucorrect_scale, n_uscale)
      call prmem_request( dpcorrect_scale,2*n_dpscale)
      call prmem_request(dpucorrect_scale,3*n_dpuscale)

      call set_to_zero1(  ucorrect_scale,    n_uscale,rec_queue)
      call set_to_zero1( dpcorrect_scale,2* n_dpscale,rec_queue)
      call set_to_zero1(dpucorrect_scale,3*n_dpuscale,rec_queue)

!$acc parallel loop gang vector private(buff,cpt,i,i1) async
!$acc&         present(ucorrect_ik,ucorrect_scale,dpcorrect_ik,
!$acc&  dpcorrect_scale,dpucorrect_ik,dpucorrect_scale,
!$acc&  scan_factorp,scan_factorn,scan_factordp)
      do ii = 1,npolelocnl
         iglob     = ipole(poleglobnl(ii))
         iscan     = scan_factorp (ii)
         iscandp   = scan_factordp(ii)
         iscandpu  = scan_factorn (ii)
         k         = 0
         !kd        = 0
         !kp        = 0
         kdp       = 0
         kdpu      = 0
         cpt       = 1
         ntot      = numscal_n(iglob)
         nnp14     = numscal_p(iglob)
         iscalbeg  = scalbeg_p(iglob)
         iscalbegp = scalbeg_n(iglob)
         xi        = x(iglob)
         yi        = y(iglob)
         zi        = z(iglob)
!$acc loop seq
         do j = 1, nnp14
            kglob  = allscal_p(iscalbeg+j)
            iscal  = typscal_p(iscalbeg+j)
            uscale = uscalevec(iscal)
            dscale = dscalevec(iscal)
            if (uscale.eq.0.0_ti_p.and.dscale.eq.0.0_ti_p) cycle

            xk   =  x(kglob)
            yk   =  y(kglob)
            zk   =  z(kglob)
            pos1 =  xi - xk
            pos2 =  yi - yk
            pos3 =  zi - zk
            call midpointimage_inl(xk,yk,zk,pos1,pos2,pos3)
            if ((zk.lt.zbeg).or.(zk.ge.zend)
     &      .or.(yk.lt.ybeg).or.(yk.ge.yend)
     &      .or.(xk.lt.xbeg).or.(xk.ge.xend)) cycle
            if (uscale.ne.0.0_ti_p) then
               k       = k+1
               kdpu    = kdpu+1
               kglob_s = kglob
               ucorrect_ik(2*(iscan+k-1)+1) = iglob
               ucorrect_ik(2*(iscan+k-1)+2) = kglob
               ucorrect_scale(iscan+k)      = uscale
               dpucorrect_ik   (2*(iscandpu+kdpu-1)+1) = iglob
               dpucorrect_ik   (2*(iscandpu+kdpu-1)+2) = kglob
               dpucorrect_scale(3*(iscandpu+kdpu-1)+3) = uscale
            end if
            if (dscale.ne.0.0_ti_p) then
               !kd  = kd+1
               kdp = kdp+1
               dpcorrect_ik   (2*(iscandp+kdp-1)+1) = iglob
               dpcorrect_ik   (2*(iscandp+kdp-1)+2) = kglob
               dpcorrect_scale(2*(iscandp+kdp-1)+1) = dscale
               if (kglob.ne.kglob_s) then !look if not found previously
                  kdpu = kdpu+1
                  dpucorrect_ik   (2*(iscandpu+kdpu-1)+1) = iglob
                  dpucorrect_ik   (2*(iscandpu+kdpu-1)+2) = kglob
                  dpucorrect_scale(3*(iscandpu+kdpu-1)+1) = dscale
               else
                  dpucorrect_scale(3*(iscandpu+kdpu-1)+1) = dscale
               end if
            end if
         end do

         ! gather (4n - 1p) interactions
         do while (cpt.le.nnp14.and.typscal_p(iscalbeg+cpt).eq.1)
            buff(cpt) = allscal_p(iscalbeg+cpt)
            cpt       = cpt + 1
         end do
!$acc loop seq
         do j = 1, ntot
            kglob  = allscal_n(iscalbegp+j)
            iscal  = typscal_n(iscalbegp+j)
            pscale = pscalevec(iscal)
            if (iscal.ne.4.and.pscale.eq.0.0_ti_p) cycle

            xk   =  x(kglob)
            yk   =  y(kglob)
            zk   =  z(kglob)
            pos1 =  xi - xk
            pos2 =  yi - yk
            pos3 =  zi - zk
            call midpointimage_inl(xk,yk,zk,pos1,pos2,pos3)
            if ((zk.lt.zbeg).or.(zk.ge.zend)
     &      .or.(yk.lt.ybeg).or.(yk.ge.yend)
     &      .or.(xk.lt.xbeg).or.(xk.ge.xend)) cycle
            if (iscal.eq.4) then ! deal with 4-1 interaction if 4n
               do i = 1,cpt-1
                  if (buff(i).eq.kglob) then
                     pscale = pscale41
                     exit
                  end if
               end do
               if (pscale.eq.0.0_ti_p) cycle
            end if
            !kp = kp+1
            do i = 1,kdpu  !look among p interactions previously gathered
               if (dpucorrect_ik(2*(iscandpu+i-1)+2).eq.kglob) exit
            end do
            do i1 = 1,kdp
               if (dpcorrect_ik(2*(iscandp+i1-1)+2).eq.kglob) exit
            end do

            ! Is p interaction already found among (d,p,u) triplet
            if (i.le.kdpu) then
              dpucorrect_scale(3*(iscandpu+i-1)+2) = pscale
            else !add pscale to the triplet
               kdpu= kdpu+1
               dpucorrect_ik   (2*(iscandpu+kdpu-1)+1) = iglob
               dpucorrect_ik   (2*(iscandpu+kdpu-1)+2) = kglob
               dpucorrect_scale(3*(iscandpu+kdpu-1)+2) = pscale
            end if

            ! Is p interaction already found among (d,p) couple
            if (i1.le.kdp) then
               dpcorrect_scale(2*(iscandp +i1-1)+2) = pscale
            else
               kdp=kdp+1
               dpcorrect_ik   (2*(iscandp+kdp-1)+1) = iglob
               dpcorrect_ik   (2*(iscandp+kdp-1)+2) = kglob
               dpcorrect_scale(2*(iscandp+kdp-1)+2) = pscale
            end if

         end do
      end do

!$acc end data

      end

      subroutine u_dp_scaling_extent
      use atoms  ,only: x,y,z
      use atmlst ,only: poleglob
      use couple
      use domdec ,only: rank,nlocnl,
     &                  xbegproc,ybegproc,zbegproc,
     &                  xendproc,yendproc,zendproc
      use inform ,only: deb_Path,tindPath
      use kscalfactor_inl
      use tinheader ,only: ti_p
      use tinMemory ,only: prmem_request
      use mpole  ,only: npolelocnl, npolebloc, ipole
      use polgrp
      use polpot
      use sizes  ,only: tinkerdebug
#ifdef _OPENACC
      use thrust
#endif
      use utils  ,only: set_to_zero1
      use utilgpu,only: rec_queue
#ifdef _OPENACC
     &           , rec_stream
#endif
      !use timestat
      implicit none
      integer i,i1,ii,j,k,kdp,kdp_s,cpt
      integer*1:: iscal
      integer iscalbeg,iscalbegp,iglob,kglob,kglob_s
      integer iscan,iscandp
      integer ntot,nnp14
      integer workl,ngangs
      logical fdp
      real(t_p) pos1,pos2,pos3
      real(t_p) xi,yi,zi,xk,yk,zk
      real(t_p) zbeg,zend,ybeg,yend,xbeg,xend
      real(t_p) uscalevec(5),uscale
      real(t_p) pscalevec(5),pscale,pscale41
      real(t_p) dscalevec(5),dscale

      if (deb_Path) print '(2X,A)','u_dp_scaling_extent'

!$acc data create(uscalevec,pscalevec,dscalevec)
!$acc&     present(poleglob,ipole,n_uscale,n_dpscale,
!$acc&  allscal_n,typscal_n,scalbeg_n,numscal_n,
!$acc&  allscal_p,typscal_p,scalbeg_p,numscal_p)
!$acc&     async(rec_queue)

      xbeg = xbegproc(rank+1)
      xend = xendproc(rank+1)
      ybeg = ybegproc(rank+1)
      yend = yendproc(rank+1)
      zbeg = zbegproc(rank+1)
      zend = zendproc(rank+1)
      pscale41 = 1.0_ti_p-p4scale*p41scale
      kglob_s=0
c
c     Construct induced dipoles scalar product scaling interactions
c
      n_uscale   = 0
      n_dpscale  = 0
!$acc serial async(rec_queue)
      n_uscale   = 0
      !n_dscale   = 0
      !n_pscale   = 0
      n_dpscale  = 0
      !n_dpuscale = 0
      uscalevec  = (/1.0_ti_p-u1scale,1.0_ti_p-u2scale,
     &               1.0_ti_p-u3scale,1.0_ti_p-u4scale,0.0_ti_p/)
      dscalevec  = (/1.0_ti_p-d1scale,1.0_ti_p-d2scale,
     &               1.0_ti_p-d3scale,1.0_ti_p-d4scale,0.0_ti_p/)
      pscalevec  = (/0.0_ti_p,1.0_ti_p-p2scale,1.0_ti_p-p3scale,
     &                        1.0_ti_p-p4scale,1.0_ti_p-p5scale/)
!$acc end serial
      ngangs = min(maxgangs,npolebloc/vecl)
      call prmem_request(workS  ,worksize1,max_factorp,async=.true.)
      call prmem_request(workS_p,worksize1,max_factorp,async=.true.)
c
c     Analyze polar scaling interactions
c
!$acc parallel loop num_gangs(ngangs) vector_length(vecl)
!$acc&         private(cpt)
!$acc&         present(n_factorp,n_factordp,workS,workS_p)
!$acc&         async(rec_queue)
      do ii = 1,npolebloc
         !Fetch and init data
         iglob = ipole(poleglob(ii))
         k     = 0
         kdp   = 0
         cpt   = 1
         workl = iand(ii,worksize1-1)+1
         ntot     = numscal_n(iglob)
         nnp14    = numscal_p(iglob)
         iscalbeg = scalbeg_p(iglob)
         iscalbegp= scalbeg_n(iglob)
!$acc loop seq
         do j = 1, nnp14  ! Loop on p interactions
            kglob  = allscal_p(iscalbeg+j)
            iscal  = typscal_p(iscalbeg+j)
            uscale = uscalevec(iscal)
            dscale = dscalevec(iscal)
            if (uscale.eq.0.0_ti_p.and.dscale.eq.0.0_ti_p) cycle

            if (uscale.ne.0.0_ti_p) then
               k  = k+1       !Find uscale insteraction
            end if
            if (dscale.ne.0.0_ti_p) then
               kdp = kdp+1    !Find new dp scale interaction
               workS_p(workl,kdp) = kglob
            end if
         end do

         kdp_s = kdp
         ! gather (4n - 1p) interactions
         do while (cpt.le.nnp14.and.typscal_p(iscalbeg+cpt).eq.1)
            workS(workl,cpt) = allscal_p(iscalbeg+cpt)
            cpt       = cpt + 1
         end do

!$acc loop seq
         do j = 1, ntot      ! Loop on n interactions
            kglob  = allscal_n(iscalbegp+j)
            iscal  = typscal_n(iscalbegp+j)
            pscale = pscalevec(iscal)
            if (iscal.ne.4.and.pscale.eq.0.0_ti_p) cycle

            if (iscal.eq.4) then ! deal with 4-1 interaction if 4n
               do i = 1,cpt-1
                  if (workS(workl,i).eq.kglob) then
                     pscale = pscale41
                     exit
                  end if
               end do
               if (pscale.eq.0.0_ti_p) cycle
            end if
            fdp = .false.
            do i = 1,kdp_s  !look among p interactions previously gathered for (d,p) couple
               if (workS_p(workl,i).eq.kglob) then
                  fdp = .true.
                  exit
               end if
            end do
            if (.not.fdp) then
               kdp = kdp+1
            end if
         end do

         n_factorp (ii)= k
         n_factordp(ii)= kdp
         n_uscale      = n_uscale   + k
         n_dpscale     = n_dpscale  + kdp
      end do

!$acc update host(n_uscale,n_dpscale) async(rec_queue)
c
c     Procede to scan
c
#ifdef _OPENACC
!$acc wait(rec_queue)
!$acc host_data use_device(n_factorp,n_factordp,
!$acc&     scan_factorp,scan_factordp)
      if (n_uscale.ne.0)
     &call thrust_exclusive_scan(n_factorp, npolebloc, scan_factorp
     &                          ,rec_stream)
      if (n_dpscale.ne.0)
     &call thrust_exclusive_scan(n_factordp, npolebloc, scan_factordp
     &                          ,rec_stream)
!$acc end host_data
#else
      scan_factorp (1) = 0
      scan_factordp(1) = 0
      do ii = 2,npolebloc
         scan_factorp (ii) = n_factorp (ii-1) + scan_factorp (ii-1)
         scan_factordp(ii) = n_factordp(ii-1) + scan_factordp(ii-1)
      end do
#endif
      if (btest(tinkerdebug,tindPath)) print '(A,4I10)',
     &   'u_dp_scaling_extent',n_uscale,n_dpscale,rank
c
      ! Allocate memory to store u_scaling,dp_scaling & dpu_scaling
      call prmem_request(  ucorrect_ik,2*n_uscale  )
      call prmem_request( dpcorrect_ik,2*n_dpscale )
      call prmem_request(  ucorrect_scale, n_uscale)
      call prmem_request( dpcorrect_scale,2*n_dpscale)

      ! Zero memory
      call set_to_zero1(  ucorrect_scale,    n_uscale,rec_queue)
      call set_to_zero1( dpcorrect_scale,2* n_dpscale,rec_queue)

!$acc parallel loop num_gangs(ngangs) vector_length(vecl)
!$acc&         private(i,i1,cpt) async(rec_queue)
!$acc&         present(ucorrect_ik,ucorrect_scale,dpcorrect_ik,
!$acc&  dpcorrect_scale,scan_factorp,scan_factordp,
!$acc&  workS)
      do ii = 1,npolebloc
         iglob     = ipole(poleglob(ii))
         iscan     = scan_factorp (ii)
         iscandp   = scan_factordp(ii)
         k         = 0
         kdp       = 0
         cpt       = 1
         workl     = iand(ii,worksize1-1)+1
         ntot      = numscal_n(iglob)
         nnp14     = numscal_p(iglob)
         iscalbeg  = scalbeg_p(iglob)
         iscalbegp = scalbeg_n(iglob)
!$acc loop seq
         do j = 1, nnp14
            kglob  = allscal_p(iscalbeg+j)
            iscal  = typscal_p(iscalbeg+j)
            uscale = uscalevec(iscal)
            dscale = dscalevec(iscal)
            if (uscale.eq.0.0_ti_p.and.dscale.eq.0.0_ti_p) cycle

            if (uscale.ne.0.0_ti_p) then
               k       = k+1
               ucorrect_ik(2*(iscan+k-1)+1) = iglob
               ucorrect_ik(2*(iscan+k-1)+2) = kglob
               ucorrect_scale(iscan+k)      = uscale
            end if
            if (dscale.ne.0.0_ti_p) then
               kdp = kdp+1
               dpcorrect_ik   (2*(iscandp+kdp-1)+1) = iglob
               dpcorrect_ik   (2*(iscandp+kdp-1)+2) = kglob
               dpcorrect_scale(2*(iscandp+kdp-1)+1) = dscale
            end if
         end do

         ! gather (4n - 1p) interactions
         do while (cpt.le.nnp14.and.typscal_p(iscalbeg+cpt).eq.1)
            workS(workl,cpt) = allscal_p(iscalbeg+cpt)
            cpt       = cpt + 1
         end do
!$acc loop seq
         do j = 1, ntot
            kglob  = allscal_n(iscalbegp+j)
            iscal  = typscal_n(iscalbegp+j)
            pscale = pscalevec(iscal)
            if (iscal.ne.4)  then
               if (pscale.eq.0.0_ti_p) cycle
            else  ! deal with 4-1 interaction if 4n
               do i = 1,cpt-1
                  if (workS(workl,i).eq.kglob) then
                     pscale = pscale41; exit;
                  end if
               end do
               if (pscale.eq.0.0_ti_p) cycle
            end if
            do i1 = 1,kdp
               if (dpcorrect_ik(2*(iscandp+i1-1)+2).eq.kglob) exit
            end do

            ! Is p interaction already found among (d,p) couple
            if (i1.le.kdp) then  ! yes
               dpcorrect_scale(2*(iscandp +i1-1)+2) = pscale
            else  ! No
               kdp=kdp+1
               dpcorrect_ik   (2*(iscandp+kdp-1)+1) = iglob
               dpcorrect_ik   (2*(iscandp+kdp-1)+2) = kglob
               dpcorrect_scale(2*(iscandp+kdp-1)+2) = pscale
            end if
         end do
      end do

!$acc end data

      end

      subroutine charge_scaling
      use atoms  ,only: x,y,z
      use atmlst ,only: chgglobnl
      use chgpot
      use charge ,only: nionlocnl,iion,pchg,chglist
      use couple
      use domdec ,only: rank,nlocnl,
     &                  xbegproc,ybegproc,zbegproc,
     &                  xendproc,yendproc,zendproc
      use inform ,only: deb_Path
      use kscalfactor_inl
      use tinheader
      use tinMemory ,only: prmem_request
#ifdef _OPENACC
      use thrust
      use utilgpu,only: rec_stream
#endif
      implicit none
      integer i,ii,j,k,cap,iichg,kkchg
      integer*1:: iscal
      integer iscalbeg,iglob,kglob,iscan
      integer nn12,nn13,nn14,nn15,ntot
      real(t_p) pos1,pos2,pos3
      real(t_p) xi,yi,zi,xk,yk,zk
      real(t_p) zbeg,zend,ybeg,yend,xbeg,xend
      real(t_p) cscalevec(5),cscale

!$acc data create(cscalevec)
!$acc&     present(chgglobnl,iion,pair_factorn,scal_factorn,
!$acc& scan_factorn,n_factorn,n_cscale,
!$acc& allscal_n,typscal_n,scalbeg_n,numscal_n)
!$acc&     async

      if (deb_Path)
     &   write(*,'(2x,a)') 'charge_scaling'

      xbeg = xbegproc(rank+1)
      xend = xendproc(rank+1)
      ybeg = ybegproc(rank+1)
      yend = yendproc(rank+1)
      zbeg = zbegproc(rank+1)
      zend = zendproc(rank+1)

!$acc serial async
      n_cscale  = 0
      cscalevec = (/0.0_ti_p,1.0_ti_p-c2scale,1.0_ti_p-c3scale,
     &                       1.0_ti_p-c4scale,1.0_ti_p-c5scale/)
!$acc end serial

!$acc parallel loop async
      do ii = 1,nlocnl
         n_factorn(ii)=0
      end do
c
c     Filter scaling factor for vscale
c
!$acc parallel loop gang vector async
      do ii = 1,nionlocnl
         iglob = iion(chgglobnl(ii))
         k     = 0
         ntot     = numscal_n(iglob)
         iscalbeg = scalbeg_n(iglob)
         xi    = x(iglob)
         yi    = y(iglob)
         zi    = z(iglob)
!$acc loop seq
         do j = 1, ntot
            kglob  = allscal_n(iscalbeg+j)
            iscal  = typscal_n(iscalbeg+j)
            cscale = cscalevec(iscal)
            if (cscale.eq.0.0_ti_p) cycle

            xk   =  x(kglob)
            yk   =  y(kglob)
            zk   =  z(kglob)
            pos1 =  xi - xk
            pos2 =  yi - yk
            pos3 =  zi - zk
            call midpointimage_inl(xk,yk,zk,pos1,pos2,pos3)
            if ((zk.lt.zbeg).or.(zk.ge.zend)
     &      .or.(yk.lt.ybeg).or.(yk.ge.yend)
     &      .or.(xk.lt.xbeg).or.(xk.ge.xend)) cycle
            k  = k+1
            pair_factorn(k,ii) = kglob
            scal_factorn(k,ii) = cscale
         end do

         n_factorn(ii) = k
         n_cscale      = n_cscale + k
      end do

!$acc update host(n_cscale) async

c
c     Procede to scan
c
#ifdef _OPENACC
!$acc wait
!$acc host_data use_device(n_factorn,scan_factorn)
      call thrust_exclusive_scan(n_factorn, nionlocnl, scan_factorn
     &                          ,rec_stream)
!$acc end host_data
#else
      scan_factorn(1) = 0
      do ii = 2,nionlocnl
         scan_factorn(ii) = n_factorn(ii-1) + scan_factorn(ii-1)
      end do
#endif
c     print*,n_cscale,nionlocnl

      ! Allocate memory to store m_scaling
      call prmem_request(ccorrect_ik,n_cscale,2,cdim=.false.)
      call prmem_request(ccorrect_scale,2*(n_cscale+1))
c
c     Fill scaling the scaling_factor intercations along with their
c     types and value
c
!$acc parallel loop present(ccorrect_ik,ccorrect_scale)
!$acc&         gang vector
!$acc&         async
      do ii = 1, nionlocnl
         iichg = chgglobnl(ii)
         iglob = iion(iichg)
         iscan = scan_factorn(ii)
         do j = 1, n_factorn(ii)
            k  = pair_factorn(j,ii)
            ccorrect_ik(iscan+j,1)  = iglob
            ccorrect_ik(iscan+j,2)  = k
            ccorrect_scale(2*(iscan+j)+1) =-scal_factorn(j,ii)
            ccorrect_scale(2*(iscan+j)+2) = pchg(iichg)*
     &                  pchg(chglist(k))
         end do
      end do

!$acc end data

      end

      subroutine scaling_postop
      use couple
      use polgrp
      use potent ,only: use_polar
      use tinMemory
      implicit none 

      if (.not.allocated(pair_factorn)) then
         print*, "ERROR !! call scaling_postop an inconvenient way"
         print*, "         Alteration of the source code "
         call fatal
      end if

!$acc exit data delete(pair_factorn,scal_factorn,n_factorn,
!$acc&    scan_factorn) async
      deallocate (pair_factorn)
      deallocate (scal_factorn)
      deallocate (   n_factorn)
      deallocate (scan_factorn)
      if (use_polar) then
         if (allocated(pair_factorp)) then
!$acc exit data delete(pair_factorp) async
            deallocate (pair_factorp)
         end if
         if (allocated(scal_factorp)) then
!$acc exit data delete(scal_factorp) async
            deallocate (scal_factorp)
         end if
!$acc exit data delete(n_factorp,scan_factorp,
!$acc&     n_factordp,scan_factordp) async
         deallocate (   n_factorp,   n_factordp)
         deallocate (scan_factorp,scan_factordp)
      end if

      end

      subroutine check_outboundnl
      use atmlst ,only: poleglobnl,vdwglobnl
      use domdec ,only: nbloc,loc
      use mpole  ,only: npolelocnl,ipole
      use potent
      use utilgpu,only: warning,sugest_vdw
      use vdw    ,only: ivdw,nvdwlocnl
      implicit none
      integer ii,i,iipole,iivdw

      if (use_vdw) then
!$acc parallel loop present(vdwglobnl,loc,ivdw,
!$acc&   warning,sugest_vdw) async
         do ii = 1, nvdwlocnl
            iivdw = vdwglobnl(ii)
            i     = loc(ivdw(iivdw))
            if (i.eq.0.or.i.gt.nbloc) then
               print*,warning
               print*,"VLST"
               cycle
            end if
         end do
      end if
      if (use_mpole) then
!$acc parallel loop present(poleglobnl,loc,ipole,
!$acc&   warning) async
         do ii = 1, npolelocnl
            iipole = poleglobnl(ii)
            i      = loc(ipole(iipole))
            if ((i.eq.0).or.(i.gt.nbloc)) then
               print*, warning
               print*, "ELST"
               cycle
            end if
         end do
      end if

      end

      subroutine build_scaling_factor(istep)
      use domdec  ,only: rank,nproc,ndir
      use inform  ,only: deb_Path
      use interfaces,only: polar_scaling_p
      use neigh   ,only: ineigup
      use potent
      use utilcomm,only: no_commdir
      use timestat,only: timer_enter,timer_exit,
     &                   timer_scalfactor
      implicit none
      integer,intent(in)::istep
      integer ndirp
      integer,save::period=4

      if (istep.eq.0) then
         if (use_pmecore) then
            ndirp = ndir
         else
            ndirp = nproc
         end if

         ! look for update period
         if (ndirp.le.2) then
            period = 32
         else if (ndirp.le.4) then
            period = 20
         else if (ndirp.le.8) then
            period = 12
         else if (ndirp.le.16) then
            period = 8
         else if (ndirp.le.32) then
            period = 6
         else if (ndirp.le.64) then
            period = 4
         else
            period = 2
         endif

         !skip pmecore for the first time
         goto 10
      end if

      ! exclude reciproqual processes
      if (use_pmecore.and.rank.ge.ndir) return

      ! Check for update
      if (mod(istep,period*ineigup).ne.0) return

10    continue

      if (use_polar.and.no_commdir) call orderpole

      if (deb_Path) print*,"build_scaling_factor"

      call check_outboundnl
      call timer_enter(timer_scalfactor)
      call scaling_preop(istep)

      if (use_vdw)    call vdw_scalingnl
      if (use_mpole)  call mpole_scaling
      if (use_charge) call charge_scaling
      if (use_polar) then
         call polar_scaling_p
         if (no_commdir) call u_dp_scaling_extent
      end if

      call scaling_postop
      call timer_exit(timer_scalfactor)

      end subroutine
