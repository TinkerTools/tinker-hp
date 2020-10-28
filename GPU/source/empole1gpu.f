
c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     "empole1" : driver for calculation of the multipole and dipole polarization
c     energy and derivatives with respect to Cartesian coordinates
c
c
#include "tinker_precision.h"
      module empole1gpu_inl
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

      subroutine empole1gpu
      use domdec
      use energi
      use potent
      use tinheader ,only:ti_p,re_p
      implicit none
c
c     choose the method for summing over multipole interactions
c
      call empole1cgpu
c
c     zero out energy and derivative terms which are not in use
c
      if (.not. use_mpole) then
         em = 0.0_ti_p
      end if
      if (.not. use_polar) then
         ep = 0.0_ti_p
      end if
      end
c
c
c     "empole1c" calculates the multipole energy and derivatives
c     with respect to Cartesian coordinates using particle mesh Ewald
c     summation and a neighbor list
c
c
      subroutine empole1cgpu
      use sizes
      use atmlst
      use atoms
      use boxes
      use chgpot
      use deriv
      use domdec
      use energi
      use ewald
      use inform     ,only: deb_Path
      use interfaces ,only: reorder_nblist
     &               ,torquegpu,commpoleglob
      use math
      use mpole
      use mpi
      use neigh
      use potent
      use shunt
      use timestat
      use utilgpu
      use virial
      implicit none
      integer i,j,ii
      integer iipole,iglob,ierr
      real(t_p) zero,one,two,three,half
      real(t_p) e,emdir,f
      real(t_p) term,fterm
      real(t_p) cii,dii,qii
      real(t_p) xd,yd,zd
      real(t_p) xq,yq,zq
      real(t_p) xv,yv,zv,vterm
      real(t_p) ci,dix,diy,diz
      real(t_p) qixx,qixy,qixz
      real(t_p) qiyy,qiyz,qizz
      real(t_p) xdfield,ydfield
      real(t_p) zdfield
      real(t_p) trq(3,npoleloc)
      real(t_p) time0,time1
      parameter(zero=0.0_ti_p,  one=1.0_ti_p,
     &           two=2.0_ti_p,three=3.0_ti_p,
     &          half=0.5_ti_p)
c
      if (npole .eq. 0)  return
c
c     zero out the atomic multipole energy and derivatives
c
      if(deb_Path) write(*,*) 'empole1cgpu'

!$acc data create(emdir,emrec)
!$acc&     present(em,dem)
!$acc&     async(rec_queue)
c
!$acc serial async
      em     = zero
      emdir  = zero
      emrec  = zero
!$acc end serial
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
      call timer_enter( timer_other )
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
c!$acc parallel loop async(dir_queue) default(present)
c              do i = 1,npolelocnl
c                 nshortelstc(i) = nshortelst(i)
c              end do
            end if
            if (mlst_enable) then
               call switch('EWALD     ')
               call reorder_nblist(elst,nelst,nelstc
     &                    ,npolelocnl,off2,ipole,poleglobnl)
c!$acc parallel loop async(dir_queue) default(present)
c              do i = 1,npolelocnl
c                 nelstc(i) = nelst(i)
c              end do
            end if
         end if
      end if
      call timer_exit ( timer_other,quiet_timers )

!!$acc update host(rpole,pole) async(rec_queue)
#ifdef _OPENACC
      if (dir_queue.ne.rec_queue) then
         call stream_wait_async(dir_stream,rec_stream,dir_event)
      end if
#endif
c
c     compute the reciprocal space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.le.ndir-1))
     &   then
         call timer_enter( timer_real )
         if (use_mreal) then
            if (use_mpoleshortreal) then
               call emrealshort1cgpu
            else if (use_mpolelong) then
               call emreallong1cgpu
            else
               call emreal1cgpu
            end if
         end if

         if (use_mself) then
!$acc data present(poleglob,rpole,ipole,x,y,z,dem,loc,vir)
c
c     compute the Ewald self-energy term over all the atoms
c
         term  = two * aewald * aewald
         fterm = -f * aewald / sqrtpi
!$acc parallel loop async(def_queue)
         do i = 1, npoleloc
            iipole = poleglob(i)
            ci     = rpole( 1,iipole)
            dix    = rpole( 2,iipole)
            diy    = rpole( 3,iipole)
            diz    = rpole( 4,iipole)
            qixx   = rpole( 5,iipole)
            qixy   = rpole( 6,iipole)
            qixz   = rpole( 7,iipole)
            qiyy   = rpole( 9,iipole)
            qiyz   = rpole(10,iipole)
            qizz   = rpole(13,iipole)
            cii    = ci*ci
            dii    =       dix*dix  +  diy*diy  + diz*diz
            qii    = two*(qixy*qixy + qixz*qixz + qiyz*qiyz)
     &                  + qixx*qixx + qiyy*qiyy + qizz*qizz
            e      = fterm*(cii + term*(dii/three +
     &                                  two*term*qii/5.0_ti_p))
            emdir  = emdir + e
         end do
c
c     compute the cell dipole boundary correction term
c
         if (boundary .eq. 'VACUUM') then
!$acc data create(trq,term,xd,yd,zd,xq,yq,zq,xdfield,
!$acc&  ydfield,zdfield)
!$acc&     async(def_queue)

!$acc serial async(def_queue)
            xd = zero
            yd = zero
            zd = zero
!$acc end serial
!$acc parallel loop async(def_queue)
            do i = 1, npoleloc
               iipole = poleglob(i)
               iglob  = ipole(iipole)
               xd = xd + rpole(2,iipole) + rpole(1,iipole)*x(iglob)
               yd = yd + rpole(3,iipole) + rpole(1,iipole)*y(iglob)
               zd = zd + rpole(4,iipole) + rpole(1,iipole)*z(iglob)
            end do
!$acc wait(def_queue)
!$acc host_data use_device(xd,yd,zd)
            call MPI_ALLREDUCE(MPI_IN_PLACE,xd,1,MPI_TPREC,MPI_SUM,
     $         COMM_TINKER,ierr)
            call MPI_ALLREDUCE(MPI_IN_PLACE,yd,1,MPI_TPREC,MPI_SUM,
     $         COMM_TINKER,ierr)
            call MPI_ALLREDUCE(MPI_IN_PLACE,zd,1,MPI_TPREC,MPI_SUM,
     $         COMM_TINKER,ierr)
!$acc end host_data
            if (rank.eq.0) then
!$acc serial async(def_queue)
              em = em + term*(xd*xd+yd*yd+zd*zd)
!$acc end serial
            end if
            term    = (two/3.0_ti_p) * f * (pi/volbox)
            emdir   = emdir + term*(xd*xd+yd*yd+zd*zd)
            xdfield = -two * term * xd
            ydfield = -two * term * yd
            zdfield = -two * term * zd
!$acc parallel loop async(def_queue)
            do ii = 1, npoleloc
               iipole   = poleglob(ii)
               iglob    = ipole(iipole)
               i        = loc(iglob)
               dem(1,i) = dem(1,i) + two*term*rpole(1,iipole)*xd
               dem(2,i) = dem(2,i) + two*term*rpole(1,iipole)*yd
               dem(3,i) = dem(3,i) + two*term*rpole(1,iipole)*zd

               trq(1,ii)=rpole(3,iipole)*zdfield-rpole(4,iipole)*ydfield
               trq(2,ii)=rpole(4,iipole)*xdfield-rpole(2,iipole)*zdfield
               trq(3,ii)=rpole(2,iipole)*ydfield-rpole(3,iipole)*xdfield
            end do
            call torquegpu(npoleloc,poleglob,loc,trq,dem,def_queue)
c
c     boundary correction to virial due to overall cell dipole
c
!$acc kernels async(def_queue)
            xd = zero
            yd = zero
            zd = zero
            xq = zero
            yq = zero
            zq = zero
!$acc end kernels
!$acc parallel loop async(def_queue)
            do i = 1, npoleloc
               iipole = poleglob(i)
               iglob  = ipole(iipole)
               xd     = xd + rpole(2,iipole)
               yd     = yd + rpole(3,iipole)
               zd     = zd + rpole(4,iipole)
               xq     = xq + rpole(1,iipole)*x(iglob)
               yq     = yq + rpole(1,iipole)*y(iglob)
               zq     = zq + rpole(1,iipole)*z(iglob)
            end do
!$acc wait(def_queue)
!$acc host_data use_device(xd,yd,zd,xq,yq,zq)
            call MPI_ALLREDUCE(MPI_IN_PLACE,xd,1,MPI_TPREC,MPI_SUM,
     $         COMM_TINKER,ierr)
            call MPI_ALLREDUCE(MPI_IN_PLACE,yd,1,MPI_TPREC,MPI_SUM,
     $         COMM_TINKER,ierr)
            call MPI_ALLREDUCE(MPI_IN_PLACE,zd,1,MPI_TPREC,MPI_SUM,
     $         COMM_TINKER,ierr)
            call MPI_ALLREDUCE(MPI_IN_PLACE,xq,1,MPI_TPREC,MPI_SUM,
     $         COMM_TINKER,ierr)
            call MPI_ALLREDUCE(MPI_IN_PLACE,yq,1,MPI_TPREC,MPI_SUM,
     $         COMM_TINKER,ierr)
            call MPI_ALLREDUCE(MPI_IN_PLACE,zq,1,MPI_TPREC,MPI_SUM,
     $         COMM_TINKER,ierr)
!$acc end host_data
            if (rank.eq.0) then
!$acc kernels async(def_queue)
            xv = xd * xq
            yv = yd * yq
            zv = zd * zq
            vterm = term * (xd*xd + yd*yd + zd*zd + two*(xv+yv+zv)
     &                    + xq*xq + yq*yq + zq*zq)
            vir(1,1) = vir(1,1) + two*term*(xq*xq+xv) + vterm
            vir(2,1) = vir(2,1) + two*term*(xq*yq+xv)
            vir(3,1) = vir(3,1) + two*term*(xq*zq+xv)
            vir(1,2) = vir(1,2) + two*term*(yq*xq+yv)
            vir(2,2) = vir(2,2) + two*term*(yq*yq+yv) + vterm
            vir(3,2) = vir(3,2) + two*term*(yq*zq+yv)
            vir(1,3) = vir(1,3) + two*term*(zq*xq+zv)
            vir(2,3) = vir(2,3) + two*term*(zq*yq+zv)
            vir(3,3) = vir(3,3) + two*term*(zq*zq+zv) + vterm
!$acc end kernels
            end if

!$acc end data
c
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
         if (use_mrec) then
         call timer_enter( timer_rec )
         call emrecip1gpu
         call timer_exit( timer_rec,quiet_timers )
         end if
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
c     Add both contribution to the energy
c
!$acc serial async(rec_queue)
      em = em + emdir + emrec
!$acc end serial
c
!$acc end data
      end
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine emreal1c  --  Ewald real space derivs via list  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "emreal1c" evaluates the real space portion of the Ewald
c     summation energy and gradient due to multipole interactions
c     via a neighbor list
c
c
      subroutine emreal1cgpu
      use atmlst  ,only: poleglobnl
      use atoms   ,only: x,y,z
      use deriv   ,only: dem
      use domdec  ,only: nbloc
      use interfaces,only: emreal1c_core_p,torquegpu
      use mpole   ,only: xaxis,yaxis,zaxis,npolelocnl,ipole
      use potent  ,only: use_mpoleshortreal,use_mpolelong
      use tinheader, only: ti_p
      use utils   ,only: set_to_zero1
      use utilgpu ,only: dir_queue,def_queue,rec_queue
#ifdef _OPENACC
     & ,rec_stream,dir_stream,def_stream,rec_event,stream_wait_async
#endif
      use virial  ,only: vxx=>g_vxx,vxy=>g_vxy,vxz=>g_vxz
     &            ,vyy=>g_vyy,vyz=>g_vyz,vzz=>g_vzz
      !use mpi
      use timestat,only:timer_enter,timer_exit,timer_emreal
     &            ,quiet_timers
      implicit none
      integer i,ii,iipole,iglob
      integer iax,iay,iaz
      real(t_p) zero
      real(t_p) r2,f,e
      real(t_p) alsq2,alsq2n
      real(t_p) xix,yix,zix
      real(t_p) xiy,yiy,ziy
      real(t_p) xiz,yiz,ziz
      real(t_p) fix(3,npolelocnl)
      real(t_p) fiy(3,npolelocnl)
      real(t_p) fiz(3,npolelocnl)
      real(t_p) tem(3,nbloc)
      character*10 mode
      logical*1,parameter::extract=.true.
      parameter(zero=0.0_ti_p)

      !if(deb_Path) write(*,'(2x,a)') 'emreal1cgpu'
      call timer_enter( timer_emreal )
      def_queue = dir_queue

!$acc data create(tem,fix,fiy,fiz) async(def_queue)

      call set_to_zero1(tem,3*nbloc,def_queue)

      call switch('EWALD')

#ifdef _OPENACC
      if (rec_queue.ne.dir_queue) then  !start async overlapping
         def_stream = dir_stream
         call stream_wait_async(rec_stream,dir_stream,rec_event)
      end if
#endif

      ! Pointer on one many subroutine
      call emreal1c_core_p(tem,vxx,vxy,vxz,vyy,vyz,vzz)
c
c     resolve site torques then increment forces and virial
c
      call torquegpu(tem,fix,fiy,fiz,dem,extract)

!$acc parallel loop vector_length(32) async(def_queue)
!$acc&         present(poleglobnl,x,y,z,ipole,
!$acc&  xaxis,yaxis,zaxis,vxx,vxy,vxz,vyy,vyz,vzz)
      do ii = 1, npolelocnl
         iipole = poleglobnl(ii)
         iglob  = ipole(iipole)
         iaz    = zaxis(iipole)
         iax    = xaxis(iipole)
         iay    = yaxis(iipole)
         if (iaz.le.0) iaz = iglob
         if (iax.le.0) iax = iglob
         if (iay.le.0) iay = iglob
         xiz    = x(iaz) - x(iglob)
         yiz    = y(iaz) - y(iglob)
         ziz    = z(iaz) - z(iglob)
         xix    = x(iax) - x(iglob)
         yix    = y(iax) - y(iglob)
         zix    = z(iax) - z(iglob)
         xiy    = x(iay) - x(iglob)
         yiy    = y(iay) - y(iglob)
         ziy    = z(iay) - z(iglob)
         vxx    = vxx + xix*fix(1,ii) + xiy*fiy(1,ii) + xiz*fiz(1,ii)
         vxy    = vxy + yix*fix(1,ii) + yiy*fiy(1,ii) + yiz*fiz(1,ii)
         vxz    = vxz + zix*fix(1,ii) + ziy*fiy(1,ii) + ziz*fiz(1,ii)
         vyy    = vyy + yix*fix(2,ii) + yiy*fiy(2,ii) + yiz*fiz(2,ii)
         vyz    = vyz + zix*fix(2,ii) + ziy*fiy(2,ii) + ziz*fiz(2,ii)
         vzz    = vzz + zix*fix(3,ii) + ziy*fiy(3,ii) + ziz*fiz(3,ii)
      end do
c
!$acc end data
c
      call timer_exit( timer_emreal )
      end

c
c     ##################################################################################
c     ##                                                                              ##
c     ##  subroutine emrealshort1c  --  Ewald short range real space derivs via list  ##
c     ##                                                                              ##
c     ##################################################################################
c
c
c     "emrealshort1cgpu" evaluates the short range real space portion of the Ewald
c     summation energy and gradient due to multipole interactions
c     via a neighbor list
c
c
      subroutine emrealshort1cgpu
      use atmlst  ,only: poleglobnl
      use atoms   ,only: x,y,z
      use deriv   ,only: dem
      use domdec  ,only: nbloc, rank
      use inform ,only: deb_Path
      use interfaces,only: emrealshortlong1c_core_p,torquegpu
      use mpole   ,only: xaxis,yaxis,zaxis,npolelocnl,ipole
      use potent  ,only: use_mpoleshortreal,use_mpolelong
      use tinheader, only: ti_p
      use utils   ,only: set_to_zero1
      use utilgpu ,only: dir_queue,def_queue,rec_queue
#ifdef _OPENACC
     & ,rec_stream,dir_stream,def_stream,rec_event,stream_wait_async
#endif
      use virial  ,only: vxx=>g_vxx,vxy=>g_vxy,vxz=>g_vxz
     &            ,vyy=>g_vyy,vyz=>g_vyz,vzz=>g_vzz
      !use mpi
      use timestat,only:timer_enter,timer_exit,timer_emreal
      implicit none
      integer i,ii,iipole,iglob
      integer iax,iay,iaz
      real(t_p) zero
      real(t_p) r2,f,e
      real(t_p) alsq2,alsq2n
      real(t_p) xix,yix,zix
      real(t_p) xiy,yiy,ziy
      real(t_p) xiz,yiz,ziz
      real(t_p) fix(3,npolelocnl)
      real(t_p) fiy(3,npolelocnl)
      real(t_p) fiz(3,npolelocnl)
      real(t_p) tem(3,nbloc)
      character*10 mode
      logical*1,parameter::extract=.true.
      parameter(zero=0.0_ti_p)

      if(deb_Path) write(*,'(2x,a)') 'emrealshort1cgpu'
      call timer_enter( timer_emreal )
      def_queue = dir_queue

!$acc data create(tem,fix,fiy,fiz) async(def_queue)

      call set_to_zero1(tem,3*nbloc,def_queue)

      call switch('SHORTEWALD')

#ifdef _OPENACC
      if (rec_queue.ne.dir_queue) then  !start async overlapping
         def_stream = dir_stream
         call stream_wait_async(rec_stream,dir_stream,rec_event)
      end if
#endif

      ! Pointer on one many subroutine
      call emrealshortlong1c_core_p(tem,vxx,vxy,vxz,vyy,vyz,vzz)
c
c     resolve site torques then increment forces and virial
c
      call torquegpu(tem,fix,fiy,fiz,dem,extract)

!$acc parallel loop vector_length(32) async(def_queue)
!$acc&         present(poleglobnl,x,y,z,ipole,
!$acc&  xaxis,yaxis,zaxis,vxx,vxy,vxz,vyy,vyz,vzz)
      do ii = 1, npolelocnl
         iipole = poleglobnl(ii)
         iglob  = ipole(iipole)
         iaz    = zaxis(iipole)
         iax    = xaxis(iipole)
         iay    = yaxis(iipole)
         if (iaz.le.0) iaz = iglob
         if (iax.le.0) iax = iglob
         if (iay.le.0) iay = iglob
         xiz    = x(iaz) - x(iglob)
         yiz    = y(iaz) - y(iglob)
         ziz    = z(iaz) - z(iglob)
         xix    = x(iax) - x(iglob)
         yix    = y(iax) - y(iglob)
         zix    = z(iax) - z(iglob)
         xiy    = x(iay) - x(iglob)
         yiy    = y(iay) - y(iglob)
         ziy    = z(iay) - z(iglob)
         vxx    = vxx + xix*fix(1,ii) + xiy*fiy(1,ii) + xiz*fiz(1,ii)
         vxy    = vxy + yix*fix(1,ii) + yiy*fiy(1,ii) + yiz*fiz(1,ii)
         vxz    = vxz + zix*fix(1,ii) + ziy*fiy(1,ii) + ziz*fiz(1,ii)
         vyy    = vyy + yix*fix(2,ii) + yiy*fiy(2,ii) + yiz*fiz(2,ii)
         vyz    = vyz + zix*fix(2,ii) + ziy*fiy(2,ii) + ziz*fiz(2,ii)
         vzz    = vzz + zix*fix(3,ii) + ziy*fiy(3,ii) + ziz*fiz(3,ii)
      end do
c
!$acc end data
c
      call timer_exit( timer_emreal )
      end

c
c     ################################################################################
c     ##                                                                            ##
c     ##  subroutine emreallong1c  --  Ewald long range real space derivs via list  ##
c     ##                                                                            ##
c     ################################################################################
c
c
c     "emreallong1c" evaluates the long range part real space portion of the Ewald
c     summation energy and gradient due to multipole interactions
c     via a neighbor list
c
c
      subroutine emreallong1cgpu
      use atmlst  ,only: poleglobnl
      use atoms   ,only: x,y,z
      use deriv   ,only: dem
      use domdec  ,only: nbloc,rank
      use inform  ,only: deb_Path
      use interfaces,only: emrealshortlong1c_core_p,torquegpu
      use mpole   ,only: xaxis,yaxis,zaxis,npolelocnl,ipole
      use potent  ,only: use_mpoleshortreal,use_mpolelong
      use tinheader, only: ti_p
      use utils   ,only: set_to_zero1
      use utilgpu ,only: dir_queue,def_queue,rec_queue
#ifdef _OPENACC
     & ,rec_stream,dir_stream,def_stream,rec_event,stream_wait_async
#endif
      use virial  ,only: vxx=>g_vxx,vxy=>g_vxy,vxz=>g_vxz
     &            ,vyy=>g_vyy,vyz=>g_vyz,vzz=>g_vzz
      !use mpi
      use timestat,only:timer_enter,timer_exit,timer_emreal
      implicit none
      integer i,ii,iipole,iglob
      integer iax,iay,iaz
      real(t_p) zero
      real(t_p) r2,f,e
      real(t_p) alsq2,alsq2n
      real(t_p) xix,yix,zix
      real(t_p) xiy,yiy,ziy
      real(t_p) xiz,yiz,ziz
      real(t_p) fix(3,npolelocnl)
      real(t_p) fiy(3,npolelocnl)
      real(t_p) fiz(3,npolelocnl)
      real(t_p) tem(3,nbloc)
      character*10 mode
      logical*1,parameter::extract=.true.
      parameter(zero=0.0_ti_p)

      if(deb_Path) write(*,'(2x,a)') 'emreallong1cgpu'
      call timer_enter( timer_emreal )
      def_queue = dir_queue

!$acc data create(tem,fix,fiy,fiz)
!$acc&     async(def_queue)

      call set_to_zero1(tem,3*nbloc,def_queue)

      call switch('EWALD')

#ifdef _OPENACC
      if (rec_queue.ne.dir_queue) then  !start async overlapping
         def_stream = dir_stream
         call stream_wait_async(rec_stream,dir_stream,rec_event)
      end if
#endif

      ! Pointer on one many subroutine
      call emrealshortlong1c_core_p(tem,vxx,vxy,vxz,vyy,vyz,vzz)
c
c     resolve site torques then increment forces and virial
c
      call torquegpu(tem,fix,fiy,fiz,dem,extract)

!$acc parallel loop vector_length(32) async(def_queue)
!$acc&         present(poleglobnl,x,y,z,ipole,
!$acc&  xaxis,yaxis,zaxis,vxx,vxy,vxz,vyy,vyz,vzz)
      do ii = 1, npolelocnl
         iipole = poleglobnl(ii)
         iglob  = ipole(iipole)
         iaz    = zaxis(iipole)
         iax    = xaxis(iipole)
         iay    = yaxis(iipole)
         if (iaz.le.0) iaz = iglob
         if (iax.le.0) iax = iglob
         if (iay.le.0) iay = iglob
         xiz    = x(iaz) - x(iglob)
         yiz    = y(iaz) - y(iglob)
         ziz    = z(iaz) - z(iglob)
         xix    = x(iax) - x(iglob)
         yix    = y(iax) - y(iglob)
         zix    = z(iax) - z(iglob)
         xiy    = x(iay) - x(iglob)
         yiy    = y(iay) - y(iglob)
         ziy    = z(iay) - z(iglob)
         vxx    = vxx + xix*fix(1,ii) + xiy*fiy(1,ii) + xiz*fiz(1,ii)
         vxy    = vxy + yix*fix(1,ii) + yiy*fiy(1,ii) + yiz*fiz(1,ii)
         vxz    = vxz + zix*fix(1,ii) + ziy*fiy(1,ii) + ziz*fiz(1,ii)
         vyy    = vyy + yix*fix(2,ii) + yiy*fiy(2,ii) + yiz*fiz(2,ii)
         vyz    = vyz + zix*fix(2,ii) + ziy*fiy(2,ii) + ziz*fiz(2,ii)
         vzz    = vzz + zix*fix(3,ii) + ziy*fiy(3,ii) + ziz*fiz(3,ii)
      end do
c
!$acc end data
c
      call timer_exit( timer_emreal )
      end


c======================================================================
c         ===================================================
c                   ===============================
c
c     SECTION : Core comput routines of emreal1c
c
c                   ===============================
c         ===================================================
c======================================================================

      ! Compute empole interactions for dynamic
      ! Reference version
      subroutine emreal1c_core1(tem,vxx,vxy,vxz,vyy,vyz,vzz)
      use atmlst ,only:poleglobnl
      use atoms  ,only:x,y,z
      !use bound  ,only:use_bounds
      use chgpot ,only:electric,dielec
      use couple
      use deriv  ,only:dem
      use domdec ,only:rank,nbloc,loc
      use empole1gpu_inl ,only: image_inl,mpole1_couple
      use energi ,only:em
      !use erf_mod
      use ewald  ,only:aewald
      use inform ,only:deb_Path
      use iounit ,only:iout
      use math   ,only:sqrtpi
      !use mpi
      use mplpot
      use mpole  ,only:rpole,ipole,polelocnl,npolelocnl
      use neigh  ,only:nelst,elst
      use shunt  ,only:off2
      use tinheader, only: ti_p
      use utilgpu,only:dir_queue,rec_queue,def_queue,warning,
     &                 real3,real3_red,rpole_elt,maxscaling
      use virial
      implicit none

      real(r_p),intent(inout):: vxx,vyy,vzz
      real(r_p),intent(inout):: vxy,vxz,vyz
      real(t_p),intent(inout):: tem(3,nbloc)

      integer i,j,k,iglob,kglob,kbis
      integer ii,kk,ki,kkk,iipole,kkpole
      integer kbstart,kbend,kb,kbq,kbr
      integer inl,nnelst
      integer nn12,nn13,nn14,ntot
      integer iax,iay,iaz
      integer iscal(maxscaling)
      real(t_p) mscalevec(5)
      real(t_p) zero
      real(t_p) r2,f,e
      real(t_p) alsq2,alsq2n
      real(t_p) xi,yi,zi
      real(t_p) xr,yr,zr
      type(rpole_elt) ip,kp
      type(real3_red) frc_r
      type(real3) ttmi,ttmk,frc
      real(t_p) mscale
      real(t_p) fscal(maxscaling)
      parameter(zero=0.0_ti_p)

      if(deb_Path) write(*,'(2x,a)') 'emreal1c_core1'

!$acc data present(tem,vxx,vxy,vxz,vyy,vyz,vzz)
!$acc&     present(dem,vir,poleglobnl,ipole,loc,x,y,z,rpole,
!$acc&  elst,nelst,i12,i13,i14,i15,allscal_n,typscal_n,scalbeg_n,
!$acc&  numscal_n)
!$acc&     present(em)
!$acc&     create(mscalevec)
!$acc&     async(def_queue)
c
c     set conversion factor, cutoff and switching coefficients
c
      f      = electric / dielec
      alsq2  = 2 * aewald**2
      alsq2n = zero
      if (aewald .gt. zero)
     &   alsq2n = 1.0_ti_p / (sqrtpi*aewald)
!$acc serial async(def_queue)
      mscalevec = [zero,m2scale,m3scale,m4scale,m5scale]
!$acc end serial
c
c     compute the real space portion of the Ewald summation
c
!$acc parallel loop gang vector_length(32)
!$acc&         private(ki,iscal,fscal)
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
c     get number of atoms  directly (1-2) or 1-3 or 1-4 bonded
c
         ki   = 0
         ntot = numscal_n(iglob)
         nn12 = scalbeg_n(iglob)
!$acc loop vector
         do j = 1,ntot
            iscal(j) = allscal_n(nn12+j)
            fscal(j) = mscalevec(typscal_n(nn12+j))
         end do
c
c     evaluate all sites within the cutoff distance
c
!$acc loop vector
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
c
c     find interaction mscale coefficients for connected atoms
c
            mscale = 1.0_ti_p
            if (ki<ntot) then
               do j=1,ntot
                  if (iscal(j)==kglob) then
                     mscale = fscal(j)
!$acc atomic update
                     ki  = ki+1
                     exit
                  end if
               end do
            end if

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
            call mpole1_couple(r2,xr,yr,zr,ip,kp,1.0_ti_p-mscale,
     &                         aewald,f,alsq2n,alsq2,
     &                         e,frc,frc_r,ttmi,ttmk,.false.)

            ! update energy
            em       = em  + e
c
c     increment force-based gradient and torque on first site
c
!$acc atomic update
            dem(1,i)    = dem(1,i) + frc_r%x
!$acc atomic update
            dem(2,i)    = dem(2,i) + frc_r%y
!$acc atomic update
            dem(3,i)    = dem(3,i) + frc_r%z
!$acc atomic update
            dem(1,kbis) = dem(1,kbis) - frc_r%x
!$acc atomic update
            dem(2,kbis) = dem(2,kbis) - frc_r%y
!$acc atomic update
            dem(3,kbis) = dem(3,kbis) - frc_r%z
c
c     increment force-based gradient and torque on second site
c
!$acc atomic update
            tem(1,i)    = tem(1,i) + ttmi%x
!$acc atomic update
            tem(2,i)    = tem(2,i) + ttmi%y
!$acc atomic update
            tem(3,i)    = tem(3,i) + ttmi%z
!$acc atomic update
            tem(1,kbis) = tem(1,kbis) + ttmk%x
!$acc atomic update
            tem(2,kbis) = tem(2,kbis) + ttmk%y
!$acc atomic update
            tem(3,kbis) = tem(3,kbis) + ttmk%z
c
c     increment the virial due to pairwise Cartesian forces
c
            vxx    = vxx  - xr * frc%x
            vxy    = vxy  - yr * frc%x
            vxz    = vxz  - zr * frc%x
            vyy    = vyy  - yr * frc%y
            vyz    = vyz  - zr * frc%y
            vzz    = vzz  - zr * frc%z
         end do
      end do
!$acc end data
c
      end

      ! Compute emreal interactions with scale set to 1 and procede to
      ! correction if necessary
      subroutine emreal1c_core2(tem,vxx,vxy,vxz,vyy,vyz,vzz)
      use atmlst ,only:poleglobnl
      use atoms  ,only:x,y,z,n
      !use bound  ,only:use_bounds
      use chgpot ,only:electric,dielec
      !use couple
      use deriv  ,only:dem
      use domdec ,only:rank,nbloc,loc
      use empole1gpu_inl ,only: image_inl,mpole1_couple
      use energi ,only:em
      !use erf_mod
      use ewald  ,only:aewald
      !use iounit ,only:iout
      use inform ,only: deb_Path
      use interfaces,only: emreal_correct_interactions
      use math   ,only:sqrtpi
      !use mpi
      use mplpot ,only:n_mscale,mcorrect_ik,mcorrect_scale
      use mpole  ,only:rpole,ipole,polelocnl,npolelocnl
      use neigh  ,only:nelst,elst
      use shunt  ,only:off2
      use utilgpu,only:dir_queue,rec_queue,def_queue,warning,
     &                 real3,real3_red,rpole_elt,maxscaling
      use virial
      implicit none

      real(r_p),intent(inout):: vxx,vyy,vzz
      real(r_p),intent(inout):: vxy,vxz,vyz
      real(t_p),intent(inout):: tem(3,nbloc)

      integer i,j,k,iglob,kglob,kbis
      integer ii,kk,iipole,kkpole
      integer nnelst
#ifdef TINKER_DEBUG
      integer inter(n)
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

      if(deb_Path) write(*,'(2x,a)') 'emreal1c_core2'

#ifdef TINKER_DEBUG
!$acc enter data create(inter) async(def_queue)
!$acc parallel loop async present(inter)
      do i = 1,n
         inter(i) = 0
      end do
#endif

c
c     set conversion factor, cutoff and switching coefficients
c
      f      = electric / dielec
      alsq2  = 2 * aewald**2
      alsq2n = zero
      if (aewald .gt. zero)
     &   alsq2n = 1.0 / (sqrtpi*aewald)
c
c     compute the real space portion of the Ewald summation
c
!$acc parallel loop gang vector_length(32)
!$acc&         present(tem,vxx,vxy,vxz,vyy,vyz,vzz)
!$acc&         present(dem,em,vir,poleglobnl,ipole,loc,x,y,z,rpole,
!$acc&  elst,nelst) private(ip)
#ifdef TINKER_DEBUG
!$acc&         present(inter)
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
!$acc loop vector private(kp,frc,frc_r,ttmi,ttmk)
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
            call mpole1_couple(r2,xr,yr,zr,ip,kp,zero,
     &                         aewald,f,alsq2n,alsq2,
     &                         e,frc,frc_r,ttmi,ttmk,.false.)

            ! update energy
            em       = em  + e

#ifdef TINKER_DEBUG
!$acc atomic
            inter(iglob) = inter(iglob)+1
#endif
c
c     increment force-based gradient and torque on first site
c
!$acc atomic update
            dem(1,i)    = dem(1,i)    + frc_r%x
!$acc atomic update
            dem(2,i)    = dem(2,i)    + frc_r%y
!$acc atomic update
            dem(3,i)    = dem(3,i)    + frc_r%z
!$acc atomic update
            dem(1,kbis) = dem(1,kbis) - frc_r%x
!$acc atomic update
            dem(2,kbis) = dem(2,kbis) - frc_r%y
!$acc atomic update
            dem(3,kbis) = dem(3,kbis) - frc_r%z
c
c     increment force-based gradient and torque on second site
c
!$acc atomic update
            tem(1,i)    = tem(1,i)    + ttmi%x
!$acc atomic update
            tem(2,i)    = tem(2,i)    + ttmi%y
!$acc atomic update
            tem(3,i)    = tem(3,i)    + ttmi%z
!$acc atomic update
            tem(1,kbis) = tem(1,kbis) + ttmk%x
!$acc atomic update
            tem(2,kbis) = tem(2,kbis) + ttmk%y
!$acc atomic update
            tem(3,kbis) = tem(3,kbis) + ttmk%z
c
c     increment the virial due to pairwise Cartesian forces
c
            vxx    = vxx  - xr * frc%x
            vxy    = vxy  - yr * frc%x
            vxz    = vxz  - zr * frc%x
            vyy    = vyy  - yr * frc%y
            vyz    = vyz  - zr * frc%y
            vzz    = vzz  - zr * frc%z
         end do
      end do

#ifdef TINKER_DEBUG
 34   format(2I10,3F12.4)
 35   format(A30,I16)
!$acc wait
!$acc exit data copyout(inter)
      print 35, 'total mpole interactions ', sum(inter)
c      do i =1,250000
c         print 34,i,inter(i)
c      end do
#endif

      call emreal_correct_interactions(tem,vxx,vxy,vxz,vyy,vyz,vzz)

      end

      !
      ! Compute emreal pair interaction into one single loop
      ! using a precompute
      !
      subroutine emreal1c_core3(tem,vxx,vxy,vxz,vyy,vyz,vzz)
      use atmlst ,only:poleglobnl
      use atoms  ,only:x,y,z
      !use bound  ,only:use_bounds
      use chgpot ,only:electric,dielec
      use couple
      use deriv  ,only:dem
      use domdec ,only:rank,nbloc,loc
      use empole1gpu_inl ,only: image_inl,mpole1_couple
      use energi ,only:em
      !use erf_mod
      use ewald  ,only:aewald
      use inform ,only: deb_Path
      use iounit ,only:iout
      use math   ,only:sqrtpi
      !use mpi
      use mplpot
      use mpole  ,only:rpole,ipole,polelocnl,npolelocnl
      use neigh  ,only:nelst,elst
      use precompute_pole ,only: iipole_precompute,kkpole_precompute,
     &                     rx_precompute,ry_precompute,rz_precompute,
     &                     ms_precompute,nprecomp,mpole_precompute
      use shunt  ,only:off2
      use utilgpu,only:dir_queue,rec_queue,def_queue,maxscaling,
     &                 warning,real3,real3_red,rpole_elt
      use virial
      implicit none

      real(r_p),intent(inout):: vxx,vyy,vzz
      real(r_p),intent(inout):: vxy,vxz,vyz
      real(t_p),intent(inout):: tem(3,nbloc)

      integer i,j,k,iglob,kglob,kbis
      integer ii,kk,ki,kkk,iipole,kkpole
      integer kbstart,kbend,kb,kbq,kbr
      integer inl,nnelst
      integer nn12,nn13,nn14,ntot
      integer iax,iay,iaz
      integer iscal(maxscaling)
      real(t_p) mscalevec(5)
      real(t_p) zero,one
      real(t_p) r2,f,e
      real(t_p) alsq2,alsq2n
      real(t_p) xi,yi,zi
      real(t_p) xr,yr,zr
      type(rpole_elt) ip,kp
      type(real3_red) frc_r
      type(real3) ttmi,ttmk,frc
      real(t_p) mscale
      real(t_p) fscal(maxscaling)
      parameter(zero=0.0,one=1.0)

      if(deb_Path) write(*,'(2x,a)') 'emreal1c_core3'

      call mpole_precompute

!$acc data present(tem,vxx,vxy,vxz,vyy,vyz,vzz)
!$acc&     present(dem,vir,poleglobnl,ipole,loc,x,y,z,rpole,
!$acc&  iipole_precompute,kkpole_precompute,rx_precompute,
!$acc&  ry_precompute,rz_precompute,ms_precompute)
!$acc&     present(em)
!$acc&     create(mscalevec)
!$acc&     async(def_queue)
c
c     set conversion factor, cutoff and switching coefficients
c
      f      = electric / dielec
      alsq2  = 2 * aewald**2
      alsq2n = zero
      if (aewald .gt. zero)
     &   alsq2n = one / (sqrtpi*aewald)
!$acc serial async(def_queue)
      mscalevec = [zero,m2scale,m3scale,m4scale,m5scale]
!$acc end serial
c
c     compute the real space portion of the Ewald summation
c
!$acc parallel loop
!$acc&         private(ip,kp,frc,frc_r,ttmi,ttmk)
!$acc&         async(def_queue)
      do ii = 1, nprecomp
         iipole = iipole_precompute(ii)
         iglob  = ipole(iipole)
         i      = loc(iglob)
         kkpole = kkpole_precompute(ii)
         kglob  = ipole(kkpole)
         kbis   = loc(kglob)

         xr     = rx_precompute(ii)
         yr     = ry_precompute(ii)
         zr     = rz_precompute(ii)
         r2     = xr*xr + yr*yr + zr*zr

         mscale = ms_precompute(ii)

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

         kp%c   = rpole( 1,kkpole)
         kp%dx  = rpole( 2,kkpole)
         kp%dy  = rpole( 3,kkpole)
         kp%dz  = rpole( 4,kkpole)
         kp%qxx = rpole( 5,kkpole)
         kp%qxy = rpole( 6,kkpole)
         kp%qxz = rpole( 7,kkpole)
         kp%qyy = rpole( 9,kkpole)
         kp%qyz = rpole(10,kkpole)
         kp%qzz = rpole(13,kkpole)

         ! compute mpole one interaction
         call mpole1_couple(r2,xr,yr,zr,ip,kp,mscale,
     &                      aewald,f,alsq2n,alsq2,
     &                      e,frc,frc_r,ttmi,ttmk,.false.)

         ! update energy
         em     = em  + e
c
c     increment force-based gradient and torque on first site
c
!$acc atomic update
         dem(1,i)    = dem(1,i) + frc_r%x
!$acc atomic update
         dem(2,i)    = dem(2,i) + frc_r%y
!$acc atomic update
         dem(3,i)    = dem(3,i) + frc_r%z
!$acc atomic update
         dem(1,kbis) = dem(1,kbis) - frc_r%x
!$acc atomic update
         dem(2,kbis) = dem(2,kbis) - frc_r%y
!$acc atomic update
         dem(3,kbis) = dem(3,kbis) - frc_r%z
c
c     increment force-based gradient and torque on second site
c
!$acc atomic update
         tem(1,i)    = tem(1,i) + ttmi%x
!$acc atomic update
         tem(2,i)    = tem(2,i) + ttmi%y
!$acc atomic update
         tem(3,i)    = tem(3,i) + ttmi%z
!$acc atomic update
         tem(1,kbis) = tem(1,kbis) + ttmk%x
!$acc atomic update
         tem(2,kbis) = tem(2,kbis) + ttmk%y
!$acc atomic update
         tem(3,kbis) = tem(3,kbis) + ttmk%z
c
c     increment the virial due to pairwise Cartesian forces
c
         vxx    = vxx  - xr * frc%x
         vxy    = vxy  - yr * frc%x
         vxz    = vxz  - zr * frc%x
         vyy    = vyy  - yr * frc%y
         vyz    = vyz  - zr * frc%y
         vzz    = vzz  - zr * frc%z
      end do
!$acc end data
c
      end


      ! CUDA C wrapper on emreal1c (check file cu_mpole1.cu)
      subroutine emreal1c_core4(tem,vxx,vxy,vxz,vyy,vyz,vzz)
      use atmlst ,only:poleglobnl
      use atoms  ,only:x,y,z,n
      use chgpot ,only:electric,dielec
      use cell
      use deriv  ,only:dem
      use domdec  ,only: xbegproc,ybegproc,zbegproc
     &            ,nproc,rank,xendproc,yendproc,zendproc
     &            ,nbloc,loc
      use energi ,only:em
      !use erf_mod
      use ewald  ,only:aewald
      use inform ,only: deb_Path
      use interfaces,only: emreal_correct_interactions
#ifdef _CUDA
     &              , cu_emreal1c
#endif
      use math   ,only:sqrtpi
      !use mpi
      use mplpot ,only:n_mscale,mcorrect_ik,mcorrect_scale
      use mpole  ,only:rpole,ipole,polelocnl,npolelocnl
     &           ,npolelocnlb,npolelocnlb_pair,npolebloc
     &           ,npolelocnlb2_pair
      use neigh  ,only:ipole_s=>celle_glob,pglob_s=>celle_pole
     &           , loc_s=>celle_loc, ieblst_s=>ieblst, eblst_s=>eblst
     &           , x_s=>celle_x, y_s=>celle_y, z_s=>celle_z
      use shunt  ,only:off2
      use utilgpu,only:dir_queue,rec_queue,def_queue,warning
     &           ,real3,real3_red,rpole_elt,RED_BUFF_SIZE
     &           ,ered_buff,vred_buff,reduce_energy_virial
     &           ,zero_evir_red_buffer
#ifdef  _OPENACC
     &           ,dir_stream,def_stream
#endif
      use virial
      implicit none

      real(r_p),intent(inout):: vxx,vyy,vzz
      real(r_p),intent(inout):: vxy,vxz,vyz
      real(t_p),intent(inout):: tem(3,nbloc)

      integer i
      real(t_p) r2,f,e
      real(t_p) alsq2,alsq2n
      real(t_p) p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend

      if(deb_Path) write(*,'(2x,a)') 'emreal1c_core4'

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

c
c     call C wrapper to compute the real space portion of the Ewald summation
c
#ifdef _CUDA
      def_stream = dir_stream
      call zero_evir_red_buffer(def_queue)

!$acc host_data use_device(ipole_s,pglob_s,loc_s,ieblst_s,eblst_s,
!$acc&    x_s,y_s,z_s,rpole,dem,tem,ered_buff,vred_buff)

      call cu_emreal1c( ipole_s,pglob_s,loc_s,ieblst_s
     &                , eblst_s(2*npolelocnlb_pair+1)
     &                , x_s,y_s,z_s,rpole
     &                , dem,tem,ered_buff,vred_buff
     &                , npolelocnlb,npolelocnlb2_pair,npolebloc,n
     &                , off2,f,alsq2,alsq2n,aewald
     &                , xcell, ycell, zcell,xcell2,ycell2,zcell2
     &                ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
     &                , def_stream)

!$acc end host_data

      call reduce_energy_virial(em,vxx,vxy,vxz,vyy,vyz,vzz,def_queue)
#else
      print 100
 100  format('emreal1c_core4 is a specific device routine !!',/,
     &       'you are not supposed to get inside with your compile ',
     &       'mode')
      call fatal
#endif

      call emreal_correct_interactions(tem,vxx,vxy,vxz,vyy,vyz,vzz)

      end

      ! CUDA Fortran version of emreal1c using C2 nblist
      subroutine emreal1c_core5(tem,vxx,vxy,vxz,vyy,vyz,vzz)
      use atmlst ,only:poleglobnl
      use atoms  ,only:x,y,z,n
      use chgpot ,only:electric,dielec
      use cell
      use deriv  ,only:dem
      use domdec  ,only: xbegproc,ybegproc,zbegproc
     &            ,nproc,rank,xendproc,yendproc,zendproc
     &            ,nbloc,loc
#ifdef _CUDA
      use empole1cu ,only: emreal1c_core_cu
#endif
      use energi ,only:em
      use ewald  ,only:aewald
      use inform ,only: deb_Path
      use interfaces,only: emreal_correct_interactions
      use math   ,only:sqrtpi
      !use mpi
      use mplpot ,only:n_mscale,mcorrect_ik,mcorrect_scale
      use mpole  ,only:rpole,ipole,polelocnl,npolelocnl
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
     &           ,ered_buff,vred_buff,reduce_energy_virial
     &           ,zero_evir_red_buffer
#ifdef  _OPENACC
     &           ,dir_stream,def_stream
#endif
      use virial
      implicit none

      real(r_p),intent(inout):: vxx,vyy,vzz
      real(r_p),intent(inout):: vxy,vxz,vyz
      real(t_p),intent(inout):: tem(3,nbloc)

      integer i
#ifdef TINKER_DEBUG
      integer inter(n)
#endif
      real(t_p) r2,f,e
      real(t_p) alsq2,alsq2n
      real(t_p) p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
      logical,save:: first_in=.true.
      integer,save:: gS

      if(deb_Path) write(*,'(2x,a)') 'emreal1c_core5'

#ifdef TINKER_DEBUG
!$acc enter data create(inter) async(def_queue)
!$acc parallel loop async present(inter)
      do i = 1,n
         inter(i) = 0
      end do
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

      if (first_in) then
#ifdef _CUDA
         call cudaMaxGridSize("emreal1c_core_cu",gS)
#endif
         first_in = .false.
      end if

#ifdef _CUDA
      def_stream = dir_stream
      call zero_evir_red_buffer(def_queue)

!$acc host_data use_device(ipole_s,pglob_s,loc_s,ieblst_s,eblst_s
!$acc&    ,x_s,y_s,z_s,rpole,dem,tem,ered_buff,vred_buff
#ifdef TINKER_DEBUG
!$acc&    ,inter
#endif
!$acc&    )

      !call CUDA kernel to compute the real space portion of the Ewald summation
      call emreal1c_core_cu<<<gS,BLOCK_DIM,0,def_stream>>>
     &        ( ipole_s, pglob_s, loc_s, ieblst_s
     &        , eblst_s(2*npolelocnlb_pair+1)
     &        , npolelocnlb, npolelocnlb2_pair, npolebloc, n
     &        , x_s, y_s, z_s, rpole
     &        , off2, f, alsq2, alsq2n, aewald
     &        , dem, tem, ered_buff, vred_buff
     &        , p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
#ifdef TINKER_DEBUG
     &        , inter
#endif
     &        )
      call check_launch_kernel(" emreal1c_core_cu")

!$acc end host_data

      call reduce_energy_virial(em,vxx,vxy,vxz,vyy,vyz,vzz,def_queue)
#else
 100  format('emreal1c_core5 is a specific device routine !!',/,
     &       'you are not supposed to get inside with your compile ',
     &       'mode.')
      print 100
      call fatal
#endif

#ifdef TINKER_DEBUG
 34   format(2I10,3F12.4)
 35   format(A30,2I16)
!$acc wait
!$acc exit data copyout(inter)
      print 35,'npole block inter',npolelocnlb_pair,npolelocnlb2_pair
      print 35,'total mpole interactions ', sum(inter)
c      do i =1,250000
c         print 34,i,inter(i)
c      end do
#endif

      call emreal_correct_interactions(tem,vxx,vxy,vxz,vyy,vyz,vzz)

      end




c======================================================================
c         ===================================================
c                   ===============================
c
c     SECTION : short and long range interactions cores
c
c                   ===============================
c         ===================================================
c======================================================================


      ! Compute emreal interactions with scale set to 1 and procede to
      ! correction if necessary
      subroutine emrealshortlong1c_core(tem,vxx,vxy,vxz,vyy,vyz,vzz)
      use atmlst ,only:poleglobnl
      use atoms  ,only:x,y,z,n
      use cutoff ,only:mpoleshortcut,shortheal
      use chgpot ,only:electric,dielec
      use deriv  ,only:dem
      use domdec ,only:rank,nbloc,loc
      use empole1gpu_inl ,only: image_inl,mpole1_couple_shortlong
      use energi ,only:em
      use ewald  ,only:aewald
      use inform ,only: deb_Path
      use interfaces,only: emreal_correct_interactions_shortlong
     &              ,long_mode,short_mode
      use math   ,only:sqrtpi
      !use mpi
      use mplpot ,only:n_mscale,mcorrect_ik,mcorrect_scale
      use mpole  ,only:rpole,ipole,polelocnl,npolelocnl
      use neigh  ,only:nelst,elst,shortelst,nshortelst
      use potent ,only:use_mpoleshortreal,use_mpolelong
      use shunt  ,only:off,off2
      use tinheader ,only:ti_p
      use utilgpu,only:dir_queue,rec_queue,def_queue,warning,
     &                 real3,real3_red,rpole_elt,maxscaling
      use virial
      implicit none

      real(r_p),intent(inout):: vxx,vyy,vzz
      real(r_p),intent(inout):: vxy,vxz,vyz
      real(t_p),intent(inout):: tem(3,nbloc)

      integer i,j,k,iglob,kglob,kbis
      integer ii,kk,iipole,kkpole
      integer nnelst,mode
      integer,pointer,save::lst(:,:),nlst(:)
#ifdef TINKER_DEBUG
      integer inter(n)
#endif
      real(t_p) zero,r_cut
      real(t_p) r2,f,e
      real(t_p) alsq2,alsq2n
      real(t_p) xi,yi,zi
      real(t_p) xr,yr,zr
      type(rpole_elt) ip,kp
      type(real3_red) frc_r
      type(real3) ttmi,ttmk,frc
      real(t_p) mscale
      real(t_p) mpoleshortcut2
      parameter(zero=0.0)

      if(deb_Path) write(*,'(3x,a)') 'emrealshortlong1c_core'

#ifdef TINKER_DEBUG
!$acc enter data create(inter) async(def_queue)
!$acc parallel loop async present(inter)
      do i = 1,n
         inter(i) = 0
      end do
#endif

c
c     set conversion factor, cutoff and switching coefficients
c
      f      = electric / dielec
      alsq2  = 2 * aewald**2
      alsq2n = zero
      if (aewald .gt. zero)
     &   alsq2n = 1.0 / (sqrtpi*aewald)

      ! Configure data to be use in next loop
      if (use_mpoleshortreal) then
         mode  = short_mode
         r_cut = off
         mpoleshortcut2 = 0.0_ti_p
          lst =>  shortelst
         nlst => nshortelst
      else
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
!$acc&         present(tem,vxx,vxy,vxz,vyy,vyz,vzz)
!$acc&         present(dem,em,vir,poleglobnl,ipole,loc,x,y,z,rpole,
!$acc&  lst,nlst) private(ip)
#ifdef TINKER_DEBUG
!$acc&         present(inter)
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
!$acc loop vector private(kp,frc,frc_r,ttmi,ttmk)
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
            call mpole1_couple_shortlong(r2,xr,yr,zr,ip,kp,zero,
     &                    r_cut,shortheal,aewald,f,alsq2n,alsq2,
     &                       e,frc,frc_r,ttmi,ttmk,.false.,mode)

            ! update energy
            em       = em  + e

#ifdef TINKER_DEBUG
!$acc atomic
            inter(iglob) = inter(iglob)+1
#endif
c
c     increment force-based gradient and torque on first site
c
!$acc atomic update
            dem(1,i)    = dem(1,i)    + frc_r%x
!$acc atomic update
            dem(2,i)    = dem(2,i)    + frc_r%y
!$acc atomic update
            dem(3,i)    = dem(3,i)    + frc_r%z
!$acc atomic update
            dem(1,kbis) = dem(1,kbis) - frc_r%x
!$acc atomic update
            dem(2,kbis) = dem(2,kbis) - frc_r%y
!$acc atomic update
            dem(3,kbis) = dem(3,kbis) - frc_r%z
c
c     increment force-based gradient and torque on second site
c
!$acc atomic update
            tem(1,i)    = tem(1,i)    + ttmi%x
!$acc atomic update
            tem(2,i)    = tem(2,i)    + ttmi%y
!$acc atomic update
            tem(3,i)    = tem(3,i)    + ttmi%z
!$acc atomic update
            tem(1,kbis) = tem(1,kbis) + ttmk%x
!$acc atomic update
            tem(2,kbis) = tem(2,kbis) + ttmk%y
!$acc atomic update
            tem(3,kbis) = tem(3,kbis) + ttmk%z
c
c     increment the virial due to pairwise Cartesian forces
c
            vxx    = vxx  - xr * frc%x
            vxy    = vxy  - yr * frc%x
            vxz    = vxz  - zr * frc%x
            vyy    = vyy  - yr * frc%y
            vyz    = vyz  - zr * frc%y
            vzz    = vzz  - zr * frc%z
         end do
      end do
!$acc exit data detach(nlst,lst) async(def_queue)

#ifdef TINKER_DEBUG
 34   format(2I10,3F12.4)
 35   format(A30,I16)
!$acc wait
!$acc exit data copyout(inter)
      print 35, 'total mpole interactions ', sum(inter)
#endif

      call emreal_correct_interactions_shortlong
     &             (tem,vxx,vxy,vxz,vyy,vyz,vzz)

      end


      ! CUDA Fortran routine
      subroutine emrealshortlong1c_core2(tem,vxx,vxy,vxz,vyy,vyz,vzz)
      use atmlst ,only:poleglobnl
      use atoms  ,only:x,y,z,n
      use cutoff ,only:mpoleshortcut,shortheal
      use chgpot ,only:electric,dielec
      use cell
      use deriv  ,only:dem
      use domdec  ,only: xbegproc,ybegproc,zbegproc
     &            ,nproc,rank,xendproc,yendproc,zendproc
     &            ,nbloc,loc
#ifdef _CUDA
      use empole1cu ,only: emrealshortlong1c_core_cu
#endif
      use energi ,only:em
      use ewald  ,only:aewald
      use inform ,only: deb_Path
      use interfaces,only: emreal_correct_interactions_shortlong
     &              ,long_mode,short_mode
      use math   ,only:sqrtpi
      !use mpi
      use mplpot ,only:n_mscale,mcorrect_ik,mcorrect_scale
      use mpole  ,only:rpole,ipole,polelocnl,npolelocnl
     &           ,npolelocnlb,npolelocnlb_pair,npolebloc
     &           ,npolelocnlb2_pair,nshortpolelocnlb2_pair
      use neigh  ,only:ipole_s=>celle_glob,pglob_s=>celle_pole
     &           , loc_s=>celle_loc, ieblst_s=>ieblst, eblst_s=>eblst
     &           , x_s=>celle_x, y_s=>celle_y, z_s=>celle_z
     &           , seblst_s=>shorteblst,iseblst_s=>ishorteblst
      use potent ,only:use_mpoleshortreal,use_mpolelong
      use shunt  ,only:off,off2
      use tinheader,only: ti_p
#ifdef _CUDA
      use utilcu ,only:BLOCK_DIM,check_launch_kernel
#endif
      use utilgpu,only:dir_queue,rec_queue,def_queue,warning
     &           ,real3,real3_red,rpole_elt,RED_BUFF_SIZE
     &           ,ered_buff,vred_buff,reduce_energy_virial
     &           ,zero_evir_red_buffer
#ifdef  _OPENACC
     &           ,dir_stream,def_stream
#endif
      use virial
      implicit none

      real(r_p),intent(inout):: vxx,vyy,vzz
      real(r_p),intent(inout):: vxy,vxz,vyz
      real(t_p),intent(inout):: tem(3,nbloc)

      integer i,mode
#ifdef TINKER_DEBUG
      integer inter(n)
#endif
      real(t_p) r2,f,e
      real(t_p) alsq2,alsq2n
      real(t_p) p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
      real(t_p) mpoleshortcut2,r_cut
      logical,save:: first_in=.true.
      integer,save:: gS

      if(deb_Path) write(*,'(3x,a)') 'emrealshortlong1c_core2'

#ifdef TINKER_DEBUG
!$acc enter data create(inter) async(def_queue)
!$acc parallel loop async present(inter)
      do i = 1,n
         inter(i) = 0
      end do
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
         mode = short_mode
         mpoleshortcut2 = 0.0_ti_p
         r_cut = off
      else if (use_mpolelong) then
         mode = long_mode
         mpoleshortcut2 = (mpoleshortcut-shortheal)**2
         r_cut = mpoleshortcut
      end if

      if (first_in) then
#ifdef _CUDA
        call cudaMaxGridSize("emrealshortlong1c_core_cu",gS)
#endif
         first_in = .false.
      end if

#ifdef _CUDA
      def_stream = dir_stream
      call zero_evir_red_buffer(def_queue)


      if (use_mpoleshortreal) then
!$acc host_data use_device(ipole_s,pglob_s,loc_s,iseblst_s,seblst_s
!$acc&    ,x_s,y_s,z_s,rpole,dem,tem,ered_buff,vred_buff
#ifdef TINKER_DEBUG
!$acc&    ,inter
#endif
!$acc&    )
      !call CUDA kernel to compute the real space portion of the Ewald summation
      call emrealshortlong1c_core_cu<<<*,BLOCK_DIM,0,def_stream>>>
     &        ( ipole_s, pglob_s, loc_s, iseblst_s
     &        , seblst_s(2*npolelocnlb_pair+1)
     &        , npolelocnlb, nshortpolelocnlb2_pair, npolebloc, n
     &        , x_s, y_s, z_s, rpole
     &        , shortheal, r_cut, mpoleshortcut2, off2, f
     &        , alsq2, alsq2n, aewald, mode
     &        , dem, tem, ered_buff, vred_buff
     &        , p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
#ifdef TINKER_DEBUG
     &        , inter
#endif
     &        )
!$acc end host_data

      else
!$acc host_data use_device(ipole_s,pglob_s,loc_s,ieblst_s,eblst_s
!$acc&    ,x_s,y_s,z_s,rpole,dem,tem,ered_buff,vred_buff
#ifdef TINKER_DEBUG
!$acc&    ,inter
#endif
!$acc&    )
      !call CUDA kernel to compute the real space portion of the Ewald summation
      call emrealshortlong1c_core_cu<<<*,BLOCK_DIM,0,def_stream>>>
     &        ( ipole_s, pglob_s, loc_s, ieblst_s
     &        , eblst_s(2*npolelocnlb_pair+1)
     &        , npolelocnlb, npolelocnlb2_pair, npolebloc, n
     &        , x_s, y_s, z_s, rpole
     &        , shortheal, r_cut, mpoleshortcut2, off2, f
     &        , alsq2, alsq2n, aewald, mode
     &        , dem, tem, ered_buff, vred_buff
     &        , p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
#ifdef TINKER_DEBUG
     &        , inter
#endif
     &        )
!$acc end host_data
      end if
      call check_launch_kernel(" emrealshortlong1c_core_cu")


      call reduce_energy_virial(em,vxx,vxy,vxz,vyy,vyz,vzz,def_queue)
#else
 100  format('emrealshortlong1c_core2 is a specific device routine !!',
     &       /,'you are not supposed to get inside with your compile ',
     &       'mode.')
      print 100
      call fatal
#endif

#ifdef TINKER_DEBUG
 34   format(2I10,3F12.4)
 35   format(A30,2I16)
!$acc wait
!$acc exit data copyout(inter)
      if (use_mpoleshortreal) then
      print 35,'npole short block inter',npolelocnlb_pair,
     &         nshortpolelocnlb2_pair
      else
      print 35,'npole long block inter',npolelocnlb_pair,
     &         npolelocnlb2_pair
      end if
      print 35,'total mpole interactions ', sum(inter)
c      do i =1,250000
c         print 34,i,inter(i)
c      end do
#endif

      call emreal_correct_interactions_shortlong
     &            (tem,vxx,vxy,vxz,vyy,vyz,vzz)

      end



c======================================================================
c         ===================================================
c                   ===============================
c
c     SECTION : Scaling factors interactions
c
c                   ===============================
c         ===================================================
c======================================================================

      ! Procede to corrections of scaling interactions
      subroutine
     &   emreal_correct_interactions(tem,vxx,vxy,vxz,vyy,vyz,vzz)
      use atoms  ,only:x,y,z
      use chgpot ,only:electric,dielec
      use deriv  ,only:dem
      use domdec ,only:rank,nbloc,loc
      use empole1gpu_inl ,only: image_inl,mpole1_couple
      use energi ,only:em
      use ewald  ,only:aewald
      use inform ,only:deb_Path
      use math   ,only:sqrtpi
      use mplpot ,only:n_mscale,mcorrect_ik,mcorrect_scale
      use mpole  ,only:rpole,ipole,polelocnl,npolelocnl
      use shunt  ,only:off2
      use utilgpu,only:real3,real3_red,rpole_elt,def_queue
      use virial
      implicit none

      real(r_p),intent(inout):: vxx,vyy,vzz
      real(r_p),intent(inout):: vxy,vxz,vyz
      real(t_p),intent(inout):: tem(:,:)

      integer   i,j,k,iglob,kglob,kbis
      integer   ii,kk,iipole,kkpole
      real(t_p) zero
      real(t_p) r2,f,e
      real(t_p) alsq2,alsq2n
      real(t_p) xr,yr,zr
      type(rpole_elt) ip,kp
      type(real3_red) frc_r
      type(real3) ttmi,ttmk,frc
      real(t_p) mscale
      logical   corre
      parameter(zero=0.0)

      if (deb_Path)
     &   write(*,'(2x,a)') 'emreal_correct_interactions'
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
!$acc parallel loop present(tem,vxx,vxy,vxz,vyy,vyz,vzz)
!$acc&         present(dem,em,vir,ipole,loc,x,y,z,rpole,
!$acc&  mcorrect_ik,mcorrect_scale)
!$acc&         private(ip,kp,frc,frc_r,ttmi,ttmk)
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

c        corre  = .false.
c        mscale = 0.0
c20      continue
         ! compute mpole one interaction
         call mpole1_couple(r2,xr,yr,zr,ip,kp,mscale,
     &                      aewald,f,alsq2n,alsq2,
     &                      e,frc,frc_r,ttmi,ttmk,.true.)

         !if (.not.corre) then
         !    frc%x=-frc%x;   frc%y=-frc%y;   frc%z=-frc%z;
         !  frc_r%x=-frc_r%x; frc_r%y=-frc_r%y; frc_r%z=-frc_r%z;
         !   ttmi%x=-ttmi%x; ttmi%y=-ttmi%y; ttmi%z=-ttmi%z;
         !   ttmk%x=-ttmk%x; ttmk%y=-ttmk%y; ttmk%z=-ttmk%z;
         !   e = -e
         !end if
         ! update energy
         em     = em  + e
c
c     increment force-based gradient and torque on first site
c
!$acc atomic update
         dem(1,i)    = dem(1,i)    + frc_r%x
!$acc atomic update
         dem(2,i)    = dem(2,i)    + frc_r%y
!$acc atomic update
         dem(3,i)    = dem(3,i)    + frc_r%z
!$acc atomic update
         dem(1,kbis) = dem(1,kbis) - frc_r%x
!$acc atomic update
         dem(2,kbis) = dem(2,kbis) - frc_r%y
!$acc atomic update
         dem(3,kbis) = dem(3,kbis) - frc_r%z
c
c     increment force-based gradient and torque on second site
c
!$acc atomic update
         tem(1,i)    = tem(1,i)    + ttmi%x
!$acc atomic update
         tem(2,i)    = tem(2,i)    + ttmi%y
!$acc atomic update
         tem(3,i)    = tem(3,i)    + ttmi%z
!$acc atomic update
         tem(1,kbis) = tem(1,kbis) + ttmk%x
!$acc atomic update
         tem(2,kbis) = tem(2,kbis) + ttmk%y
!$acc atomic update
         tem(3,kbis) = tem(3,kbis) + ttmk%z
c
c     increment the virial due to pairwise Cartesian forces c
         vxx    = vxx  - xr * frc%x
         vxy    = vxy  - yr * frc%x
         vxz    = vxz  - zr * frc%x
         vyy    = vyy  - yr * frc%y
         vyz    = vyz  - zr * frc%y
         vzz    = vzz  - zr * frc%z

         !if (mscale.eq.0.0) then
         !   corre = .true.
         !   mscale = mcorrect_scale(ii)
         !   goto 20
         !end if
      end do

      end

      ! Procede to corrections of scaling interactions
      subroutine
     &   emreal_correct_interactions_shortlong
     &           (tem,vxx,vxy,vxz,vyy,vyz,vzz)
      use atoms  ,only:x,y,z
      use cutoff ,only:mpoleshortcut,shortheal
      use chgpot ,only:electric,dielec
      use deriv  ,only:dem
      use domdec ,only:rank,nbloc,loc
      use empole1gpu_inl ,only: image_inl,mpole1_couple_shortlong
      use energi ,only:em
      use ewald  ,only:aewald
      use inform ,only:deb_Path
      use interfaces ,only:long_mode,short_mode
      use math   ,only:sqrtpi
      use mplpot ,only:n_mscale,mcorrect_ik,mcorrect_scale
      use mpole  ,only:rpole,ipole,polelocnl,npolelocnl
      use potent ,only:use_mpoleshortreal,use_mpolelong
      use shunt  ,only:off2,off
      use tinheader ,only:ti_p
      use utilgpu,only:real3,real3_red,rpole_elt,def_queue
      use virial
      implicit none

      real(r_p),intent(inout):: vxx,vyy,vzz
      real(r_p),intent(inout):: vxy,vxz,vyz
      real(t_p),intent(inout):: tem(:,:)

      integer   i,j,k,iglob,kglob,kbis
      integer   ii,kk,iipole,kkpole
      integer   mode
      real(t_p) zero
      real(t_p) r2,f,e
      real(t_p) alsq2,alsq2n
      real(t_p) xr,yr,zr
      type(rpole_elt) ip,kp
      type(real3_red) frc_r
      type(real3) ttmi,ttmk,frc
      real(t_p) mpoleshortcut2,r_cut
      real(t_p) mscale
      logical   corre
      parameter(zero=0.0)

      if (deb_Path)
     &   write(*,'(2x,a)') 'emreal_correct_interactions_shortlong'
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
         mode = short_mode
         r_cut = off
      else if (use_mpolelong) then
         mpoleshortcut2 = (mpoleshortcut-shortheal)**2
         mode = long_mode
         r_cut = mpoleshortcut
      else
         print*,'unknown mode for emreal_correct_interactions_shortlong'
         call fatal
      end if

c
c     Apply correction to scaling factor interactions
c
!$acc parallel loop present(tem,vxx,vxy,vxz,vyy,vyz,vzz)
!$acc&         present(dem,em,vir,ipole,loc,x,y,z,rpole,
!$acc&  mcorrect_ik,mcorrect_scale)
!$acc&         private(ip,kp,frc,frc_r,ttmi,ttmk)
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
         call mpole1_couple_shortlong(r2,xr,yr,zr,ip,kp,mscale,
     &                   r_cut,shortheal,aewald,f,alsq2n,alsq2,
     &                       e,frc,frc_r,ttmi,ttmk,.true.,mode)

         ! update energy
         em     = em  + e
c
c     increment force-based gradient and torque on first site
c
!$acc atomic update
         dem(1,i)    = dem(1,i)    + frc_r%x
!$acc atomic update
         dem(2,i)    = dem(2,i)    + frc_r%y
!$acc atomic update
         dem(3,i)    = dem(3,i)    + frc_r%z
!$acc atomic update
         dem(1,kbis) = dem(1,kbis) - frc_r%x
!$acc atomic update
         dem(2,kbis) = dem(2,kbis) - frc_r%y
!$acc atomic update
         dem(3,kbis) = dem(3,kbis) - frc_r%z
c
c     increment force-based gradient and torque on second site
c
!$acc atomic update
         tem(1,i)    = tem(1,i)    + ttmi%x
!$acc atomic update
         tem(2,i)    = tem(2,i)    + ttmi%y
!$acc atomic update
         tem(3,i)    = tem(3,i)    + ttmi%z
!$acc atomic update
         tem(1,kbis) = tem(1,kbis) + ttmk%x
!$acc atomic update
         tem(2,kbis) = tem(2,kbis) + ttmk%y
!$acc atomic update
         tem(3,kbis) = tem(3,kbis) + ttmk%z
c
c     increment the virial due to pairwise Cartesian forces c
         vxx    = vxx  - xr * frc%x
         vxy    = vxy  - yr * frc%x
         vxz    = vxz  - zr * frc%x
         vyy    = vyy  - yr * frc%y
         vyz    = vyz  - zr * frc%y
         vzz    = vzz  - zr * frc%z
      end do

      end




c======================================================================
c         ===================================================
c                   ===============================
c
c     SECTION : Reciproqual part
c
c                   ===============================
c         ===================================================
c======================================================================

c
c
c     #######################################################################
c     ##                                                                   ##
c     ##  subroutine emrecip1gpu  --  PME recip multipole energy & derivs  ##
c     ##                                                                   ##
c     #######################################################################
c
c
c     "emrecip1" evaluates the reciprocal space portion of the particle
c     mesh Ewald summation energy and gradient due to multipoles
c
      subroutine emrecip1gpu
      use atmlst
      use atoms
      use bound
      use boxes
      use chgpot
      use deriv
      use domdec
      use energi
      use ewald
      use fft
      use inform    ,only: deb_Path
      use interfaces,only: torquegpu,fphi_mpole_site_p
     &              ,grid_mpole_site_p,bspline_fill_sitegpu
      use polar_temp,only: cmp=>fphid,fmp=>fphip !Use register pool
     &              , trqrec=>fuind
      use pme
      use pme1
      use math
      use mpole
      use potent
      use timestat
      use virial
      use utilgpu
      use utils   ,only:set_to_zero1,rpole_scale,
     &                  rpole_ind_extract,comput_norm
      use mpi
      implicit none
      integer status(MPI_STATUS_SIZE),tag,ierr
      integer i,j,k,ii,iipole,iglob
      integer k1,k2,k3
      integer m1,m2,m3
      integer ntot,nff
      integer nf1,nf2,nf3
      integer deriv1(10)
      integer deriv2(10)
      integer deriv3(10)
      integer gridin_size,gridout_size
      real(t_p) zero,one,two,half
      real(r_p) e
      real(t_p) eterm,f
      real(t_p) r1,r2,r3
      real(t_p) h1,h2,h3
      real(t_p) f1,f2,f3
      real(r_p) vxx,vyy,vzz
      real(r_p) vxy,vxz,vyz
      real(t_p) volterm,denom
      real(t_p) hsq,expterm
      real(t_p) term,pterm
      real(t_p) vterm,struc2
      real(t_p) trq(3),fix(3)
      real(t_p) fiy(3),fiz(3)
c     real(t_p) trqrec(3,npolerecloc)
c     real(t_p) cmp(10,max(npolerecloc,1)),
c    &          fmp(10,max(npolerecloc,1))
c     real(t_p) qgridmpi(2,n1mpimax,n2mpimax,n3mpimax,nrec_recep)
      integer reqsend(nproc),reqrec(nproc)
      integer req2send(nproc),req2rec(nproc)
      integer nprocloc,commloc,rankloc,proc
      real(8) time0,time1
      parameter(zero=0.0_ti_p,  one=1.0_ti_p,
     &           two=2.0_ti_p, half=0.5_ti_p)
c
c     indices into the electrostatic field array
c
      parameter(
     &  deriv1=(/ 2, 5,  8,  9, 11, 16, 18, 14, 15, 20 /),
     &  deriv2=(/ 3, 8,  6, 10, 14, 12, 19, 16, 20, 17 /),
     &  deriv3=(/ 4, 9, 10,  7, 15, 17, 13, 20, 18, 19 /))
c
      if (deb_Path) write(*,'(2x,a)') 'emrecip1gpu'
      call timer_enter( timer_emrecip )
      if (use_pmecore) then
        nprocloc = nrec
        commloc  = comm_rec
        rankloc  = rank_bis
      else
        nprocloc = nproc
        commloc  = MPI_COMM_WORLD
        rankloc  = rank
      end if
c
c     return if the Ewald coefficient is zero
c
c     if (rank.eq.0) then
c     print*,'emrecip'
c     print*,'in grid',n1mpimax,n2mpimax,n3mpimax,nrec_recep,nrec_send
c     print*,'out grid',isize2(rankloc+1),jsize2(rankloc+1),
c    $                  ksize2(rankloc+1)

c     gridin_size  = 2*(nrec_send+1)*n1mpimax*n2mpimax*n3mpimax
c     gridout_size = 2*isize2(rankloc+1)*jsize2(rankloc+1)*
c    $                     ksize2(rankloc+1)
c     print*,gridin_size,gridout_size
c     end if

      if (aewald .lt. 1.0d-6)  return
      f = electric / dielec
c
c     perform dynamic allocation of some global arrays
c
      j = max(npolerecloc,1)
      call mallocMpiGrid
      call prmem_request(cphirec,10,j,async=.true.)
      call prmem_request(fphirec,20,j,async=.true.)
      call prmem_request(cmp    ,10,j,async=.true.)
      call prmem_request(fmp    ,10,j,async=.true.)
      call prmem_request(trqrec , 3,j,async=.true.)

!$acc data create(e,vxx,vxy,vxz,vyy,vyz,vzz)
!$acc&     present(polerecglob,ipole,rpole,istart2,jstart2,
!$acc&  kstart2,iend2,jend2,kend2,qfac_2d,cphirec,fphirec,
!$acc&  qgridin_2d,vir,demrec,use_bounds,recip,cmp,fmp,trqrec)
!$acc&     present(emrec)
!$acc&     async(rec_queue)
c
c     zero out the temporary virial accumulation variables
c
!$acc serial async
      vxx = zero
      vxy = zero
      vxz = zero
      vyy = zero
      vyz = zero
      vzz = zero
      e   = zero
!$acc end serial
c
      call timer_enter( timer_other )
!$acc parallel loop collapse(2) async(rec_queue)
      do i=1,npolerecloc
         do j=1,20
            if (j.lt.10)
     &      cphirec(j,i) = zero
            fphirec(j,i) = zero
         end do
      end do
      call timer_exit ( timer_other,quiet_timers )
c
c     zero out pme grid
c
      call timer_enter( timer_grid1 )
      call set_to_zero1(qgridin_2d(1,1,1,1,1),
     &     2*n1mpimax*n2mpimax*n3mpimax*(nrec_send+1),rec_queue)
c
c     copy multipole moments and coordinates to local storage
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
        tag = nprocloc*prec_send(i) + rankloc + 1
!$acc host_data use_device(qgridin_2d)
        call MPI_ISEND(qgridin_2d(1,1,1,1,i+1),
     $   2*n1mpimax*n2mpimax*n3mpimax,MPI_TPREC,proc,tag,
     $   commloc,reqsend(i),ierr)
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
     $        qgridin_2d(1,1,1,1,1),qgridmpi(1,1,1,1,i),
     $        qgridin_2d(1,1,1,1,1))
      end do
      call timer_exit( timer_recreccomm,quiet_timers )

c     call MPI_ALLREDUCE(MPI_in_place,time0,1,MPI_TPREC,MPI_MAX,
c    &     MPI_COMM_WORLD,i)
c     call MPI_ALLREDUCE(MPI_in_place,time1,1,MPI_TPREC,MPI_SUM,
c    &     MPI_COMM_WORLD,i)
c     if (rank.eq.0) print*,' gridin Li ',time0
c     if (rank.eq.0) print*,' gridin L1 ',time1
c
c     Perform 3-D FFT forward transform
c
#ifdef _OPENACC
      call cufft2d_frontmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
     $                      n3mpimax)
#else
      call   fft2d_frontmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
     $                      n3mpimax)
#endif
c
c     initialize variables required for the scalar summation
c
      pterm   = (pi/aewald)**2
      volterm = pi * volbox
      nff     = nfft1 * nfft2
      nf1     = (nfft1+1) / 2
      nf2     = (nfft2+1) / 2
      nf3     = (nfft3+1) / 2
c
      call timer_enter( timer_scalar )
!$acc serial async(rec_queue)
      if ((istart2(rankloc+1).eq.1).and.
     &    (jstart2(rankloc+1).eq.1).and.
     &    (kstart2(rankloc+1).eq.1)) then
         qfac_2d(1,1,1) = zero
      end if
!$acc end serial
c
c     make the scalar summation over reciprocal lattice
c
!$acc parallel loop collapse(3)
!!$acc&         reduction(+:vxx,vxy,vxz,vyy,vyz,vzz)
!$acc&         async(rec_queue)
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
            expterm = zero
            if (term .gt. -50.0_ti_p) then
               denom = volterm*hsq*bsmod1(k1)*bsmod2(k2)*bsmod3(k3)
               expterm = exp(term) / denom
               if (.not. use_bounds) then
                  expterm = expterm * (1.0_ti_p-cos(pi*xbox*sqrt(hsq)))
               else if (octahedron) then
                  if (mod(m1+m2+m3,2) .ne. 0)  expterm = zero
               end if
               struc2 = qgridout_2d(1,k1-istart2(rankloc+1)+1,
     &          k2-jstart2(rankloc+1)+1,k3-kstart2(rankloc+1)+1)**2
     &                + qgridout_2d(2,k1-istart2(rankloc+1)+1,
     &          k2-jstart2(rankloc+1)+1,k3-kstart2(rankloc+1)+1)**2
               eterm = half * electric * expterm * struc2
               vterm = (two/hsq) * (1.0_ti_p-term) * eterm
               vxx = vxx + h1*h1*vterm - eterm
               vxy = vxy + h1*h2*vterm
               vxz = vxz + h1*h3*vterm
               vyy = vyy + h2*h2*vterm - eterm
               vyz = vyz + h2*h3*vterm
               vzz = vzz + h3*h3*vterm - eterm
             qfac_2d(k1-istart2(rankloc+1)+1,k2-jstart2(rankloc+1)+1,k3-
     $       kstart2(rankloc+1)+1) = expterm
            end if
 10         continue
          end do
        end do
      end do
c
c     save the virial for use in polarization computation
c
!$acc serial async(rec_queue)
      vmxx = vxx
      vmxy = vxy
      vmxz = vxz
      vmyy = vyy
      vmyz = vyz
      vmzz = vzz
c
c     account for zeroth grid point for nonperiodic system
c
      if ((istart2(rankloc+1).eq.1).and.
     &    (jstart2(rankloc+1).eq.1).and.
     &    (kstart2(rankloc+1).eq.1)) then
        if (.not. use_bounds) then
           expterm = half * pi / xbox
           struc2  = qgrid2in_2d(1,1,1,1,1)**2 +
     &               qgrid2in_2d(2,1,1,1,1)**2
           e       = f * expterm * struc2
           emrec   = emrec + e
        end if
      end if
!$acc end serial
!!$acc update host(vmxx,vmxy,vmxz,vmyy,vmyz,vmzz) async(rec_queue)
c
c     complete the transformation of the PME grid
c
!$acc parallel loop collapse(3) async(rec_queue)
      do k = 1, ksize2(rankloc+1)
         do j = 1, jsize2(rankloc+1)
            do i = 1, isize2(rankloc+1)
              term = qfac_2d(i,j,k)
              qgridout_2d(1,i,j,k) = term*qgridout_2d(1,i,j,k)
              qgridout_2d(2,i,j,k) = term*qgridout_2d(2,i,j,k)
            end do
         end do
      end do
      call timer_exit( timer_scalar,quiet_timers )

c     nff   = 2*isize2(rankloc+1)*jsize2(rankloc+1)*ksize2(rankloc+1)
c     time0 = comput_norm(qgridout_2d,nff,0)
c     time1 = comput_norm(qgridout_2d,nff,1)
c     call MPI_ALLREDUCE(MPI_in_place,time0,1,MPI_real8,MPI_MAX,
c    &     MPI_COMM_WORLD,i)
c     call MPI_ALLREDUCE(MPI_in_place,time1,1,MPI_real8,MPI_SUM,
c    &     MPI_COMM_WORLD,i)
c     if (rank.eq.0) print*,' gridout Li ',time0
c     if (rank.eq.0) print*,' gridout L1 ',time1
c
c     perform 3-D FFT backward transform and get potential
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
      call timer_enter( timer_recreccomm )
      do i = 1, nrec_send
         proc = prec_send(i)
         tag  = nprocloc*rankloc + prec_send(i) + 1
!$acc host_data use_device(qgridin_2d)
         call MPI_IRECV(qgridin_2d(1,1,1,1,i+1),
     $                  2*n1mpimax*n2mpimax*n3mpimax,MPI_TPREC,
     $                  prec_send(i),tag,commloc,req2rec(i),ierr)
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
     $                  2*n1mpimax*n2mpimax*n3mpimax,MPI_TPREC,
     $                  prec_recep(i),tag,commloc,req2send(i),ierr)
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

      call timer_enter ( timer_grid2 )
      call fphi_mpole_site_p
!$acc parallel loop collapse(2) async(rec_queue)
      do i = 1, npolerecloc
         do j = 1, 20
            fphirec(j,i) = electric * fphirec(j,i)
         end do
      end do
      call fphi_to_cphi_sitegpu(fphirec,cphirec)
      call timer_exit( timer_grid2,quiet_timers )
c
c     increment the permanent multipole energy and gradient
c
      call timer_enter( timer_fmanage )
!$acc parallel loop async(rec_queue)
      do i = 1, npolerecloc
         f1 = zero
         f2 = zero
         f3 = zero
!$acc loop seq
         do k = 1, 10
            e  = e  + fmp(k,i)*fphirec(       k ,i)
            f1 = f1 + fmp(k,i)*fphirec(deriv1(k),i)
            f2 = f2 + fmp(k,i)*fphirec(deriv2(k),i)
            f3 = f3 + fmp(k,i)*fphirec(deriv3(k),i)
         end do
         f1 = real(nfft1,t_p) * f1
         f2 = real(nfft2,t_p) * f2
         f3 = real(nfft3,t_p) * f3
         h1 = recip(1,1)*f1 + recip(1,2)*f2 + recip(1,3)*f3
         h2 = recip(2,1)*f1 + recip(2,2)*f2 + recip(2,3)*f3
         h3 = recip(3,1)*f1 + recip(3,2)*f2 + recip(3,3)*f3
         iipole = polerecglob(i)
         iglob  = ipole(iipole)
         ii     = locrec1(iglob)
         demrec(1,ii) = demrec(1,ii) + h1
         demrec(2,ii) = demrec(2,ii) + h2
         demrec(3,ii) = demrec(3,ii) + h3
      end do
c
c     distribute torques into the permanent multipole gradient
c
!$acc parallel loop async(rec_queue)
      do i = 1, npolerecloc
         iipole = polerecglob(i)
         trqrec(1,i) = cmp(4,i)*cphirec(3,i) - cmp(3,i)*cphirec(4,i)
     &               + two*(cmp(7,i)-cmp(6,i))*cphirec(10,i)
     &               + cmp(9,i)*cphirec(8,i) + cmp(10,i)*cphirec(6,i)
     &               - cmp(8,i)*cphirec(9,i) - cmp(10,i)*cphirec(7,i)
         trqrec(2,i) = cmp(2,i)*cphirec(4,i) - cmp(4,i)*cphirec(2,i)
     &               + two*(cmp(5,i)-cmp(7,i))*cphirec(9,i)
     &               + cmp(8,i)*cphirec(10,i) + cmp(9,i)*cphirec(7,i)
     &               - cmp(9,i)*cphirec(5,i) - cmp(10,i)*cphirec(8,i)
         trqrec(3,i) = cmp(3,i)*cphirec(2,i) - cmp(2,i)*cphirec(3,i)
     &               + two*(cmp(6,i)-cmp(5,i))*cphirec(8,i)
     &               + cmp(8,i)*cphirec(5,i) + cmp(10,i)*cphirec(9,i)
     &               - cmp(8,i)*cphirec(6,i) - cmp(9,i)*cphirec(10,i)
c
c     permanent multipole contribution to the internal virial
c
         vxx =vxx - cmp(2,i)*cphirec(2,i) - two*cmp(5,i)*cphirec(5,i)
     &            - cmp(8,i)*cphirec(8,i) - cmp(9,i)*cphirec(9,i)
         vxy =vxy - half*(cmp(3,i)*cphirec(2,i)+cmp(2,i)*cphirec(3,i))
     &            -(cmp(5,i)+cmp(6,i))*cphirec(8,i)
     &            -half*cmp(8,i)*(cphirec(5,i)+cphirec(6,i))
     &            -half*(cmp(9,i)*cphirec(10,i)+cmp(10,i)*cphirec(9,i))
         vxz =vxz - half*(cmp(4,i)*cphirec(2,i)+cmp(2,i)*cphirec(4,i))
     &            -(cmp(5,i)+cmp(7,i))*cphirec(9,i)
     &            -half*cmp(9,i)*(cphirec(5,i)+cphirec(7,i))
     &            -half*(cmp(8,i)*cphirec(10,i)+cmp(10,i)*cphirec(8,i))
         vyy =vyy - cmp(3,i)*cphirec(3,i) - two*cmp(6,i)*cphirec(6,i)
     &            - cmp(8,i)*cphirec(8,i) - cmp(10,i)*cphirec(10,i)
         vyz =vyz - half*(cmp(4,i)*cphirec(3,i)+cmp(3,i)*cphirec(4,i))
     &            - (cmp(6,i)+cmp(7,i))*cphirec(10,i)
     &            - half*cmp(10,i)*(cphirec(6,i)+cphirec(7,i))
     &            - half*(cmp(8,i)*cphirec(9,i)+cmp(9,i)*cphirec(8,i))
         vzz =vzz - cmp(4,i)*cphirec(4,i) - two*cmp(7,i)*cphirec(7,i)
     &            - cmp(9,i)*cphirec(9,i) - cmp(10,i)*cphirec(10,i)
      end do

      call torquegpu(npolerecloc,polerecglob,locrec1,
     &                trqrec,demrec,rec_queue)

      call timer_exit(timer_fmanage,quiet_timers )
c
c     increment the internal virial tensor components
c
c
c     Proceed to atomic update to avoid collision with direct queue
c     even if it's highly unlikely
c
c#ifdef _OPENACC
c      if (dir_queue.ne.rec_queue) then
c         call stream_wait_async(rec_stream,dir_stream,rec_event)
c      end if
c#endif
!$acc serial async(rec_queue)
      emrec    = emrec + half*e
      vir(1,1) = vir(1,1) + vxx
      vir(2,1) = vir(2,1) + vxy
      vir(3,1) = vir(3,1) + vxz
      vir(1,2) = vir(1,2) + vxy
      vir(2,2) = vir(2,2) + vyy
      vir(3,2) = vir(3,2) + vyz
      vir(1,3) = vir(1,3) + vxz
      vir(2,3) = vir(2,3) + vyz
      vir(3,3) = vir(3,3) + vzz
!$acc end serial

!$acc end data
c
      call timer_exit( timer_emrecip )
      end

c
c      subroutine emrecip1bgpu
c      use atmlst
c      use atoms
c      use bound
c      use boxes
c      use chgpot
c      use deriv
c      use domdec
c      use energi
c      use ewald
c      use fft
c      use inform ,only: deb_Path
c      use pme
c      use math
c      use mpole
c      use potent
c      use timestat
c      use virial
c      use utilgpu
c      use mpi
c      implicit none
c      integer status(MPI_STATUS_SIZE),tag,ierr
c      integer i,j,k,ii,iipole,iglob
c      integer k1,k2,k3
c      integer m1,m2,m3
c      integer ntot,nff
c      integer nf1,nf2,nf3
c      integer deriv1(10)
c      integer deriv2(10)
c      integer deriv3(10)
c      real(t_p) zero,one,two,half
c      real(t_p) e,eterm,f
c      real(t_p) r1,r2,r3
c      real(t_p) h1,h2,h3
c      real(t_p) f1,f2,f3
c      real(t_p) vxx,vyy,vzz
c      real(t_p) vxy,vxz,vyz
c      real(t_p) volterm,denom
c      real(t_p) hsq,expterm
c      real(t_p) term,pterm
c      real(t_p) vterm,struc2
c      real(t_p) trq(3),fix(3)
c      real(t_p) fiy(3),fiz(3)
c      real(t_p) trqrec(3)
c      real(t_p) cmp(10,max(npolerecloc,1)),
c     &          fmp(10,max(npolerecloc,1))
c      real(t_p) qgridmpi(2,n1mpimax,n2mpimax,n3mpimax,
c     &          nrec_recep)
c      integer reqsend(nproc),reqrec(nproc)
c      integer req2send(nproc),req2rec(nproc)
c      integer nprocloc,commloc,rankloc,proc
c      real(t_p) time0,time1
c      parameter(zero=0.0_ti_p,  one=1.0_ti_p,
c     &           two=2.0_ti_p, half=0.5_ti_p)
c!$acc routine(cmp_to_fmp_site)
c!$acc routine(bspline_fill_site)
c!$acc routine(grid_mpole_site)
c!$acc routine(fphi_mpole_site)
c!$acc routine(fphi_to_cphi_site)
c!$acc routine(torque_rec_acc)
cc
c      if(deb_Path) write(*,'(2x,a)') 'emrecip1bgpu'
c      call timer_enter( timer_emrecip )
c      if (use_pmecore) then
c        nprocloc = nrec
c        commloc  = comm_rec
c        rankloc  = rank_bis
c      else
c        nprocloc = nproc
c        commloc  = MPI_COMM_WORLD
c        rankloc  = rank
c      end if
cc
cc     indices into the electrostatic field array
cc
c      data deriv1  / 2, 5,  8,  9, 11, 16, 18, 14, 15, 20 /
c      data deriv2  / 3, 8,  6, 10, 14, 12, 19, 16, 20, 17 /
c      data deriv3  / 4, 9, 10,  7, 15, 17, 13, 20, 18, 19 /
cc
cc     return if the Ewald coefficient is zero
cc
c      if (aewald .lt. 1.0d-6)  return
c      f = electric / dielec
cc
cc     perform dynamic allocation of some global arrays
cc
c      if (allocated(cphirec)) then
c!$acc exit data delete(cphirec) async(rec_queue)
c         deallocate (cphirec)
c      end if
c      if (allocated(fphirec)) then
c!$acc exit data delete(fphirec) async(rec_queue)
c         deallocate(fphirec)
c      end if
c      allocate (cphirec(10,max(npolerecloc,1)))
c      allocate (fphirec(20,max(npolerecloc,1)))
c
cc
cc     zero out the temporary virial accumulation variables
cc
c      vxx = zero
c      vxy = zero
c      vxz = zero
c      vyy = zero
c      vyz = zero
c      vzz = zero
c      e   = zero
c
c!$acc enter data create(cphirec,fphirec,qgridmpi,cmp,fmp,trqrec,
c!$acc& vmxx,vmxy,vmxz,vmyy,vmyz,vmzz)
c!$acc&           copyin(vxx,vxy,vxz,vyy,vyz,vzz,e,
c!$acc& deriv1,deriv2,deriv3)
c!$acc&           async(rec_queue)
c
c!$acc data present(polerecglob,ipole,rpole,istart2,jstart2,
c!$acc&  kstart2,iend2,jend2,kend2,cphirec,fphirec,qfac_2d,
c!$acc&  qgridin_2d,qgridmpi,cmp,fmp,vir,demrec,use_bounds)
c!$acc&     present(vmxx,vmxy,vmxz,vmyy,vmyz,vmzz,vxx,vxy,vxz,vyy,
c!$acc&  vyz,vzz,emrec,e,deriv1,deriv2,deriv3,trqrec)
c
cc
c!$acc parallel loop collapse(2) async(rec_queue)
c      do i=1,npolerecloc
c         do j=1,20
c            if (j.lt.10)
c     &      cphirec(j,i) = zero
c            fphirec(j,i) = zero
c         end do
c      end do
cc
cc     zero out pme grid
cc
c!$acc parallel loop collapse(5) async(rec_queue)
c      do i=1,nrec_send+1
c         do j = 1,n3mpimax
c            do k = 1,n2mpimax
c               do k1= 1,n1mpimax
c                  do k2 = 1,2
c             qgridin_2d(k2,k1,k,j,i) = zero
c                  end do
c               end do
c            end do
c         end do
c      end do
cc
cc     copy multipole moments and coordinates to local storage
cc
c!$acc parallel loop async(rec_queue)
c      do i = 1, npolerecloc
c         iipole   = polerecglob(i)
c         iglob    = ipole(iipole)
c         cmp(1,i) = rpole(1,iipole)
c         cmp(2,i) = rpole(2,iipole)
c         cmp(3,i) = rpole(3,iipole)
c         cmp(4,i) = rpole(4,iipole)
c         cmp(5,i) = rpole(5,iipole)
c         cmp(6,i) = rpole(9,iipole)
c         cmp(7,i) = rpole(13,iipole)
c         cmp(8,i) = two * rpole(6,iipole)
c         cmp(9,i) = two * rpole(7,iipole)
c         cmp(10,i) = two * rpole(10,iipole)
c         call cmp_to_fmp_site(cmp(1,i),fmp(1,i))
c         call bspline_fill_site(iglob,i)
c      end do
cc
cc     assign permanent multipoles to PME grid and perform
cc     the 3-D FFT forward transformation
cc
c      time0 = mpi_wtime()
c!$acc parallel loop num_gangs(npolerecloc/(128*10))
c!$acc&              async(rec_queue)
c      do i = 1, npolerecloc
c        iipole = polerecglob(i)
c        iglob  = ipole(iipole)
c        call grid_mpole_site(iglob,i,fmp(1,i))
c      end do
c
c      time1 = mpi_wtime()
c      timegrid1 = timegrid1 + time1-time0
cc
cc     MPI : Begin reception
cc
c      time0 = mpi_wtime()
c      do i = 1, nrec_recep
c        tag = nprocloc*rankloc + prec_recep(i) + 1
c!$acc host_data use_device(qgridmpi)
c        call MPI_IRECV(qgridmpi(1,1,1,1,i),2*n1mpimax*n2mpimax*
c     $   n3mpimax,MPI_TPREC,prec_recep(i),tag,
c     $   commloc,reqrec(i),ierr)
c!$acc end host_data
c      end do
cc
cc     MPI : begin sending
cc
c      do i = 1, nrec_send
c!$acc wait(rec_queue)
c        proc = prec_send(i)
c        tag = nprocloc*prec_send(i) + rankloc + 1
c!$acc host_data use_device(qgridin_2d)
c        call MPI_ISEND(qgridin_2d(1,1,1,1,i+1),
c     $   2*n1mpimax*n2mpimax*n3mpimax,MPI_TPREC,proc,tag,
c     $   commloc,reqsend(i),ierr)
c!$acc end host_data
c      end do
cc
c      do i = 1, nrec_recep
c         call MPI_WAIT(reqrec(i),status,ierr)
c      end do
c      do i = 1, nrec_send
c         call MPI_WAIT(reqsend(i),status,ierr)
c      end do
cc
cc     do the reduction 'by hand'
cc
c      do i = 1, nrec_recep
c         call aaddgpuAsync(2*n1mpimax*n2mpimax*n3mpimax,
c     $        qgridin_2d(1,1,1,1,1),
c     $        qgridmpi(1,1,1,1,i),qgridin_2d(1,1,1,1,1))
c      end do
c      time1 = mpi_wtime()
c      timerecreccomm = timerecreccomm + time1 - time0
cc
cc     Perform 3-D FFT forward transform
cc
c      time0 = mpi_wtime()
c#ifdef _OPENACC
c      call cufft2d_frontmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
c     $ n3mpimax)
c#else
c      call fft2d_frontmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
c     $ n3mpimax)
c#endif
c      time1 = mpi_wtime()
c      timeffts = timeffts + time1-time0
cc
cc     initialize variables required for the scalar summation
cc
c      pterm = (pi/aewald)**2
c      volterm = pi * volbox
c      nff = nfft1 * nfft2
c      nf1 = (nfft1+1) / 2
c      nf2 = (nfft2+1) / 2
c      nf3 = (nfft3+1) / 2
cc
c      time0 = mpi_wtime()
c!$acc kernels async(rec_queue)
c      if ((istart2(rankloc+1).eq.1).and.
c     &    (jstart2(rankloc+1).eq.1).and.
c     $    (kstart2(rankloc+1).eq.1)) then
c         qfac_2d(1,1,1) = zero
c      end if
cc
cc     make the scalar summation over reciprocal lattice
cc
c!$acc loop collapse(3)
c      do k3 = kstart2(rankloc+1),kend2(rankloc+1)
c        do k2 = jstart2(rankloc+1),jend2(rankloc+1)
c          do k1 = istart2(rankloc+1),iend2(rankloc+1)
c            m1 = k1 - 1
c            m2 = k2 - 1
c            m3 = k3 - 1
c            if (k1 .gt. nf1)  m1 = m1 - nfft1
c            if (k2 .gt. nf2)  m2 = m2 - nfft2
c            if (k3 .gt. nf3)  m3 = m3 - nfft3
c            if ((m1.eq.0).and.(m2.eq.0).and.(m3.eq.0)) goto 10
c            r1 = real(m1,t_p)
c            r2 = real(m2,t_p)
c            r3 = real(m3,t_p)
c            h1 = recip(1,1)*r1 + recip(1,2)*r2 + recip(1,3)*r3
c            h2 = recip(2,1)*r1 + recip(2,2)*r2 + recip(2,3)*r3
c            h3 = recip(3,1)*r1 + recip(3,2)*r2 + recip(3,3)*r3
c            hsq = h1*h1 + h2*h2 + h3*h3
c            term = -pterm * hsq
c            expterm = zero
c            if (term .gt. -50.0_ti_p) then
c               denom = volterm*hsq*bsmod1(k1)*bsmod2(k2)*bsmod3(k3)
c               expterm = exp(term) / denom
c               if (.not. use_bounds) then
c                  expterm = expterm * (1.0_ti_p-cos(pi*xbox*sqrt(hsq)))
c               else if (octahedron) then
c                  if (mod(m1+m2+m3,2) .ne. 0)  expterm = zero
c               end if
c               struc2 = qgridout_2d(1,k1-istart2(rankloc+1)+1,
c     &          k2-jstart2(rankloc+1)+1,k3-kstart2(rankloc+1)+1)**2
c     &          + qgridout_2d(2,k1-istart2(rankloc+1)+1,k2-
c     &          jstart2(rankloc+1)+1,k3-kstart2(rankloc+1)+1)**2
c               eterm = half * electric * expterm * struc2
c               vterm = (two/hsq) * (1.0_ti_p-term) * eterm
c               vxx = vxx + h1*h1*vterm - eterm
c               vxy = vxy + h1*h2*vterm
c               vxz = vxz + h1*h3*vterm
c               vyy = vyy + h2*h2*vterm - eterm
c               vyz = vyz + h2*h3*vterm
c               vzz = vzz + h3*h3*vterm - eterm
c             qfac_2d(k1-istart2(rankloc+1)+1,k2-jstart2(rankloc+1)+1,k3-
c     $       kstart2(rankloc+1)+1) = expterm
c           end if
c 10        continue
c          end do
c        end do
c      end do
cc
cc     save the virial for use in polarization computation
cc
c      vmxx = vxx
c      vmxy = vxy
c      vmxz = vxz
c      vmyy = vyy
c      vmyz = vyz
c      vmzz = vzz
cc
cc     account for zeroth grid point for nonperiodic system
cc
c      if ((istart2(rankloc+1).eq.1).and.
c     &    (jstart2(rankloc+1).eq.1).and.
c     $    (kstart2(rankloc+1).eq.1)) then
c        if (.not. use_bounds) then
c           expterm = half * pi / xbox
c           struc2 = qgrid2in_2d(1,1,1,1,1)**2 +
c     $       qgrid2in_2d(2,1,1,1,1)**2
c           e = f * expterm * struc2
c           emrec = emrec + e
c        end if
c      end if
cc
cc     complete the transformation of the PME grid
cc
c!$acc loop collapse(3)
c      do k = 1, ksize2(rankloc+1)
c         do j = 1, jsize2(rankloc+1)
c           do i = 1, isize2(rankloc+1)
c              term = qfac_2d(i,j,k)
c              qgridout_2d(1,i,j,k) = term*qgridout_2d(1,i,j,k)
c              qgridout_2d(2,i,j,k) = term*qgridout_2d(2,i,j,k)
c            end do
c         end do
c      end do
c!$acc end kernels
c      time1 = mpi_wtime()
c      timescalar = timescalar + time1-time0
cc
cc     perform 3-D FFT backward transform and get potential
cc
c      time0 = mpi_wtime()
c#ifdef _OPENACC
c      call cufft2d_backmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
c     $ n3mpimax)
c#else
c      call fft2d_backmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
c     $ n3mpimax)
c#endif
c      time1 = mpi_wtime()
c      timeffts = timeffts + time1-time0
cc
cc     MPI : Begin reception
cc
c      time0 = mpi_wtime()
c      do i = 1, nrec_send
c        proc = prec_send(i)
c        tag  = nprocloc*rankloc + prec_send(i) + 1
c!$acc host_data use_device(qgridin_2d)
c        call MPI_IRECV(qgridin_2d(1,1,1,1,i+1),
c     $   2*n1mpimax*n2mpimax*n3mpimax,MPI_TPREC,
c     $   prec_send(i),tag,commloc,req2rec(i),ierr)
c!$acc end host_data
c      end do
cc
cc     MPI : begin sending
cc
c      do i = 1, nrec_recep
c!$acc wait(rec_queue)
c        tag = nprocloc*prec_recep(i) + rankloc + 1
c!$acc host_data use_device(qgridin_2d)
c        call MPI_ISEND(qgridin_2d,
c     $   2*n1mpimax*n2mpimax*n3mpimax,
c     $   MPI_TPREC,prec_recep(i),tag,commloc,req2send(i),
c     $   ierr)
c!$acc end host_data
c      end do
cc
c      do i = 1, nrec_send
c        call MPI_WAIT(req2rec(i),status,ierr)
c      end do
c      do i = 1, nrec_recep
c        call MPI_WAIT(req2send(i),status,ierr)
c      end do
c      time1 = mpi_wtime()
c      timerecreccomm = timerecreccomm + time1-time0
c
c      time0 = mpi_wtime()
c!$acc parallel loop async(rec_queue)
c      do i = 1, npolerecloc
c        iipole = polerecglob(i)
c        iglob = ipole(iipole)
c        call fphi_mpole_site(iglob,i)
c      end do
c      time1 = mpi_wtime()
c      timegrid2 = timegrid2 + time1-time0
c!$acc parallel loop async(rec_queue)
c      do i = 1, npolerecloc
c!$acc loop seq
c         do j = 1, 20
c            fphirec(j,i) = electric * fphirec(j,i)
c         end do
c         call fphi_to_cphi_site (fphirec(1,i),cphirec(1,i))
c      end do
cc
cc     increment the permanent multipole energy and gradient
cc
c!$acc parallel loop async(rec_queue)
c      do i = 1, npolerecloc
c         f1 = zero
c         f2 = zero
c         f3 = zero
c!$acc loop seq
c         do k = 1, 10
c            e  = e  + fmp(k,i)*fphirec(       k ,i)
c            f1 = f1 + fmp(k,i)*fphirec(deriv1(k),i)
c            f2 = f2 + fmp(k,i)*fphirec(deriv2(k),i)
c            f3 = f3 + fmp(k,i)*fphirec(deriv3(k),i)
c         end do
c         f1 = real(nfft1,t_p) * f1
c         f2 = real(nfft2,t_p) * f2
c         f3 = real(nfft3,t_p) * f3
c         h1 = recip(1,1)*f1 + recip(1,2)*f2 + recip(1,3)*f3
c         h2 = recip(2,1)*f1 + recip(2,2)*f2 + recip(2,3)*f3
c         h3 = recip(3,1)*f1 + recip(3,2)*f2 + recip(3,3)*f3
c         iipole = polerecglob(i)
c         iglob  = ipole(iipole)
c         ii     = locrec1(iglob)
c         demrec(1,ii) = demrec(1,ii) + h1
c         demrec(2,ii) = demrec(2,ii) + h2
c         demrec(3,ii) = demrec(3,ii) + h3
c      end do
cc
cc     distribute torques into the permanent multipole gradient
cc
c!$acc parallel loop async(rec_queue)
c!$acc&         private(trq,fix,fiy,fiz)
c      do i = 1, npolerecloc
c         iipole = polerecglob(i)
c         trq(1) = cmp(4,i)*cphirec(3,i) - cmp(3,i)*cphirec(4,i)
c     &               + two*(cmp(7,i)-cmp(6,i))*cphirec(10,i)
c     &               + cmp(9,i)*cphirec(8,i) + cmp(10,i)*cphirec(6,i)
c     &               - cmp(8,i)*cphirec(9,i) - cmp(10,i)*cphirec(7,i)
c         trq(2) = cmp(2,i)*cphirec(4,i) - cmp(4,i)*cphirec(2,i)
c     &               + two*(cmp(5,i)-cmp(7,i))*cphirec(9,i)
c     &               + cmp(8,i)*cphirec(10,i) + cmp(9,i)*cphirec(7,i)
c     &               - cmp(9,i)*cphirec(5,i) - cmp(10,i)*cphirec(8,i)
c         trq(3) = cmp(3,i)*cphirec(2,i) - cmp(2,i)*cphirec(3,i)
c     &               + two*(cmp(6,i)-cmp(5,i))*cphirec(8,i)
c     &               + cmp(8,i)*cphirec(5,i) + cmp(10,i)*cphirec(9,i)
c     &               - cmp(8,i)*cphirec(6,i) - cmp(9,i)*cphirec(10,i)
c         call torque_rec_acc (iipole,trq,fix,fiy,fiz,demrec)
c      end do
c
c!$acc parallel loop async(rec_queue)
c      do i = 1, npolerecloc
c         iipole = polerecglob(i)
cc        trqrec(1,i) = cmp(4,i)*cphirec(3,i) - cmp(3,i)*cphirec(4,i)
cc    &               + two*(cmp(7,i)-cmp(6,i))*cphirec(10,i)
cc    &               + cmp(9,i)*cphirec(8,i) + cmp(10,i)*cphirec(6,i)
cc    &               - cmp(8,i)*cphirec(9,i) - cmp(10,i)*cphirec(7,i)
cc        trqrec(2,i) = cmp(2,i)*cphirec(4,i) - cmp(4,i)*cphirec(2,i)
cc    &               + two*(cmp(5,i)-cmp(7,i))*cphirec(9,i)
cc    &               + cmp(8,i)*cphirec(10,i) + cmp(9,i)*cphirec(7,i)
cc    &               - cmp(9,i)*cphirec(5,i) - cmp(10,i)*cphirec(8,i)
cc        trqrec(3,i) = cmp(3,i)*cphirec(2,i) - cmp(2,i)*cphirec(3,i)
cc    &               + two*(cmp(6,i)-cmp(5,i))*cphirec(8,i)
cc    &               + cmp(8,i)*cphirec(5,i) + cmp(10,i)*cphirec(9,i)
cc    &               - cmp(8,i)*cphirec(6,i) - cmp(9,i)*cphirec(10,i)
cc
cc     permanent multipole contribution to the internal virial
cc
c         vxx =vxx - cmp(2,i)*cphirec(2,i) - two*cmp(5,i)*cphirec(5,i)
c     &            - cmp(8,i)*cphirec(8,i) - cmp(9,i)*cphirec(9,i)
c         vxy =vxy - half*(cmp(3,i)*cphirec(2,i)+cmp(2,i)*cphirec(3,i))
c     &            -(cmp(5,i)+cmp(6,i))*cphirec(8,i)
c     &            -half*cmp(8,i)*(cphirec(5,i)+cphirec(6,i))
c     &            -half*(cmp(9,i)*cphirec(10,i)+cmp(10,i)*cphirec(9,i))
c         vxz =vxz - half*(cmp(4,i)*cphirec(2,i)+cmp(2,i)*cphirec(4,i))
c     &            -(cmp(5,i)+cmp(7,i))*cphirec(9,i)
c     &            -half*cmp(9,i)*(cphirec(5,i)+cphirec(7,i))
c     &            -half*(cmp(8,i)*cphirec(10,i)+cmp(10,i)*cphirec(8,i))
c         vyy =vyy - cmp(3,i)*cphirec(3,i) - two*cmp(6,i)*cphirec(6,i)
c     &            - cmp(8,i)*cphirec(8,i) - cmp(10,i)*cphirec(10,i)
c         vyz =vyz - half*(cmp(4,i)*cphirec(3,i)+cmp(3,i)*cphirec(4,i))
c     &            - (cmp(6,i)+cmp(7,i))*cphirec(10,i)
c     &            - half*cmp(10,i)*(cphirec(6,i)+cphirec(7,i))
c     &            - half*(cmp(8,i)*cphirec(9,i)+cmp(9,i)*cphirec(8,i))
c         vzz =vzz - cmp(4,i)*cphirec(4,i) - two*cmp(7,i)*cphirec(7,i)
c     &            - cmp(9,i)*cphirec(9,i) - cmp(10,i)*cphirec(10,i)
c      end do
c
cc
cc     increment the internal virial tensor components
cc
c!$acc kernels async(rec_queue)
c      emrec    = emrec + half*e
c      vir(1,1) = vir(1,1) + vxx
c      vir(2,1) = vir(2,1) + vxy
c      vir(3,1) = vir(3,1) + vxz
c      vir(1,2) = vir(1,2) + vxy
c      vir(2,2) = vir(2,2) + vyy
c      vir(3,2) = vir(3,2) + vyz
c      vir(1,3) = vir(1,3) + vxz
c      vir(2,3) = vir(2,3) + vyz
c      vir(3,3) = vir(3,3) + vzz
c!$acc end kernels
c
c!$acc end data
c!$acc exit data delete(qgridmpi,cmp,fmp,trqrec)
c!$acc&          delete(vxx,vxy,vxz,vyy,vyz,vzz,e,deriv1,
c!$acc&  deriv2,deriv3)
c!$acc&          copyout(vmxx,vmxy,vmxz,vmyy,vmyz,vmzz)
c!$acc&          async(rec_queue)
cc
c      call timer_exit( timer_emrecip )
c      end
