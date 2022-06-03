c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine epolar1  --  polarization energy & derivs  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "epolar1" calculates the induced dipole polarization energy
c     and derivatives with respect to Cartesian coordinates
c
c
#include "tinker_precision.h"
#include "tinker_types.h"
      module epolar1gpu_inl
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

      subroutine epolar1gpu
      use polpot
      use potent
      implicit none
c
c     choose the method for summing over polarization interactions
c
      if (use_lambdadyn) then
        call elambdapolar1cgpu
      else
        if (use_polarshortreal) then
          if (polalgshort.eq.3) then
            !call epolar1tcggpu !FIXME
          else
            call epolar1cgpu
          end if
        else
          if (polalg.eq.3) then
            !call epolar1tcggpu
          else
            call epolar1cgpu
          end if
        end if
      end if
      end
c
c
c     ###################################################################
c     ##                                                               ##
c     ##  subroutine epolar1c  --  Ewald polarization derivs via list  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     "epolar1c" calculates the dipole polarization energy and
c     derivatives with respect to Cartesian coordinates using
c     particle mesh Ewald summation and a neighbor list
c
c
      subroutine epolar1cgpu
      use atmlst
      use atoms
      use boxes
      use chgpot
      use deriv
      use domdec
      use energi
      use epolar1gpu_inl
      use ewald
      use inform ,only: deb_Path
      use interfaces,only:torquegpu,epreal1c_p
     &              ,epreal1c_core3
      use iounit
      use math
      use mpole
      use elec_wspace   ,only: trq=>r2Work4
      use pme
      use polar
      use polpot
      use potent
      use precompute_pole,only:precompute_tmat,
     &                         precompute_solvpole
      use virial
      use mpi
      use sizes
      use timestat
      use utils
      use utilgpu
      use vec

      implicit none
      integer i,j,ii,iglob,iipole,ierr
      ener_rtyp,save:: epself
      real(t_p) f
      real(t_p) term,term1,fterm
      real(t_p) dix,diy,diz
      real(t_p) uix,uiy,uiz,uii
      real(t_p) xd,yd,zd
      real(t_p) xq,yq,zq
      real(t_p) xu,yu,zu
      real(t_p) xup,yup,zup
      real(t_p) xv,yv,zv,vterm
      real(t_p) xufield,yufield
      real(t_p) zufield
      real(t_p) time0,time1,time2
      real(t_p) fix(3),fiy(3),fiz(3)
      logical,save:: f_in=.true.
c
      if (npole .eq. 0)  return
      if(deb_Path)write(*,*) 'epolar1cgpu'
      if (f_in) then
         f_in=.false.
!$acc enter data create(epself)
      end if

c
!$acc serial async(rec_queue) present(epself,ep,ep_r,eprec)
      epself= 0
      ep    = 0.0_re_p
      ep_r  = 0
      eprec = 0.0_re_p
!$acc end serial
c
c     set the energy unit conversion factor
c
      f = electric / dielec

c
c     compute the induced dipoles at each polarizable atom
c
      call timer_enter( timer_polarsolve )
      if (use_polarshortreal) then
        if (polalg.eq.5) then
#ifdef _OPENACC
          ! FIXME
          write(0,*) 'short dcinduce solver unavalaible on device'
          __TINKER_FATAL__
#else
          call dcinduce_shortrealgpu
#endif
        else
          call newinduce_shortrealgpu
        end if
      else if (use_pmecore) then
        if (polalg.eq.5) then
#ifdef _OPENACC
          ! FIXME
          write(0,*) 'dcinduce solver (pme-core) unavalaible on device'
          __TINKER_FATAL__
#else
          call dcinduce_pme
#endif
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
      call timer_exit( timer_polarsolve )
c
c     Reset precompute switch if necessary
c
      if (precompute_solvpole) precompute_tmat=.true.

#ifdef _OPENACC
      if (dir_queue.ne.rec_queue) call start_dir_stream_cover
#endif
c
c     compute the real space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.le.ndir-1))
     &   then
         call timer_enter( timer_real )
         if (use_polarshortreal) then
            call epreal1cgpu
         else
            call epreal1c_p
         end if

         if (use_pself) then
            def_queue = dir_queue
#ifdef _OPENACC
           if (dir_queue.ne.rec_queue) call start_dir_stream_cover
#endif
           call prmem_request(trq,3,npoleloc,queue=def_queue)
c
c     compute the Ewald self-energy term over all the atoms
c
           term  = 2.0_ti_p * aewald * aewald
           fterm = -f * aewald / sqrtpi
           term1 = (4.0_ti_p/3.0_ti_p) * f * aewald**3 / sqrtpi
!$acc parallel loop async(def_queue) default(present) present(epself)
!$acc&         reduction(+:epself)
           do ii = 1, npoleloc
              iipole   = poleglob(ii)
              dix      = rpole(2,iipole)
              diy      = rpole(3,iipole)
              diz      = rpole(4,iipole)
              uix      = uind (1,iipole)
              uiy      = uind (2,iipole)
              uiz      = uind (3,iipole)
              uii      = dix*uix + diy*uiy + diz*uiz
              epself   = epself + tp2enr(fterm*term*uii / 3.0_ti_p)
c
c     compute the self-energy torque term due to induced dipole
c
              uix       = 0.5_ti_p * (uix     + uinp(1,iipole))
              uiy       = 0.5_ti_p * (uiy     + uinp(2,iipole))
              uiz       = 0.5_ti_p * (uiz     + uinp(3,iipole))
              trq(1,ii) = term1 * (diy*uiz - diz*uiy)
              trq(2,ii) = term1 * (diz*uix - dix*uiz)
              trq(3,ii) = term1 * (dix*uiy - diy*uix)
           end do

           call torquegpu(npoleloc,nbloc,poleglob,loc,trq,dep,def_queue)
c
c     compute the cell dipole boundary correction term
c
           if (boundary .eq. 'VACUUM') then
c
            !TODO Offload this part on device
              xd  = 0.0_ti_p
              yd  = 0.0_ti_p
              zd  = 0.0_ti_p
              xu  = 0.0_ti_p
              yu  = 0.0_ti_p
              zu  = 0.0_ti_p
              xup = 0.0_ti_p
              yup = 0.0_ti_p
              zup = 0.0_ti_p
              do i = 1, npoleloc
                 iipole = poleglob(i)
                 iglob  = ipole(iipole)
               xd     = xd  + rpole(2,iipole) + rpole(1,iipole)*x(iglob)
               yd     = yd  + rpole(3,iipole) + rpole(1,iipole)*y(iglob)
               zd     = zd  + rpole(4,iipole) + rpole(1,iipole)*z(iglob)
                 xu     = xu  + uind(1,iipole)
                 yu     = yu  + uind(2,iipole)
                 zu     = zu  + uind(3,iipole)
                 xup    = xup + uinp(1,iipole)
                 yup    = yup + uinp(2,iipole)
                 zup    = zup + uinp(3,iipole)
              end do
              call MPI_ALLREDUCE(MPI_IN_PLACE,xd,1,MPI_TPREC,MPI_SUM,
     $           COMM_TINKER,ierr)
              call MPI_ALLREDUCE(MPI_IN_PLACE,yd,1,MPI_TPREC,MPI_SUM,
     $           COMM_TINKER,ierr)
              call MPI_ALLREDUCE(MPI_IN_PLACE,zd,1,MPI_TPREC,MPI_SUM,
     $           COMM_TINKER,ierr)
              call MPI_ALLREDUCE(MPI_IN_PLACE,xu,1,MPI_TPREC,MPI_SUM,
     $           COMM_TINKER,ierr)
              call MPI_ALLREDUCE(MPI_IN_PLACE,yu,1,MPI_TPREC,MPI_SUM,
     $           COMM_TINKER,ierr)
              call MPI_ALLREDUCE(MPI_IN_PLACE,zu,1,MPI_TPREC,MPI_SUM,
     $           COMM_TINKER,ierr)
              call MPI_ALLREDUCE(MPI_IN_PLACE,xup,1,MPI_TPREC,MPI_SUM,
     $           COMM_TINKER,ierr)
              call MPI_ALLREDUCE(MPI_IN_PLACE,yup,1,MPI_TPREC,MPI_SUM,
     $           COMM_TINKER,ierr)
              call MPI_ALLREDUCE(MPI_IN_PLACE,zup,1,MPI_TPREC,MPI_SUM,
     $           COMM_TINKER,ierr)
              term = (2.0_ti_p/3.0_ti_p) * f * (pi/volbox)
              if (rank.eq.0) then
                 epself = epself  + tp2enr(term*(xd*xu+yd*yu+zd*zu))
              end if
              do ii = 1, npoleloc
                 iipole   = poleglob(ii)
                 iglob    = ipole(iipole)
                 i        = loc(iglob)
                 dep(1,i) = dep(1,i) + term*rpole(1,iipole)*(xu+xup)
                 dep(2,i) = dep(2,i) + term*rpole(1,iipole)*(yu+yup)
                 dep(3,i) = dep(3,i) + term*rpole(1,iipole)*(zu+zup)
              end do
              xufield = -term * (xu+xup)
              yufield = -term * (yu+yup)
              zufield = -term * (zu+zup)
              do i = 1, npoleloc
                 iipole   = poleglob(i)
                 trq(1,i) =  rpole(3,iipole)*zufield
     &                     - rpole(4,iipole)*yufield
                 trq(2,i) =  rpole(4,iipole)*xufield
     &                     - rpole(2,iipole)*zufield
                 trq(3,i) =  rpole(2,iipole)*yufield
     &                     - rpole(3,iipole)*xufield
              end do
              call torquegpu(npoleloc,nbloc
     &                      ,poleglob,loc,trq,dep,def_queue)
c
c     boundary correction to virial due to overall cell dipole
c
              xd = 0.0_ti_p
              yd = 0.0_ti_p
              zd = 0.0_ti_p
              xq = 0.0_ti_p
              yq = 0.0_ti_p
              zq = 0.0_ti_p
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
              call MPI_ALLREDUCE(MPI_IN_PLACE,xd,1,MPI_TPREC,MPI_SUM,
     $           COMM_TINKER,ierr)
              call MPI_ALLREDUCE(MPI_IN_PLACE,yd,1,MPI_TPREC,MPI_SUM,
     $           COMM_TINKER,ierr)
              call MPI_ALLREDUCE(MPI_IN_PLACE,zd,1,MPI_TPREC,MPI_SUM,
     $           COMM_TINKER,ierr)
              call MPI_ALLREDUCE(MPI_IN_PLACE,xq,1,MPI_TPREC,MPI_SUM,
     $           COMM_TINKER,ierr)
              call MPI_ALLREDUCE(MPI_IN_PLACE,yq,1,MPI_TPREC,MPI_SUM,
     $           COMM_TINKER,ierr)
              call MPI_ALLREDUCE(MPI_IN_PLACE,zq,1,MPI_TPREC,MPI_SUM,
     $           COMM_TINKER,ierr)
              if (rank.eq.0) then
                xv    = xq * (xu+xup)
                yv    = yq * (yu+yup)
                zv    = zq * (zu+zup)
                vterm = xv + yv + zv + xu*xup + yu*yup + zu*zup
     &                     + xd*(xu+xup) + yd*(yu+yup) + zd*(zu+zup)
                vterm    = term * vterm
                vir(1,1) = vir(1,1) + term*xv + vterm
                vir(2,1) = vir(2,1) + term*xv
                vir(3,1) = vir(3,1) + term*xv
                vir(1,2) = vir(1,2) + term*yv
                vir(2,2) = vir(2,2) + term*yv + vterm
                vir(3,2) = vir(3,2) + term*yv
                vir(1,3) = vir(1,3) + term*zv
                vir(2,3) = vir(2,3) + term*zv
                vir(3,3) = vir(3,3) + term*zv + vterm
              end if
c
           end if  !( boundary .eq. VACUUM )
c
#ifdef _OPENACC
           if (dir_queue.ne.rec_queue) call end_dir_stream_cover
#endif
         end if  !( use_pself )
         call timer_exit( timer_real,quiet_timers )
      end if  !( use_pmecore direct test )
c
c     compute the reciprocal space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.gt.ndir-1))
     &   then
         if (use_prec) then
            call timer_enter( timer_rec )
            call eprecip1gpu
            call timer_exit( timer_rec,quiet_timers )
         end if
      end if
c
#ifdef _OPENACC
      !End async compute overlapping if necessary
      if (dir_queue.ne.rec_queue) then
         call end_dir_stream_cover
      end if
#endif
c
c     Sum contribution of all energy
c
!$acc serial async(rec_queue) present(ep,epself,ep_r,eprec)
      ep = ep + enr2en( epself+ep_r ) + eprec
!$acc end serial
c
      end
c
c
c
      subroutine elambdapolar1cgpu
      use atmlst
      use atoms
      use boxes
      use chgpot
      use deriv
      use domdec
      use energi
      use ewald
      use epolar1gpu_inl
      use iounit
      use math
      use mpi
      use mpole
      use mutant
      use polar
      use polpot
      use potent
      use tinheader,only: zerom,zeromd
      use uprior
      use utilgpu
      use virial
      implicit none
      integer i,iipole,j,k,ierr,altopt
      integer(mipk) sizd8, sizr8
      real(r_p) elambdatemp,plambda,temp0
      mdyn_rtyp, allocatable :: delambdap0(:,:),delambdap1(:,:)
      real(r_p), allocatable :: delambdaprec0(:,:),delambdaprec1(:,:)
      real(r_p) elambdap0,elambdap1
      real(r_p) dplambdadelambdae,d2plambdad2elambdae
      parameter(
#ifdef _OPENACC
     &          altopt = 0
#else
     &          altopt = 1
#endif
     &         )
c
      if (npole .eq. 0)  return
      if (.not.(use_mpole)) then
!$acc serial async(rec_queue) present(delambdae)
        delambdae = 0.0
!$acc end serial
      end if
c
      allocate (delambdaprec0(3,nlocrec2))
      allocate (delambdaprec1(3,nlocrec2))
      allocate (delambdap0(3,nbloc))
      allocate (delambdap1(3,nbloc))
!$acc enter data create(delambdaprec0,delambdaprec1,delambdap0
!$acc&          ,delambdap1,elambdap0,elambdap1) async(rec_queue)
      elambdatemp = elambda  
      sizd8  = 3*nbloc
      sizr8  = 3*nlocrec2
c
c     polarization is interpolated between elambda=1 and elambda=0, for lambda.gt.plambda,
c     otherwise the value taken is for elambda=0
c
      if (elambda.gt.bplambda) then
         elambda = 1.0
         call altelec(altopt)
         call rotpolegpu
         call epolar1cgpu

!$acc serial async(rec_queue) present(elambdap1,ep)
         elambdap1  = ep
         ep         = 0
!$acc end serial
         call mem_move(delambdap1,dep,sizd8,rec_stream)
         call mem_move(delambdaprec1,deprec,sizr8,rec_stream)
      else
!$acc serial async(rec_queue) present(elambdap1,ep)
         elambdap1 = 0.0
         ep        = 0
!$acc end serial
         call mem_set(delambdap1,zeromd,sizd8,rec_stream)
         call mem_set(delambdaprec1,zerom,sizr8,rec_stream)
      end if

      elambda = 0.0
      call altelec(altopt)
      call rotpolegpu
      call mem_set(dep,zeromd,sizd8,rec_stream)
      call mem_set(deprec,zerom,sizr8,rec_stream)
      call epolar1cgpu
!$acc serial async(rec_queue) present(elambdap0,ep)
      elambdap0  = ep
!$acc end serial
      call mem_move(delambdap0,dep,sizd8,rec_stream)
      call mem_move(delambdaprec0,deprec,sizr8,rec_stream)
 
      elambda = elambdatemp 
c
c     interpolation of "plambda" between bplambda and 1 as a function of
c     elambda: 
c       plambda = 0 for elambda.le.bplambda
c       u = (elambda-bplambda)/(1-bplambda)
c       plambda = u**3 for elambda.gt.plambda
c       ep = (1-plambda)*ep0 +  plambda*ep1
c
      if (elambda.le.bplambda) then
           plambda          = 0.0
          dplambdadelambdae = 0.0
        d2plambdad2elambdae = 0.0
      else
           plambda          =     ((elambda-bplambda)/(1-bplambda))**3
          dplambdadelambdae = 3.0*((elambda-bplambda)/(1-bplambda))**2
        d2plambdad2elambdae = 6.0*((elambda-bplambda)/(1-bplambda))
      end if

!$acc serial async(rec_queue) present(elambdap0,elambdap1,ep
!$acc&      ,delambdae)
      ep        =  plambda*elambdap1 + (1-plambda)*elambdap0
      delambdae = delambdae + (elambdap1-elambdap0)*dplambdadelambdae
!$acc end serial
!$acc parallel loop async(rec_queue) collapse(2) default(present)
      do i = 1,nlocrec2; do j = 1,3
      deprec(j,i) = (1-plambda)*delambdaprec0(j,i)
     &            +    plambda *delambdaprec1(j,i)
      end do; end do
!$acc parallel loop async(rec_queue) collapse(2) default(present)
      do i = 1,nbloc; do j = 1,3
      temp0    = (1-plambda)*mdr2md(delambdap0(j,i))
     &         +    plambda *mdr2md(delambdap1(j,i))
      dep(j,i) = rp2mdr(temp0)      
      end do; end do
!$acc update host(delambdae) async(rec_queue)
c
c     reset lambda to initial value
c
      call altelec(altopt)
      call rotpolegpu
!$acc exit data delete(delambdaprec0,delambdaprec1,delambdap0
!$acc&         ,delambdap1,elambdap0,elambdap1) async(rec_queue)
      end
c
c
c
c     #################################################################
c     ##                                                             ##
c     ## subroutine epreal1cpgu -- Ewald real space derivs via list  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "epreal1d" evaluates the real space portion of the Ewald
c     summation energy and gradient due to dipole polarization
c     via a neighbor list
c
c
      subroutine epreal1cgpu
      use atmlst
      use atoms
      use bound
      use chgpot
      use couple
      use deriv
      use domdec
      use energi
      use elec_wspace,only: fix=>r2Work1,fiy=>r2Work2,fiz=>r2Work3
     &               ,trqvec=>r2Work4
      !use erf_mod
      use epolar1gpu_inl
      use ewald
      use inform     ,only: deb_Path
      use inter
      use interfaces ,only: epreal1c_core_p,torquegpu
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
      use virial
      use mpi
      use tinMemory
      use utils
      use utilgpu
      use timestat

      implicit none
      integer i,iglob,j,k
      integer ii,iipole,ivec
      integer kpole,kglob,kbis
      integer iax,iay,iaz
      real(t_p) alsq2,alsq2n,f
      real(r_p),save:: vxx,vxy,vxz,vyy,vyz,vzz
      real(t_p) xix,xiy,xiz
      real(t_p) yix,yiy,yiz
      real(t_p) zix,ziy,ziz
      logical*1,parameter::extract=.false.
      logical,save:: f_in=.true.
      character*10 mode

      if(deb_Path)write(*,'(2x,a)') 'epreal1cgpu'

      call timer_enter( timer_epreal )
#ifdef _OPENACC
      def_queue  = dir_queue
      def_stream = dir_stream
#endif

      if (f_in) then
         f_in=.false.
c
c     zero out temporary accumulation virial components
c
         vxx = 0.0_re_p
         vxy = 0.0_re_p
         vxz = 0.0_re_p
         vyy = 0.0_re_p
         vyz = 0.0_re_p
         vzz = 0.0_re_p
!$acc enter data copyin(vxx,vxy,vxz,vyy,vyz,vzz)
      end if

      call prmem_request(fix,3,npolelocnl)
      call prmem_request(fiy,3,npolelocnl)
      call prmem_request(fiz,3,npolelocnl)
      call prmem_request(trqvec,3,npolelocnl)
c
c     set arrays to store fields
c
      call set_to_zero1(trqvec,3*npolelocnl,def_queue)
c
c     set conversion factor, cutoff and switching coefficients
c
      if (use_polarshortreal) then
         mode = 'SHORTEWALD'
      else
         mode = 'EWALD'
      end if
      call switch (mode)

      f      = 0.5_ti_p * electric / dielec
      alsq2  = 2.0_ti_p * aewald**2
      alsq2n = 0.0_ti_p
      if (aewald .gt. 0.0_ti_p) alsq2n = 1.0_ti_p / (sqrtpi * aewald)

      call epreal1c_core_p(trqvec,vxx,vxy,vxz,vyy,vyz,vzz)
c
c     torque is induced field and gradient cross permanent moments
c
      call torquegpu(trqvec,fix,fiy,fiz,dep,extract)

      if (use_virial) then

!$acc parallel loop async(def_queue) default(present)
!$acc&         present(vxx,vxy,vxz,vyy,vyz,vzz)
      do k = 1, npolelocnl
         iipole = poleglobnl(k)
         iglob  = ipole(iipole)
         iax    = xaxis(iipole)
         iay    = yaxis(iipole)
         iaz    = zaxis(iipole)
         if (iax.gt.0)  then
            xix = x(iax) - x(iglob) ! xix
            yix = y(iax) - y(iglob) ! yix
            zix = z(iax) - z(iglob) ! zix
         else
            xix = 0.0_ti_p ! xix
            yix = 0.0_ti_p ! yix
            zix = 0.0_ti_p ! zix
         endif
         if (iay.gt.0)  then
            xiy = x(iay) - x(iglob) ! xiy
            yiy = y(iay) - y(iglob) ! yiy
            ziy = z(iay) - z(iglob) ! ziy
         else
            xiy = 0.0_ti_p ! xiy
            yiy = 0.0_ti_p ! yiy
            ziy = 0.0_ti_p ! ziy
         endif
         if (iaz.gt.0) then
            xiz = x(iaz) - x(iglob) ! xiz
            yiz = y(iaz) - y(iglob) ! yiz
            ziz = z(iaz) - z(iglob) ! ziz
         else
            xiz = 0.0_ti_p ! xiz
            yiz = 0.0_ti_p ! yiz
            ziz = 0.0_ti_p ! ziz
         endif

         vxx    = vxx +  xix*fix(1,k)  +  xiy*fiy(1,k)  + xiz*fiz(1,k)
         vxy    = vxy +  yix*fix(1,k)  +  yiy*fiy(1,k)  + yiz*fiz(1,k)
         vxz    = vxz +  zix*fix(1,k)  +  ziy*fiy(1,k)  + ziz*fiz(1,k)
         vyy    = vyy +  yix*fix(2,k)  +  yiy*fiy(2,k)  + yiz*fiz(2,k)
         vyz    = vyz +  zix*fix(2,k)  +  ziy*fiy(2,k)  + ziz*fiz(2,k)
         vzz    = vzz +  zix*fix(3,k)  +  ziy*fiy(3,k)  + ziz*fiz(3,k)
      end do
c
c     increment the virial due to pairwise Cartesian forces
c
      if (use_polarshortreal) then
!$acc serial present(virsave) async(def_queue)
!$acc&       present(vxx,vxy,vxz,vyy,vyz,vzz)
!$acc&       present(g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
         g_vxx  = g_vxx + vxx
         g_vxy  = g_vxy + vxy
         g_vxz  = g_vxz + vxz
         g_vyy  = g_vyy + vyy
         g_vyz  = g_vyz + vyz
         g_vzz  = g_vzz + vzz
         virsave(1,1)  = virsave(1,1) + vxx
         virsave(2,1)  = virsave(2,1) + vxy
         virsave(3,1)  = virsave(3,1) + vxz
         virsave(1,2)  = virsave(1,2) + vxy
         virsave(2,2)  = virsave(2,2) + vyy
         virsave(3,2)  = virsave(3,2) + vyz
         virsave(1,3)  = virsave(1,3) + vxz
         virsave(2,3)  = virsave(2,3) + vyz
         virsave(3,3)  = virsave(3,3) + vzz
         vxx = 0.0_re_p
         vxy = 0.0_re_p
         vxz = 0.0_re_p
         vyy = 0.0_re_p
         vyz = 0.0_re_p
         vzz = 0.0_re_p
!$acc end serial
      else
!$acc serial async(def_queue)
!$acc&       present(vxx,vxy,vxz,vyy,vyz,vzz)
!$acc&       present(g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
         g_vxx  = g_vxx + vxx
         g_vxy  = g_vxy + vxy
         g_vxz  = g_vxz + vxz
         g_vyy  = g_vyy + vyy
         g_vyz  = g_vyz + vyz
         g_vzz  = g_vzz + vzz
         vxx = 0.0_re_p
         vxy = 0.0_re_p
         vxz = 0.0_re_p
         vyy = 0.0_re_p
         vyz = 0.0_re_p
         vzz = 0.0_re_p
!$acc end serial
      end if !( use_polarshortreal )

      end if !( use_virial )
c
      call timer_exit( timer_epreal )
      end

      subroutine epreal1c_core1(trqvec,vxx,vxy,vxz,vyy,vyz,vzz)
      use atmlst
      use atoms
      use bound
      use chgpot
      use couple
      use deriv
      use domdec
      use energi
      !use erf_mod
      use epolar1gpu_inl
      use ewald
      use inter
      use inform ,only: deb_Path
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
      use virial
      use mpi
      use utils
      use utilgpu
      use timestat

      implicit none

      real(t_p) trqvec(3,npolelocnl)
      real(r_p) vxx,vxy,vxz,vyy,vyz,vzz

      integer i,iglob,j,k,kk,ksp,kd,iploc,kploc
      integer nnelst,nnelst1
      integer ii,iii,iipole,ivec
      integer kpole,kglob,kbis
      integer nn12,nn13,nn14,ntot,nn4,nn15
      integer nnp11,nnp12,nnp13,nnp14
      real(t_p) f,e,r2
      real(t_p) one,two,half
      real(t_p) alsq2,alsq2n
      real(t_p) xi,yi,zi
      real(t_p) pdi,pti
      real(t_p) pgamma,damp
      real(t_p) pscale,dscale,uscale

      integer   ipscal(maxscaling),idscal(92)
      real(t_p) fpscal(maxscaling),fdscal(92),fuscal(92)
      type(rpole_elt):: ip,kp
      type(real6) :: dpui,dpuk
      type(real3) :: pos,posi,ufli,uflk,trqi,trqk,frc
      type(mdyn3_r):: frc_r
      parameter(half=0.5_ti_p)
      parameter(one=1.0_ti_p, two=2.0_ti_p)
      real(t_p),save:: uscalevec(5),dscalevec(5),pscalevec(5)
      logical,save::f_in=.true.

      if(deb_Path) write(*,'(2x,a)') 'epreal1c_core1'
c
c     set conversion factor, cutoff and switching coefficients
c
      if (use_polarshortreal) then
         print*, 'epreal1c_core1 is unefficient for short interactions'
      end if

      f      = 0.5_ti_p * electric / dielec
      alsq2  = 2.0_ti_p * aewald**2
      alsq2n = 0.0_ti_p
      if (aewald .gt. 0.0_ti_p) alsq2n = 1.0_ti_p / (sqrtpi * aewald)

      if (f_in) then
         uscalevec = [u1scale,u2scale,u3scale,u4scale,0.0_ti_p]
         dscalevec = [d1scale,d2scale,d3scale,d4scale,0.0_ti_p]
         pscalevec = [0.0_ti_p,p2scale,p3scale,p4scale,p5scale]
!$acc enter data copyin(uscalevec,dscalevec,pscalevec)
         f_in = .false.
      end if

!$acc parallel loop gang vector_length(32)
!$acc&         private(ksp,kd,ip,dpui,ipscal,idscal,fpscal,fdscal,
!$acc&   fuscal,nn4,nn15,nnp11)
!$acc&         present(trqvec,vxx,vxy,vxz,vyy,vyz,vzz,ep)
!$acc&         present(poleglobnl,ipole,rpole,thole,pdamp,loc,x,y,z,
!$acc&   uind,uinp,allscal_n,allscal_p,numscal_n,numscal_p,
!$acc&   scalbeg_n,scalbeg_p,typscal_n,typscal_p,
!$acc&   elst,nelstc,dep,vir,polelocnl,
!$acc&   uscalevec,pscalevec,dscalevec)
!$acc&         async(def_queue)
      MAINLOOP:
     &do ii = 1, npolelocnl
         iipole  = poleglobnl(ii)
         iglob   = ipole (iipole)
         i       = loc(iglob)
         if (i.eq.0.or.i.gt.nbloc) then
            print*, warning
            cycle MAINLOOP
         endif
         nnelst  = nelstc (ii)
c
c        No neighbours
c
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
         ip%qxy  = rpole( 6, iipole)
         ip%qxz  = rpole( 7, iipole)
         ip%qyy  = rpole( 9, iipole)
         ip%qyz  = rpole(10, iipole)
         ip%qzz  = rpole(13, iipole)

         dpui%x  = uind ( 1, iipole)
         dpui%y  = uind ( 2, iipole)
         dpui%z  = uind ( 3, iipole)
         dpui%xx = uinp ( 1, iipole)
         dpui%yy = uinp ( 2, iipole)
         dpui%zz = uinp ( 3, iipole)
         ksp     = 0
         kd      = 0
c
c       fill scaling factor along kglob and interaction type
c
         nn4   = 0
         nn15  = 0
         nnp11 = 0
         ntot  = numscal_n(iglob)
         nn12  = scalbeg_n(iglob)
         nnp14 = numscal_p(iglob)
         nnp12 = scalbeg_p(iglob)
         if (nnp14>92) print*,"epreal_core1: not enough space for",
     &      "scaling factor"
!$acc loop vector
         do j = 1,ntot
            ipscal(j) = allscal_n(nn12+j)
            fpscal(j) = pscalevec(typscal_n(nn12+j))
            if (typscal_n(nn12+j).eq.4) then
!$acc atomic
               nn4 = nn4 +1
            end if
            if (typscal_n(nn12+j).eq.5) then
!$acc atomic
               nn15 = nn15 +1
            end if
         end do

!$acc loop vector
         do j = 1,nnp14
            idscal(j) = allscal_p(nnp12+j)
            fdscal(j) = dscalevec(typscal_p(nnp12+j))
            fuscal(j) = uscalevec(typscal_p(nnp12+j))
            if (typscal_p(nnp12+j).eq.1) then
!$acc atomic
               nnp11 = nnp11 +1
            end if
         end do
         nn13 = ntot - nn4 - nn15
         nn14 = ntot - nn15
c
c     loop on the neighbors
c
!$acc loop vector private(frc,frc_r,trqi,trqk,dpuk,kp,pos)
         do k = 1, nnelst
            kpole  = elst(k,ii)
            kglob  = ipole(kpole)
            kbis   = loc  (kglob)
            kploc  = polelocnl(kpole)
            pos%x  = x(kglob) - posi%x
            pos%y  = y(kglob) - posi%y
            pos%z  = z(kglob) - posi%z

            call image_inl(pos%x,pos%y,pos%z)
            ! cutoff
            r2     = pos%x**2 + pos%y**2 + pos%z**2
            if (r2>off2 .or. kbis.eq.0) cycle
c
c      set exclusion coefficients for connected atoms
c
            pscale = 1.0_ti_p
            uscale = 1.0_ti_p
            dscale = 1.0_ti_p
            if (ksp<ntot) then
!$acc          loop seq
               do j=1,ntot
                  if (ipscal(j)==kglob) then
                     pscale  = fpscal(j)
                     ! deal with 4-1 interaction
                     if (nn13.lt.j.and.j.le.nn14) then
                        do kk=1,nnp11
                           if (idscal(kk).eq.kglob) then
                              pscale = pscale*p41scale
                              goto 20
                           end if
                        end do
                     end if
  20                 continue
!$acc                atomic update
                     ksp = ksp+1
                     goto 30
                  end if
               end do
            end if
  30        continue
            if (kd<nnp14) then
!$acc          loop seq
               do j=1,nnp14
                  if (idscal(j)==kglob) then
                     uscale  = fuscal(j)
                     dscale  = fdscal(j)
!$acc                atomic update
                     kd = kd+1
                     exit
                  end if
               end do
            end if
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
              frc%x=0.0;  frc%y=0.0;  frc%z=0.0;
            frc_r%x=0.0;frc_r%y=0.0;frc_r%z=0.0;
             trqi%x=0.0; trqi%y=0.0; trqi%z=0.0;
             trqk%x=0.0; trqk%y=0.0; trqk%z=0.0;
            call epolar1_couple(dpui,ip,dpuk,kp,r2,pos,
     &                  aewald,alsq2,alsq2n,pgamma,damp,f,
     &                  dscale,pscale,uscale,
     &                  e,frc,frc_r,trqi,trqk,.false.)
c
c     increment energy
c
            ep       = ep + e
c
c     increment gradient and virial due to Cartesian forces
c
!$acc atomic
            dep(1,i) = dep(1,i) - frc_r%x
!$acc atomic
            dep(2,i) = dep(2,i) - frc_r%y
!$acc atomic
            dep(3,i) = dep(3,i) - frc_r%z
!$acc atomic
            dep(1,kbis) = dep(1,kbis) + frc_r%x
!$acc atomic
            dep(2,kbis) = dep(2,kbis) + frc_r%y
!$acc atomic
            dep(3,kbis) = dep(3,kbis) + frc_r%z

            vxx     = vxx + pos%x * frc%x
            vxy     = vxy + pos%y * frc%x
            vxz     = vxz + pos%z * frc%x
            vyy     = vyy + pos%y * frc%y
            vyz     = vyz + pos%z * frc%y
            vzz     = vzz + pos%z * frc%z
c
c     increment torque
c
!$acc atomic
            trqvec(1,ii) = trqvec(1,ii) + trqi%x
!$acc atomic
            trqvec(2,ii) = trqvec(2,ii) + trqi%y
!$acc atomic
            trqvec(3,ii) = trqvec(3,ii) + trqi%z
!$acc atomic
            trqvec(1,kploc) = trqvec(1,kploc) + trqk%x
!$acc atomic
            trqvec(2,kploc) = trqvec(2,kploc) + trqk%y
!$acc atomic
            trqvec(3,kploc) = trqvec(3,kploc) + trqk%z
         enddo

      end do  MAINLOOP
c
      end

      subroutine epreal1c_core2(trqvec,vxx,vxy,vxz,vyy,vyz,vzz)
      use atmlst  ,only: poleglobnl
      use atoms   ,only: x,y,z
      use chgpot  ,only: dielec,electric
      use deriv   ,only: dep
      use domdec  ,only: rank,loc
      use energi  ,only: ep=>ep_r
      use epolar1gpu_inl
      use ewald   ,only: aewald
      use inform  ,only: deb_Path
      use interfaces ,only: epreal1c_correct_scale
      use math    ,only: sqrtpi
      use mpole   ,only: npolelocnl,ipole,rpole,polelocnl
      use neigh   ,only: elst,nelst,nelstc,shortelst,nshortelst
      use polar   ,only: uind,uinp,thole,pdamp
      use potent  ,only: use_polarshortreal
      use shunt   ,only: off2
      use tinheader  ,only: ti_p
      use virial
      use tinTypes,only: real3,real6,rpole_elt,mdyn3_r
      use utilgpu ,only: def_queue

      implicit none

      real(t_p),intent(inout):: trqvec(3,npolelocnl)
      real(r_p),intent(inout):: vxx,vxy,vxz,vyy,vyz,vzz

      integer i,iglob,j,k,iploc,kploc
      integer nnelst
      integer ii,iipole
      integer kpole,kglob,kbis
      integer,pointer,save :: lst(:,:),nlst(:)
      real(t_p) alsq2,alsq2n
      real(t_p) r2,pgamma,damp
      real(t_p) f,e
      real(t_p) pdi,pti
      real(t_p) one
      real(t_p) pscale,dscale,uscale

      type(rpole_elt):: ip,kp
      type(real6) :: dpui,dpuk
      type(real3) :: pos,posi,ufli,uflk,trqi,trqk,frc
      type(mdyn3_r):: frc_r
      parameter(one=1.0_ti_p)

      if(deb_Path)write(*,'(2x,a)') 'epreal1c_core2'
c
c     set conversion factor, cutoff and switching coefficients
c
      dscale = 1.0_ti_p
      pscale = 1.0_ti_p
      uscale = 1.0_ti_p
      f      = 0.5_ti_p * electric / dielec
      alsq2  = 2.0_ti_p * aewald**2
      alsq2n = 0.0_ti_p
      if (aewald .gt. 0.0_ti_p) alsq2n = 1.0_ti_p / (sqrtpi * aewald)

      if (use_polarshortreal) then
          lst =>  shortelst
         nlst => nshortelst
      else
          lst =>  elst
         nlst => nelstc
      end if

!$acc parallel loop gang vector_length(32)
!$acc&         present(trqvec,vxx,vxy,vxz,vyy,vyz,vzz,ep)
!$acc&         present(poleglobnl,ipole,rpole,thole,pdamp,
!$acc&  loc,x,y,z,uind,uinp,lst,nlst,dep,vir,polelocnl)
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
!$acc loop vector private(frc,frc_r,trqi,trqk,dpuk,kp,pos)
         do k = 1, nnelst
            kpole    = lst(k,ii)
            kglob    = ipole(kpole)
            kbis     = loc  (kglob)
            kploc    = polelocnl(kpole)
            pos%x    = x(kglob) - posi%x
            pos%y    = y(kglob) - posi%y
            pos%z    = z(kglob) - posi%z

            call image_inl(pos%x,pos%y,pos%z)
            ! cutoff
            r2       = pos%x**2 + pos%y**2 + pos%z**2
            if (r2>off2 .or. kbis.eq.0) cycle
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
              frc%x=0.0;  frc%y=0.0;  frc%z=0.0;
            frc_r%x=0.0;frc_r%y=0.0;frc_r%z=0.0;
             trqi%x=0.0; trqi%y=0.0; trqi%z=0.0;
             trqk%x=0.0; trqk%y=0.0; trqk%z=0.0;
            call epolar1_couple(dpui,ip,dpuk,kp,r2,pos,
     &                  aewald,alsq2,alsq2n,pgamma,damp,f,
     &                  1.0_ti_p,1.0_ti_p,1.0_ti_p,
     &                  e,frc,frc_r,trqi,trqk,.false.)
c
c     increment energy
c
            ep          = ep + tp2enr(e)
c
c     increment gradient and virial due to Cartesian forces
c
!$acc atomic
            dep(1,i)    = dep(1,i)    - frc_r%x
!$acc atomic
            dep(2,i)    = dep(2,i)    - frc_r%y
!$acc atomic
            dep(3,i)    = dep(3,i)    - frc_r%z
!$acc atomic
            dep(1,kbis) = dep(1,kbis) + frc_r%x
!$acc atomic
            dep(2,kbis) = dep(2,kbis) + frc_r%y
!$acc atomic
            dep(3,kbis) = dep(3,kbis) + frc_r%z

            vxx         = vxx + pos%x * frc%x
            vxy         = vxy + pos%y * frc%x
            vxz         = vxz + pos%z * frc%x
            vyy         = vyy + pos%y * frc%y
            vyz         = vyz + pos%z * frc%y
            vzz         = vzz + pos%z * frc%z
c
c     increment torque
c
!$acc atomic
            trqvec(1,ii)    = trqvec(1,ii)    + trqi%x
!$acc atomic
            trqvec(2,ii)    = trqvec(2,ii)    + trqi%y
!$acc atomic
            trqvec(3,ii)    = trqvec(3,ii)    + trqi%z
!$acc atomic
            trqvec(1,kploc) = trqvec(1,kploc) + trqk%x
!$acc atomic
            trqvec(2,kploc) = trqvec(2,kploc) + trqk%y
!$acc atomic
            trqvec(3,kploc) = trqvec(3,kploc) + trqk%z
         enddo

      end do  MAINLOOP
c
      call epreal1c_correct_scale(trqvec,vxx,vxy,vxz,vyy,vyz,vzz)

      end

      subroutine epreal1c_core3(trqvec,vxx,vxy,vxz,vyy,vyz,vzz)
      use atmlst  ,only: poleglobnl
      use atoms   ,only: x,y,z,n
      use chgpot  ,only: dielec,electric
      use deriv   ,only: dep
      use domdec  ,only: xbegproc,ybegproc,zbegproc
     &            ,nproc,rank,xendproc,yendproc,zendproc
     &            ,nbloc,loc
      use energi  ,only: ep=>ep_r
#ifdef _CUDA
      use epolar1cu ,only: epreal1c_core_cu
#endif
      use ewald   ,only: aewald
      use inform  ,only: deb_Path
      use interfaces ,only: epreal1c_correct_scale
      use math    ,only: sqrtpi
      use mpole   ,only: npolelocnl,ipole,rpole,polelocnl
     &            ,ipole,npolelocnlb,npolebloc
     &            ,npolelocnlb_pair,npolelocnlb2_pair
     &            ,nspnlb2=>nshortpolelocnlb2_pair
      use neigh   ,only:ipole_s=>celle_glob,pglob_s=>celle_pole
     &            , loc_s=>celle_loc, plocnl_s=>celle_plocnl
     &            , ieblst_s=>ieblst, eblst_s=>eblst
     &            , iseblst_s=>ishorteblst, seblst_s=>shorteblst
     &            , x_s=>celle_x, y_s=>celle_y, z_s=>celle_z
      use polar   ,only: uind,uinp,thole,pdamp
      use potent  ,only: use_polarshortreal
      use polpot  ,only: n_dpuscale,dpucorrect_ik,dpucorrect_scale
      use shunt   ,only: off2
      use tinheader,only: ti_p
      use virial
#ifdef _CUDA
      use cudafor
      use utilcu  ,only: BLOCK_DIM,check_launch_kernel
#endif
      use utilgpu ,only: def_queue
     &            ,real3,real6,real3_red,rpole_elt
     &            ,ered_buff=>ered_buf1,vred_buff
     &            ,reduce_energy_virial
     &            ,RED_BUFF_SIZE,zero_evir_red_buffer
     &            ,maxBlock,BLOCK_SIZE
#ifdef  _OPENACC
     &            ,dir_stream,def_stream,rec_stream,nSMP
#endif

      implicit none

      real(t_p),intent(inout):: trqvec(3,npolelocnl)
      real(r_p),intent(inout):: vxx,vxy,vxz,vyy,vyz,vzz

      integer i
      integer start,start1,sized
      integer,save::ndec
      integer lst_beg,gS1
      integer,save:: gS=0
      real(t_p) alsq2,alsq2n,f
      real(t_p) p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
      logical,save::first_in=.true.
      logical,parameter::dyn_gS=.false.

      if(deb_Path)write(*,'(2x,a)') 'epreal1c_core3'
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

#ifdef _CUDA
      if (first_in) then
         ! Compute though occupancy the right gridSize to launch the kernel with
         call cudaMaxGridSize("epreal1c_core_cu",gS)
         if (deb_Path) print*, 'epreal1c blockSize ',gS
         ndec=1
         if (def_stream.ne.rec_stream) ndec=4
         first_in = .false.
      end if
      !print*,'pair block ',npolelocnlb2_pair

      call zero_evir_red_buffer(def_queue)

!$acc host_data use_device(ipole_s,pglob_s,loc_s,plocnl_s,ieblst_s,
!$acc&    iseblst_s,eblst_s,seblst_s,x_s,y_s,z_s,rpole,pdamp,thole,
!$acc&    loc,x,y,z,ipole,polelocnl,dpucorrect_ik,dpucorrect_scale,
!$acc&    uind,uinp,dep,trqvec,ered_buff,vred_buff)

      if (use_polarshortreal) then

      if (dyn_gS) gS = min( nspnlb2/8,maxBlock )
      call epreal1c_core_cu<<<gS,BLOCK_DIM,0,def_stream>>>
     &     (ipole_s,pglob_s,loc_s,plocnl_s
     &     ,iseblst_s,seblst_s(lst_beg)
     &     ,x_s,y_s,z_s,rpole,pdamp,thole,uind,uinp
     &     ,dep,trqvec,ered_buff,vred_buff
     &     ,npolelocnlb,nspnlb2,npolebloc,n
     &     ,off2,f,alsq2,alsq2n,aewald
     &     ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend)

      else

      sized = npolelocnlb2_pair/ndec
      if (dyn_gS) gS = min( sized/8,maxBlock )

      ! Split electrostatic kernel to ease recovering process in MPI
      do i = 1,ndec
         start  = (i-1)*sized + 1
         start1 = lst_beg+(start-1)*BLOCK_SIZE
         if (i.eq.ndec) sized = npolelocnlb2_pair-start+1
         call epreal1c_core_cu<<<gS,BLOCK_DIM,0,def_stream>>>
     &        (ipole_s,pglob_s,loc_s,plocnl_s
     &        ,ieblst_s(start),eblst_s(start1)
     &        ,x_s,y_s,z_s,rpole,pdamp,thole,uind,uinp
     &        ,dep,trqvec,ered_buff,vred_buff
     &        ,npolelocnlb,sized,npolebloc,n
     &        ,off2,f,alsq2,alsq2n,aewald
     &        ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend)
         call check_launch_kernel(" epreal1c_core_cu")
      end do

      end if

      gS1 = n_dpuscale/(2*BLOCK_DIM)
      call epreal1c_scaling_cu<<<gS1,BLOCK_DIM,0,def_stream>>>
     &     ( dpucorrect_ik,dpucorrect_scale,ipole,loc
     &     , polelocnl,x,y,z,pdamp,thole,rpole,uind,uinp
     &     , dep,trqvec,ered_buff,vred_buff
     &     , n,n_dpuscale,nbloc
     &     , off2,aewald,alsq2,alsq2n,f)
      call check_launch_kernel(" epreal1c_scaling_cu")

!$acc end host_data

      call reduce_energy_virial(ep,vxx,vxy,vxz,vyy,vyz,vzz
     &                         ,ered_buff,def_queue)
c     call epreal1c_correct_scale(trqvec,vxx,vxy,vxz,vyy,vyz,vzz)
#else
      print 100
 100  format('epreal1c_core3 is a specific device routine',/,
     &       'you are not supposed to get inside with your compile',
     &       'type.')
      __TINKER_FATAL__
#endif

      end

#ifdef _CUDA
      subroutine mpreal1c_core(trqvec,vxx,vxy,vxz,vyy,vyz,vzz)
      use atmlst  ,only: poleglobnl
      use atoms   ,only: x,y,z,n
      use chgpot  ,only: dielec,electric
      use cutoff  ,only: mpoleshortcut,shortheal
      use deriv   ,only: dep,dem
      use domdec  ,only: xbegproc,ybegproc,zbegproc
     &            ,nproc,rank,xendproc,yendproc,zendproc
     &            ,nbloc,loc
      use energi  ,only: ep=>ep_r
      use empole1cu ,only: emreal_scaling_cu
      use epolar1cu ,only: mpreal1c_core_cu
      use ewald   ,only: aewald
      use inform  ,only: deb_Path
      use interfaces ,only: epreal1c_correct_scale
     &               , m_normal,m_short,m_long
      use math    ,only: sqrtpi
      use mpole   ,only: npolelocnl,ipole,rpole,polelocnl
     &            ,npolelocnlb,npolebloc
     &            ,npolelocnlb_pair,npolelocnlb2_pair
     &            ,nspnlb2=>nshortpolelocnlb2_pair
      use mplpot  ,only:n_mscale,mcorrect_ik,mcorrect_scale
      use neigh   ,only:ipole_s=>celle_glob,pglob_s=>celle_pole
     &            ,loc_s=>celle_loc, plocnl_s=>celle_plocnl
     &            ,ieblst_s=>ieblst, eblst_s=>eblst
     &            ,iseblst_s=>ishorteblst, seblst_s=>shorteblst
     &            ,x_s=>celle_x, y_s=>celle_y, z_s=>celle_z
      use polar   ,only: uind,uinp,thole,pdamp
      use polpot  ,only: n_dpuscale,dpucorrect_ik,dpucorrect_scale
      use potent  ,only: use_polarshortreal,use_mpoleshortreal
     &            ,use_mpolelong
      use shunt   ,only: off,off2
      use tinheader  ,only: ti_p
      use virial
      use cudafor
      use utilcu  ,only: BLOCK_DIM,check_launch_kernel
      use utilgpu ,only: def_queue
     &            ,real3,real6,real3_red,rpole_elt
     &            ,ered_buff=>ered_buf1,vred_buff
     &            ,reduce_energy_virial
     &            ,RED_BUFF_SIZE,zero_evir_red_buffer
     &            ,maxBlock,BLOCK_SIZE
     &            ,dir_stream,def_stream,rec_stream,nSMP
      implicit none

      real(t_p),intent(inout):: trqvec(3,npolelocnl)
      real(r_p),intent(inout):: vxx,vxy,vxz,vyy,vyz,vzz

      integer i
      integer start,start1,sized
      integer,save::ndec
      integer lst_beg,gS1
      integer,save:: gS=0
      real(t_p) alsq2,alsq22,alsq24,alsq2n,f,fem,r_cut,sh_cut2
      real(t_p) p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
      logical,save::first_in=.true.
      logical,parameter::dyn_gS=.false.

      if(deb_Path)write(*,'(2x,a)') 'mpreal1c_core'
c
c     set conversion factor, cutoff and switching coefficients
c
      f      = 0.5_ti_p * electric / dielec
      fem    = 1.0_ti_p * electric / dielec
      alsq2  = 2.0_ti_p * aewald**2
      alsq22 = alsq2**2
      alsq24 = alsq22**2
      alsq2n = 0.0_ti_p
      if (aewald .gt. 0.0_ti_p) alsq2n = 1.0_ti_p / (sqrtpi * aewald)

      p_xbeg = xbegproc(rank+1)
      p_xend = xendproc(rank+1)
      p_ybeg = ybegproc(rank+1)
      p_yend = yendproc(rank+1)
      p_zbeg = zbegproc(rank+1)
      p_zend = zendproc(rank+1)
      lst_beg= 2*npolelocnlb_pair+1
      sh_cut2= 0.0_ti_p
      r_cut  = 0.0_ti_p

      if (use_mpoleshortreal) then
         sh_cut2 = 0.0_ti_p
         r_cut   = off
      else if (use_mpolelong) then
         sh_cut2 = (mpoleshortcut-shortheal)**2
         r_cut   = mpoleshortcut
      end if

      if (first_in) then
         ! Compute though occupancy the right gridSize to launch the kernel with
         call cudaMaxGridSize("mpreal1c_core_cu",gS)
         if (deb_Path) print*, ' mpreal1c blockSize ',gS
         ndec=1
         if (def_stream.ne.rec_stream) ndec=4
         first_in = .false.
      end if
      !print*,'pair block ',npolelocnlb2_pair

      call zero_evir_red_buffer(def_queue)

!$acc host_data use_device(ipole_s,pglob_s,loc_s,plocnl_s,ieblst_s,
!$acc&    iseblst_s,eblst_s,seblst_s,x_s,y_s,z_s,rpole,pdamp,thole,
!$acc&    mcorrect_ik,mcorrect_scale,dpucorrect_ik,dpucorrect_scale,
!$acc&    polelocnl,loc,ipole,x,y,z,
!$acc&    uind,uinp,dep,trqvec,ered_buff,vred_buff)

      if (use_polarshortreal.or.use_mpoleshortreal) then

      if (dyn_gS) gS = min(max(nspnlb2,8)/8,maxBlock)
      call mpreal1c_core_cu<<<gS,BLOCK_DIM,0,def_stream>>>
     &     (ipole_s,pglob_s,loc_s,plocnl_s
     &     ,iseblst_s,seblst_s(lst_beg)
     &     ,x_s,y_s,z_s,rpole,pdamp,thole,uind,uinp
     &     ,dep,trqvec,ered_buff,vred_buff
     &     ,npolelocnlb,nspnlb2,npolebloc,n,m_short
     &     ,off2,f,alsq2,alsq2n,aewald,r_cut,sh_cut2,shortheal
     &     ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend)

      else if (use_mpolelong) then

         write(0,*) 'ERROR mpreal1c_core !!'
         write(0,*) '      cannot compute long range interactions'
         __TINKER_FATAL__

      if (dyn_gS) gS = max(npolelocnlb2_pair,12)/12
      call mpreal1c_core_cu<<<gS,BLOCK_DIM,0,def_stream>>>
     &     (ipole_s,pglob_s,loc_s,plocnl_s
     &     ,ieblst_s,eblst_s(lst_beg)
     &     ,x_s,y_s,z_s,rpole,pdamp,thole,uind,uinp
     &     ,dep,trqvec,ered_buff,vred_buff
     &     ,npolelocnlb,npolelocnlb2_pair,npolebloc,n,m_long
     &     ,off2,f,alsq2,alsq2n,aewald,r_cut,sh_cut2,shortheal
     &     ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend)

      else

      sized = npolelocnlb2_pair/ndec
      if (dyn_gS) gS = min(max(sized,8)/8,maxBlock)

      ! Split electrostatic kernel to ease recovering process in MPI
      do i = 1,ndec
         start  = (i-1)*sized + 1
         start1 = lst_beg+(start-1)*BLOCK_SIZE
         if (i.eq.ndec) sized = npolelocnlb2_pair-start+1
         call mpreal1c_core_cu<<<gS,BLOCK_DIM,0,def_stream>>>
     &        (ipole_s,pglob_s,loc_s,plocnl_s
     &        ,ieblst_s(start),eblst_s(start1)
     &        ,x_s,y_s,z_s,rpole,pdamp,thole,uind,uinp
     &        ,dep,trqvec,ered_buff,vred_buff
     &        ,npolelocnlb,sized,npolebloc,n,m_normal
     &        ,off2,f,alsq2,alsq2n,aewald,r_cut,sh_cut2,shortheal
     &        ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend)
         call check_launch_kernel(" mpreal1c_core_cu")
      end do

      end if

      if (n_mscale.gt.0) then
         gS1 = max(n_mscale/(BLOCK_DIM),1)
         if (use_mpoleshortreal) then
         call emrealShLg_scaling_cu<<<gS1,BLOCK_DIM,0,def_stream>>>
     &        ( mcorrect_ik,mcorrect_scale,ipole,loc,polelocnl
     &        , x,y,z,rpole
     &        , dep,trqvec,ered_buff,vred_buff
     &        , n,nbloc,n_mscale,size(mcorrect_ik,1),.true.
     &        , r_cut,sh_cut2,off2
     &        , shortheal,fem,aewald,alsq2,alsq2n )
         call check_launch_kernel(" emrealShLg_scaling_cu")
         else
         call emreal_scaling_cu<<<gS1,BLOCK_DIM,0,def_stream>>>
     &        ( mcorrect_ik,mcorrect_scale,ipole,loc,polelocnl
     &        , x,y,z,rpole
     &        , dep,trqvec,ered_buff,vred_buff
     &        , n,nbloc,n_mscale,size(mcorrect_ik,1),.true.
     &        , off2,fem,aewald,alsq2,alsq2n )
         call check_launch_kernel(" emreal_scaling_cu")
         end if
      end if

      if (n_dpuscale.gt.0) then
         gS1 = max(n_dpuscale/(BLOCK_DIM),1)
         call epreal1c_scaling_cu<<<gS1,BLOCK_DIM,0,def_stream>>>
     &        ( dpucorrect_ik,dpucorrect_scale,ipole,loc
     &        , polelocnl,x,y,z,pdamp,thole,rpole,uind,uinp
     &        , dep,trqvec,ered_buff,vred_buff
     &        , n,n_dpuscale,nbloc
     &        , off2,aewald,alsq2,alsq2n,f )
         call check_launch_kernel(" epreal1c_scaling_cu")
      end if

!$acc end host_data

      call reduce_energy_virial(ep,vxx,vxy,vxz,vyy,vyz,vzz
     &                         ,ered_buff,def_queue)

      end
#endif

      ! Loop on scale interaction for correction
      subroutine epreal1c_correct_scale
     &           (trqvec,vxx,vxy,vxz,vyy,vyz,vzz)
      use atmlst  ,only: poleglobnl
      use atoms   ,only: x,y,z
      use chgpot  ,only: dielec,electric
      use deriv   ,only: dep
      use domdec  ,only: rank,loc
      use energi  ,only: ep=>ep_r
      use epolar1gpu_inl
      use ewald   ,only: aewald
      use inform  ,only: deb_Path
      use math    ,only: sqrtpi
      use mpole   ,only: npolelocnl,ipole,rpole,polelocnl
      use neigh   ,only: elst,nelst
      use polar   ,only: uind,uinp,thole,pdamp
      use polpot  ,only: n_dpuscale,dpucorrect_ik,dpucorrect_scale
      use shunt   ,only: off2
      use tinheader ,only: ti_p
      use tinTypes,only: real3,real6,mdyn3_r,rpole_elt
      use utilgpu ,only: def_queue
      use virial
      implicit none

      real(t_p),intent(inout):: trqvec(:,:)
      real(r_p),intent(inout):: vxx,vxy,vxz,vyy,vyz,vzz

      integer i,k,iglob,kglob,iploc,kploc
      integer ii,iipole,kpole,j,kbis
      integer nnelst
      real(t_p) alsq2,alsq2n
      real(t_p) r2,pgamma,damp
      real(t_p) f,e
      real(t_p) pdi,pti
      real(t_p) one
      real(t_p) pscale,dscale,uscale
      type(rpole_elt):: ip,kp
      type(real6)    :: dpui,dpuk
      type(real3)    :: pos,posi,ufli,uflk,trqi,trqk,frc
      type(mdyn3_r)  :: frc_r

      parameter(one=1.0_ti_p)

      if(deb_Path)
     &   write(*,'(2x,a)') 'epreal1c_correct_scale'
c
c     set conversion factor, cutoff and switching coefficients
c
      f      = 0.5_ti_p * electric / dielec
      alsq2  = 2.0_ti_p * aewald**2
      alsq2n = 0.0_ti_p
      if (aewald .gt. 0.0_ti_p) alsq2n = 1.0_ti_p / (sqrtpi * aewald)

!$acc parallel loop present(trqvec,vxx,vxy,vxz,vyy,vyz,vzz)
!$acc&     present(poleglobnl,ipole,rpole,thole,pdamp,loc,x,y,z,
!$acc&  uind,uinp,elst,nelst,dep,ep,vir,polelocnl,dpucorrect_ik,
!$acc&  dpucorrect_scale)
!$acc&         private(pos,ip,kp,dpui,dpuk,trqk,trqi,frc_r,frc)
!$acc&         async(def_queue)
      do ii = 1, n_dpuscale
         iipole   = dpucorrect_ik(2*(ii-1)+1)
         kpole    = dpucorrect_ik(2*(ii-1)+2)

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

         dscale   = dpucorrect_scale(3*(ii-1)+1)
         pscale   = dpucorrect_scale(3*(ii-1)+2)
         uscale   = dpucorrect_scale(3*(ii-1)+3)

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
           frc%x=0.0;  frc%y=0.0;  frc%z=0.0;
         frc_r%x=0.0;frc_r%y=0.0;frc_r%z=0.0;
          trqi%x=0.0; trqi%y=0.0; trqi%z=0.0;
          trqk%x=0.0; trqk%y=0.0; trqk%z=0.0;
         call epolar1_couple(dpui,ip,dpuk,kp,r2,pos,
     &               aewald,alsq2,alsq2n,pgamma,damp,f,
     &               dscale,pscale,uscale,
     &               e,frc,frc_r,trqi,trqk,.true.)
c
c     increment energy
c
         ep       = ep + tp2enr(e)
c
c     increment gradient and virial due to Cartesian forces
c
!$acc atomic
         dep(1,i) = dep(1,i) - frc_r%x
!$acc atomic
         dep(2,i) = dep(2,i) - frc_r%y
!$acc atomic
         dep(3,i) = dep(3,i) - frc_r%z
!$acc atomic
         dep(1,k) = dep(1,k) + frc_r%x
!$acc atomic
         dep(2,k) = dep(2,k) + frc_r%y
!$acc atomic
         dep(3,k) = dep(3,k) + frc_r%z

         vxx      = vxx + pos%x*frc%x
         vxy      = vxy + pos%y*frc%x
         vxz      = vxz + pos%z*frc%x
         vyy      = vyy + pos%y*frc%y
         vyz      = vyz + pos%z*frc%y
         vzz      = vzz + pos%z*frc%z
c
c     increment torque
c
!$acc atomic
         trqvec(1,iploc) = trqvec(1,iploc) + trqi%x
!$acc atomic
         trqvec(2,iploc) = trqvec(2,iploc) + trqi%y
!$acc atomic
         trqvec(3,iploc) = trqvec(3,iploc) + trqi%z
!$acc atomic
         trqvec(1,kploc) = trqvec(1,kploc) + trqk%x
!$acc atomic
         trqvec(2,kploc) = trqvec(2,kploc) + trqk%y
!$acc atomic
         trqvec(3,kploc) = trqvec(3,kploc) + trqk%z
      end do
c
      end
c
c
c     ###################################################################
c     ##                                                               ##
c     ##  subroutine eprecip1  --  PME recip polarize energy & derivs  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     "eprecip1" evaluates the reciprocal space portion of the particle
c     mesh Ewald summation energy and gradient due to dipole polarization
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
      module eprecip1gpu_mod
      logical  :: f_in=.true.
      real(r_p):: er,vxx,vyy,vzz,vxy,vxz,vyz
      real(t_p):: a(3,3)
      real(t_p),allocatable,target:: qgrip_b(:)
      real(t_p),pointer::            qgrip(:,:,:,:)
      end module
      subroutine eprecip1gpu
      use atmlst
      use atoms
      use bound
      use boxes
      use chgpot
      use deriv
      use domdec
      use energi    ,only: eprec
      use eprecip1gpu_mod
      use ewald
      use fft
      use inform    ,only: deb_Path
      use interfaces,only: fphi_uind_site1_p
     &              ,torquegpu,grid_uind_site_p
     &              ,grid_mpole_site_p,epreal1c_cp
      use math
      use mpole
      use pme
      use pme1
      use polar
      use polar_temp
      use polpot
      use potent
      use timestat
      use utils
      use utilgpu
      use virial
      use mpi
      implicit none
      integer status(MPI_STATUS_SIZE),tag,ierr,siz
      integer i,j,k,m,ii,iipole,iipol,iglob,iloc
      integer j1,j2,j3
      integer k1,k2,k3
      integer m1,m2,m3
      integer ntot,nff
      integer nf1,nf2,nf3
      integer kstat2,ked2,istat2,ied2,jstat2,jed2
      integer nprocloc,commloc,rankloc,proc
      integer(mipk) siz8
      real(t_p) eterm,f
      real(t_p) r1,r2,r3
      real(t_p) h1,h2,h3
      real(t_p) f1,f2,f3
      real(t_p) volterm,denom
      real(t_p) hsq,expterm
      real(t_p) term,pterm
      real(t_p) vterm,struc2
      real(t_p) time0,time1
      real(t_p) fiy(3),fiz(3),fix(3)
      real(t_p) cphim(4),cphid(4)
      real(t_p) cphip(4)
      integer  ,dimension(nproc)::reqsend,reqrec,req2send,req2rec
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
c
      if (deb_Path) write(*,'(2x,a)') 'eprecip1gpu'
      call timer_enter( timer_eprecip )

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
c     perform dynamic allocation of some global arrays
c
c     call mallocMpiGrid
      call prmem_request(cmp,   10,npolerecloc,async=.true.)
      call prmem_request(fmp,   10,npolerecloc,async=.true.)
      call prmem_request(fphidp,20,npolerecloc,async=.true.)
      call prmem_request(trqrec,03,npolerecloc,async=.true.)

      call prmem_request(fuind, 3,npolerecloc,async=.true.)
      call prmem_request(fuinp, 3,npolerecloc,async=.true.)
      call prmem_request(fphid,10,npolerecloc,async=.true.)
      call prmem_request(fphip,10,npolerecloc,async=.true.)

      if (f_in) then
!$acc enter data create(er,a)
!$acc&      create(vxx,vxy,vxz,vyy,vyz,vzz)
         f_in=.false.
      end if
c
c
c
c     zero out the temporary virial accumulation variables
c
c     vxx = 0.0_ti_p
c     vxy = 0.0_ti_p
c     vxz = 0.0_ti_p
c     vyy = 0.0_ti_p
c     vyz = 0.0_ti_p
c     vzz = 0.0_ti_p
c
c     get the fractional to Cartesian transformation matrix
c
      !call frac_to_cartgpu
c
c     initialize variables required for the scalar summation
c
      f       = electric / dielec
      ntot    = nfft1 * nfft2 * nfft3
      pterm   = (pi/aewald)**2
      volterm = pi * volbox
      nff     = nfft1 * nfft2
      nf1     = (nfft1+1) / 2
      nf2     = (nfft2+1) / 2
      nf3     = (nfft3+1) / 2
c
c     remove scalar sum virial from prior multipole 3-D FFT
c
!$acc serial async(rec_queue)
!$acc&       present(vxx,vxy,vxz,vyy,vyz,vzz,eprec)
!$acc&       present(vmxx,vmxy,vmxz,vmyy,vmyz,vmzz)
      eprec   = 0.0_re_p
      vxx     = -vmxx
      vxy     = -vmxy
      vxz     = -vmxz
      vyy     = -vmyy
      vyz     = -vmyz
      vzz     = -vmzz
!$acc end serial

      call timer_enter( timer_grid1 )
!$acc parallel loop collapse(2) async(rec_queue) default(present)
      do i = 1, npolerecloc
         do j = 1, 10
            iipole   = polerecglob(i)
            cmp(j,i) = rpole_scale(j)*rpole(rpole_ind_extract(j),iipole)
         end do
      end do

      call cmp_to_fmp_sitegpu(cmp,fmp)
c
c     zero out the PME grid
c
      call set_to_zero1(qgrid2in_2d(1,1,1,1,1),
     &          2*n1mpimax*n2mpimax*n3mpimax*(nrec_send+1),rec_queue)
c
c     convert Cartesian induced dipoles to fractional coordinates
c
!$acc parallel loop async(rec_queue)
      do i = 1, 3
         a(1,i) = real(nfft1,t_p) * recip(i,1)
         a(2,i) = real(nfft2,t_p) * recip(i,2)
         a(3,i) = real(nfft3,t_p) * recip(i,3)
      end do
!$acc parallel loop collapse(2) async(rec_queue) default(present)
      do i = 1, npolerecloc
         do j = 1, 3
            iipole = polerecglob(i)
            !iglob  = ipole(iipole)
            fuind(j,i)  =  a(j,1)*uind(1,iipole)
     &                   + a(j,2)*uind(2,iipole)
     &                   + a(j,3)*uind(3,iipole)
            fuinp(j,i)  =  a(j,1)*uinp(1,iipole)
     &                   + a(j,2)*uinp(2,iipole)
     &                   + a(j,3)*uinp(3,iipole)
         end do
      end do

      call grid_uind_site_p(fuind,fuinp,qgrid2in_2d)
      call timer_exit( timer_grid1,quiet_timers )

      !Exchange Grid Between reciprocal process
      call commGridFront( qgrid2in_2d,r_comm )

#ifdef _OPENACC
      ! Recover MPI communication with real space computations
      if (dir_queue.ne.rec_queue) then
         call start_dir_stream_cover
         call epreal1c_cp
      end if
#endif

      !Exchange Grid Between reciprocal process
      call commGridFront( qgrid2in_2d,r_wait )
c
#ifdef _OPENACC
      call cufft2d_frontmpi(qgrid2in_2d,qgrid2out_2d,n1mpimax,n2mpimax,
     &                      n3mpimax)
#else
      call   fft2d_frontmpi(qgrid2in_2d,qgrid2out_2d,n1mpimax,n2mpimax,
     &                      n3mpimax)
#endif
c
c     account for zeroth grid point for nonperiodic system
c
      call timer_enter( timer_scalar )
!$acc serial async(rec_queue) default(present) present(eprec,er)
      if (.not. use_bounds) then
         expterm = 0.5_ti_p * pi / xbox
         struc2  = qgrid2in_2d(1,1,1,1,1)**2 +
     &             qgrid2in_2d(2,1,1,1,1)**2
         er       = f * expterm * struc2
         eprec   = eprec + er
      end if
      if ((istart2(rankloc+1).eq.1) .and.
     &    (jstart2(rankloc+1).eq.1) .and.
     &    (kstart2(rankloc+1).eq.1)) then
         if (.not. use_bounds) then
            expterm = 0.5_ti_p * pi / xbox
            struc2  = qgrid2in_2d(1,1,1,1,1)**2 +
     &                qgrid2in_2d(2,1,1,1,1)**2
            er       = f * expterm * struc2
            eprec   = eprec + er
        end if
      end if
      er = 0.0_ti_p
!$acc end serial
c
c     complete the transformation of the PME grid
c
!$acc parallel loop collapse(3) async(rec_queue) default(present)
      do k = 1, ksize2(rankloc+1)
         do j = 1, jsize2(rankloc+1)
           do i = 1, isize2(rankloc+1)
             term                  = qfac_2d(i,j,k)
             qgrid2out_2d(1,i,j,k) = term*qgrid2out_2d(1,i,j,k)
             qgrid2out_2d(2,i,j,k) = term*qgrid2out_2d(2,i,j,k)
           end do
         end do
      end do
      call timer_exit( timer_scalar,quiet_timers )
c
c     perform 3-D FFT backward transform and get potential
c
#ifdef _OPENACC
      call cufft2d_backmpi(qgrid2in_2d,qgrid2out_2d,n1mpimax,n2mpimax,
     &                     n3mpimax)
#else
      call   fft2d_backmpi(qgrid2in_2d,qgrid2out_2d,n1mpimax,n2mpimax,
     &                     n3mpimax)
#endif

      !Exchange Grid Between reciprocal process
      call commGridBack( qgrid2in_2d,r_comm )

      ! Wait for Grid communication
      call commGridBack( qgrid2in_2d,r_wait )

#ifdef _OPENACC
      ! sync streams
      if (dir_queue.ne.rec_queue) call end_dir_stream_cover
#endif

c     fphip(1,:) = 0.0_ti_p !flush first line to zero? never use after
c     fphid(1,:) = 0.0_ti_p !flush first line to zero? same
      call timer_enter ( timer_grid2 )
      call fphi_uind_site1_p(fphid,fphip,fphidp)

!$acc parallel loop collapse(2) async(rec_queue) default(present)
      do ii = 1, npolerecloc
         do j = 1, 10
            fphid(j,ii) = electric * fphid(j,ii)
            fphip(j,ii) = electric * fphip(j,ii)
         end do
      end do
!$acc parallel loop collapse(2) async(rec_queue)
      do ii = 1, npolerecloc
         do j = 1, 20
            fphidp (j,ii) = electric * fphidp(j,ii)
            fphirec(j,ii) = electric * fphirec(j,ii)
         end do
      end do
      call timer_exit( timer_grid2,quiet_timers )
c
c     increment the induced dipole energy and gradient
c
      call timer_enter( timer_fmanage )
!$acc parallel loop reduction(+:er)
!$acc&         async(rec_queue)
!$acc&         default(present) present(er)
      do i = 1, npolerecloc
         iipole  = polerecglob(i)
         iglob   = ipole(iipole)
         ii      = locrec(iglob)
         f1      = 0.0_ti_p
         f2      = 0.0_ti_p
         f3      = 0.0_ti_p

!$acc loop seq reduction(+:f1,f2,f3)
         do k = 1, 3
            j1 = deriv1(k+1)
            j2 = deriv2(k+1)
            j3 = deriv3(k+1)
            er = er +  fuind(k,i)             *fphirec(k+1,i)
            f1 = f1 + (fuind(k,i)+fuinp(k ,i))*fphirec(j1,i)
     &              +  fuind(k,i)*fphip(j1,i)
     &              +  fuinp(k,i)*fphid(j1,i)
            f2 = f2 + (fuind(k,i)+fuinp(k ,i))*fphirec(j2,i)
     &              +  fuind(k,i)*fphip(j2,i)
     &              +  fuinp(k,i)*fphid(j2,i)
            f3 = f3 + (fuind(k,i)+fuinp(k ,i))*fphirec(j3,i)
     &              +  fuind(k,i)*fphip(j3,i)
     &              +  fuinp(k,i)*fphid(j3,i)
         end do
!$acc loop seq reduction(+:f1,f2,f3)
         do k = 1, 10
            f1 = f1 + fmp(k,i)*fphidp(deriv1(k),i)
            f2 = f2 + fmp(k,i)*fphidp(deriv2(k),i)
            f3 = f3 + fmp(k,i)*fphidp(deriv3(k),i)
         end do

         f1           = 0.5_ti_p * real(nfft1,t_p) * f1
         f2           = 0.5_ti_p * real(nfft2,t_p) * f2
         f3           = 0.5_ti_p * real(nfft3,t_p) * f3
         h1           = recip(1,1)*f1 + recip(1,2)*f2 + recip(1,3)*f3
         h2           = recip(2,1)*f1 + recip(2,2)*f2 + recip(2,3)*f3
         h3           = recip(3,1)*f1 + recip(3,2)*f2 + recip(3,3)*f3
         deprec(1,ii) = deprec(1,ii) + h1
         deprec(2,ii) = deprec(2,ii) + h2
         deprec(3,ii) = deprec(3,ii) + h3
      end do
!$acc serial async(rec_queue) present(eprec,er)
      eprec = eprec + 0.5*er
!$acc end serial
c
c     set the potential to be the induced dipole average
c
!$acc parallel loop collapse(2) async(rec_queue) default(present)
      do i = 1, npolerecloc
         do k = 1, 10
            fphidp(k,i) = 0.5_ti_p * fphidp(k,i)
         end do
      end do
      call fphi_to_cphi_sitegpu(fphidp,cphirec)
c
c     distribute torques into the induced dipole gradient
c
!$acc parallel loop async(rec_queue) default(present)
      do i = 1, npolerecloc
         iipole      = polerecglob(i)
         trqrec(1,i) = cmp(4,i)*cphirec(3,i) - cmp(3,i)*cphirec(4,i)
     &               + 2.0_ti_p*(cmp(7,i)-cmp(6,i))*cphirec(10,i)
     &               + cmp(9,i)*cphirec(8,i) + cmp(10,i)*cphirec(6,i)
     &               - cmp(8,i)*cphirec(9,i) - cmp(10,i)*cphirec(7,i)
         trqrec(2,i) = cmp(2,i)*cphirec(4,i) - cmp(4,i)*cphirec(2,i)
     &               + 2.0_ti_p*(cmp(5,i)-cmp(7,i))*cphirec(9,i)
     &               + cmp(8,i)*cphirec(10,i) + cmp(9,i)*cphirec(7,i)
     &               - cmp(9,i)*cphirec(5,i) - cmp(10,i)*cphirec(8,i)
         trqrec(3,i) = cmp(3,i)*cphirec(2,i) - cmp(2,i)*cphirec(3,i)
     &               + 2.0_ti_p*(cmp(6,i)-cmp(5,i))*cphirec(8,i)
     &               + cmp(8,i)*cphirec(5,i) + cmp(10,i)*cphirec(9,i)
     &               - cmp(8,i)*cphirec(6,i) - cmp(9,i)*cphirec(10,i)
      end do

      call torquegpu(npolerecloc,nlocrec2,polerecglob,locrec
     &              ,trqrec,deprec,rec_queue)

      call timer_exit ( timer_fmanage,quiet_timers )

      if (.not.use_virial) goto 100

      siz = 2*isize2(rankloc+1)*jsize2(rankloc+1)*ksize2(rankloc+1)
      call prmem_request(qgrip_b,siz,async=.false.)
      qgrip(1:2,1:isize2(rankloc+1),1:jsize2(rankloc+1)
     &      ,1:ksize2(rankloc+1)) => qgrip_b(1:)

c
c     induced dipole contribution to the internal virial
c
      call timer_enter(timer_other)
!$acc parallel loop gang vector
!$acc&         private(cphim,cphid,cphip)
!$acc&         reduction(+:vxx,vxy,vxz,vyy,vyz,vzz)
!$acc&         present(polerecglob,fphirec,fphid,fphip,uind,uinp
!$acc&   ,cphirec,cmp,vxx,vxy,vxz,vyy,vyz,vzz)
!$acc&         async(rec_queue)
      do i = 1, npolerecloc
         iipole = polerecglob(i)
         do j = 2, 4
            cphim(j) = 0.0_ti_p
            cphid(j) = 0.0_ti_p
            cphip(j) = 0.0_ti_p
            do k = 2, 4
               cphim(j) = cphim(j) + ftc(j,k)*fphirec(k,i)
               cphid(j) = cphid(j) + ftc(j,k)*fphid  (k,i)
               cphip(j) = cphip(j) + ftc(j,k)*fphip  (k,i)
            end do
         end do
         vxx = vxx - cphirec(2,i)*cmp(2,i)
     &         - 0.5*(cphim(2)*(uind(1,iipole)+uinp(1,iipole))
     &         +cphid(2)*uinp(1,iipole)+cphip(2)*uind(1,iipole))
         vxy = vxy - 0.5*(cphirec(2,i)*cmp(3,i)+cphirec(3,i)*
     &         cmp(2,i))
     &         - 0.25*(cphim(2)*(uind(2,iipole)+uinp(2,iipole))
     &         +cphim(3)*(uind(1,iipole)+uinp(1,iipole))
     &         +cphid(2)*uinp(2,iipole)+cphip(2)*uind(2,iipole)
     &         +cphid(3)*uinp(1,iipole)+cphip(3)*uind(1,iipole))
         vxz = vxz - 0.5*(cphirec(2,i)*cmp(4,i)+cphirec(4,i)*
     &         cmp(2,i))
     &         - 0.25*(cphim(2)*(uind(3,iipole)+uinp(3,iipole))
     &         +cphim(4)*(uind(1,iipole)+uinp(1,iipole))
     &         +cphid(2)*uinp(3,iipole)+cphip(2)*uind(3,iipole)
     &         +cphid(4)*uinp(1,iipole)+cphip(4)*uind(1,iipole))
         vyy = vyy - cphirec(3,i)*cmp(3,i)
     &         - 0.5*(cphim(3)*(uind(2,iipole)+uinp(2,iipole))
     &         +cphid(3)*uinp(2,iipole)+cphip(3)*uind(2,iipole))
         vyz = vyz - 0.5*(cphirec(3,i)*cmp(4,i)+cphirec(4,i)*cmp(3,i))
     &         - 0.25*(cphim(3)*(uind(3,iipole)+uinp(3,iipole))
     &         +cphim(4)*(uind(2,iipole)+uinp(2,iipole))
     &         +cphid(3)*uinp(3,iipole)+cphip(3)*uind(3,iipole)
     &         +cphid(4)*uinp(2,iipole)+cphip(4)*uind(2,iipole))
         vzz = vzz - cphirec(4,i)*cmp(4,i)
     &         - 0.5*(cphim(4)*(uind(3,iipole)+uinp(3,iipole))
     &         +cphid(4)*uinp(3,iipole)+cphip(4)*uind(3,iipole))
         vxx = vxx - 2.0*cmp(5,i)*cphirec(5,i) - cmp(8,i)*cphirec(8,i)
     &         - cmp(9,i)*cphirec(9,i)
         vxy = vxy - (cmp(5,i)+cmp(6,i))*cphirec(8,i)
     &         - 0.5*(cmp(8,i)*(cphirec(6,i)+cphirec(5,i))
     &         +cmp(9,i)*cphirec(10,i)+cmp(10,i)*cphirec(9,i))
         vxz = vxz - (cmp(5,i)+cmp(7,i))*cphirec(9,i)
     &         - 0.5*(cmp(9,i)*(cphirec(5,i)+cphirec(7,i))
     &          +cmp(8,i)*cphirec(10,i)+cmp(10,i)*cphirec(8,i))
         vyy = vyy - 2.0*cmp(6,i)*cphirec(6,i) - cmp(8,i)*
     &         cphirec(8,i)
     &         - cmp(10,i)*cphirec(10,i)
         vyz = vyz - (cmp(6,i)+cmp(7,i))*cphirec(10,i)
     &         - 0.5*(cmp(10,i)*(cphirec(6,i)+cphirec(7,i))
     &         +cmp(8,i)*cphirec(9,i)+cmp(9,i)*cphirec(8,i))
         vzz = vzz - 2.0*cmp(7,i)*cphirec(7,i) -
     &             cmp(9,i)*cphirec(9,i)
     &             - cmp(10,i)*cphirec(10,i)
      end do
      call timer_exit( timer_other,quiet_timers )

c
c     assign permanent and induced multipoles to the PME grid
c     and perform the 3-D FFT forward transformation
c
c
c    zero out the grid
c
      call timer_enter( timer_grid1 )
      call set_to_zero1(qgridin_2d(1,1,1,1,1),
     &            2*n1mpimax*n2mpimax*n3mpimax*(nrec_send+1),rec_queue)
c
!$acc parallel loop collapse(2) async(rec_queue) default(present)
      do i = 1, npolerecloc
         do j = 2, 4
            iipole = polerecglob(i)
            cmp(j,i) = cmp(j,i) + uinp(j-1,iipole)
         end do
      end do

      call cmp_to_fmp_sitegpu(cmp,fmp)
      call grid_mpole_site_p(fmp)
      call timer_exit( timer_grid1,quiet_timers )
c
      !Exchange Grid Between reciprocal process
      call commGridFront ( qgridin_2d,r_comm )

      ! Wait for Grid communication
      call commGridFront ( qgridin_2d,r_wait )
c
#ifdef _OPENACC
      call cufft2d_frontmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
     $                    n3mpimax)
#else
      call   fft2d_frontmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
     $                    n3mpimax)
#endif

      ! save qgridout_2d inside qgriq
      call timer_enter(timer_other)
      call utils_amove(2*isize2(rankloc+1)*jsize2(rankloc+1)
     &                  *ksize2(rankloc+1),qgridout_2d,qgrip,rec_queue)
      call timer_exit ( timer_other,quiet_timers )
c
c     zero out the PME grid
c
      call timer_enter( timer_grid1 )
      call set_to_zero1(qgridin_2d(1,1,1,1,1),
     &      2*n1mpimax*n2mpimax*n3mpimax*(nrec_send+1),rec_queue)

!$acc parallel loop collapse(2) async(rec_queue) default(present)
      do i = 1, npolerecloc
         do j = 2, 4
            iipole = polerecglob(i)
            cmp(j,i) = cmp(j,i) + uind(j-1,iipole) - uinp(j-1,iipole)
         end do
      end do
c
      call cmp_to_fmp_sitegpu(cmp,fmp)
      call grid_mpole_site_p(fmp)
      call timer_exit( timer_grid1,quiet_timers )
c
      !Exchange Grid Between reciprocal process
      call commGridFront ( qgridin_2d,r_comm )

      ! Wait for Grid communication
      call commGridFront ( qgridin_2d,r_wait )
c
#ifdef _OPENACC
      call cufft2d_frontmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
     &                      n3mpimax)
#else
      call   fft2d_frontmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
     &                      n3mpimax)
#endif
      istat2 = istart2(rankloc+1)
      ied2   = iend2  (rankloc+1)
      jstat2 = jstart2(rankloc+1)
      jed2   = jend2  (rankloc+1)
      kstat2 = kstart2(rankloc+1)
      ked2   = kend2  (rankloc+1)
c
c     make the scalar summation over reciprocal lattice
c
      call timer_enter(timer_other)
!$acc serial async(rec_queue) default(present)
      if ((istat2.eq.1).and.(jstat2.eq.1).and.(kstat2.eq.1)) then
           qfac_2d(1,1,1) = 0.0_ti_p
      end if
!$acc end serial
!$acc parallel loop collapse(3) async(rec_queue) default(present)
!$acc&         present(vxx,vxy,vxz,vyy,vyz,vzz)
!$acc&         reduction(+:vxx,vxy,vxz,vyy,vyz,vzz)
      do k3 = kstat2,ked2
        do k2 = jstat2,jed2
          do k1 = istat2,ied2
             m1   = k1 - 1
             m2   = k2 - 1
             m3   = k3 - 1
             if (k1 .gt. nf1)  m1 = m1 - nfft1
             if (k2 .gt. nf2)  m2 = m2 - nfft2
             if (k3 .gt. nf3)  m3 = m3 - nfft3
             if ((m1.eq.0).and.(m2.eq.0).and.(m3.eq.0)) cycle
             r1   = real(m1,t_p)
             r2   = real(m2,t_p)
             r3   = real(m3,t_p)
             h1   = recip(1,1)*r1 + recip(1,2)*r2 + recip(1,3)*r3
             h2   = recip(2,1)*r1 + recip(2,2)*r2 + recip(2,3)*r3
             h3   = recip(3,1)*r1 + recip(3,2)*r2 + recip(3,3)*r3
             hsq  = h1*h1 + h2*h2 + h3*h3
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
             struc2 = qgridout_2d(1,k1-istat2+1,k2-jstat2+1,k3-kstat2+1)
     &                     *qgrip(1,k1-istat2+1,k2-jstat2+1,k3-kstat2+1)
     &              + qgridout_2d(2,k1-istat2+1,k2-jstat2+1,k3-kstat2+1)
     &                     *qgrip(2,k1-istat2+1,k2-jstat2+1,k3-kstat2+1)
                eterm = 0.5_ti_p * f * expterm * struc2
                vterm = (2.0_ti_p/hsq) * (1.0_ti_p-term) * eterm
                vxx   = vxx + h1*h1*vterm - eterm
                vxy   = vxy + h1*h2*vterm
                vxz   = vxz + h1*h3*vterm
                vyy   = vyy + h2*h2*vterm - eterm
                vyz   = vyz + h2*h3*vterm
                vzz   = vzz + h3*h3*vterm - eterm
             end if
             qfac_2d(k1-istat2+1,k2-jstat2+1,k3-kstat2+1) = expterm
          end do
        end do
      end do
c
c     increment the internal virial tensor components
c
c
c     Proceed to atomic update to avoid collision with direct queue
c     even if it's highly unlikely
c
!$acc serial async(rec_queue) default(present)
!$acc&       present(vir,vxx,vxy,vxz,vyy,vyz,vzz)
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
      call timer_exit(timer_other,quiet_timers)
c

 100  continue

c!$acc update host(vxx,vxy,vxz,eprec) async(rec_queue)
c!$acc update host(deprec,fmp) async
c!$acc update host(fphip,fphid,fphidp) async(rec_queue)
c      do i = 0,nproc-1
c         if (rank.eq.i) then
c         print*,comput_norm( uind,3*n,2),
c     &          comput_norm( uinp,3*n,2),
c     &          comput_norm(fuind,3*npolerecloc,2)
c         print*,comput_normr(deprec,3*nlocrec2,2)
c     &         ,sum(abs(a))
c     &         ,comput_norm(fuinp,3*npolerecloc,2)
c          print*,vxx,vxy,eprec,ep
c          print*,comput_norm(qgrid2in_2d,
c     &       2*n1mpimax*n2mpimax*n3mpimax*(nrec_send+1),2)
c         end if
c         call mpi_barrier(MPI_COMM_WORLD,j)
c      end do

      call timer_exit( timer_eprecip )
      end
