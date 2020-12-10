c
c     sorbonne university
c     washington university in saint louis
c     university of texas at austin
c
c     ###################################################################
c     ##                                                               ##
c     ##  subroutine ehal1gpu  --  buffered 14-7 energy & derivatives  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     "ehal1" calculates the buffered 14-7 van der waals energy and
c     its first derivatives with respect to cartesian coordinates
c
c
#include "tinker_precision.h"
#include "tinker_types.h"

      module ehal1gpu_inl
        integer(1) one1,two1
        integer onei
        real(t_p)  half,one
        parameter( half=0.5, one=1.0, one1=1, two1=2,
     &             onei=1 )
        contains
#include "convert.f.inc"
#include "image.f.inc"
#include "pair_ehal1.f.inc"
      end module

      subroutine ehal1gpu

      use domdec     ,only: ndir,rank
      use energi     ,only: ev
      use interfaces ,only: ehal1c_p,ehalshort1c_p
     &               ,ehallong1c_p
      use utilgpu    ,only: def_queue,dir_queue
      use potent
      use virial     ,only: vir
      use vdwpot     ,only: use_vcorr
      implicit none
      real(r_p) elrc,vlrc
c
      ! PME-Core case
      if (use_pmecore.and.rank.ge.ndir) return
c
c     choose the method for summing over pairwise interactions
c
      def_queue = dir_queue
      if (use_vdwshort) then
        call ehalshort1c_p
      else if (use_vdwlong) then
        call ehallong1c_p
      else
        call ehal1c_p
      end if
c
c     apply long range van der Waals correction if desired
c
      if (use_vcorr) then
         print*, 'vcorr'
!$acc data  create(elrc,vlrc) async(def_queue)
!$acc&      present(ev,vir)

         call evcorr1gpu (elrc,vlrc)
!$acc kernels async(def_queue)
         ev = ev + elrc
         vir(1,1) = vir(1,1) + vlrc
         vir(2,2) = vir(2,2) + vlrc
         vir(3,3) = vir(3,3) + vlrc
!$acc end kernels

!$acc end data
      end if
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine ehal1c  --  buffered 14-7 vdw derivs via list  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ehal1c" calculates the buffered 14-7 van der Waals energy and
c     its first derivatives with respect to Cartesian coordinates
c     using a pairwise neighbor list
c
c
      subroutine ehal1cgpu1
      use atmlst
      use atoms
      use bound
      use couple
      use deriv
      use domdec
      use ehal1gpu_inl
      use energi
      use group
      use iounit
      use inform ,only: deb_Path
      use inter
      use molcul
      use mutant ,only: mut=>mutInt
     &           , vcouple,vlambda,scalpha,scexp
      use neigh
      use shunt
      use usage
      use virial
      use vdw
      use vdwpot
      use vec
      use vec_vdw
      use utilgpu

      implicit none
      integer i,j,k,kk,ksave
      integer kt,kglob,kbis,kvloc,kv,ki
      integer iglob,iivdw
      integer ii,iv,it,ivloc
      integer nnvlst,nnvlst2
      integer nn12,nn13,nn14,ntot
      integer interac
      integer iscal(maxscaling), inter_(maxscaling)
      integer(1) muti,mutik

      real(r_p)  vxx_vdw,vxy_vdw,vxz_vdw
      real(r_p)  vyy_vdw,vyz_vdw,vzz_vdw
      real(t_p)  xi,yi,zi,redi,e,de
      real(t_p)  rdn,rdn1,redk
      real(t_p)  invrik,rik,rik2,rik3,rik4,rik5,rik6,rik7
      real(t_p)  dedx,dedy,dedz
      real(r_p)  devx,devy,devz,devt
      real(t_p)  invrho,rv7orho
      real(t_p)  dtau,gtau,tau,tau7,rv7
      real(t_p)  rv2,eps2
      real(t_p)  xpos,ypos,zpos
      real(t_p)  dtaper,taper
      real(t_p)  vscale
      real(t_p)  fscal(maxscaling)
      real(t_p)  xred(nbloc)
      real(t_p)  yred(nbloc)
      real(t_p)  zred(nbloc)
      character*10 mode
c
      if(deb_Path) write (*,*) 'ehal1cgpu1'

#ifdef _OPENACC
      if (dir_queue.ne.rec_queue) then
!!$acc wait(rec_queue) async(dir_queue)
         call stream_wait_async(rec_stream,dir_stream,rec_event)
      end if
#endif

!$acc data create(xred,yred,zred,vxx_vdw,vxy_vdw,vxz_vdw,vyy_vdw,
!$acc&            vyz_vdw,vzz_vdw)
!$acc&     present(loc,ired,kred,x,y,z,vdwglobnl,ivdw,loc,jvdw,
!$acc&  vir,dev,vdwglob,vlst,nvlst,radmin,radmin4,epsilon,mut,
!$acc&  epsilon4,warning,vir,dev)
!$acc&     present(ev)
!$acc&     async(def_queue)

!$acc serial async(def_queue)
c
c     zero out the van der Waals energy and first derivatives
c
      ev      = 0.0_re_p
c
c     increment the internal virial tensor components
c
      vxx_vdw = 0.0_ti_p
      vxy_vdw = 0.0_ti_p
      vxz_vdw = 0.0_ti_p
      vyy_vdw = 0.0_ti_p
      vyz_vdw = 0.0_ti_p
      vzz_vdw = 0.0_ti_p
!$acc end serial

c
c     apply any reduction factor to the atomic coordinates
c
!$acc parallel loop default(present)
!$acc&         async(def_queue)
      do k = 1,nvdwbloc
         iglob   = ivdw (vdwglob (k))
         i       = loc  (iglob)
         iv      = ired (iglob)
         rdn     = kred (iglob)
         rdn1    = 1.0_ti_p - rdn
         xred(i) = rdn * x(iglob) + rdn1 * x(iv)
         yred(i) = rdn * y(iglob) + rdn1 * y(iv)
         zred(i) = rdn * z(iglob) + rdn1 * z(iv)
      enddo
c
c     set the coefficients for the switching function
c
      mode = 'VDW'
      call switch (mode)
c
c     find van der Waals energy and derivatives via neighbor list
c
!$acc parallel loop gang vector_length(32)
!$acc&         private(ki,iscal,fscal,inter_)
!$acc&         async(def_queue)
      MAINLOOP:
     &do ii = 1, nvdwlocnl
         iivdw = vdwglobnl(ii)
         iglob = ivdw(iivdw)
         i     = loc(iglob)
         if (i.eq.0) then
            print*,warning
            print*,sugest_vdw
            cycle MAINLOOP
         end if

         nnvlst = nvlst(ii)
         if(nnvlst.eq.0) cycle MAINLOOP
         iv    = ired(iglob)
         redi  = merge (1.0_ti_p,kred(iglob),(iglob.eq.iv))
         ivloc = loc(iv)
         it    = jvdw(iglob)
         muti  = mut(iglob)
         xi    = xred(i)
         yi    = yred(i)
         zi    = zred(i)
         ki    = 0
c
c     get number of atoms  directly (1-2) or 1-3 or 1-4 bonded
c
         nn12  = n12(iglob)
         nn13  = n13(iglob) + nn12 
         nn14  = n14(iglob) + nn13 
         ntot  = n15(iglob) + nn14 
         if (ntot.gt.maxscaling) 
     &      print*,'too many scaling factor, ehal'
c
c     fill scaling factor along kglob and interaction type
c
!$acc loop vector
         do j=1,nn13
            if      (j.le.nn12) then
               iscal(j) = i12 (j,iglob)
               fscal(j) = v2scale
               inter_(j) = 4
            else
               iscal(j) = i13 (j-nn12,iglob)
                fscal(j) = v3scale
               inter_(j) = 8
            end if
         end do
!$acc loop vector
         do j=nn13+1,ntot
            if      (j.le.nn14) then
               iscal(j) = i14 (j-nn13,iglob)
               fscal(j) = v4scale
               inter_(j) = 16
            else
               iscal(j) = i15 (j-nn14,iglob)
               fscal(j) = v5scale
               inter_(j) = 32
            end if
         end do

!$acc loop vector 
         do k = 1, nnvlst
            kglob  = vlst(k,ii)
            kbis   = loc (kglob)
            kvloc  = loc (ired(kglob))
            kt     =    jvdw (kglob)
            mutik  = muti + mut(kglob)
            if ((kvloc.eq.0).or.(kvloc.gt.nbloc).or.(kbis.gt.nbloc))
     &      then
               print*, warning
               print*, sugest_vdw
               cycle
            end if
            if (kbis.ne.kvloc) then
               redk = kred(kglob)
            else
               redk = 1.0_ti_p
            endif
            vscale  = 1.0_ti_p
            interac = 0
c
c     compute the energy contribution for this interaction
c
            xpos   = xi - xred(kbis)
            ypos   = yi - yred(kbis)
            zpos   = zi - zred(kbis)
            call image_inl(xpos,ypos,zpos)
c
c     decide whether to compute the current interaction
c     and check for an interaction distance less than the cutoff
c
            rik2   = xpos**2 + ypos**2 + zpos**2
            if (rik2>off2) cycle
            ! Annihilation
            if (vcouple.eq.onei.and.mutik.eq.two1) mutik=one1
c
c     set interaction vscale coefficients for connected atoms
c
            if (ki<ntot) then
!$acc loop seq
               do j=1,ntot
                  if (iscal(j)==kglob) then
                     vscale  = fscal(j)
                     interac = inter_(j)
!$acc atomic update
                     ki  = ki+1
                     exit
                  end if
               end do
            end if
c
c     replace 1-4 interactions
c
            if (btest(interac,4)) then
               rv2  =  radmin4(kt,it)
               eps2 = epsilon4(kt,it)
            else
               rv2  =  radmin (kt,it)
               eps2 = epsilon (kt,it)
            endif

            call ehal1_couple(xpos,ypos,zpos,rik2,rv2,eps2,vscale
     &                       ,cut2,cut,off,ghal,dhal
     &                       ,scexp,vlambda,scalpha,mutik
     &                       ,e,dedx,dedy,dedz)
            ev    =   ev + e
c
c          increment the van der Waals derivatives
c
            devx = redk*dedx
            devy = redk*dedy
            devz = redk*dedz
!$acc atomic update
            dev(1,kbis)  = dev(1,kbis)  - devx
!$acc atomic update
            dev(2,kbis)  = dev(2,kbis)  - devy
!$acc atomic update
            dev(3,kbis)  = dev(3,kbis)  - devz

            if (kbis.ne.kvloc) then
               devx = (1.0_ti_p - redk)* dedx
               devy = (1.0_ti_p - redk)* dedy
               devz = (1.0_ti_p - redk)* dedz
!$acc atomic update
               dev(1,kvloc) = dev(1,kvloc) - devx
!$acc atomic update
               dev(2,kvloc) = dev(2,kvloc) - devy
!$acc atomic update
               dev(3,kvloc) = dev(3,kvloc) - devz
            end if

            devx  = redi * dedx
            devy  = redi * dedy
            devz  = redi * dedz
!$acc atomic update
            dev(1,i) = dev(1,i) + devx
!$acc atomic update
            dev(2,i) = dev(2,i) + devy
!$acc atomic update
            dev(3,i) = dev(3,i) + devz

            if (iglob.ne.iv) then
               devx  = (1.0_ti_p - redi)* dedx
               devy  = (1.0_ti_p - redi)* dedy
               devz  = (1.0_ti_p - redi)* dedz
!$acc atomic update
               dev(1,ivloc) = dev(1,ivloc) + devx
!$acc atomic update
               dev(2,ivloc) = dev(2,ivloc) + devy
!$acc atomic update
               dev(3,ivloc) = dev(3,ivloc) + devz
            end if
c
c     increment the total van der Waals energy 
c
            vxx_vdw = vxx_vdw + xpos * dedx
            vxy_vdw = vxy_vdw + ypos * dedx
            vxz_vdw = vxz_vdw + zpos * dedx
            vyy_vdw = vyy_vdw + ypos * dedy
            vyz_vdw = vyz_vdw + zpos * dedy
            vzz_vdw = vzz_vdw + zpos * dedz

         end do
      end do MAINLOOP

!$acc kernels async(def_queue)
      vir(1,1) = vir(1,1) + vxx_vdw
      vir(2,1) = vir(2,1) + vxy_vdw
      vir(3,1) = vir(3,1) + vxz_vdw
      vir(1,2) = vir(1,2) + vxy_vdw
      vir(2,2) = vir(2,2) + vyy_vdw
      vir(3,2) = vir(3,2) + vyz_vdw
      vir(1,3) = vir(1,3) + vxz_vdw
      vir(2,3) = vir(2,3) + vyz_vdw
      vir(3,3) = vir(3,3) + vzz_vdw
!$acc end kernels

!$acc end data
      end

      subroutine ehal1cgpu2
      use atmlst    ,only: vdwglobnl,vdwglob
      use atoms     ,only: x,y,z,n
      use couple    ,only: i12,n12
      use deriv     ,only: dev=>debond
      use domdec    ,only: loc,rank,nbloc
      use ehal1gpu_inl
      use energi    ,only: ev=>ev_r
      use inform    ,only: deb_Path
      use interfaces,only: ehal1c_correct_scaling
      use mutant    ,only: scalpha,scexp,vlambda,vcouple,mut=>mutInt
      use neigh     ,only: vlst,nvlst
      use tinheader ,only: ti_p,re_p
      use tinMemory ,only: prmem_request
      use sizes     ,only: maxvalue,tinkerdebug
      use shunt     ,only: c0,c1,c2,c3,c4,c5,off2,off,cut2,cut
      use vdw       ,only: ired,kred,jvdw,ivdw,radmin,
     &                     epsilon,nvdwbloc,nvdwlocnl,
     &                     nvdwclass,skipvdw12
      use vdwpot    ,only: vcorrect_ik,vcorrect_scale,n_vscale,dhal,ghal
      use vdw_locArray
      use utilgpu   ,only: def_queue,dir_queue,rec_queue
#ifdef _OPENACC
     &                    ,dir_stream
     &                    ,rec_stream,rec_event,stream_wait_async
#endif
      use virial
      implicit none
      integer i,j,k,kk,ksave
      integer kt,kglob,kbis,kvloc,kv,ki
      integer iglob,iivdw
      integer ii,iv,it,ivloc
      integer,save:: ncall=0
      integer nnvlst,nnvlst2
      integer nn12,nn13,nn14,ntot
      integer in12,ai12(maxvalue)
#ifdef TINKER_DEBUG
      integer interac(n)
#endif
      integer(1) muti,mutik
      logical ik12
      real(t_p) xi,yi,zi,redi,e,de
      real(t_p) rdn,rdn1,redk
      real(t_p) invrik,rik,rik2,rik3,rik4,rik5,rik6,rik7
      real(t_p) dedx,dedy,dedz
      mdyn_rtyp devx,devy,devz
      real(t_p) invrho,rv7orho
      real(t_p) dtau,gtau,tau,tau7,rv7
      real(t_p) rv2,eps2
      real(t_p) xpos,ypos,zpos
      real(t_p) dtaper,taper
      real(t_p) vscale,vscale4
      logical   do_scale4
      character*10 mode

c
      if(deb_Path) write (*,*) 'ehal1cgpu2'
      ncall = ncall + 1

#ifdef _OPENACC
      if (dir_queue.ne.rec_queue) then
!!$acc wait(rec_queue) async(dir_queue)
         call stream_wait_async(rec_stream,dir_stream,rec_event)
      end if
#endif

      call prmem_request(xred,nbloc,queue=def_queue)
      call prmem_request(yred,nbloc,queue=def_queue)
      call prmem_request(zred,nbloc,queue=def_queue)
#ifdef TINKER_DEBUG
!$acc enter data create(interac) async(def_queue)
!$acc kernels present(interac) async(def_queue)
      interac(:) = 0
!$acc end kernels
#endif
      if (ncall.eq.1) then
!$acc wait
      end if

!$acc data present(xred,yred,zred)
!$acc&     present(loc,ired,kred,x,y,z,vdwglobnl,ivdw,loc,jvdw,
!$acc&  vir,dev,vdwglob,vlst,nvlst,radmin,epsilon,mut,i12,n12)
!$acc&     present(ev,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)

c
c     apply any reduction factor to the atomic coordinates
c
!$acc parallel loop default(present)
!$acc&         async(def_queue)
      do k = 1,nvdwbloc
         iglob   = ivdw (vdwglob (k))
         i       = loc  (iglob)
         iv      = ired (iglob)
         rdn     = kred (iglob)
         rdn1    = 1.0_ti_p - rdn
         xred(i) = rdn * x(iglob) + rdn1 * x(iv)
         yred(i) = rdn * y(iglob) + rdn1 * y(iv)
         zred(i) = rdn * z(iglob) + rdn1 * z(iv)
      enddo

      call resetForces_buff(def_queue)
c
c     set the coefficients for the switching function
c
      mode = 'VDW'
      call switch (mode)
c
c     find van der Waals energy and derivatives via neighbor list
c
!$acc parallel loop gang vector_length(32)
!$acc&         private(ai12)
#ifdef TINKER_DEBUG
!$acc&         present(interac)
#endif
!$acc&         async(def_queue)
      MAINLOOP:
     &do ii = 1, nvdwlocnl
         iivdw = vdwglobnl(ii)
         iglob = ivdw(iivdw)
         i     = loc(iglob)
         it    = jvdw(iglob)
         xi    = xred(i)
         yi    = yred(i)
         zi    = zred(i)
         muti  = mut(iglob)
         if (skipvdw12) then
            in12 = n12(iglob)
!$acc loop vector
            do j = 1,in12
               ai12(j) = i12(j,iglob)
            end do
         end if


!$acc loop vector 
         do k = 1, nvlst(ii)
            kglob  = vlst(k,ii)
            kbis   = loc  (kglob)
            kt     = jvdw (kglob)
            mutik  = muti + mut(kglob)

            if (skipvdw12) then
               ik12 = .false.
!$acc loop seq
               do j = 1, in12
                  if (ai12(j).eq.kglob) ik12=.true.
               end do
               if (ik12) cycle
            end if
            !vscale = 1.0_ti_p
c
c     compute the energy contribution for this interaction
c
            xpos   = xi - xred(kbis)
            ypos   = yi - yred(kbis)
            zpos   = zi - zred(kbis)
            call image_inl(xpos,ypos,zpos)
c
c     decide whether to compute the current interaction
c     and check for an interaction distance less than the cutoff
c
            rik2   = xpos**2 + ypos**2 + zpos**2
            if (rik2>off2) cycle

            ! Annihilate
            if (vcouple.eq.1.and.mutik.eq.two1) mutik=one1

            rv2  =  radmin (kt,it)
            eps2 = epsilon (kt,it)

            call ehal1_couple(xpos,ypos,zpos,rik2,rv2,eps2,1.0_ti_p
     &                       ,cut2,cut,off,ghal,dhal
     &                       ,scexp,vlambda,scalpha,mutik
     &                       ,e,dedx,dedy,dedz)

#ifdef TINKER_DEBUG
!$acc atomic
            interac(iglob) = interac(iglob) +1
c           if (rank.eq.0) then
c           if (iglob.eq.12)
c    &         print*, kglob,real(xpos,4),real(ypos,4),zpos
c           if (kglob.eq.12)
c    &         print*, iglob,real(xpos,4),real(ypos,4),zpos
c           end if
#endif

            !Increment the total van der Waals energy 
            ev   =   ev + tp2enr(e)

            devx = tp2mdr(dedx)
            devy = tp2mdr(dedy)
            devz = tp2mdr(dedz)
            !Increment the van der Waals derivatives
!$acc atomic
            dev(1,kbis)  = dev(1,kbis)  - devx
!$acc atomic
            dev(2,kbis)  = dev(2,kbis)  - devy
!$acc atomic
            dev(3,kbis)  = dev(3,kbis)  - devz

!$acc atomic
            dev(1,i) = dev(1,i) + devx
!$acc atomic
            dev(2,i) = dev(2,i) + devy
!$acc atomic
            dev(3,i) = dev(3,i) + devz

            ! Increment van der waals virial contribution
            g_vxx = g_vxx + xpos * dedx
            g_vxy = g_vxy + ypos * dedx
            g_vxz = g_vxz + zpos * dedx
            g_vyy = g_vyy + ypos * dedy
            g_vyz = g_vyz + zpos * dedy
            g_vzz = g_vzz + zpos * dedz

         end do
      end do MAINLOOP


#ifdef TINKER_DEBUG
 34   format(2I10,3F12.4)
 35   format(A30,I16,3x,F16.6,I3)
!$acc wait
!$acc exit data copyout(interac)
!$acc update host(dev,ev)
      write(*,35)'nev & ev & rank ',sum(interac),enr2en(ev),rank
c     if (rank.eq.0) then
c     do i =1,3000
c        print 34,i,interac(i)
c     end do
c     end if
#endif

      call ehal1c_correct_scaling(xred,yred,zred,
     &           g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)

      call vdw_gradient_reduce
!$acc end data

      end subroutine

      subroutine ehalshort1cgpu
      use atmlst    ,only: vdwglobnl,vdwglob
      use atoms     ,only: x,y,z,n
      use couple    ,only: i12,n12
      use cutoff    ,only: vdwshortcut,shortheal
      use deriv     ,only: dev=>debond
      use domdec    ,only: loc,rank,nbloc
      use ehal1gpu_inl
      use energi    ,only: ev=>ev_r
      use inform    ,only: deb_Path
      use interfaces,only: ehal1c_correct_scaling_shortlong
      use mutant    ,only: scalpha,scexp,vlambda,vcouple,mut=>mutInt
      use neigh     ,only: shortvlst,nshortvlst
      use tinheader ,only: ti_p,re_p
      use tinMemory ,only: prmem_request
      use shunt     ,only: c0,c1,c2,c3,c4,c5,off2,off,cut2
      use sizes     ,only: maxvalue,tinkerdebug
      use vdw       ,only: ired,kred,jvdw,ivdw,radmin,
     &                     epsilon,nvdwbloc,nvdwlocnl,
     &                     nvdwclass,skipvdw12
      use vdw_locArray
      use vdwpot    ,only: vcorrect_ik,vcorrect_scale,n_vscale,dhal,ghal
      use utilgpu   ,only: def_queue,dir_queue,rec_queue
#ifdef _OPENACC
     &                    ,dir_stream
     &                    ,rec_stream,rec_event,stream_wait_async
#endif
      use virial
      implicit none
      integer i,j,k,kk,ksave
      integer kt,kglob,kbis,kvloc,kv,ki
      integer iglob,iivdw
      integer ii,iv,it,ivloc
      integer nnvlst,nnvlst2
      integer nn12,nn13,nn14,ntot
      integer,save:: ncall=0
      integer in12,ai12(maxvalue)
#ifdef TINKER_DEBUG
      integer interac(n)
#endif
      integer(1) muti,mutik
      real(t_p)  xi,yi,zi,redi,e,de
      real(t_p)  rdn,rdn1,redk
      real(t_p)  invrik,rik,rik2,rik3,rik4,rik5,rik6,rik7
      real(t_p)  dedx,dedy,dedz
      mdyn_rtyp  devx,devy,devz
      real(t_p)  invrho,rv7orho
      real(t_p)  dtau,gtau,tau,tau7,rv7
      real(t_p)  rv2,eps2
      real(t_p)  xpos,ypos,zpos
      real(t_p)  dtaper,taper
      real(t_p)  vscale,vscale4
      logical    do_scale4,ik12
      character*10 mode

      if(deb_Path) write (*,*) 'ehalshort1cgpu'
      ncall = ncall + 1

#ifdef _OPENACC
      if (dir_queue.ne.rec_queue) then
!!$acc wait(rec_queue) async(dir_queue)
         call stream_wait_async(rec_stream,dir_stream,rec_event)
      end if
#endif

      call prmem_request(xred,nbloc,queue=def_queue)
      call prmem_request(yred,nbloc,queue=def_queue)
      call prmem_request(zred,nbloc,queue=def_queue)
#ifdef TINKER_DEBUG
!$acc enter data create(interac) async(def_queue)
!$acc kernels present(interac) async(def_queue)
      interac(:) = 0
!$acc end kernels
#endif

      if (ncall.eq.1) then
!$acc wait
      end if

!$acc data present(xred,yred,zred)
!$acc&     present(loc,ired,kred,x,y,z,vdwglobnl,ivdw,loc,jvdw,
!$acc&  vir,dev,vdwglob,shortvlst,nshortvlst,radmin,epsilon,mut,
!$acc&  i12,n12)
!$acc&     present(ev,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
#ifdef TINKER_DEBUG
!$acc&     present(interac)
#endif

c
c     apply any reduction factor to the atomic coordinates
c
!$acc parallel loop async(def_queue)
      do k = 1,nvdwbloc
         iglob   = ivdw (vdwglob (k))
         i       = loc  (iglob)
         iv      = ired (iglob)
         rdn     = kred (iglob)
         rdn1    = 1.0_ti_p - rdn
         xred(i) = rdn * x(iglob) + rdn1 * x(iv)
         yred(i) = rdn * y(iglob) + rdn1 * y(iv)
         zred(i) = rdn * z(iglob) + rdn1 * z(iv)
      enddo

      call resetForces_buff(def_queue)
c
c     set the coefficients for the switching function
c
      mode = 'SHORTVDW'
      call switch (mode)
c
c     find van der Waals energy and derivatives via neighbor list
c
!$acc parallel loop gang vector_length(32)
!$acc&         private(ai12)
!$acc&         async(def_queue)
      MAINLOOP:
     &do ii = 1, nvdwlocnl
         iivdw = vdwglobnl(ii)
         iglob = ivdw(iivdw)
         i     = loc(iglob)
         muti  = mut(iglob)
         if (skipvdw12) then
            in12 = n12(iglob)
!$acc loop vector
            do j = 1,in12
               ai12(j) = i12(j,iglob)
            end do
         end if

         it    = jvdw(iglob)
         xi    = xred(i)
         yi    = yred(i)
         zi    = zred(i)

!$acc loop vector 
         do k = 1, nshortvlst(ii)
            kglob  = shortvlst(k,ii)
            kbis   = loc (kglob)
            kt     = jvdw (kglob)
            mutik  = muti + mut(kglob)

            if (skipvdw12) then
               ik12 = .false.
!$acc loop seq
               do j = 1, in12
                  if (ai12(j).eq.kglob) ik12=.true.
               end do
               if (ik12) cycle
            end if
            xpos   = xi - xred(kbis)
            ypos   = yi - yred(kbis)
            zpos   = zi - zred(kbis)
            !vscale = 1.0_ti_p
            call image_inl(xpos,ypos,zpos)
c
c     decide whether to compute the current interaction
c     and check for an interaction distance less than the cutoff
c
            rik2   = xpos**2 + ypos**2 + zpos**2
            if (rik2>off2) cycle

            ! Annihilate
            if (vcouple.eq.1.and.mutik.eq.two1) mutik=one1

            rv2  =  radmin (kt,it)
            eps2 = epsilon (kt,it)

            call ehal1_couple_short(xpos,ypos,zpos,rik2,rv2,eps2
     &                 ,1.0_ti_p,cut2,off
     &                 ,scexp,vlambda,scalpha,mutik
     &                 ,shortheal,ghal,dhal,e,dedx,dedy,dedz)

#ifdef TINKER_DEBUG
!$acc atomic
            interac(iglob) = interac(iglob) +1
c           if (rank.eq.0) then
c           if (iglob.eq.12)
c    &         print*, kglob,real(xpos,4),real(ypos,4),zpos
c           if (kglob.eq.12)
c    &         print*, iglob,real(xpos,4),real(ypos,4),zpos
c           end if
#endif
            !Increment the total van der Waals energy 
            ev   =   ev + tp2enr(e)

            !Increment the van der Waals derivatives
            devx = tp2mdr(dedx)
            devy = tp2mdr(dedy)
            devz = tp2mdr(dedz)
!$acc atomic
            dev(1,kbis)  = dev(1,kbis)  - devx
!$acc atomic
            dev(2,kbis)  = dev(2,kbis)  - devy
!$acc atomic
            dev(3,kbis)  = dev(3,kbis)  - devz

!$acc atomic
            dev(1,i) = dev(1,i) + devx
!$acc atomic
            dev(2,i) = dev(2,i) + devy
!$acc atomic
            dev(3,i) = dev(3,i) + devz

            !Increment van der waals virial contribution
            g_vxx = g_vxx + xpos * dedx
            g_vxy = g_vxy + ypos * dedx
            g_vxz = g_vxz + zpos * dedx
            g_vyy = g_vyy + ypos * dedy
            g_vyz = g_vyz + zpos * dedy
            g_vzz = g_vzz + zpos * dedz

         end do
      end do MAINLOOP


#ifdef TINKER_DEBUG
 34   format(2I10,3F12.4)
 35   format(A40,I16,3x,F16.6,I3)
!$acc wait
!$acc exit data copyout(interac)
!$acc update host(dev,ev)
      write(*,35)'nev & ev short & rank ',sum(interac),enr2en(ev),rank
c     if (rank.eq.0) then
c     do i =1,3000
c        print 34,i,interac(i)
c     end do
c     end if
#endif

      mode = 'SHORT '
      call ehal1c_correct_scaling_shortlong(xred,yred,zred,
     &            g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz,mode)

      call vdw_gradient_reduce
!$acc end data

      end subroutine

      subroutine ehallong1cgpu
      use atmlst    ,only: vdwglobnl,vdwglob
      use atoms     ,only: x,y,z,n
      use couple    ,only: i12,n12
      use cutoff    ,only: vdwshortcut,shortheal
      use deriv     ,only: dev=>debond
      use domdec    ,only: loc,rank,nbloc
      use ehal1gpu_inl
      use energi    ,only: ev=>ev_r
      use inform    ,only: deb_Path
      use interfaces,only: ehal1c_correct_scaling_shortlong
      use mutant    ,only: scalpha,scexp,vlambda,vcouple,mut=>mutInt
      use neigh     ,only: vlst,nvlst
      use tinheader ,only: ti_p,re_p
      use tinMemory ,only: prmem_request
      use sizes     ,only: maxvalue,tinkerdebug
      use shunt     ,only: c0,c1,c2,c3,c4,c5,off2,off,cut2,cut
      use vdw       ,only: ired,kred,jvdw,ivdw,radmin,
     &                     epsilon,nvdwbloc,nvdwlocnl,
     &                     nvdwclass,skipvdw12
      use vdwpot    ,only: vcorrect_ik,vcorrect_scale,n_vscale,dhal,ghal
      use vdw_locArray
      use utilgpu   ,only: def_queue,dir_queue,rec_queue
#ifdef _OPENACC
     &                    ,dir_stream
     &                    ,rec_stream,rec_event,stream_wait_async
#endif
      use virial
      implicit none
      integer i,j,k,kk,ksave
      integer kt,kglob,kbis,kvloc,kv,ki
      integer iglob,iivdw
      integer ii,iv,it,ivloc
      integer nnvlst,nnvlst2
      integer nn12,nn13,nn14,ntot
      integer in12,ai12(maxvalue)
      integer,save:: ncall=0
#ifdef TINKER_DEBUG
      integer interac(n)
#endif
      real(t_p)  xi,yi,zi,redi,e,de
      real(t_p)  rdn,rdn1,redk
      real(t_p)  invrik,rik,rik2,rik3,rik4,rik5,rik6,rik7
      real(t_p)  dedx,dedy,dedz
      mdyn_rtyp  devx,devy,devz
      real(t_p)  invrho,rv7orho
      real(t_p)  dtau,gtau,tau,tau7,rv7
      real(t_p)  rv2,eps2
      real(t_p)  xpos,ypos,zpos
      real(t_p)  dtaper,taper
      real(t_p)  vscale,vscale4
      real(t_p)  vdwshortcut2
      logical    do_scale4,ik12
      integer(1) muti,mutik
      character*10 mode

c
      if(deb_Path) write (*,*) 'ehallong1cgpu'
      ncall = ncall + 1

#ifdef _OPENACC
      if (dir_queue.ne.rec_queue) then
!!$acc wait(rec_queue) async(dir_queue)
         call stream_wait_async(rec_stream,dir_stream,rec_event)
      end if
#endif

      call prmem_request(xred,nbloc,queue=def_queue)
      call prmem_request(yred,nbloc,queue=def_queue)
      call prmem_request(zred,nbloc,queue=def_queue)
#ifdef TINKER_DEBUG
!$acc enter data create(interac) async(def_queue)
!$acc kernels present(interac) async(def_queue)
      interac(:) = 0
!$acc end kernels
#endif

!$acc data present(xred,yred,zred)
!$acc&     present(loc,ired,kred,x,y,z,vdwglobnl,ivdw,loc,jvdw,
!$acc&  vir,dev,vdwglob,vlst,nvlst,radmin,epsilon,i12,n12)
!$acc&     present(ev,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)

c     ev  = 0

c
c     apply any reduction factor to the atomic coordinates
c
!$acc parallel loop default(present)
!$acc&         async(def_queue)
      do k = 1,nvdwbloc
         iglob   = ivdw (vdwglob(k))
         i       = loc  (iglob)
         iv      = ired (iglob)
         rdn     = kred (iglob)
         rdn1    = 1.0_ti_p - rdn
         xred(i) = rdn * x(iglob) + rdn1 * x(iv)
         yred(i) = rdn * y(iglob) + rdn1 * y(iv)
         zred(i) = rdn * z(iglob) + rdn1 * z(iv)
      enddo
      call resetForces_buff(def_queue)
c
c     set the coefficients for the switching function
c
      mode = 'VDW'
      call switch (mode)
      vdwshortcut2 = (vdwshortcut-shortheal)**2
c
c     find van der Waals energy and derivatives via neighbor list
c
!$acc parallel loop gang vector_length(32)
!$acc&         private(ai12)
#ifdef TINKER_DEBUG
!$acc&         present(interac)
#endif
!$acc&         async(def_queue)
      MAINLOOP:
     &do ii = 1, nvdwlocnl
         iivdw = vdwglobnl(ii)
         iglob = ivdw(iivdw)
         i     = loc(iglob)
         it    = jvdw(iglob)
         xi    = xred(i)
         yi    = yred(i)
         zi    = zred(i)
         muti  = mut(iglob)

         if (skipvdw12) then
            in12 = n12(iglob)
!$acc loop vector
            do j = 1,in12
               ai12(j) = i12(j,iglob)
            end do
         end if

!$acc loop vector 
         do k = 1, nvlst(ii)
            kglob  = vlst(k,ii)
            kbis   = loc (kglob)
            kt     = jvdw(kglob)
            mutik  = muti + mut (kglob)
            if (skipvdw12) then
               ik12 = .false.
!$acc loop seq
               do j = 1, in12
                  if (ai12(j).eq.kglob) ik12=.true.
               end do
               if (ik12) cycle
            end if
            !vscale  = 1.0_ti_p
            xpos   = xi - xred(kbis)
            ypos   = yi - yred(kbis)
            zpos   = zi - zred(kbis)
            call image_inl(xpos,ypos,zpos)
c
c     decide whether to compute the current interaction
c     and check for an interaction distance less than the cutoff
c
            rik2   = xpos**2 + ypos**2 + zpos**2
            if (rik2<vdwshortcut2.or.rik2>off2) cycle

            ! Annihilate
            if (vcouple.eq.1.and.mutik.eq.two1) mutik=one1

            rv2  =  radmin (kt,it)
            eps2 = epsilon (kt,it)

            call ehal1_couple_long(xpos,ypos,zpos,rik2,rv2,eps2
     &                 ,1.0_ti_p,cut2,cut,off,vdwshortcut
     &                 ,scexp,vlambda,scalpha,mutik
     &                 ,shortheal,ghal,dhal,e,dedx,dedy,dedz)

#ifdef TINKER_DEBUG
!$acc atomic
            interac(iglob) = interac(iglob) +1
#endif

            !Increment interaction energy
            ev    =   ev + tp2enr(e)

            !Increment the van der Waals derivatives
            devx = tp2mdr(dedx)
            devy = tp2mdr(dedy)
            devz = tp2mdr(dedz)
!$acc atomic
            dev(1,kbis) = dev(1,kbis) - devx
!$acc atomic
            dev(2,kbis) = dev(2,kbis) - devy
!$acc atomic
            dev(3,kbis) = dev(3,kbis) - devz

!$acc atomic
            dev(1,i) = dev(1,i) + devx
!$acc atomic
            dev(2,i) = dev(2,i) + devy
!$acc atomic
            dev(3,i) = dev(3,i) + devz

            !Increment the total van der Waals virial
            g_vxx = g_vxx + xpos * dedx
            g_vxy = g_vxy + ypos * dedx
            g_vxz = g_vxz + zpos * dedx
            g_vyy = g_vyy + ypos * dedy
            g_vyz = g_vyz + zpos * dedy
            g_vzz = g_vzz + zpos * dedz

         end do
      end do MAINLOOP


#ifdef TINKER_DEBUG
 34   format(2I10,3F12.4)
 35   format(A40,I16,3x,F16.6,I3)
!$acc wait
!$acc exit data copyout(interac)
!$acc update host(dev,ev)
      write(*,35)'nev & ev long & rank ',sum(interac),enr2en(ev),rank
c     if (rank.eq.0) then
c     do i =1,3000
c        print 34,i,interac(i)
c     end do
c     end if
#endif

      mode = 'LONG  '
      call ehal1c_correct_scaling_shortlong(xred,yred,zred,
     &            g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz,mode)

      call vdw_gradient_reduce

!$acc end data

      end subroutine

#ifdef _CUDA
      subroutine ehal1c_cu
      use atmlst    ,only: vdwglobnl,vdwglob
      use atoms     ,only: x,y,z,n
      use deriv     ,only: dev=>debond
      use domdec    ,only: loc,rank,nbloc,nproc
     &              ,xbegproc,xendproc,ybegproc,yendproc,zbegproc
     &              ,zendproc,glob
      use ehal1cu
      use energi    ,only: ev=>ev_r
      use inform    ,only: deb_Path
      use interfaces,only: ehal1c_correct_scaling
      use mutant    ,only: scalpha,scexp,vlambda,mut=>mutInt
      use neigh     ,only: cellv_glob,cellv_loc,cellv_jvdw
     &              ,vblst,ivblst
      use tinheader ,only: ti_p
      use shunt     ,only: c0,c1,c2,c3,c4,c5,off2,off,cut2,cut
      use utilcu    ,only: check_launch_kernel
      use utilgpu   ,only: def_queue,dir_queue,rec_queue,dir_stream
     &              ,rec_stream,rec_event,stream_wait_async
     &              ,warp_size,def_stream,inf
     &              ,ered_buff=>ered_buf1,vred_buff,reduce_energy_virial
     &              ,zero_evir_red_buffer,prmem_request
      use vdw       ,only: ired,kred,jvdw,ivdw,radmin_c
     &              ,epsilon_c,nvdwbloc,nvdwlocnl
     &              ,nvdwlocnlb,nvdwclass
     &              ,nvdwlocnlb_pair,nvdwlocnlb2_pair
      use vdwpot    ,only: vcorrect_ik,vcorrect_scale,n_vscale,dhal,ghal
      use vdw_locArray
      use virial    ,only: vir,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz
      implicit none
      integer i,k
      integer iglob,iivdw,iv,hal_Gs
      integer ierrSync,lst_start
#ifdef TINKER_DEBUG
      integer inter(n)
#endif
      real(t_p)  xbeg,xend,ybeg,yend,zbeg,zend
      real(r_p)  vxx,vxy,vxz
      real(r_p)  vyy,vyz,vzz
      real(t_p)  rdn,rdn1
      character*10 mode

      call prmem_request(xred    ,nvdwlocnlb,queue=def_queue)
      call prmem_request(yred    ,nvdwlocnlb,queue=def_queue)
      call prmem_request(zred    ,nvdwlocnlb,queue=def_queue)
      call prmem_request(xredc   ,nbloc     ,queue=def_queue)
      call prmem_request(yredc   ,nbloc     ,queue=def_queue)
      call prmem_request(zredc   ,nbloc     ,queue=def_queue)
      call prmem_request(loc_ired,nvdwlocnlb,queue=def_queue)
      call prmem_request(loc_kred,nvdwlocnlb,queue=def_queue)

c
      if(deb_Path) write (*,*) 'ehal1c_cu'
      def_stream = dir_stream
      xbeg = xbegproc(rank+1)
      xend = xendproc(rank+1)
      ybeg = ybegproc(rank+1)
      yend = yendproc(rank+1)
      zbeg = zbegproc(rank+1)
      zend = zendproc(rank+1)
      lst_start = 2*nvdwlocnlb_pair+1

#ifdef _OPENACC
      if (dir_queue.ne.rec_queue)
     &   call stream_wait_async(rec_stream,dir_stream,rec_event)
#endif

#ifdef TINKER_DEBUG
      inter(:) = 0
!$acc enter data copyin(inter)
#endif

c!$acc data present(xred,yred,zred,xredc,yredc,zredc,loc_ired,
c!$acc&  loc_kred)
c!$acc&     present(loc,ired,kred,x,y,z,ivdw,ivblst,vblst,
c!$acc&  cellv_glob,cellv_loc,cellv_jvdw,
c!$acc&  dev,vdwglob,radmin_c,epsilon_c)
c!$acc&     present(ev)
c!$acc&     async(def_queue)

c
c     zero out the van der Waals energy and virial
c     ! No need it has already been done at gradient entry
c!$acc serial async(def_queue)
c      ev  = 0
c      vxx = 0.0_ti_p
c      vxy = 0.0_ti_p
c      vxz = 0.0_ti_p
c      vyy = 0.0_ti_p
c      vyz = 0.0_ti_p
c      vzz = 0.0_ti_p
c!$acc end serial

c
c     apply any reduction factor to the atomic coordinates
c
!$acc parallel loop default(present) async(def_queue)
      do k = 1,nvdwlocnlb
         if (k.le.nvdwlocnl) then
            iglob    = cellv_glob(k)
            iv       = ired (iglob)
            rdn      = kred (iglob)
            rdn1     = 1.0_ti_p - rdn
            cellv_loc(k) = loc(iglob)
            loc_ired(k)  = loc(iv)
            if (iglob.eq.iv) then
               loc_kred(k) = rdn
            else
               loc_kred(k) = 1.0_ti_p
            end if
            xred(k)  = rdn * x(iglob) + rdn1 * x(iv)
            yred(k)  = rdn * y(iglob) + rdn1 * y(iv)
            zred(k)  = rdn * z(iglob) + rdn1 * z(iv)
         else
            ! Exclusion buffer to prevent interaction compute
            cellv_loc(k) = nbloc
            loc_ired(k)  = nbloc
            xred(k) = inf
            yred(k) = inf
            zred(k) = inf
         end if
      end do

!$acc parallel loop default(present) async(def_queue)
      do k = 1,nvdwbloc
         iglob    = ivdw(vdwglob(k))
         i        = loc  (iglob)
         iv       = ired (iglob)
         rdn      = kred (iglob)
         rdn1     = 1.0_ti_p - rdn
         xredc(i)  = rdn * x(iglob) + rdn1 * x(iv)
         yredc(i)  = rdn * y(iglob) + rdn1 * y(iv)
         zredc(i)  = rdn * z(iglob) + rdn1 * z(iv)
      end do

      call resetForces_buff(def_queue)
      call zero_evir_red_buffer(def_queue)
c
c     set the coefficients for the switching function
c
      !print*, nvdwlocnlb_pair
      mode = 'VDW'
      hal_Gs = nvdwlocnlb_pair/4
      call switch (mode)
c!$acc serial loop present(ivblst,vblst)
c      do i =1,32
c         print*,ivblst(872970),vblst(lst_start+872969*32+i-1)
c      end do

c
c     Call Vdw kernel in CUDA using C1 nblist
c
c!$acc host_data use_device(xred,yred,zred,cellv_glob,cellv_loc,dev
c!$acc&   ,loc_ired,vblst,cellv_jvdw,epsilon_c,radmin_c,ired
c!$acc&   ,kred,ered_buff,vred_buff,mut
c#ifdef TINKER_DEBUG
c!$acc&   ,inter
c#endif
c!$acc&   )
c
c      call ehal1_cu<<<hal_Gs,VDW_BLOCK_DIM,0,def_stream>>>
c     &             (xred,yred,zred,cellv_glob,cellv_loc,loc_ired
c     &             ,vblst,cellv_jvdw
c     &             ,epsilon_c,radmin_c,ired,kred,dev
c     &             ,ered_buff,vred_buff
c     &             ,nvdwlocnlb_pair,n,nbloc,nvdwlocnl,nvdwlocnlb
c     &             ,nvdwclass
c     &             ,c0,c1,c2,c3,c4,c5,cut2,cut,off2,off,ghal,dhal
c     &             ,scexp,vlambda,scalpha,mut
c#ifdef TINKER_DEBUG
c     &             ,inter
c#endif
c     &             )
c      call check_launch_kernel(" ehal1_cu1 ")
c
c!$acc end host_data


c
c     Call Vdw kernel in CUDA using C2 nblist
c
!$acc host_data use_device(xred,yred,zred,cellv_glob,cellv_loc
!$acc&    ,loc_ired,ivblst,vblst,cellv_jvdw,epsilon_c,mut
!$acc&    ,radmin_c,ired,kred,dev,ered_buff,vred_buff
#ifdef TINKER_DEBUG
!$acc&    ,inter
#endif
!$acc&    )

      call ehal1_cu2<<<hal_Gs,VDW_BLOCK_DIM,0,def_stream>>>
     &             (xred,yred,zred,cellv_glob,cellv_loc,loc_ired
     &             ,ivblst,vblst(lst_start),cellv_jvdw
     &             ,epsilon_c,radmin_c,ired,kred,dev
     &             ,ered_buff,vred_buff
     &             ,nvdwlocnlb2_pair,n,nbloc,nvdwlocnl,nvdwlocnlb
     &             ,nvdwclass
     &             ,c0,c1,c2,c3,c4,c5,cut2,cut,off2,off,ghal,dhal
     &             ,scexp,vlambda,scalpha,mut
     &             ,xbeg,xend,ybeg,yend,zbeg,zend
#ifdef TINKER_DEBUG
     &             ,inter,rank
#endif
     &             )
      call check_launch_kernel(" ehal1_cu2 ")

!$acc end host_data

      call reduce_energy_virial(ev,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz
     &                         ,ered_buff,def_queue)

#ifdef TINKER_DEBUG
 34   format(2I10,3F12.4)
 36   format(A30,2I10)
 35   format(A30,I16,3x,F16.6,I16)
!$acc wait
!$acc exit data copyout(inter)
!$acc update host(dev,ev)
      write(*,36)'nvdw pair block ',nvdwlocnlb_pair,nvdwlocnlb2_pair
      write(*,35)'nev & ev & rank ',sum(inter),enr2en(ev),rank
c     if (rank.eq.0) then
c      do i =1,3000
c         print 34,i,inter(i)
c      end do
c     end if
#endif

      call ehal1c_correct_scaling(xredc,yredc,zredc,
     &            g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)

      call vdw_gradient_reduce

c!$acc end data
      end subroutine

      subroutine ehalshortlong1c_cu
      use atmlst    ,only: vdwglobnl,vdwglob
      use atoms     ,only: x,y,z,n
      use cutoff    ,only: vdwshortcut,shortheal
      use deriv     ,only: dev=>debond
      use domdec    ,only: loc,rank,nbloc,nproc
     &              ,xbegproc,xendproc,ybegproc,yendproc,zbegproc
     &              ,zendproc,glob
      use ehal1cu
      use energi    ,only: ev=>ev_r
      use inform    ,only: deb_Path
      use interfaces,only: ehal1c_correct_scaling_shortlong
      use mutant    ,only: scalpha,scexp,vlambda,mut=>mutInt
      use neigh     ,only: cellv_glob,cellv_loc,cellv_jvdw
     &              ,vblst,ivblst,shortvblst,ishortvblst
      use potent    ,only: use_vdwshort,use_vdwlong
      use tinheader ,only: ti_p
      use shunt     ,only: c0,c1,c2,c3,c4,c5,off2,off,cut2,cut
      use utilcu    ,only: check_launch_kernel
      use utilgpu   ,only: def_queue,dir_queue,rec_queue,dir_stream
     &              ,rec_stream,rec_event,stream_wait_async
     &              ,warp_size,def_stream,inf
     &              ,ered_buff=>ered_buf1,vred_buff,reduce_energy_virial
     &              ,zero_evir_red_buffer,prmem_request
      use vdw       ,only: ired,kred,jvdw,ivdw,radmin_c
     &              ,epsilon_c,nvdwbloc,nvdwlocnl
     &              ,nvdwlocnlb,nvdwclass
     &              ,nvdwlocnlb_pair,nvdwlocnlb2_pair
     &              ,nshortvdwlocnlb2_pair
      use vdwpot    ,only: vcorrect_ik,vcorrect_scale,n_vscale,dhal,ghal
      use vdw_locArray
      use virial    ,only: vir,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz
      implicit none
      integer i,k
      integer iglob,iivdw,iv,hal_Gs
      integer ierrSync,lst_start
#ifdef TINKER_DEBUG
      integer inter(n)
#endif
      real(t_p)  vdwshortcut2
      real(t_p)  xbeg,xend,ybeg,yend,zbeg,zend
      real(t_p)  rdn,rdn1
      character*10 mode

      call prmem_request(xred    ,nvdwlocnlb,queue=def_queue)
      call prmem_request(yred    ,nvdwlocnlb,queue=def_queue)
      call prmem_request(zred    ,nvdwlocnlb,queue=def_queue)
      call prmem_request(xredc   ,nbloc     ,queue=def_queue)
      call prmem_request(yredc   ,nbloc     ,queue=def_queue)
      call prmem_request(zredc   ,nbloc     ,queue=def_queue)
      call prmem_request(loc_ired,nvdwlocnlb,queue=def_queue)
      call prmem_request(loc_kred,nvdwlocnlb,queue=def_queue)

c
      if(deb_Path) write (*,*) 'ehalshortlong1c_cu'
      def_stream = dir_stream
      xbeg = xbegproc(rank+1)
      xend = xendproc(rank+1)
      ybeg = ybegproc(rank+1)
      yend = yendproc(rank+1)
      zbeg = zbegproc(rank+1)
      zend = zendproc(rank+1)
      lst_start = 2*nvdwlocnlb_pair+1

#ifdef _OPENACC
      if (dir_queue.ne.rec_queue)
     &   call stream_wait_async(rec_stream,dir_stream,rec_event)
#endif


#ifdef TINKER_DEBUG
      inter(:) = 0
!$acc enter data copyin(inter)
#endif

c!$acc data present(xred,yred,zred,xredc,yredc,zredc,loc_ired,
c!$acc&  loc_kred)
c!$acc&     present(loc,ired,kred,x,y,z,ivdw,ivblst,vblst,
c!$acc&  cellv_glob,cellv_loc,cellv_jvdw,
c!$acc&  dev,vdwglob,radmin_c,epsilon_c)
c!$acc&     present(ev)
c!$acc&     async(def_queue)

c
c     apply any reduction factor to the atomic coordinates
c
!$acc parallel loop default(present) async(def_queue)
      do k = 1,nvdwlocnlb
         if (k.le.nvdwlocnl) then
            iglob    = cellv_glob(k)
            iv       = ired (iglob)
            rdn      = kred (iglob)
            rdn1     = 1.0_ti_p - rdn
            cellv_loc(k) = loc(iglob)
            loc_ired(k)  = loc(iv)
            if (iglob.eq.iv) then
               loc_kred(k) = rdn
            else
               loc_kred(k) = 1.0_ti_p
            end if
            xred(k)  = rdn * x(iglob) + rdn1 * x(iv)
            yred(k)  = rdn * y(iglob) + rdn1 * y(iv)
            zred(k)  = rdn * z(iglob) + rdn1 * z(iv)
         else
            ! Exclusion buffer to prevent interaction compute
            cellv_loc(k) = nbloc
            loc_ired(k)  = nbloc
            xred(k) = inf
            yred(k) = inf
            zred(k) = inf
         end if
      end do

!$acc parallel loop default(present) async(def_queue)
      do k = 1,nvdwbloc
         iglob     = ivdw(vdwglob(k))
         i         = loc  (iglob)
         iv        = ired (iglob)
         rdn       = kred (iglob)
         rdn1      = 1.0_ti_p - rdn
         xredc(i)  = rdn * x(iglob) + rdn1 * x(iv)
         yredc(i)  = rdn * y(iglob) + rdn1 * y(iv)
         zredc(i)  = rdn * z(iglob) + rdn1 * z(iv)
      end do

      call zero_evir_red_buffer(def_queue)
      call resetForces_buff(def_queue)
c
c     set the coefficients for the switching function
c
      if (use_vdwshort) then
         mode = 'SHORTVDW'
         call switch (mode)
         mode = 'SHORT '
      else if (use_vdwlong) then
         mode = 'VDW'
         call switch(mode)
         mode = 'LONG  '
      end if

      vdwshortcut2 = (vdwshortcut-shortheal)**2
c
c     Call Vdw kernel in CUDA using C2 nblist
c
!$acc host_data use_device(xred,yred,zred,cellv_glob,cellv_loc
!$acc&    ,loc_ired,ivblst,vblst,shortvblst,ishortvblst,cellv_jvdw
!$acc&    ,epsilon_c,radmin_c,mut,ired,kred,dev,ered_buff,vred_buff
#ifdef TINKER_DEBUG
!$acc&    ,inter
#endif
!$acc&    )

      if (use_vdwshort) then

      hal_Gs = nshortvdwlocnlb2_pair/4
      call ehal1short_cu <<<hal_Gs,VDW_BLOCK_DIM,0,def_stream>>>
     &           (xred,yred,zred,cellv_glob,cellv_loc,loc_ired
     &           ,ishortvblst,shortvblst(lst_start),cellv_jvdw
     &           ,epsilon_c,radmin_c
     &           ,ired,kred,dev,ered_buff,vred_buff
     &           ,nshortvdwlocnlb2_pair,n,nbloc,nvdwlocnl,nvdwlocnlb
     &           ,nvdwclass
     &           ,c0,c1,c2,c3,c4,c5,cut2,off2,off
     &           ,scexp,vlambda,scalpha,mut
     &           ,shortheal,ghal,dhal,use_vdwshort
     &           ,xbeg,xend,ybeg,yend,zbeg,zend
#ifdef TINKER_DEBUG
     &           ,inter,rank
#endif
     &           )
      call check_launch_kernel(" ehal1short_cu ")

      else if (use_vdwlong) then

      hal_Gs = nvdwlocnlb2_pair/4
      call ehal1long_cu <<<hal_Gs,VDW_BLOCK_DIM,0,def_stream>>>
     &           (xred,yred,zred,cellv_glob,cellv_loc,loc_ired
     &           ,ivblst,vblst(lst_start),cellv_jvdw,epsilon_c,radmin_c
     &           ,ired,kred,dev,ered_buff,vred_buff
     &           ,nvdwlocnlb2_pair,n,nbloc,nvdwlocnl,nvdwlocnlb
     &           ,nvdwclass
     &           ,c0,c1,c2,c3,c4,c5,cut2,cut,off2,off,vdwshortcut2
     &           ,scexp,vlambda,scalpha,mut
     &           ,vdwshortcut,shortheal,ghal,dhal,use_vdwshort
     &           ,xbeg,xend,ybeg,yend,zbeg,zend
#ifdef TINKER_DEBUG
     &           ,inter,rank
#endif
     &           )
      call check_launch_kernel(" ehal1long_cu ")
      end if

!$acc end host_data

      call reduce_energy_virial(ev,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz
     &                         ,ered_buff,def_queue)

#ifdef TINKER_DEBUG
 34   format(2I10,3F12.4)
 36   format(A30,2I10)
 35   format(A30,I16,3x,F16.6,I16)
!$acc wait
!$acc exit data copyout(inter)
!$acc update host(dev,ev)
      write(*,36)'nvdw pair block ',nvdwlocnlb_pair,nvdwlocnlb2_pair
      write(*,35)'nev & ev & rank ',sum(inter),enr2en(ev),rank
#endif

      call ehal1c_correct_scaling_shortlong(xredc,yredc,zredc,
     &               g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz,mode)

      call vdw_gradient_reduce

c!$acc end data
      end subroutine
#endif

      subroutine vdw_gradient_reduce
      use atoms  , only: n
      use atmlst , only: vdwglob
#ifdef USE_DETERMINISTIC_REDUCTION
      use ehal1gpu_inl ,only: mdr2md
#endif
      use vdw    , only: ivdw,ired,kred,nvdwbloc
      use domdec , only: loc,rank,nbloc
      use deriv  , debuff=>debond
      use inform , only: deb_Path
      use utilgpu, only: def_queue
      implicit none
      integer ii,j,i,il,iglob
      real(md_p) ded,de_red
      real(t_p) redi

12    format(2x,A)
      if (deb_Path) write(*,12) 'ehal1_gradient_reduce'

!$acc parallel loop collapse(2) default(present) async(def_queue)
      do ii = 1, nvdwbloc
         do j = 1, 3
            iglob  = ivdw(vdwglob(ii))
            i      = loc (iglob)
            ded    = mdr2md(debuff(j,i))
            ! Eject if not computed
            ! Necessary when (nbloc<n)
            if (ded.eq.0.0_re_p) cycle
            redi   = kred(iglob)
            il     = loc (ired(iglob))
            if ( i.ne.il ) then
               de_red = redi*ded
!$acc atomic
               dev(j,i)  = dev(j,i)  + de_red
               de_red    = real((1-redi),r_p)*ded
!$acc atomic
               dev(j,il) = dev(j,il) + de_red
            else
!$acc atomic
               dev(j,i) = dev(j,i) + ded
            end if
         end do
      end do

      end subroutine

      subroutine resetForces_buff(queue)
      use domdec , only: nbloc
      use deriv  , dbuff=>debond
      implicit none
      integer,intent(in)::queue
      integer i,j

!$acc parallel loop collapse(2) present(dbuff) async(queue)
      do i = 1,nbloc
         do j = 1, 3
            dbuff(j,i) = 0
         end do
      end do
      end subroutine

      subroutine ehal1c_correct_scaling(xred,yred,zred,
     &           vxx,vxy,vxz,vyy,vyz,vzz)

      use atmlst    ,only: vdwglobnl
      use deriv     ,only: dev=>debond
      use domdec    ,only: loc,rank
      use ehal1gpu_inl
      use energi    ,only: ev=>ev_r
      use inform    ,only: deb_Path
      use mutant    ,only: scexp,scalpha,vlambda,vcouple,mut=>mutInt
      use tinheader ,only: ti_p
      use shunt     ,only: c0,c1,c2,c3,c4,c5,off2,off,cut2,cut
      use vdw       ,only: ired,kred,jvdw,ivdw,radmin,radmin4,
     &                     epsilon,epsilon4
      use vdwpot    ,only: vcorrect_ik,vcorrect_scale,n_vscale,dhal,ghal
      use utilgpu   ,only: def_queue
      use virial
      implicit none
      integer i,j,k,kk,ksave
      integer kt,kglob,kbis,kvloc,kv,ki
      integer iglob,iivdw
      integer ii,iv,it,ivloc
      integer nnvlst,nnvlst2
      integer nn12,nn13,nn14,ntot
      integer interac
      integer(1) muti,mutik
      real(t_p)  xi,yi,zi,redi,e,de
      real(t_p)  rdn,rdn1,redk
      real(t_p)  invrik,rik,rik2,rik3,rik4,rik5,rik6,rik7
      real(t_p)  dedx,dedy,dedz
      mdyn_rtyp  devx,devy,devz
      real(t_p)  invrho,rv7orho
      real(t_p)  dtau,gtau,tau,tau7,rv7
      real(t_p)  rv2,eps2
      real(t_p)  xpos,ypos,zpos
      real(t_p)  dtaper,taper
      real(t_p)  vscale,vscale4
      logical    do_scale4
      character*10 mode

      real(t_p),intent(in):: xred(:)
      real(t_p),intent(in):: yred(:)
      real(t_p),intent(in):: zred(:)
      real(r_p)  vxx,vxy,vxz
      real(r_p)  vyy,vyz,vzz

      ! Scaling factor correction loop
      if (deb_Path) write(*,*) "ehal1c_correct_scaling"

!$acc parallel loop async(def_queue)
!$acc&     gang vector
!$acc&     present(xred,yred,zred,vxx,vxy,vxz,vyy,vyz,vzz)
!$acc&     present(loc,ivdw,loc,jvdw,vir,dev,radmin,mut,
!$acc&  radmin4,epsilon,epsilon4,vcorrect_ik,vcorrect_scale)
!$acc&     present(ev)
      do ii = 1,n_vscale
         iglob  = vcorrect_ik(ii,1)
         kglob  = vcorrect_ik(ii,2)
         vscale = vcorrect_scale(ii)
         i      = loc(iglob)
         kbis   = loc(kglob)

         it     = jvdw(iglob)
         kt     = jvdw(kglob)
         muti   = mut(iglob)
         mutik  = muti + mut(kglob)

         do_scale4 = .false.
         vscale4   = 0

         if (vscale.lt.0) then 
            vscale4 = -vscale
            vscale  = one
         end if
c
c     compute the energy contribution for this interaction
c
         xpos   = xred(i) - xred(kbis)
         ypos   = yred(i) - yred(kbis)
         zpos   = zred(i) - zred(kbis)
         call image_inl(xpos,ypos,zpos)
c
c     decide whether to compute the current interaction
c     and check for an interaction distance less than the cutoff
c
         rik2   = xpos**2 + ypos**2 + zpos**2
         if (rik2>off2) cycle
         ! Annihilate
         if (vcouple.eq.onei.and.mutik.eq.two1) mutik=one1
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

         call ehal1_couple(xpos,ypos,zpos,rik2,rv2,eps2,vscale
     &                    ,cut2,cut,off,ghal,dhal
     &                    ,scexp,vlambda,scalpha,mutik
     &                    ,e,dedx,dedy,dedz)

         if (.not.do_scale4) then
         e    = -e
         dedx = -dedx; dedy = -dedy; dedz = -dedz;
         end if

         ev           =   ev + tp2enr(e)
c        if (rank.eq.0.and.mod(ii,1).eq.0)
c    &      print*,iglob,kglob,vscale,enr2en(ev)

         devx = tp2mdr(dedx)
         devy = tp2mdr(dedy)
         devz = tp2mdr(dedz)
!$acc atomic
         dev(1,kbis)  = dev(1,kbis)  - devx
!$acc atomic
         dev(2,kbis)  = dev(2,kbis)  - devy
!$acc atomic
         dev(3,kbis)  = dev(3,kbis)  - devz

!$acc atomic
         dev(1,i)     = dev(1,i) + devx
!$acc atomic
         dev(2,i)     = dev(2,i) + devy
!$acc atomic
         dev(3,i)     = dev(3,i) + devz
c
c     increment the total van der Waals energy 
c
         vxx          = vxx + xpos * dedx
         vxy          = vxy + ypos * dedx
         vxz          = vxz + zpos * dedx
         vyy          = vyy + ypos * dedy
         vyz          = vyz + zpos * dedy
         vzz          = vzz + zpos * dedz

         ! deal with 1-4 Interactions
         if (vscale4.gt.0) then
            vscale    =  vscale4
            do_scale4 = .true.
            vscale4   = 0
            goto 20
         end if
      end do
      end

      subroutine ehal1c_correct_scaling_shortlong(xred,yred,zred,
     &           vxx,vxy,vxz,vyy,vyz,vzz,mode)

      use atmlst    ,only: vdwglobnl
      use cutoff    ,only: shortheal,vdwshortcut
      use deriv     ,only: dev=>debond
      use domdec    ,only: loc,rank
      use ehal1gpu_inl
      use energi    ,only: ev=>ev_r
      use inform    ,only: deb_Path
      use mutant    ,only: scalpha,scexp,vlambda,vcouple,mut=>mutInt
      use tinheader ,only: ti_p
      use shunt     ,only: c0,c1,c2,c3,c4,c5,off2,off,cut2,cut
      use vdw       ,only: ired,kred,jvdw,ivdw,radmin,radmin4,
     &                     epsilon,epsilon4
      use vdwpot    ,only: vcorrect_ik,vcorrect_scale,n_vscale,dhal,ghal
      use utilgpu   ,only: def_queue
      use virial
      implicit none
      integer i,j,k,kk,ksave
      integer kt,kglob,kbis,kvloc,kv,ki
      integer iglob,iivdw
      integer ii,iv,it,ivloc
      integer nnvlst,nnvlst2
      integer nn12,nn13,nn14,ntot
      integer interac
      real(t_p) xi,yi,zi,redi,e,de
      real(t_p) rdn,rdn1,redk
      real(t_p) invrik,rik,rik2,rik3,rik4,rik5,rik6,rik7
      real(t_p) dedx,dedy,dedz
      mdyn_rtyp devx,devy,devz
      real(t_p) invrho,rv7orho
      real(t_p) dtau,gtau,tau,tau7,rv7
      real(t_p) rv2,eps2,vdwshortcut2
      real(t_p) xpos,ypos,zpos
      real(t_p) dtaper,taper
      real(t_p) vscale,vscale4
      logical   do_scale4,short
      integer(1) muti,mutik

      real(t_p),intent(in):: xred(:)
      real(t_p),intent(in):: yred(:)
      real(t_p),intent(in):: zred(:)
      real(r_p)  vxx,vxy,vxz
      real(r_p)  vyy,vyz,vzz
      character*10,intent(in):: mode

      ! Scaling factor correction loop
      if      (mode.eq.'SHORT ') then
         vdwshortcut2 = 0.0_ti_p
         short = .true.
      else if (mode.eq.'LONG  ') then
         vdwshortcut2 = (vdwshortcut-shortheal)**2
         short = .false.
      else
         print*,'mode ',mode, ' is unrecognized for',
     &          'ehal1_correct_scaling_shortlong'
         call fatal
      end if

!$acc parallel loop async(def_queue)
!$acc&     gang vector
!$acc&     present(xred,yred,zred,vxx,vxy,vxz,vyy,vyz,vzz)
!$acc&     present(loc,ivdw,loc,jvdw,vir,dev,radmin,
!$acc&  radmin4,epsilon,epsilon4,vcorrect_ik,vcorrect_scale)
!$acc&     present(ev)
      do ii = 1,n_vscale
         iglob  = vcorrect_ik(ii,1)
         kglob  = vcorrect_ik(ii,2)
         vscale = vcorrect_scale(ii)
         i      = loc(iglob)
         kbis   = loc(kglob)

         it     = jvdw(iglob)
         kt     = jvdw(kglob)
         mutik  = mut(iglob) + mut(kglob)

         do_scale4 = .false.
         vscale4   = 0

         if ( vscale.lt.0 ) then 
            vscale4 = -vscale
            vscale  = 1
         end if
c
c     compute the energy contribution for this interaction
c
         xpos   = xred(i) - xred(kbis)
         ypos   = yred(i) - yred(kbis)
         zpos   = zred(i) - zred(kbis)
         call image_inl(xpos,ypos,zpos)
c
c     decide whether to compute the current interaction
c     and check for an interaction distance less than the cutoff
c
         rik2   = xpos**2 + ypos**2 + zpos**2
         if (rik2<vdwshortcut2.or.rik2>off2) cycle

         ! Annihilate
         if (vcouple.eq.onei.and.mutik.eq.two1) mutik=one1
 
         !replace 1-4 interactions
 20      continue
         if (do_scale4) then
            rv2  = radmin4 (kt,it)
            eps2 = epsilon4(kt,it)
         else
            rv2  =  radmin (kt,it)
            eps2 = epsilon (kt,it)
         end if

         ! Compute pairwaise (short/long)
         if (short) then
            call ehal1_couple_short(xpos,ypos,zpos,rik2,rv2,eps2,vscale
     &                       ,cut2,off
     &                       ,scexp,vlambda,scalpha,mutik
     &                       ,shortheal,ghal,dhal,e,dedx,dedy,dedz)
         else
            call ehal1_couple_long(xpos,ypos,zpos,rik2,rv2,eps2,vscale
     &                       ,cut2,cut,off,vdwshortcut
     &                       ,scexp,vlambda,scalpha,mutik
     &                       ,shortheal,ghal,dhal,e,dedx,dedy,dedz)
         end if

         if (.not.do_scale4) then
         e    = -e
         dedx = -dedx; dedy = -dedy; dedz = -dedz;
         end if

         ev           =   ev + tp2enr(e)

         devx         = tp2mdr(dedx)
         devy         = tp2mdr(dedy)
         devz         = tp2mdr(dedz)
!$acc atomic
         dev(1,kbis)  = dev(1,kbis)  - devx
!$acc atomic
         dev(2,kbis)  = dev(2,kbis)  - devy
!$acc atomic
         dev(3,kbis)  = dev(3,kbis)  - devz

!$acc atomic
         dev(1,i)     = dev(1,i) + devx
!$acc atomic
         dev(2,i)     = dev(2,i) + devy
!$acc atomic
         dev(3,i)     = dev(3,i) + devz
c
c     increment the total van der Waals energy 
c
         vxx          = vxx + real(xpos * dedx,r_p)
         vxy          = vxy + real(ypos * dedx,r_p)
         vxz          = vxz + real(zpos * dedx,r_p)
         vyy          = vyy + real(ypos * dedy,r_p)
         vyz          = vyz + real(zpos * dedy,r_p)
         vzz          = vzz + real(zpos * dedz,r_p)

         ! deal with 1-4 Interactions
         if (vscale4.gt.0) then
            vscale    =  vscale4
            do_scale4 = .true.
            vscale4   = 0
            goto 20
         end if
      end do
      end


      subroutine searchpair(nlst,lst,maxlst,int1,int2)
      use atmlst    ,only: vdwglobnl,vdwglob
      use atoms
      use deriv     ,only: dev=>debond
      use cell
      use domdec
      use ehal1gpu_inl
      use energi    ,only: ev
      use inform    ,only: deb_Path,dibuff
      use interfaces,only: ehal1c_correct_scaling
      use mutant    ,only: scalpha,scexp,vlambda,vcouple,mut=>mutInt
      use neigh 
      use tinheader ,only: ti_p,re_p
      use tinMemory ,only: prmem_request
      use shunt     ,only: c0,c1,c2,c3,c4,c5,off2,off,cut2,cut
      use vdw       ,only: ired,kred,jvdw,ivdw,radmin,
     &                     epsilon,nvdwbloc,nvdwlocnl,
     &                     nvdwclass
      use vdwpot    ,only: vcorrect_ik,vcorrect_scale,n_vscale,dhal,ghal
      use vdw_locArray
      use utilgpu   ,only: def_queue,dir_queue,rec_queue
#ifdef _OPENACC
     &                    ,dir_stream
     &                    ,rec_stream,rec_event,stream_wait_async
#endif
      use virial
      implicit none
      integer,intent(in)::maxlst,int1,int2
      integer nlst(nvdwlocnl)
      integer lst(maxlst,nvdwlocnl)

      integer i,j,k,kk,ksave
      integer kt,kglob,kbis,kvloc,kv,ki
      integer iglob,iivdw
      integer ii,iv,it,ivloc
      integer,save:: ncall=0
      integer nnvlst,nnvlst2
      integer nn12,nn13,nn14,ntot
      integer(1) muti,mutik
      integer xw,yw,zw,proc
      integer icell_len,icell,imageCell,distImage2Cell
      integer,parameter:: ncell2buff=2
      integer ncell2buffb
      real(t_p) xi,yi,zi,redi,e,de
      real(t_p) rdn,rdn1,redk
      real(t_p) invrik,rik,rik2,rik3,rik4,rik5,rik6,rik7
      real(t_p) dedx,dedy,dedz
      real(r_p) devx,devy,devz,devt
      real(t_p) invrho,rv7orho
      real(t_p) dtau,gtau,tau,tau7,rv7
      real(t_p) rv2,eps2
      real(t_p) xpos,ypos,zpos
      real(t_p) dtaper,taper
      real(t_p) vscale,vscale4
      real(t_p) xr,yr,zr,mbuf,vbuf,bigbuf
      real(t_p) lenx,leny,lenz
      real(t_p) lenx_cell,leny_cell,lenz_cell
      real(t_p) xmin,xmax,ymin,ymax,zmin,zmax
      logical   do_scale4
      character*10 mode
c
      write (*,*) 'searchpair',int1,int2
      ncall = ncall + 1

#ifdef _OPENACC
      if (dir_queue.ne.rec_queue) then
!!$acc wait(rec_queue) async(dir_queue)
         call stream_wait_async(rec_stream,dir_stream,rec_event)
      end if
#endif

      call prmem_request(xred,nbloc,queue=def_queue)
      call prmem_request(yred,nbloc,queue=def_queue)
      call prmem_request(zred,nbloc,queue=def_queue)
      if (ncall.eq.1) then
!$acc wait
      end if

!$acc data present(xred,yred,zred)
!$acc&     present(loc,ired,kred,x,y,z,vdwglobnl,ivdw,loc,jvdw,
!$acc&  vir,dev,vdwglob,lst,nlst,radmin,epsilon,mut)
!$acc&     present(g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)

c
c     set the coefficients for the switching function
c
      mode = 'SHORTVDW'
      call switch (mode)
c
c     apply any reduction factor to the atomic coordinates
c
!$acc parallel loop default(present)
!$acc&         async(def_queue)
      do k = 1,nvdwbloc
         iglob   = ivdw (vdwglob (k))
         i       = loc  (iglob)
         iv      = ired (iglob)
         rdn     = kred (iglob)
         rdn1    = 1.0_ti_p - rdn
         xred(i) = rdn * x(iglob) + rdn1 * x(iv)
         yred(i) = rdn * y(iglob) + rdn1 * y(iv)
         zred(i) = rdn * z(iglob) + rdn1 * z(iv)
      end do

      call resetForces_buff(def_queue)

!$acc parallel loop
      do ii = 1,nvdwlocnl
         iivdw = vdwglobnl(ii)
         iglob = ivdw(iivdw)
         if (iglob.eq.int1) print*, "found vdwglobnl",int1,ii,rank
         if (iglob.eq.int2) print*, "found vdwglobnl",int2,ii,rank
      end do

      ncell2buffb = ncell2buff
      xmin = xbegproc(rank+1)
      xmax = xendproc(rank+1)
      ymin = ybegproc(rank+1)
      ymax = yendproc(rank+1)
      zmin = zbegproc(rank+1)
      zmax = zendproc(rank+1)
      do i = 1, nbig_recep
        proc = pbig_recep(i)
        if (xbegproc(proc+1).le.xmin) xmin = xbegproc(proc+1)
        if (xendproc(proc+1).ge.xmax) xmax = xendproc(proc+1)
        if (ybegproc(proc+1).le.ymin) ymin = ybegproc(proc+1)
        if (yendproc(proc+1).ge.ymax) ymax = yendproc(proc+1)
        if (zbegproc(proc+1).le.zmin) zmin = zbegproc(proc+1)
        if (zendproc(proc+1).ge.zmax) zmax = zendproc(proc+1)
      end do

      if ((mbuf2).gt.vbuf2) then
        bigbuf = sqrt(mbuf2)/ncell2buffb
      else
        bigbuf = sqrt(vbuf2)/ncell2buffb
      end if
      lenx       = abs(xmax-xmin)
      nx_cell    = max(2*ncell2buffb+1,int(lenx/(bigbuf)))
      lenx_cell  = lenx/nx_cell
      leny       = abs(ymax-ymin)
      ny_cell    = max(2*ncell2buffb+1,int(leny/(bigbuf)))
      leny_cell  = leny/ny_cell
      lenz       = abs(zmax-zmin)
      nz_cell    = max(2*ncell2buffb+1,int(lenz/(bigbuf)))
      lenz_cell  = lenz/nz_cell

!$acc parallel loop async
      do i = 1,ncell_tot
         cell_len(i) = 0
      end do
c!$acc parallel loop async
c      do i = 1,n
c         indcelltemp(i) = 0
c      end do
       if (rank.eq.0) then
          print*,'nx',nx_cell,ny_cell,nz_cell
          print*,'lx',lenx,leny,lenz,bigbuf
       end if
!$acc serial async present(xcell2,ycell2) copyin(zcell2)
       print*,'cell',xcell2,ycell2,zcell2
!$acc end serial

!$acc parallel loop async
      do i = 1,nlocnl
         iglob = ineignl(i)
         xr    = x(iglob)
         yr    = y(iglob)
         zr    = z(iglob)
         call image_inl(xr,yr,zr)
         if ((xcell2-xr).lt.eps_cell) xr = xr-0.05*lenx_cell
         if ((ycell2-yr).lt.eps_cell) yr = yr-0.05*leny_cell
         if ((zcell2-zr).lt.eps_cell) zr = zr-0.05*lenz_cell
         xw  = int((xr-xmin)/lenx_cell)
         yw  = int((yr-ymin)/leny_cell)
         zw  = int((zr-zmin)/lenz_cell)
         if (iglob.eq.int1) print*,'box1',rank,xw,yw,zw
         if (iglob.eq.int1) print*,'box1',rank,x(iglob),xr
     &                     ,y(iglob),yr,z(iglob),zr
         if (iglob.eq.int2) print*,'box2',rank,xw,yw,zw
         if (iglob.eq.int2) print*,'box2',rank,x(iglob),xr
     &                     ,y(iglob),yr,z(iglob),zr
         icell = (xw + nx_cell*yw + nx_cell*ny_cell*zw) + 1
         repartcell (iglob) = icell
!$acc atomic
         cell_len (icell) = cell_len(icell)+1
c         icell_len        = cell_len(icell)
c!$acc end atomic
      end do
c
c     find van der Waals energy and derivatives via neighbor list
c
!$acc parallel loop num_gangs(2) vector_length(32)
!$acc&         async(def_queue)
      MAINLOOP:
     &do ii = 1, nvdwlocnl
         iivdw = vdwglobnl(ii)
         iglob = ivdw(iivdw)
         i     = loc(iglob)
         it    = jvdw(iglob)
         xi    = xred(i)
         yi    = yred(i)
         zi    = zred(i)
         muti  = mut(iglob)

!$acc loop vector 
         do k = 1, nlst(ii)
            kglob  = lst(k,ii)
            kbis   = loc (kglob)
            kt     = jvdw (kglob)
            mutik  = muti + mut(kglob)
            !vscale = 1.0_ti_p
c
c     compute the energy contribution for this interaction
c
            xpos   = xi - xred(kbis)
            ypos   = yi - yred(kbis)
            zpos   = zi - zred(kbis)
            call image_inl(xpos,ypos,zpos)
c
c     decide whether to compute the current interaction
c     and check for an interaction distance less than the cutoff
c
            rik2   = xpos**2 + ypos**2 + zpos**2
         if (iglob.eq.int1.and.kglob.lt.int1+10) then
             print*,"comp0",kglob,rank,nlst(ii),rik2
             print*,"comp0",kglob,rank,xpos,ypos,zpos
         end if
            if (rik2>off2) cycle
            ! Annihilate
            if (vcouple.eq.1.and.mutik.eq.two1) mutik=one1

            rv2  =  radmin (kt,it)
            eps2 = epsilon (kt,it)

            call ehal1_couple(xpos,ypos,zpos,rik2,rv2,eps2,1.0_ti_p
     &                       ,cut2,cut,off,ghal,dhal
     &                       ,scexp,vlambda,scalpha,mutik
     &                       ,e,dedx,dedy,dedz)

         if (iglob.eq.int1.and.kglob.lt.int1+10) then
             print*,"comp",kglob,rank,dedx,dedy,dedz
         end if

         end do
      end do MAINLOOP

      ! Save positions before removing comments on this kernel
c!$acc parallel loop async
c!$acc&         present(dibuff,xold_nl,yold_nl,zold_nl)
c      do i = 1,nbloc
c         iglob    = dibuff(i)
c         x(iglob) = xold_nl(iglob)
c         y(iglob) = yold_nl(iglob)
c         z(iglob) = zold_nl(iglob)
c         loc(iglob) = i
c      end do
      call reassignrespa(2,2)
      call kvdw(.false.,-1)
      call nblist(0)

      mode = 'SHORTVDW'
      call switch (mode)
c
c     apply any reduction factor to the atomic coordinates
c
!$acc parallel loop default(present)
!$acc&         async(def_queue)
      do k = 1,nvdwbloc
         iglob   = ivdw (vdwglob (k))
         i       = loc  (iglob)
         iv      = ired (iglob)
         rdn     = kred (iglob)
         rdn1    = 1.0_ti_p - rdn
         xred(i) = rdn * x(iglob) + rdn1 * x(iv)
         yred(i) = rdn * y(iglob) + rdn1 * y(iv)
         zred(i) = rdn * z(iglob) + rdn1 * z(iv)
      end do
c
c     find van der Waals energy and derivatives via neighbor list
c
!$acc parallel loop num_gangs(2) vector_length(32)
!$acc&         async(def_queue)
      do ii = 1, nvdwlocnl
         iivdw = vdwglobnl(ii)
         iglob = ivdw(iivdw)
         i     = loc(iglob)
         it    = jvdw(iglob)
         xi    = xred(i)
         yi    = yred(i)
         zi    = zred(i)
         muti  = mut(iglob)

!$acc loop vector 
         do k = 1, nlst(ii)
            kglob  = lst(k,ii)
            kbis   = loc (kglob)
            kt     = jvdw (kglob)
            mutik  = muti + mut(kglob)
            !vscale = 1.0_ti_p
c
c     compute the energy contribution for this interaction
c
            xpos   = xi - xred(kbis)
            ypos   = yi - yred(kbis)
            zpos   = zi - zred(kbis)
            call image_inl(xpos,ypos,zpos)
c
c     decide whether to compute the current interaction
c     and check for an interaction distance less than the cutoff
c
            rik2   = xpos**2 + ypos**2 + zpos**2
         if (iglob.eq.int1.and.kglob.le.int1+10) then
             print*,"comb0",kglob,rank,nlst(ii),rik2
             print*,"comp0",kglob,rank,xpos,ypos,zpos
         end if
            if (rik2>off2) cycle
            ! Annihilate
            if (vcouple.eq.1.and.mutik.eq.two1) mutik=one1

            rv2  =  radmin (kt,it)
            eps2 = epsilon (kt,it)

            call ehal1_couple(xpos,ypos,zpos,rik2,rv2,eps2,1.0_ti_p
     &                       ,cut2,cut,off,ghal,dhal
     &                       ,scexp,vlambda,scalpha,mutik
     &                       ,e,dedx,dedy,dedz)

         if (iglob.eq.int1.and.kglob.le.int1+10) then
             print*,"comb",kglob,rank,dedx,dedy,dedz
         end if
         end do
      end do

      call searchVdwScaled(xred,yred,zred,
     &     g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz,int1,int2)

c      call vdw_gradient_reduce
!$acc end data
!$acc wait
      end subroutine

      subroutine searchVdwScaled(xred,yred,zred,
     &           vxx,vxy,vxz,vyy,vyz,vzz,int1,int2)

      use atmlst    ,only: vdwglobnl
      use deriv     ,only: dev=>debond
      use domdec    ,only: loc,rank,nbloc
      use ehal1gpu_inl
      use energi    ,only: ev
      use inform    ,only: deb_Path
      use mutant    ,only: scexp,scalpha,vlambda,vcouple,mut=>mutInt
      use tinheader ,only: ti_p
      use shunt     ,only: c0,c1,c2,c3,c4,c5,off2,off,cut2,cut
      use vdw       ,only: ired,kred,jvdw,ivdw,radmin,radmin4,
     &                     epsilon,epsilon4,nvdwbloc
      use vdwpot    ,only: vcorrect_ik,vcorrect_scale,n_vscale,dhal,ghal
      use utilgpu   ,only: def_queue
      use virial
      implicit none
      integer i,j,k,kk,ksave
      integer kt,kglob,kbis,kvloc,kv,ki
      integer iglob,iivdw
      integer ii,iv,it,ivloc
      integer nnvlst,nnvlst2
      integer nn12,nn13,nn14,ntot
      integer interac
      integer(1) muti,mutik
      real(t_p)  xi,yi,zi,redi,e,de
      real(t_p)  rdn,rdn1,redk
      real(t_p)  invrik,rik,rik2,rik3,rik4,rik5,rik6,rik7
      real(t_p)  dedx,dedy,dedz
      real(r_p)  devx,devy,devz,devt
      real(t_p)  invrho,rv7orho
      real(t_p)  dtau,gtau,tau,tau7,rv7
      real(t_p)  rv2,eps2
      real(t_p)  xpos,ypos,zpos
      real(t_p)  dtaper,taper
      real(t_p)  vscale,vscale4
      logical    do_scale4
      character*10 mode

      integer int1,int2
      real(t_p),intent(in):: xred(nbloc)
      real(t_p),intent(in):: yred(nbloc)
      real(t_p),intent(in):: zred(nbloc)
      real(r_p)  vxx,vxy,vxz
      real(r_p)  vyy,vyz,vzz

      ! Scaling factor correction loop
      write(*,*) "searchVdwScaled",n_vscale

!$acc parallel loop async(def_queue)
!$acc&     gang vector
!$acc&     present(xred,yred,zred,vxx,vxy,vxz,vyy,vyz,vzz)
!$acc&     present(loc,ivdw,loc,jvdw,vir,dev,radmin,mut,
!$acc&  radmin4,epsilon,epsilon4,vcorrect_ik,vcorrect_scale)
      do ii = 1,n_vscale
         iglob  = vcorrect_ik(ii,1)
         kglob  = vcorrect_ik(ii,2)
         vscale = vcorrect_scale(ii)
         i      = loc(iglob)
         kbis   = loc(kglob)

         it     = jvdw(iglob)
         kt     = jvdw(kglob)
         muti   = mut(iglob)
         mutik  = muti + mut(kglob)

         do_scale4 = .false.
         vscale4   = 0

         if (vscale.lt.0) then 
            vscale4 = -vscale
            vscale  = 1
         end if
c
c     compute the energy contribution for this interaction
c
         xpos   = xred(i) - xred(kbis)
         ypos   = yred(i) - yred(kbis)
         zpos   = zred(i) - zred(kbis)
         call image_inl(xpos,ypos,zpos)
c
c     decide whether to compute the current interaction
c     and check for an interaction distance less than the cutoff
c
         rik2   = xpos**2 + ypos**2 + zpos**2
         if (rik2>off2) cycle

         if ((iglob.eq.int1.and.kglob.eq.int2).or.
     &       (kglob.eq.int1.and.iglob.eq.int2)) then
             print*,"corr",rank,xpos,ypos,zpos,e
         end if
         ! Annihilate
         if (vcouple.eq.1.and.mutik.eq.two1) mutik=one1
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

         call ehal1_couple(xpos,ypos,zpos,rik2,rv2,eps2,vscale
     &                    ,cut2,cut,off,ghal,dhal
     &                    ,scexp,vlambda,scalpha,mutik
     &                    ,e,dedx,dedy,dedz)

         if (.not.do_scale4) then
         e    = -e
         dedx = -dedx; dedy = -dedy; dedz = -dedz;
         end if

         if ((iglob.eq.int1.and.kglob.eq.int2).or.
     &       (kglob.eq.int1.and.iglob.eq.int2)) then
             print*,"corr",rank,iglob,kglob,rik2,dedx,dedy,dedz
         end if

         ! deal with 1-4 Interactions
         if (vscale4.gt.0) then
            vscale    =  vscale4
            do_scale4 = .true.
            vscale4   = 0
            goto 20
         end if
      end do
      end

