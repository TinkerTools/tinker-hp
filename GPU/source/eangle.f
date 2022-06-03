c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine eangle  --  angle bending potential energy  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "eangle" calculates the angle bending potential energy;
c     projected in-plane angles at trigonal centers, special
c     linear or Fourier angle bending terms are optionally used
c
c
#include "tinker_precision.h"
      module eangle_inl
        contains
#include "image.f.inc"
      end module

      subroutine eangle
      use angle
      use angpot
      use atmlst
      use atoms
      use bound
      use energi
      use eangle_inl
      use group
      use math
      use tinheader
      use timestat,only:timer_enter,timer_exit,timer_eangle
      use usage
      implicit none
      integer i,ia,ib,ic,id,iangle
#ifdef USE_NVSHMEM_CUDA
      integer ipe,ind
#endif
      real(t_p) e,ideal,force
      real(t_p) fold,factor
      real(t_p) dot,cosine
      real(t_p) angle1
      real(t_p) dt,dt2,dt3,dt4
      real(t_p) xia,yia,zia
      real(t_p) xib,yib,zib
      real(t_p) xic,yic,zic
      real(t_p) xid,yid,zid
      real(t_p) xab,yab,zab
      real(t_p) xcb,ycb,zcb
      real(t_p) xad,yad,zad
      real(t_p) xbd,ybd,zbd
      real(t_p) xcd,ycd,zcd
      real(t_p) xip,yip,zip
      real(t_p) xap,yap,zap
      real(t_p) xcp,ycp,zcp
      real(t_p) rab2,rcb2
      real(t_p) rap2,rcp2
      real(t_p) xt,yt,zt
      real(t_p) rt2,delta
      logical proceed
c
c
c     zero out the angle bending energy component
c
      call timer_enter( timer_eangle )
      ea = 0.0_re_p
c
c     calculate the bond angle bending energy term
c
!$acc parallel loop
#ifdef USE_NVSHMEM_CUDA
!$acc&         present(x,y,z,use,iang,angleglob
#else
!$acc&         present(x,y,z,use,iang,angleglob
#endif
!$acc&    ,anat,ak,afld,angtypI) present(ea) async
!$acc&         reduction(+:ea)
      do iangle = 1, nangleloc
         i     = angleglob(iangle)
#ifdef USE_NVSHMEM_CUDA
         ipe   =     (i-1)/nangle_pe
         ind   = mod((i-1),nangle_pe) +1
         ia    = d_iang(ipe)%pel(1,ind)
         ib    = d_iang(ipe)%pel(2,ind)
         ic    = d_iang(ipe)%pel(3,ind)
#else
         ia    = iang(1,i)
         ib    = iang(2,i)
         ic    = iang(3,i)
#endif
         ideal = anat(i)
         force = ak(i)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (angtypI(i) .eq. ANG_IN_PLANE) then
            id = iang(4,i)
            if (proceed)  proceed = (use(ia) .or. use(ib) .or.
     &                                 use(ic) .or. use(id))
         else
            if (proceed)  proceed = (use(ia) .or. use(ib) .or. use(ic))
         end if
c
c     get the coordinates of the atoms in the angle
c
         if (proceed) then
            xia = x(ia)
            yia = y(ia)
            zia = z(ia)
            xib = x(ib)
            yib = y(ib)
            zib = z(ib)
            xic = x(ic)
            yic = y(ic)
            zic = z(ic)
c
c     compute the bond angle bending energy
c
            if (angtypI(i) .ne. ANG_IN_PLANE) then
               xab = xia - xib
               yab = yia - yib
               zab = zia - zib
               xcb = xic - xib
               ycb = yic - yib
               zcb = zic - zib
               if (use_polymer) then
                  call image_inl (xab,yab,zab)
                  call image_inl (xcb,ycb,zcb)
               end if
               rab2 = xab*xab + yab*yab + zab*zab
               rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
               if (rab2.ne.0.0_ti_p .and. rcb2.ne.0.0_ti_p) then
                  dot = xab*xcb + yab*ycb + zab*zcb
                  cosine = dot / sqrt(rab2*rcb2)
                  cosine = min(1.0_ti_p,max(-1.0_ti_p,cosine))
                  angle1 = radian * acos(cosine)
                  if (angtypI(i) .eq. ANG_HARMONIC) then
                     dt = angle1 - ideal
                     dt2 = dt * dt
                     dt3 = dt2 * dt
                     dt4 = dt2 * dt2
                     e = angunit * force * dt2
     &                   * (1.0_ti_p+cang*dt+qang*dt2+pang*dt3+sang*dt4)
                  else if (angtypI(i) .eq. ANG_LINEAR) then
                     factor = 2.0_ti_p * angunit * radian**2
                     e = factor * force * (1.0_ti_p+cosine)
                  else if (angtypI(i) .eq. ANG_FOURIER) then
                     fold = afld(i)
                     factor = 2.0_ti_p * angunit * (radian/fold)**2
                     cosine = cos((fold*angle1-ideal)/radian)
                     e = factor * force * (1.0_ti_p+cosine)
                  end if
c
c     increment the total bond angle bending energy
c
                  ea = ea + e
               end if
c
c     compute the projected in-plane angle bend energy
c
            else
               xid = x(id)
               yid = y(id)
               zid = z(id)
               xad = xia - xid
               yad = yia - yid
               zad = zia - zid
               xbd = xib - xid
               ybd = yib - yid
               zbd = zib - zid
               xcd = xic - xid
               ycd = yic - yid
               zcd = zic - zid
               if (use_polymer) then
                  call image_inl (xad,yad,zad)
                  call image_inl (xbd,ybd,zbd)
                  call image_inl (xcd,ycd,zcd)
               end if
               xt = yad*zcd - zad*ycd
               yt = zad*xcd - xad*zcd
               zt = xad*ycd - yad*xcd
               rt2 = xt*xt + yt*yt + zt*zt
               delta = -(xt*xbd + yt*ybd + zt*zbd) / rt2
               xip = xib + xt*delta
               yip = yib + yt*delta
               zip = zib + zt*delta
               xap = xia - xip
               yap = yia - yip
               zap = zia - zip
               xcp = xic - xip
               ycp = yic - yip
               zcp = zic - zip
               if (use_polymer) then
                  call image_inl (xap,yap,zap)
                  call image_inl (xcp,ycp,zcp)
               end if
               rap2 = xap*xap + yap*yap + zap*zap
               rcp2 = xcp*xcp + ycp*ycp + zcp*zcp
               if (rap2.ne.0.0_ti_p .and. rcp2.ne.0.0_ti_p) then
                  dot = xap*xcp + yap*ycp + zap*zcp
                  cosine = dot / sqrt(rap2*rcp2)
                  cosine = min(1.0_ti_p,max(-1.0_ti_p,cosine))
                  angle1 = radian * acos(cosine)
                  dt = angle1 - ideal
                  dt2 = dt * dt
                  dt3 = dt2 * dt
                  dt4 = dt2 * dt2
                  e = angunit * force * dt2
     &                   * (1.0_ti_p+cang*dt+qang*dt2+pang*dt3+sang*dt4)
c
c     increment the total bond angle bending energy
c
                  ea = ea + e
               end if
            end if
         end if
      end do
      call timer_exit ( timer_eangle )
      end
