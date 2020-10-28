c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine estrbnd1   --  stretch-bend energy and derivs  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "estrbnd1" calculates the stretch-bend potential energy and
c     first derivatives with respect to Cartesian coordinates
c
c
#include "tinker_precision.h"
      module estrbnd1gpu_inl
        contains
#include "image.f.inc"
      end module

      subroutine estrbnd1gpu
      use angle
      use angpot
      use atmlst
      use atoms
      use bond
      use bound
      use deriv
      use domdec
      use energi
      use estrbnd1gpu_inl
      use group
      use math
      use nvshmem
      use potent,only:use_amd_wat1
      use strbnd
      use tinheader ,only:ti_p,re_p
      use usage
      use virial
      use mamd
      implicit none
      integer i,j,k,istrbnd,istrbndloc
      integer ia,ib,ic
      integer ialoc,ibloc,icloc
      integer ipe,ind
      real(t_p) e,dr1,dr2,dt
      real(t_p) angle1
      real(t_p) force1,force2
      real(t_p) dot,cosine
      real(t_p) xia,yia,zia
      real(t_p) xib,yib,zib
      real(t_p) xic,yic,zic
      real(t_p) xab,yab,zab
      real(t_p) xcb,ycb,zcb
      real(t_p) rab,rab2
      real(t_p) rcb,rcb2
      real(t_p) xp,yp,zp,rp
      real(t_p) term1,term2
      real(t_p) termr,term1t,term2t
      real(t_p) ddtdxia,ddtdyia,ddtdzia
      real(t_p) ddtdxic,ddtdyic,ddtdzic
      real(t_p) ddrdxia,ddrdyia,ddrdzia
      real(t_p) ddrdxic,ddrdyic,ddrdzic
      real(r_p) dedxia,dedyia,dedzia
      real(r_p) dedxib,dedyib,dedzib
      real(r_p) dedxic,dedyic,dedzic
      real(t_p) vxx,vyy,vzz
      real(t_p) vyx,vzx,vzy
      logical proceed

      if(rank.eq.0.and.tinkerdebug) write(*,*) 'estrbnd1gpu'
c
c     calculate the stretch-bend energy and first derivatives
c
!$acc parallel loop present(eba,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz,
!$acc&         eW1aMD)
!$acc&         default(present) async
      do istrbndloc = 1, nstrbndloc
         istrbnd = strbndglob(istrbndloc)
         i      = isb(1,istrbnd)
#ifdef USE_NVSHMEM_CUDA
         ipe     =     (i-1)/nangle_pe
         ind     = mod((i-1),nangle_pe) +1
         ia      = d_iang(ipe)%pel(1,ind)
         ib      = d_iang(ipe)%pel(2,ind)
         ic      = d_iang(ipe)%pel(3,ind)
#else
         ia     = iang(1,i)
         ib     = iang(2,i)
         ic     = iang(3,i)
#endif
         ialoc  = loc(ia)
         ibloc  = loc(ib)
         icloc  = loc(ic)
         force1 = sbk(1,istrbnd)
         force2 = sbk(2,istrbnd)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (proceed)  proceed = (use(ia) .or. use(ib) .or. use(ic))
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
c     compute the value of the bond angle
c
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
               rab = sqrt(rab2)
               rcb = sqrt(rcb2)
               xp = ycb*zab - zcb*yab
               yp = zcb*xab - xcb*zab
               zp = xcb*yab - ycb*xab
               rp = sqrt(xp*xp + yp*yp + zp*zp)
               rp = max(rp,0.001_ti_p)
               dot = xab*xcb + yab*ycb + zab*zcb
               cosine = dot / (rab*rcb)
               cosine = min(1.0_ti_p,max(-1.0_ti_p,cosine))
               angle1 = radian * acos(cosine)
c
c     find chain rule terms for the bond angle deviation
c
               dt = angle1 - anat(i)
               term1 = -radian / (rab2*rp)
               term2 = radian / (rcb2*rp)
               ddtdxia = term1 * (yab*zp-zab*yp)
               ddtdyia = term1 * (zab*xp-xab*zp)
               ddtdzia = term1 * (xab*yp-yab*xp)
               ddtdxic = term2 * (ycb*zp-zcb*yp)
               ddtdyic = term2 * (zcb*xp-xcb*zp)
               ddtdzic = term2 * (xcb*yp-ycb*xp)
c
c     find chain rule terms for the bond length deviations
c
               j = isb(2,istrbnd)
               k = isb(3,istrbnd)
#ifdef USE_NVSHMEM_CUDA
               ipe = (j-1)/nbond_pe
               ind = mod((j-1),nbond_pe) +1
               dr1 = rab - d_bl(ipe)%pel(ind)
               ipe = (k-1)/nbond_pe
               ind = mod((k-1),nbond_pe) +1
               dr2 = rcb - d_bl(ipe)%pel(ind)
#else
               dr1 = rab - bl(j)
               dr2 = rcb - bl(k)
#endif
               term1 = 1.0_ti_p / rab
               term2 = 1.0_ti_p / rcb
               ddrdxia = term1 * xab
               ddrdyia = term1 * yab
               ddrdzia = term1 * zab
               ddrdxic = term2 * xcb
               ddrdyic = term2 * ycb
               ddrdzic = term2 * zcb
c
c     abbreviations used in defining chain rule terms
c
               term1 = stbnunit * force1
               term2 = stbnunit * force2
               termr = term1*dr1 + term2*dr2
               term1t = term1 * dt
               term2t = term2 * dt
c
c     get the energy and master chain rule terms for derivatives
c
               e = termr * dt
               dedxia = term1t*ddrdxia + termr*ddtdxia
               dedyia = term1t*ddrdyia + termr*ddtdyia
               dedzia = term1t*ddrdzia + termr*ddtdzia
               dedxic = term2t*ddrdxic + termr*ddtdxic
               dedyic = term2t*ddrdyic + termr*ddtdyic
               dedzic = term2t*ddrdzic + termr*ddtdzic
               dedxib = -dedxia - dedxic
               dedyib = -dedyia - dedyic
               dedzib = -dedzia - dedzic
c
c     increment the total stretch-bend energy and derivatives
c
               eba = eba + e
!$acc atomic update
               deba(1,ibloc) = deba(1,ibloc) + dedxib
!$acc atomic update
               deba(2,ibloc) = deba(2,ibloc) + dedyib
!$acc atomic update
               deba(3,ibloc) = deba(3,ibloc) + dedzib

!$acc atomic update
               deba(1,ialoc) = deba(1,ialoc) + dedxia
!$acc atomic update
               deba(2,ialoc) = deba(2,ialoc) + dedyia
!$acc atomic update
               deba(3,ialoc) = deba(3,ialoc) + dedzia
!$acc atomic update

               deba(1,icloc) = deba(1,icloc) + dedxic
!$acc atomic update
               deba(2,icloc) = deba(2,icloc) + dedyic
!$acc atomic update
               deba(3,icloc) = deba(3,icloc) + dedzic
c
c     aMD storage if waters are considered
c
               if (use_amd_wat1) then
               if (type(ia) == aMDwattype(1) .or. type(ib)
     $         == aMDwattype(1) .or. type(ic) == aMDwattype(1)) then
                  eW1aMD = eW1aMD + e
!$acc atomic
                  deW1aMD(1,ialoc) = deW1aMD(1,ialoc) + dedxia
!$acc atomic
                  deW1aMD(2,ialoc) = deW1aMD(2,ialoc) + dedyia
!$acc atomic
                  deW1aMD(3,ialoc) = deW1aMD(3,ialoc) + dedzia
!$acc atomic
                  deW1aMD(1,ibloc) = deW1aMD(1,ibloc) + dedxib
!$acc atomic
                  deW1aMD(2,ibloc) = deW1aMD(2,ibloc) + dedyib
!$acc atomic
                  deW1aMD(3,ibloc) = deW1aMD(3,ibloc) + dedzib
!$acc atomic
                  deW1aMD(1,icloc) = deW1aMD(1,icloc) + dedxic
!$acc atomic
                  deW1aMD(2,icloc) = deW1aMD(2,icloc) + dedyic
!$acc atomic
                  deW1aMD(3,icloc) = deW1aMD(3,icloc) + dedzic
               end if
               end if
c
c     increment the internal virial tensor components
c
             g_vxx = g_vxx + xab*real(dedxia,t_p) + xcb*real(dedxic,t_p)
             g_vxy = g_vxy + yab*real(dedxia,t_p) + ycb*real(dedxic,t_p)
             g_vxz = g_vxz + zab*real(dedxia,t_p) + zcb*real(dedxic,t_p)
             g_vyy = g_vyy + yab*real(dedyia,t_p) + ycb*real(dedyic,t_p)
             g_vyz = g_vyz + zab*real(dedyia,t_p) + zcb*real(dedyic,t_p)
             g_vzz = g_vzz + zab*real(dedzia,t_p) + zcb*real(dedzic,t_p)

            end if
         end if
      end do

      end
