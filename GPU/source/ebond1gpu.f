c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine ebond1  --  bond stretch energy & derivatives  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ebond1" calculates the bond stretching energy and
c     first derivatives with respect to Cartesian coordinates
c
c
#include "tinker_precision.h"
      module ebond1gpu_inl
        contains
#include "image.f.inc"
      end module

      subroutine ebond1gpu
      use atmlst
      use atoms
      use bndpot
      use bond
      use bound
      use deriv
      use domdec
      use energi
      use ebond1gpu_inl
      use group
      use inform   ,only: deb_Path
      use nvshmem
      use usage
      use virial
      use timestat ,only: timer_enter,timer_exit,timer_ebond
     &             ,quiet_timers
      use tinheader
      use mamd
      use potent,only:use_amd_wat1
      implicit none
      integer i,ia,ib,ialoc,ibloc
      integer iaglob,ibglob
      integer ibond
      integer ia1,ib1,ind,ipe
      real(t_p) e,de,ideal,force
      real(t_p) expterm,bde
      real(t_p) dt,dt2,deddt
      real(r_p) dedx,dedy,dedz
      real(t_p) xab,yab,zab,rab
      real(t_p) time0,time1
      logical proceed
!$acc routine(image_acc) seq  

      if(deb_Path) write(*,*) 'ebond1gpu'
      call timer_enter( timer_ebond )
c
c     calculate the bond stretch energy and first derivatives
c
!$acc parallel loop present(x,y,z,use,loc,bndglob,deb
!$acc&    ,eb,eW1aMD,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz) async
#ifdef USE_NVSHMEM_CUDA
!$acc&         deviceptr(d_ibnd,d_bl,d_bk)
#else
!$acc&         present(bl,bk,ibnd)
#endif
      do ibond = 1, nbondloc
         i     = bndglob(ibond)
#ifdef USE_NVSHMEM_CUDA
         ipe   = (i-1)/nbond_pe
         ind   = mod((i-1),nbond_pe) +1
         ia    = d_ibnd(ipe)%pel(1,ind)
         ib    = d_ibnd(ipe)%pel(2,ind)
         ideal = d_bl  (ipe)%pel(ind)
         force = d_bk  (ipe)%pel(ind)
#else
         ia    = ibnd(1,i)
         ib    = ibnd(2,i)
         ideal = bl(i)
         force = bk(i)
#endif
         ialoc = loc(ia)
         ibloc = loc(ib)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (proceed)  proceed = (use(ia) .or. use(ib))
c
c     compute the value of the bond length deviation
c
         if (proceed) then
            xab = x(ia) - x(ib)
            yab = y(ia) - y(ib)
            zab = z(ia) - z(ib)
            if (use_polymer)  call image_inl (xab,yab,zab)
            rab = sqrt(xab*xab + yab*yab + zab*zab)
            dt = rab - ideal
            if (dt.eq.0.0) cycle
c
c     harmonic potential uses Taylor expansion of Morse potential
c     through the fourth power of the bond length deviation
c
            if (bndtyp .eq. 'HARMONIC') then
               dt2 = dt * dt
               e = bndunit * force * dt2 * (1.0_ti_p+cbnd*dt+qbnd*dt2)
               deddt = 2.0_ti_p * bndunit * force * dt
     &                 * (1.0_ti_p+1.5_ti_p*cbnd*dt+2.0_ti_p*qbnd*dt2)
c
c     Morse potential uses energy = BDE * (1 - e**(-alpha*dt))**2)
c     with the approximations alpha = sqrt(ForceConst/BDE) = -2
c     and BDE = Bond Dissociation Energy = ForceConst/alpha**2
c
            else if (bndtyp .eq. 'MORSE') then
               expterm = exp(-2.0_ti_p*dt)
               bde   = 0.25_ti_p * bndunit * force
               e     = bde * (1.0_ti_p-expterm)**2
               deddt = 4.0_ti_p * bde * (1.0_ti_p-expterm) * expterm
            end if
c
c     compute chain rule terms needed for derivatives
c
            if (rab .eq. 0.0_ti_p) then
               de = 0.0_ti_p
            else
               de = deddt / rab
            end if
            dedx = de * xab
            dedy = de * yab
            dedz = de * zab
c
c     increment the total bond energy and first derivatives
c
            eb = eb + e
!$acc atomic update      
            deb(1,ialoc) = deb(1,ialoc) + dedx
!$acc atomic update  
            deb(2,ialoc) = deb(2,ialoc) + dedy
!$acc atomic update          
            deb(3,ialoc) = deb(3,ialoc) + dedz
c
!$acc atomic update
            deb(1,ibloc) = deb(1,ibloc) - dedx
!$acc atomic update
            deb(2,ibloc) = deb(2,ibloc) - dedy
!$acc atomic update
            deb(3,ibloc) = deb(3,ibloc) - dedz
c
c     aMD storage if waters are considered
c
            if (use_amd_wat1) then
            if (type(ia) == aMDwattype(1) .or. type(ib)
     $      == aMDwattype(1)) then
               eW1aMD = eW1aMD + e
!$acc atomic
               deW1aMD(1,ialoc) = deW1aMD(1,ialoc) + dedx
!$acc atomic
               deW1aMD(2,ialoc) = deW1aMD(2,ialoc) + dedy
!$acc atomic
               deW1aMD(3,ialoc) = deW1aMD(3,ialoc) + dedz
!$acc atomic
               deW1aMD(1,ibloc) = deW1aMD(1,ibloc) - dedx
!$acc atomic
               deW1aMD(2,ibloc) = deW1aMD(2,ibloc) - dedy
!$acc atomic
               deW1aMD(3,ibloc) = deW1aMD(3,ibloc) - dedz
            end if
            end if
c
c     increment the internal virial tensor components
c
            g_vxx = g_vxx + xab * real(dedx,t_p)
            g_vxy = g_vxy + yab * real(dedx,t_p)
            g_vxz = g_vxz + zab * real(dedx,t_p)
            g_vyy = g_vyy + yab * real(dedy,t_p)
            g_vyz = g_vyz + zab * real(dedy,t_p)
            g_vzz = g_vzz + zab * real(dedz,t_p)
         end if
      end do

      call timer_exit( timer_ebond )
      end
