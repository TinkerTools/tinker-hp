c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine eurey1  --  bond stretch energy & derivatives  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "eurey1" calculates the Urey-Bradley interaction energy and
c     its first derivatives with respect to Cartesian coordinates
c
c
#include "tinker_precision.h"
      module eurey1gpu_inl
      contains
#include "image.f.inc"
      end module

      subroutine eurey1gpu
      use atmlst
      use atoms
      use bound
      use deriv
      use domdec
      use energi
      use eurey1gpu_inl
      use group
      use tinheader ,only:ti_p,re_p
      use urey
      use urypot
      use usage
      use virial
      use timestat
      implicit none
      integer i,ia,ic,iurey
      integer ialoc,icloc
      real(t_p) e,de,ideal,force
      real(t_p) dt,dt2,deddt
      real(r_p) dedx,dedy,dedz
      real(t_p) vxx,vyy,vzz
      real(t_p) vyx,vzx,vzy
      real(t_p) xac,yac,zac,rac
      logical proceed
!$acc routine(image_acc) seq 
         
      if(rank.eq.0.and.tinkerdebug) write(*,*) 'eurey1gpu'
      call timer_enter( timer_eurey1 )
c
c     calculate the Urey-Bradley 1-3 energy and first derivatives
c
!$acc parallel loop present(ureyglob,iury,loc,x,y,z,ul,uk
!$acc&         ,use,deub,eub,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz) async
      do iurey = 1, nureyloc
         i     = ureyglob(iurey)
         ia    = iury(1,i)
         ic    = iury(3,i)
         ialoc = loc(ia)
         icloc = loc(ic)
         ideal = ul(i)
         force = uk(i)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (proceed)  proceed = (use(ia) .or. use(ic))
c
c     compute the value of the 1-3 distance deviation
c
         if (proceed) then
            xac = x(ia) - x(ic)
            yac = y(ia) - y(ic)
            zac = z(ia) - z(ic)
            if (use_polymer)  call image_inl (xac,yac,zac)
            rac = sqrt(xac*xac + yac*yac + zac*zac)
            dt  = rac - ideal
            dt2 = dt * dt
            e   = ureyunit * force * dt2 * (1.0_ti_p+cury*dt+qury*dt2)
            deddt = 2.0_ti_p * ureyunit * force * dt
     &                 * (1.0_ti_p+1.5_ti_p*cury*dt+2.0_ti_p*qury*dt2)
c
c     compute chain rule terms needed for derivatives
c
            de   = deddt / rac
            dedx = de * xac
            dedy = de * yac
            dedz = de * zac
c
c     increment the total Urey-Bradley energy and first derivatives
c
            eub  = eub + e
!$acc atomic update    
            deub(1,ialoc) = deub(1,ialoc) + dedx
!$acc atomic update    
            deub(2,ialoc) = deub(2,ialoc) + dedy
!$acc atomic update    
            deub(3,ialoc) = deub(3,ialoc) + dedz
c
!$acc atomic update    
            deub(1,icloc) = deub(1,icloc) - dedx
!$acc atomic update    
            deub(2,icloc) = deub(2,icloc) - dedy
!$acc atomic update    
            deub(3,icloc) = deub(3,icloc) - dedz
c
c     increment the internal virial tensor components
c
            g_vxx = g_vxx + xac * real(dedx,t_p)
            g_vxy = g_vxy + yac * real(dedx,t_p)
            g_vxz = g_vxz + zac * real(dedx,t_p)
            g_vyy = g_vyy + yac * real(dedy,t_p)
            g_vyz = g_vyz + zac * real(dedy,t_p)
            g_vzz = g_vzz + zac * real(dedz,t_p)
         end if
      end do
c
      call timer_exit( timer_eurey1 )
      end
