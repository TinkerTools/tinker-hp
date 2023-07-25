c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine eangle3  --  angle bending energy & analysis  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "eangle3" calculates the angle bending potential energy, also
c     partitions the energy among the atoms; projected in-plane
c     angles at trigonal centers, spceial linear or Fourier angle
c     bending terms are optionally used
c
c
#include "tinker_macro.h"
      module eangle3gpu_inl
#include "atomicOp.h.f"      
      contains

#include "ker_angle.inc.f"

      subroutine eangle3gpu_(dea,deW1aMD,atmType,afld)
      use action
      use analyz
      use angle    ,only: nangleloc,iang,anat,ak
      use angpot
      use atmlst
      use atmtyp
      use atoms    ,only: n,x,y,z
      use bound
      use domdec
      use energi
      use group
      use inform
      use iounit
      use math
      use mamd
      use potent  ,only: use_amd_wat1
      use usage
      use tinheader
      use timestat,only:timer_enter,timer_exit,timer_eangle
      use virial
      implicit none
      integer  ,intent(in)   :: atmType(:)
      real(t_p),intent(inout):: afld(:)
      real(r_p),intent(inout):: dea(1),deW1aMD(1)
      integer   i,ia,ib,ic,id,iangle,angtypii,tfea,tver
      logical   proceed,header,huge
#ifdef USE_NVSHMEM_CUDA
      integer   ipe,ind
#endif
      real(t_p) ideal,force,e,fgrp,fold
      character*9 label
      parameter(
     &    tfea=__use_groups__+__use_polymer__,
     &    tver=__use_ene__+__use_act__
     &    )

      if(deb_Path) write(*,*) 'eangle3gpu'
      call timer_enter( timer_eangle )
c
c     zero out the angle bending energy and partitioning terms
c
      nea = 0
      ea  = 0.0_re_p
c     aea = 0.0_ti_p
      header = rank.eq.0
c
c     calculate the bond angle bending energy term
c
!$acc parallel loop 
#ifdef USE_NVSHMEM_CUDA
!$acc&         present(x,y,z,loc,use,aea,angleglob
#else
!$acc&         present(x,y,z,loc,use,iang,angleglob
#endif
!$acc&    ,anat,ak,afld,angtypI,grplist,wgrp,atmType,aMDwattype
!$acc&    ,ea,eW1aMD,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz
!$acc&    ,fmat_ps,dfmat_ps) async
!$acc&         reduction(+:ea,nea)
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
         angtypii = angtypI(i)
c
c     decide whether to compute the current interaction
c
         if (angtypii .eq. ANG_IN_PLANE) then
            id   = iang(4,i)
            if(use_group)
     &         call groups4_inl(fgrp,ia,ib,ic,id,ngrp,grplist,wgrp)
            proceed = (use(ia).or. use(ib).or. use(ic).or. use(id))
         else
            if(use_group)
     &         call groups3_inl(fgrp,ia,ib,ic,ngrp,grplist,wgrp)
            proceed = (use(ia).or. use(ib).or. use(ic))
         end if
c
c     get the coordinates of the atoms in the angle
c
         if (proceed) then
c
c     compute the bond angle bending energy
c
            if (angtypii .ne. ANG_IN_PLANE) then
 
               call ker_angle(i,ia,ib,ic,loc,ideal,force
     &                 ,angunit,cang,pang,sang,qang,fgrp
     &                 ,use_group,use_polymer,angtypii,x,y,z,afld
     &                 ,ea,e,dea,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz
     &                 ,g_vzz,tver,tfea,fmat_ps,dfmat_ps)
 
              !increment the total bond angle bending energy
               nea = nea + 1
#if 0
               ib = loc(ib)
               aea(ib) = aea(ib) + e
               fold  = afld(i)
c
c     print a message if the energy of this interaction is large
c
               huge = (e .gt. 5.0_ti_p)
               if (debug .or. (verbose.and.huge)) then
                  if (header) then
                     header = .false.
                     write (iout,10)
   10                format (/,' Individual Angle Bending',
     &                          ' Interactions :',
     &                       //,' Type',18x,'Atom Names',18x,
     &                          'Ideal',4x,'Actual',6x,'Energy',/)
                  end if
                  label = 'Angle    '
                  if (angtypii .eq. ANG_LINEAR) then
                     label = 'Angle-Lin'
                  else if (angtypii .eq. ANG_FOURIER) then
                     label = 'Angle-Cos'
                     ideal = (ideal+180.0_ti_p) / fold
                     if (angle1-ideal .gt. 180.0_ti_p/fold)
     &                  ideal = ideal + 360.0_ti_p/fold
                  end if
                  write (iout,20)  label,ia,name(ia),ib,name(ib),
     &                             ic,name(ic),ideal,angle1,e
   20             format (1x,a9,1x,i7,'-',a3,i7,'-',a3,i7,
     &                       '-',a3,2x,2f10.4,f12.4)
               end if
#endif
c
c     compute the projected in-plane angle bend energy
c
            else
               call ker_angle_plan(i,ia,ib,ic,id,loc,ideal,force
     &                 ,angunit,cang,pang,sang,qang,fgrp
     &                 ,use_group,use_polymer,angtypii,x,y,z,afld
     &                 ,use_amd_wat1,atmType,aMDwattype,eW1aMD,ea,e
     &                 ,deW1aMD,dea,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz
     &                 ,tver,tfea)
               !Increment the total bond angle bending energy
               if(e.ne.0.0) nea = nea + 1
#if 0
               ib = loc(ib)
               aea(ib) = aea(ib) + e
c
c     print a message if the energy of this interaction is large
c
               huge = (e .gt. 5.0_ti_p)
               if (debug .or. (verbose.and.huge)) then
                  if (header) then
                     header = .false.
                     write (iout,30)
   30                format (/,' Individual Angle Bending',
     &                         ' Interactions :',
     &                      //,' Type',18x,'Atom Names',18x,
     &                         'Ideal',4x,'Actual',6x,'Energy',/)
                  end if
                  write (iout,40) ia,name(ia),ib,name(ib),ic,name(ic)
     &                           ,ideal,angle1,e
   40             format (' Angle-IP',2x,i7,'-',a3,i7,'-',a3,i7,
     &                    '-',a3,2x,2f10.4,f12.4)
               end if
#endif
            end if
         end if
      end do
      call timer_exit( timer_eangle )
      end

      end module

      subroutine eangle3gpu
      use angle
      use atoms
      use deriv
      use eangle3gpu_inl,only: eangle3gpu_
      use utilgpu       ,only: lam_buff
      implicit none
      ! lam_buff acts as an unsued but allocated variable to this call
      call eangle3gpu_(lam_buff,lam_buff,type,afld)
      end subroutine
