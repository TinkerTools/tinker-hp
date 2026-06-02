!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###############################################################
!     ##                                                           ##
!     ##  subroutine eangtor  --  angle-torsion cross term energy  ##
!     ##                                                           ##
!     ###############################################################
!
!
!     "eangtor" calculates the angle-torsion potential energy
!
!
#include "tinker_macro.h"
module eangtor_inl
#include "atomicOp.h.f90"
contains
#include "ker_angtor.inc.f90"

   subroutine eangtor_(iat,anat,kant,tors1,tors2,tors3,deat)
      use atmlst
      use angtor   ,only: nangtor,nangtorloc
      use atoms
      use bound
      use domdec
      use energi
      use group
      use inform   ,only: deb_Path
      use math
      use tinheader,only: ti_p
      use torpot
      use tors     ,only: itors
      use usage
      use virial
      implicit none
      integer  ,intent(in):: iat(:,:)
      real(t_p),intent(in)::anat(:),kant(:,:)&
         &,tors1(:,:),tors2(:,:),tors3(:,:)
      real(r_p):: deat(1)

      integer grp
      integer i,iiangtor,iangtor,tver,tfea
      integer ia,ib,ic,id
      logical proceed
      real(t_p) fgrp,e
      parameter(&
         &grp=__use_groups__,&
         &tver=__use_ene__,&
         &tfea=__use_polymer__+__use_groups__&
         &)
!
      if (deb_Path) print*,'eangtor'
#ifdef USE_NVSHMEM_CUDA
      ! TODO Remove this check
      ! Implement NVSHMEM Access to anat & itors, tors[123]
      print*, '  FATAL ERROR  '
      print*, 'NVSHMEM feature not implemented inside eangtor3'
      __TINKER_FATAL__
#endif
!
!     zero out the energy due to extra potential terms
!
      eat = 0
!
!     calculate the angle-torsion interaction energy term
!
!$acc parallel loop async &
!$acc         present(eat,deat,angtorglob,iat,anat,grplist,wgrp &
!$acc   ,itors,tors1,tors2,tors3,x,y,z,use,kant,loc) &
!$acc         reduction(+:eat) &
!$acc         private(fgrp)
      do iangtor = 1, nangtorloc
         iiangtor= angtorglob(iangtor)
         i       = iat(1,iiangtor)
         ia      = itors(1,i)
         ib      = itors(2,i)
         ic      = itors(3,i)
         id      = itors(4,i)
         if (IAND(tfea,grp).NE.0.and.use_group)&
            &call groups4_inl(fgrp,ia,ib,ic,id,ngrp,grplist,wgrp)
!
!     decide whether to compute the current interaction
!
         proceed = (use(ia).or. use(ib).or. use(ic).or. use(id))
!
!     compute the value of the torsional angle
!
         if (proceed) then
            call ker_angtor(iiangtor,i,ia,ib,ic,id,loc,iat,radian&
               &,atorunit,fgrp,x,y,z,anat,kant,tors1,tors2,tors3&
               &,use_polymer,use_group,use_virial&
               &,eat,e,deat,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz&
               &,tver,tfea)
         end if
      end do
   end
end module

subroutine eangtor
   use angle
   use angtor
   use eangtor_inl
   use deriv
   use tors
   use utilgpu ,only: lam_buff
   implicit none

   call eangtor_(iat,anat,kant,tors1,tors2,tors3,lam_buff)
end subroutine
