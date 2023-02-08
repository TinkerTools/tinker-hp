c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine ani  --  energy & gradient components            ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "gradient" calls subroutines to calculate the potential energy
c     and first derivatives with respect to Cartesian coordinates
c
c
#include "tinker_macro.h"
#ifdef NN_SUPPORT
      module mlpot_inl
      contains
#include "image.f.inc"
      end module

      subroutine set_embedding_weights
      use ani
      use group
      use tinheader
      implicit none

      SELECT case(ml_embedding_mode)
         case(0)
            return
         case(1)
         ! full calculation with embedding
            mlpotscale  = wgrp(2,2)
            wgrp(2,2)   = 0.0_ti_p
         case(2)
         ! only compute intra-group difference
            mlpotscale  = wgrp(2,2)
            wgrp(:,:)   = 0.0_ti_p
            wgrp(2,2)   = -mlpotscale
         case defaULT
            write(0,*) 'unknown embedding mode',ml_embedding_mode
            call fatal
      END SELECT
!$acc update device(wgrp) async
      end subroutine set_embedding_weights

      subroutine ml_potential(dograd)
      use atoms
      use action, only: nemlpot,action_data_ondevice
      use ani
      use mlpot_inl
      use bath
      use domdec
      use atmtyp
      use inform    ,only: deb_Path,minmaxone
      use cell
      use deriv
      use energi
      use potent
      use sizes     ,only: tinkerdebug
      use tinMemory ,only: prmem_request
      use utilgpu   ,only: devicenum
      use neigh     ,only: list_ani, ineignl
      use timestat
      use virial
      implicit none
      logical, intent(in):: dograd
      integer iglob, iloc, idx1, idx2, i, j
      logical,save::f_in=.true.
      integer pairs, natmnl, nstrict,icapt
      integer nadd
      integer(c_int32_t) :: ierr
      real(t_p) :: dx,dy,dz
      real(t_p) :: dx0,dy0,dz0
      integer, allocatable :: istrict(:)
      integer  ,pointer::c_glob(:),list1(:),list2(:)
      real(t_p),pointer::cell_x(:),cell_y(:),cell_z(:)

      if (.not. use_mlpot) return
      if ( mlpotscale .eq. 0.0_ti_p) then
        ! no need to calculate the ML potential
        ! but set a number of interactions 
        ! so that it shows up in analyze
        nemlpot = list_ani%natmnl
        if(action_data_ondevice) then
!$acc update device(nemlpot) async
        endif
        return
      endif

      if (deb_Path) write(*,*) '>> ml_potential'
      if (f_in) then
        if (deb_Path) write(*,*) 'init_ml_ressources'
        ierr =  init_ml_ressources(devicenum,trim(MLpot)//char(0)
     &     ,trim(model_file)//char(0),tinkerdebug)
        if (ierr /= 0) then
          write(0,*) "Error ",ierr," in init_ml_ressources"
          call fatal
        endif
        ml_resources_initialized = .true.
c        if (use_mpole.or.use_polar.or.use_bond.or.use_angle
c     &     .or.use_tors.or.use_charge.or.use_urey) then
c            use_ani_only = .false.
c         end if
      end if

      if (f_in.or.isobaric) then
         cell_a=0.0d0
         cell_a(1,1)=xcell
         cell_a(2,2)=ycell
         cell_a(3,3)=zcell
         f_in=.false.
      end if

      if (deb_Path) write(*,*) 'ML POTENTIAL starting'

      !TODO Make atomic_ani an integer(8)
      nadd   = 0
      pairs  = list_ani%npairs
      natmnl = list_ani%natmnl

      call realloc_position_buffer(natmnl,pairs)
      call set_pairlist_Cellorder(amloc,list_ani,.false.)

      if (deb_Path) write(*,*) 'ML POTENTIAL set_pairlist'

      c_glob => list_ani%c_glob
      list1(1:pairs) => list_ani%list(1:pairs)
      list2(1:pairs) => list_ani%list(1+pairs:2*pairs)
      cell_x => list_ani%cell_x
      cell_y => list_ani%cell_y
      cell_z => list_ani%cell_z

      allocate(istrict(natmnl))
!$acc data create(istrict)

      nstrict=0
!$acc parallel loop async default(present) copy(nstrict)
      do i=1,natmnl
         iglob             = c_glob(i)
         xyz_c(i         ) = x(iglob)
         xyz_c(i+  natmnl) = y(iglob)
         xyz_c(i+2*natmnl) = z(iglob)
         atomic_ani(i)     = atomic(iglob)
         if(loc(iglob)<=nloc) then
!$acc atomic capture
            nstrict = nstrict + 1
            icapt = nstrict
!$acc end atomic
            istrict(icapt) = i-1
         endif
      enddo

!$acc parallel loop async default(present)
      do j=1,pairs
         idx1 = list1(j)+1
         idx2 = list2(j)+1
         dx0        = cell_x(idx1) - cell_x(idx2)
         dy0        = cell_y(idx1) - cell_y(idx2)
         dz0        = cell_z(idx1) - cell_z(idx2)
         dx=dx0; dy=dy0; dz=dz0
         call image_inl(dx,dy,dz)
         dp(j        ) = dx - dx0 
         dp(j+  pairs) = dy - dy0
         dp(j+2*pairs) = dz - dz0
         d2(j)     = (dx*dx + dy*dy + dz*dz)**0.5
      enddo

      if (deb_Path) write(*,*) 'ML POTENTIAL call external ml_models'
!$acc wait
      call timer_enter(timer_b2)
!$acc host_data use_device(xyz_c,d_ani,atomic_ani,istrict
!$acc&                    ,list1,list2,d2, dp,aeani)
        ierr=ml_models(xyz_c,aeani,d_ani,cell_a,atomic_ani
     &                ,list1,list2,d2,dp,istrict
     &                 ,natmnl,pairs,nstrict
     &                ,merge(1,0,dograd))
!$acc end host_data
        if (ierr /= 0) then
          write(0,*) "Error ",ierr," in ml_models"
          call fatal
        else if(deb_Path) then
           write(*,*) 'ML POTENTIAL successfully got result ' 
     &                  ,'from external routine!'
        endif
!$acc wait
      call timer_exit(timer_b2)
      if (deb_Path) write(*,*) 'ML POTENTIAL save energy/forces'
      if(dograd) then
!$acc parallel loop
!$acc&  default(present) reduction(+:nadd) 
!$acc&  present(emlpot,g_vxx,g_vyy,g_vzz,g_vxy,g_vyz,g_vxz)
        do i = 1,natmnl
          iglob= c_glob(i)
          iloc = loc(iglob)
!$acc loop seq
          do j=1,3
          dmlpot(j,iloc) = dmlpot(j,iloc) + d_ani(j+3*(i-1))*mlpotscale
          enddo
          g_vxx = g_vxx + dmlpot(1,iloc)*x(iglob)
          g_vyy = g_vyy + dmlpot(2,iloc)*y(iglob)
          g_vzz = g_vzz + dmlpot(3,iloc)*z(iglob)
          g_vxy = g_vxy + 0.5d0*(dmlpot(1,iloc)*y(iglob)
     &                     +dmlpot(2,iloc)*x(iglob))
          g_vyz = g_vyz + 0.5d0*(dmlpot(2,iloc)*z(iglob)
     &                     +dmlpot(3,iloc)*y(iglob))
          g_vxz = g_vxz + 0.5d0*(dmlpot(3,iloc)*x(iglob)
     &                     +dmlpot(1,iloc)*z(iglob))
          if (iloc.le.nloc) then
            emlpot = emlpot + aeani(i)
            nadd = nadd +1
          else
            aeani(i)=0.0d0 ! remove wrong energies
          end if
        end do
      else
!$acc parallel loop 
!$acc&         default(present) present(emlpot) reduction(+:nadd)
        do i = 1,natmnl
          iloc = loc(c_glob(i))
          if (iloc.le.nloc) then
            emlpot = emlpot + aeani(i)
            nadd = nadd +1
          else
            aeani(i)=0.0d0 ! remove wrong energies
          end if
        end do
      endif

      if(mlpotscale.ne.1.0_ti_p) then
!$acc serial async present(emlpot)
        emlpot = emlpot*mlpotscale
!$acc end serial
      endif
      nemlpot = nadd
      if(action_data_ondevice) then
!$acc update device(nemlpot) async
      endif


!$acc end data
      deallocate(istrict)

      if (rank.lt.2.and.tinkerdebug.gt.0) then
!$acc wait
!$acc update host(emlpot)
 24   format(A,F14.4,' loc nl pairs',2I8,I10,F6.2,' rank',I2)
         print 24, 'e ani',emlpot
     &           ,nloc,natmnl,pairs,n/(1.0d0*natmnl),rank
         end if
         if (nadd.ne.nloc.and.naml.eq.n) then
 34   format("Found an issue during Force and energy
     & reduction !! nloc .ne. nreduction",/,2I7)
            write(0,34) nloc,nadd
      end if

      if (deb_Path) write(*,*) '<< ml_potential'


      end subroutine ml_potential
#else
      subroutine ml_potential
      use domdec
      implicit none

      if(ranktot==1) then
        write(0,*) 'Error: ML potential not activated.' 
     &     //' Must compile with NN_SUPPORT=1'
        __TINKER_FATAL__
      endif
      end subroutine ml_potential
      subroutine set_embedding_weights
      end subroutine
#endif
