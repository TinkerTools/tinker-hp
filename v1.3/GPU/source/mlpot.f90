!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!

#include "tinker_macro.h"
#ifdef NN_SUPPORT
module mlpot_inl
contains
#include "image.inc.f90"
end module

subroutine set_embedding_weights()
   use mlpot
   use domdec
   use group
   use tinheader
   implicit none

   SELECT CASE(ml_embedding_mode)
    CASE(0)
      return
    CASE(1)
      ! full calculation with embedding
      mlpotscale  = wgrp(2,2)
      wgrp(2,2)   = 0.0_ti_p
    CASE(2)
      ! only compute intra-group difference
      mlpotscale  = wgrp(2,2)
      wgrp(:,:)   = 0.0_ti_p
      wgrp(2,2)   = -mlpotscale
    CASE(3)
      ! only compute chosen intra terms-group difference
      mlpotscale  = 0.0_ti_p
      wgrp(:,:)   = 0.0_ti_p
      wgrp(2,2)   = 1.0_ti_p
    CASE DEFAULT
      write(0,*) 'unknown embedding mode',ml_embedding_mode
      call fatal
   END SELECT

!$acc update device(wgrp) async

end subroutine set_embedding_weights

subroutine ml_potential(dograd)
   use atoms
   use action, only: nemlpot,action_data_ondevice
   use mlpot
   use mlpot_inl
   use bath
   use domdec
   use atmtyp
   use inform    ,only: deb_Path,minmaxone
   use cell
   use deriv
   use energi
   use mutant
   use potent
   use sizes     ,only: tinkerdebug
   use tinMemory ,only: prmem_request
   use utilgpu   ,only: devicenum
   use neigh     ,only: list_mlpot, ineignl
   use timestat
   use virial
   implicit none
   logical, intent(in):: dograd
   integer iglob, iloc, idx1, idx2, i, j
   integer :: use_lambda_int
   logical,save::f_in=.true.
   integer pairs, natmnl,icapt, nstrict,iglob1,iglob2,nat_all
   integer nadd
   integer(c_int32_t) :: ierr
   real(t_p) :: dx,dy,dz,mlpotcut2
   real(t_p) :: dx0,dy0,dz0
   integer  ,pointer::c_glob(:),list1(:),list2(:)
   real(t_p),pointer::cell_x(:),cell_y(:),cell_z(:)
   integer :: max_pairs
   real(c_float), save :: qtot
   real(c_float) :: vlambda_c,elambda_c
   integer,save :: ipass = 0

   if (.not. use_mlpot) return
   if ( mlpotscale .eq. 0.0_ti_p) then
      ! no need to calculate the ML potential
      ! but set a number of interactions
      ! so that it shows up in analyze
      nemlpot = list_mlpot%natmnl
      if(action_data_ondevice) then
!$acc update device(nemlpot) async
      endif
      return
   endif

   if (deb_Path) write(*,*) '>> ml_potential'

   use_lambda_int = 0
   if (nmut > 0) then
      use_lambda_int = 1
      vlambda_c = real(vlambda,c_float)
      elambda_c = real(elambda,c_float)
   endif

   if (f_in) then
!$acc enter data create(dedle_ml,dedlv_ml,vir_ml)
      if (deb_Path) write(*,*) 'init_ml_ressources'
      ierr =  init_ml_ressources(ranktot,devicenum,trim(model_file)//char(0) &
      &,tinkerdebug, ml_port,use_lambda_int,mlqtot,mlqligand &
      &,gc_stride,n,atomic, nmut, imut)
      if (ierr /= 0) then
22       format("Error ",I0," detected from rank ",I0," in&
            & init_ml_ressources")
         write(0,*) ierr,ranktot
         call fatal
      endif
      ml_resources_initialized = .true.
   end if

   if (f_in.or.isobaric) then
      cell_a=0.0d0
      cell_a(1,1)=xcell
      cell_a(2,2)=ycell
      cell_a(3,3)=zcell
   end if

   if (deb_Path) write(*,*) 'ML POTENTIAL starting'

   nadd   = 0
   natmnl = list_mlpot%natmnl
   pairs  = list_mlpot%npairs

   call realloc_position_buffer(natmnl,pairs)
   call set_pairlist_Cellorder(amloc,list_mlpot,.false.)

   if (deb_Path) write(*,*) 'ML POTENTIAL set_pairlist'

   list1(1:pairs) => list_mlpot%list(1:pairs)
   list2(1:pairs) => list_mlpot%list(1+pairs:2*pairs)
   c_glob => list_mlpot%c_glob
   cell_x => list_mlpot%cell_x
   cell_y => list_mlpot%cell_y
   cell_z => list_mlpot%cell_z


   nstrict=namloc_strict
   nat_all = size(atomic_mlpot)
!$acc parallel loop async default(present) !copy(nstrict)
   do i=1,nat_all
      if (i<= natmnl) then
         iglob             = c_glob(i)
         j = 3*(i-1)
         xyz_c(1+j) = x(iglob)
         xyz_c(2+j) = y(iglob)
         xyz_c(3+j) = z(iglob)
         atomic_mlpot(i)     = atomic(iglob)
         if(loc(iglob)<=nloc) then
            trueat_mlpot(i) = 1
         else
            trueat_mlpot(i) = 0
         endif
      else
         j = 3*(i-1)
         xyz_c(1+j) = 0.
         xyz_c(2+j) = 0.
         xyz_c(3+j) = 0.
         atomic_mlpot(i) = -1
         trueat_mlpot(i) = 0
      endif
   enddo

   if (use_lambda_int > 0) then
!$acc parallel loop async default(present)
      do i=1,nat_all
         if (i<= natmnl) then
            iglob = c_glob(i)
            alch_group(i) = int(mutInt(iglob), c_int32_t)
         else
            alch_group(i) = -1
         endif
      enddo
   endif

   mlpotcut2 = mlpotcut**2
   max_pairs=size(edge_src)/2
   !$acc parallel loop async default(present)
   do i=1,max_pairs
      j=3*(i-1)
      if (i<=pairs) then
         idx1 = list1(i)+1
         idx2 = list2(i)+1
         edge_src(i) = list1(i)
         edge_src(i+max_pairs) = list2(i)
         edge_dst(i) = list2(i)
         edge_dst(i+max_pairs) = list1(i)
         dx0        = cell_x(idx2) - cell_x(idx1)
         dy0        = cell_y(idx2) - cell_y(idx1)
         dz0        = cell_z(idx2) - cell_z(idx1)
         dx=dx0; dy=dy0; dz=dz0
         call image_inl(dx,dy,dz)
         dp(1+j) = (dx - dx0)*i_xcell
         dp(2+j) = (dy - dy0)*i_ycell
         dp(3+j) = (dz - dz0)*i_zcell
         dp(1+j+3*max_pairs) = -dp(1+j)
         dp(2+j+3*max_pairs) = -dp(2+j)
         dp(3+j+3*max_pairs) = -dp(3+j)
         d2(i) = dx*dx + dy*dy + dz*dz
         d2(i+max_pairs) = d2(i)
      else
         dp(1+j) = 0.0
         dp(2+j) = 0.0
         dp(3+j) = 0.0
         dp(1+j+3*max_pairs) = 0.0
         dp(2+j+3*max_pairs) = 0.0
         dp(3+j+3*max_pairs) = 0.0
         d2(i) = mlpotcut2
         d2(i+max_pairs) = mlpotcut2
         edge_src(i) = nat_all
         edge_src(i+max_pairs) = nat_all
         edge_dst(i) = nat_all
         edge_dst(i+max_pairs) = nat_all
      endif
   enddo

   if (deb_Path) write(*,*) 'ML POTENTIAL call external ml_models'
!$acc wait
   call timer_enter(timer_b2)
!$acc host_data use_device(xyz_c,d_mlpot,atomic_mlpot,trueat_mlpot &
!$acc                    ,edge_src,edge_dst,d2, dp,aemlp &
!$acc                    ,alch_group, dedle_ml,dedlv_ml, vir_ml)
   ierr=ml_models(xyz_c,aemlp,d_mlpot, vir_ml, cell_a,atomic_mlpot&
      &,edge_src,edge_dst,d2,dp,trueat_mlpot&
      &,natmnl,size(edge_src),nat_all&
      &,merge(1,0,dograd) &
      &,use_lambda_int,elambda_c,vlambda_c,alch_group &
      &, dedle_ml,dedlv_ml)
!$acc end host_data
   if (ierr /= 0) then
24    format("Error ",I0," detected from rank ",I0," in ml_models")
      write(0,24) ierr,rank
      call fatal
   else if(deb_Path) then
      write(*,*) 'ML POTENTIAL successfully got result '&
         &,'from external routine!'
   endif
!$acc wait


   call timer_exit(timer_b2)
   if (deb_Path) write(*,*) 'ML POTENTIAL save energy/forces'
   if(dograd) then
!$acc parallel loop &
!$acc  default(present) reduction(+:nadd) &
!$acc  present(emlpot)
      do i = 1,natmnl
         iglob = c_glob(i)
         iloc  = loc(iglob)
!$acc loop seq
         do j=1,3
            dmlpot(j,iloc) = dmlpot(j,iloc) + d_mlpot(j+3*(i-1))*mlpotscale
         enddo
         if (iloc.le.nloc) then
            emlpot = emlpot + aemlp(i)
            nadd = nadd +1
         else
            aemlp(i)=0.0d0 ! remove wrong energies
         end if
      end do

!$acc serial async present(vir_ml,vir)
      vir(1,1) = vir(1,1) + vir_ml(1)  ! xx
      vir(2,1) = vir(2,1) + vir_ml(2)  ! xy
      vir(3,1) = vir(3,1) + vir_ml(3)  ! xz
      vir(1,2) = vir(1,2) + vir_ml(4)  ! yx
      vir(2,2) = vir(2,2) + vir_ml(5)  ! yy
      vir(3,2) = vir(3,2) + vir_ml(6)  ! yz
      vir(1,3) = vir(1,3) + vir_ml(7)  ! zx
      vir(2,3) = vir(2,3) + vir_ml(8)  ! zy
      vir(3,3) = vir(3,3) + vir_ml(9)  ! zz
!$acc end serial

      if (use_lambda_int > 0) then
!$acc serial async present(dedle_ml,dedlv_ml,delambdae,delambdav)
         delambdae = delambdae + dedle_ml(1)
         delambdav = delambdav + dedlv_ml(1)
!$acc end serial
!$acc update host(delambdae,delambdav) async
      endif


   else
!$acc parallel loop &
!$acc         default(present) present(emlpot) reduction(+:nadd)
      do i = 1,natmnl
         iloc = loc(c_glob(i))
         if (iloc.le.nloc) then
            emlpot = emlpot + aemlp(i)
            nadd = nadd +1
         else
            aemlp(i)=0.0d0 ! remove wrong energies
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


   if (rank.lt.2.and.tinkerdebug.gt.0) then
!$acc wait
!$acc update host(emlpot)
26    format(A,F14.4,' loc nl pairs',2I8,I10,F6.2,' rank',I2)
      print 26, 'emlpot',emlpot&
         &,nloc,natmnl,pairs,n/(1.0d0*natmnl),rank
   end if
   if (nadd.ne.nloc.and.naml.eq.n) then
34    format("Found an issue during Force and energy&
         & reduction !! nloc .ne. nreduction",/,2I7)
      write(0,34) nloc,nadd
   end if

   f_in=.false.

   if (deb_Path) write(*,*) '<< ml_potential'


end subroutine ml_potential

#else
subroutine ml_potential
   use domdec
   implicit none

   if(ranktot==1) then
      write(0,*) 'Error: ML potential not activated.'&
         &//' Must compile with NN_SUPPORT=1'
      __TINKER_FATAL__
   endif
end subroutine ml_potential
subroutine set_embedding_weights
   implicit none
end subroutine
#endif
