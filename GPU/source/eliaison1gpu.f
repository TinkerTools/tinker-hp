c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine eliaison1                                      ##
c     ##                                                            ##
c     ################################################################
c
      ! off_* defines the offset device thread in charge of *
      ! num_* defines the total number of threads in charge of *
#include "tinker_precision.h"
      module eliaison1gpu_inl
      integer,save:: eliaison1cu_call=0
      integer off_bond,num_bond,off_urey,num_urey,off_angl,num_angl
     &       ,off_tors,num_tors,off_pitr,num_pitr
     &       ,off_totr,num_totr,off_antr,num_antr
     &       ,off_opbd,num_opbd,off_stbd,num_stbd
     &       ,off_impr,num_impr,off_impt,num_impt
     &       ,off_def ,num_def
      contains
      integer function f_mul_warp( iin ) result( oout )
      use utilgpu ,only: WARP_SHFT
      integer,intent(in):: iin
      oout = merge( ishft(ishft((iin-1),-WARP_SHFT)+1,WARP_SHFT) , 0
     &            , iin.ne.0)
      end function
#include "convert.f.inc"
      end module
c
c     "eliaison1" calculates all bonded terms in a single routine
c
      subroutine eliaison1cu
#ifdef _OPENACC
      use angle
      use angpot
      use angtor
      use atmlst
      use atmtyp
      use atoms
      use atomsMirror ,only: x_=>x,y_=>y,z_=>z
      use bitor
      use bndpot
      use bond
      use bound
      use couple
      use deriv
      use domdec
      use energi
      use eliaison1gpu_inl
      use cudafor   ,only: dim3
      use eliaisoncu,only: eliaison1_kcu,lia_Gs,lia_Bd
      use eliaisoncu_dat,only: grx,gry,grz,grtx,grty,grtz
     &                  ,eliaisoncu_set_pointers=>attach_p
      use utilcu    ,only: check_launch_kernel
     &              ,transpose_z3m,TRP_BLOCK_DIM
      use group
      use improp
      use imptor
      use inform   ,only: deb_Path
      use ktrtor   ,ktt_dalloc=>malloc_device_data
      use mamd
      use math
      use nvshmem
      use opbend
      use pitors
      use potent   ,only: use_amd_wat1,use_amd_dih
      use strbnd
      use timestat ,only: timer_enter,timer_exit,timer_b1
     &             ,quiet_timers,print_timers
      use tinheader
      use torpot
      use tortor
      use tors
      use urey
      use urypot
      use usage
      use utilgpu
      use virial
      implicit none
      integer i,j,siz
      integer,parameter:: marge=15
      mdyn_rtyp ssum

      if(deb_Path) write(*,*) 'eliaison1gpu'

      if (nproc.eq.1) then
         if (eliaison1cu_call.ne.0) goto 31
      end if

      lia_Gs  = min(max(nbondloc+nureyloc+nangleloc+ntorsloc
     &              +niproploc,lia_Bd)/lia_Bd
     &             ,maxBlock)
      off_def = 0
      if (1.and.lia_Gs.le.maxBlock-marge) then
         off_bond = 0
         off_urey = off_bond + f_mul_warp(nbondloc)
         off_angl = off_urey + f_mul_warp(nureyloc)
         off_tors = off_angl + f_mul_warp(nangleloc)
         off_pitr = off_tors + f_mul_warp(ntorsloc)
         off_totr = off_pitr + f_mul_warp(npitorsloc)
         off_antr = off_totr + f_mul_warp(ntortorloc)
         off_impr = off_antr + f_mul_warp(nangtorloc)
         off_impt = off_impr + f_mul_warp(niproploc)
         off_opbd = off_impt + f_mul_warp(nitorsloc)
         off_stbd = off_opbd + f_mul_warp(nopbendloc)
         num_def  = off_stbd + f_mul_warp(nstrbndloc)
         lia_Gs   = (num_def-1)/lia_Bd + 1

         num_bond = num_def
         num_urey = num_def  - off_urey
         num_angl = num_def  - off_angl
         num_tors = num_def  - off_tors
         num_impr = num_def  - off_impr
         num_impt = num_def  - off_impt
         num_pitr = num_def  - off_pitr
         num_totr = num_def  - off_totr
         num_antr = num_def  - off_antr
         num_opbd = num_def  - off_opbd
         num_stbd = num_def  - off_stbd
      else
         lia_Gs  = max(nbondloc,nureyloc,nangleloc,ntorsloc)
     &                /lia_Bd+1
         num_def= lia_Gs*lia_Bd
         off_bond=0; num_bond=num_def
         off_urey=0; num_urey=num_def
         off_angl=0; num_angl=num_def
         off_tors=0; num_tors=num_def
         off_impr=0; num_impr=num_def
         off_impt=0; num_impt=num_def
         off_pitr=0; num_pitr=num_def
         off_totr=0; num_totr=num_def
         off_antr=0; num_antr=num_def
         off_opbd=0; num_opbd=num_def
         off_stbd=0; num_stbd=num_def
      end if

!$acc host_data use_device(de_buff1)
      call eliaisoncu_set_pointers(de_buff1,dr_stride,dr_obws,2)
!$acc end host_data
 31   continue

      siz = dr_stride

      if (eliaison1cu_call.eq.0) then
         ! allocate ktrtor data when .not.use_tortor
         call ktt_dalloc
c20   format('   ',13A9)
c30   format(A3,13I9)
c     write(*,20) "bond","urey","angle","tors","pitors","tortor"
c    &           ,"angtor","improp","imptor","opbend","stdbnd"
c    &           ,"def"
c     write(*,30) "n",nbondloc,nureyloc,nangleloc,ntorsloc,npitorsloc
c    &           ,ntortorloc,nangtorloc,niproploc,nitorsloc,nopbendloc
c    &           ,nstrbndloc
c     write(*,30) "off",off_bond,off_urey,off_angl,off_tors,off_pitr
c    &           ,off_totr,off_antr,off_impr,off_impt,off_opbd,off_stbd
c     write(*,30) "num",num_bond,num_urey,num_angl,num_tors,num_pitr
c    &           ,num_totr,num_antr,num_impr,num_impt,num_opbd,num_stbd
c    &           ,num_def
      end if

!$acc host_data use_device(bndglob,ibnd,bl,bk,type,deb,aMDwattype
!$acc&         ,ureyglob,iury,ul,uk,deub
!$acc&         ,angleglob,angtypI,iang,anat,ak,afld,dea
!$acc&         ,opbendglob,iopb,opbk,deopb
!$acc&         ,strbndglob,isb,sbk,deba
!$acc&         ,torsglob,itors,tors1,tors2,tors3,tors4,tors5,tors6
!$acc&         ,et,det
!$acc&         ,pitorsglob,ipit,kpit,ept,dept
!$acc&         ,tortorglob,itt,ibitor,tnx,tny,ttx,tty,tbf,tbx,tby
!$acc&         ,tbxy,i12,n12,atomic,ett,dett
!$acc&         ,angtorglob,iat,kant,deat
!$acc&         ,impropglob,iiprop,vprop,kprop,eid,deid
!$acc&         ,imptorglob,iitors,itors1,itors2,itors3,eit,deit
!$acc&         ,eDaMD,eW1aMD,deW1aMD
!$acc&         ,x_,y_,z_,x,y,z,loc,iuse
!$acc&         ,ered_buf1,vred_buff)
      call eliaison1_kcu<<<lia_Gs,lia_Bd,0,rec_stream>>>
     &   (!bond
     &    bndglob,ibnd,bl,bk,type,nbondloc,bndunit,cbnd,qbnd,bndtyp_i
     &   ,use_polymer,deb
     &   ,off_bond,num_bond
     &    !urey
     &   ,ureyglob,nureyloc,iury,ul,uk,ureyunit,cury,qury,deub
     &   ,off_urey,num_urey
     &    !angle
     &   ,angleglob,nangleloc,angtypI,iang,anat,ak,afld,angunit
     &   ,cang,qang,pang,sang,dea
     &   ,off_angl,num_angl
     &    !opbend
     &   ,opbendglob,nopbendloc,iopb,opbk,opbunit,opbtypInt
     &   ,copb,qopb,popb,sopb,deopb
     &   ,off_opbd,num_opbd
     &    !strbnd
     &   ,strbndglob,nstrbndloc,isb,sbk,stbnunit,deba
     &   ,off_stbd,num_stbd
     &    !tors
     &   ,torsglob,ntorsloc,itors,tors1,tors2,tors3,tors4,tors5,tors6
     &   ,torsunit,et,det
     &   ,off_tors,num_tors
     &    !pitors
     &   ,pitorsglob,npitorsloc,ipit,kpit,ptorunit,ept,dept
     &   ,off_pitr,num_pitr
     &    !tortor
     &   ,tortorglob,ntortorloc,itt,ibitor,tnx,tny,ttx,tty,tbf,tbx,tby
     &   ,tbxy,i12,n12,atomic,ttorunit,ett,dett
     &   ,off_totr,num_totr
     &    !angtor
     &   ,angtorglob,nangtorloc,iat,kant,atorunit,deat
     &   ,off_antr,num_antr
     &    !improp
     &   ,impropglob,niproploc,iiprop,vprop,kprop,idihunit,eid,deid
     &   ,off_impr,num_impr
     &    !imptor
     &   ,imptorglob,nitorsloc,iitors,itors1,itors2,itors3
     &   ,itorunit,eit,deit
     &   ,off_impt,num_impt
     &  
     &   ,aMDwattype,use_amd_wat1,use_amd_dih,eDaMD,eW1aMD,deW1aMD
     &   ,x_,y_,z_,x,y,z,loc,iuse
     &   ,n,nbloc,siz,useAll,grx,gry,grz,grtx,grty,grtz
     &   ,ered_buf1,vred_buff
     &   )
      call check_launch_kernel("<eliaison1_kcu>")
!$acc end host_data

c      block;  type(dim3) gS,bS;
c      gS  = dim3((siz-1)/TRP_BLOCK_DIM+1,1,1)
c      bS  = dim3(TRP_BLOCK_DIM,3,1)
c!$acc host_data use_device(deb)
c      call transpose_z3m<<<gS,bS,0,rec_stream>>>( grx,deb,siz )
c!$acc end host_data
c      call check_launch_kernel(" transpose_z3f ")
c      end block

      if (use_amd_dih) then
!$acc parallel loop collapse(2) async default(present)
         do i = 1,nloc; do j=1,3
            ssum        = d_x( i+(j-1)*siz )
            deamdD(j,i) = mdr2md( ssum )
         end do; end do
      end if

      if (ftot_l.and.fdebs_l.and.tdes_l) then
      else
!$acc parallel loop collapse(2) async default(present)
         do j=1,3; do i=1,siz;
            ssum = de1x(i+(j-1)*siz) + d_x(i+(j-1)*siz)
            deb ( j,i )         = mdr2md( ssum )
            de1x( i+(j-1)*siz ) = 0
            d_x ( i+(j-1)*siz ) = 0
         end do; end do;
      end if

      if (calc_e.or.use_virial)
     &   call reduce_energy_virial(eb_r,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz
     &                            ,g_vzz,ered_buf1,rec_queue)

      eliaison1cu_call = eliaison1cu_call + 1
#else
 63   format(/,' --- ERROR --- '
             /,' eliaison1gpu routine is not made for host execution ')
      write(0,63); call fatal;
#endif
      end subroutine
