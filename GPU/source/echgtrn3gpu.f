c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine echgtrn3  --  charge transfer energy & derivs  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "echgtrn1" calculates the charge transfer energy and first
c     derivatives with respect to Cartesian coordinates
c
c
#include "tinker_macro.h"
      subroutine echgtrn3gpu
      implicit none
#ifdef _OPENACC
      call echgtrn3c_cu
#else
      call echgtrn3
#endif
      end subroutine
c
c
c     ####################################################################
c     ##                                                                ##
c     ##  subroutine echgtrn1cgpu  --  charge transfer derivs via list  ##
c     ##                                                                ##
c     ####################################################################
c
c
c     "echgtrn1cgpu" calculates the charge transfer energy and first
c     derivatives using a pairwise neighbor list on device
c
c
      subroutine echgtrn3c_cu
#ifdef _OPENACC
      use action ,only: nect
      use atoms
      use atmlst
      use bound
      use cell
      use chgpot
      use chgtrn
      use ctrpot
      use couple
      use ctrpot
      use cutoff
      use deriv
      use domdec
      use echgtrncu
      use energi
      use group
      use inform ,only: deb_Path, minmaxone
      use mplpot
      use mpole
      use neigh  ,only:ipole_s=>celle_glob,pglob_s=>celle_pole
     &           , loc_s=>celle_loc, ieblst, eblst
     &           , x_s=>celle_x, y_s=>celle_y, z_s=>celle_z
     &           , seblst_s=>shorteblst,iseblst_s=>ishorteblst
      use potent
      use shunt
      use tinheader
      use usage
      use utilcu ,only: BLOCK_DIM,check_launch_kernel
      use utilgpu
      use virial
      implicit none
      integer   i,st,gS,elsI
      real(t_p) p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend,rinv
     &         ,ctrnshortcut2,f
      character*11 mode
      character*64 srami
c
      if (nct.eq.0) return
      if (deb_Path) then
         srami = 'echgtrn3c_cu'
         if (use_chgtrnshort) then
            srami=trim(srami)//' SHORT'
         else if (use_chgtrnlong) then
            srami=trim(srami)//' LONG'
         end if
         write(*,*) trim(srami)
      end if

      mode = merge('SHORTCHGTRN','CHGTRN',use_chgtrnshort)
c
c     set the coefficients for the switching function
c
      call switch (mode)
      ctrnshortcut2 = merge
     &((ctrnshortcut-shortheal)**2,zeror,use_chgtrnlong)
      p_xbeg = xbegproc(rank+1)
      p_xend = xendproc(rank+1)
      p_ybeg = ybegproc(rank+1)
      p_yend = yendproc(rank+1)
      p_zbeg = zbegproc(rank+1)
      p_zend = zendproc(rank+1)
      rinv   = 1.0d0/real(cut-off,r_p)
      st     = 2*npolelocnlb_pair+1
 
      !set conversion factor
      f      = electric / dielec

      !set kernel ressources
      elsI   = merge(nshortpolelocnlb2_pair,npolelocnlb2_pair
     &              ,use_chgtrnshort)
      gs     = get_GridDim(elsI,BLOCK_DIM)
      def_stream = rec_stream
      def_queue  = rec_queue

      if      (use_chgtrnshort) then
         __TINKER_FATAL__
      else if (use_chgtrnlong) then
c!$acc host_data use_device(ipole_s,pglob_s,loc_s,ieblst,eblst
c!$acc&    ,x_s,y_s,z_s
c!$acc&    ,chgct,dmpct,grplist,wgrp,dect,ered_buf1,vred_buff,nred_buff
c!$acc&    ,mcorrect_ik,mcorrect_scale,ipole,loc,x,y,z)
c      call echgtrn3l_kcu<<<gS,BLOCK_DIM,0,def_stream>>>
c     &     ( ipole_s,pglob_s,loc_s,ieblst,eblst,x_s,y_s,z_s,chgct
c     &     , dmpct,grplist,wgrp
c     &     , dect,ered_buf1,vred_buff,nred_buff
c     &     , npolelocnl, nbloc, npolelocnlb, npolelocnlb2_pair
c     &     , ctrntyp_ID, use_group
c     &     , f,ctrnshortcut,ctrnshortcut2,off2,off,cut2,shortheal,rinv
c     &     , p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
c     &     , mcorrect_ik, mcorrect_scale, ipole, loc, x, y, z
c     &     , n_mscale
c     &     )
c      call check_launch_kernel("echgtrn3l_kcu")
c!$acc end host_data
      else
!$acc host_data use_device(ipole_s,pglob_s,loc_s,ieblst,eblst
!$acc&    ,x_s,y_s,z_s
!$acc&    ,chgct,dmpct,grplist,wgrp,dect,ered_buf1,vred_buff,nred_buff
!$acc&    ,mcorrect_ik,mcorrect_scale,ipole,loc,x,y,z)
      call echgtrn3_kcu<<<gS,BLOCK_DIM,0,def_stream>>>
     &     ( ipole_s,pglob_s,loc_s,ieblst,eblst(st),x_s,y_s,z_s
     &     , chgct,dmpct,grplist,wgrp
     &     , dect,ered_buf1,vred_buff,nred_buff
     &     , npolelocnl, nbloc, npolelocnlb, npolelocnlb2_pair
     &     , ctrntyp_ID, use_group
     &     , f,ctrnshortcut,ctrnshortcut2,off2,off,cut2,shortheal,rinv
     &     , p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
     &     , mcorrect_ik, mcorrect_scale, ipole, loc, x, y, z
     &     , n_mscale
     &     )
      call check_launch_kernel("echgtrn3_kcu")
!$acc end host_data
      end if

      call reduce_energy_action(ect,nect,ered_buf1,def_queue)
#else
      __TINKER_FATAL__
#endif
      end
