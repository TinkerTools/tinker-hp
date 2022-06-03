c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine cutoffs  --  set distance                      ##
c     ##                                                            ##
c     ################################################################
c
c
c     "cutoffs" initializes and stores spherical energy cutoff
c     distance windows and Ewald sum cutoffs,
c     and the pairwise neighbor generation method
c
c
#include "tinker_precision.h"
      subroutine cutoffs
      use atoms  ,only: n
      use bound
      use cutoff
      use domdec
      use keys
      use iounit
      use mdstuf ,only: integrate
      use neigh
      use potent ,only: use_dewald
      use tinheader
      implicit none
      integer i,next
      real(t_p) big,value
      real(t_p) readBuf
      logical truncate
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     set defaults for spherical energy cutoff distances
c
      integrate = 'VERLET'
      big       = 1.0d12
      if (use_bounds) then
         vdwcut   = 9.0_ti_p
         dispcut  = 9.0_ti_p
         repcut   = 9.0_ti_p
         ctrncut  = 9.0_ti_p
         chgcut   = 9.0_ti_p
         mpolecut = 9.0_ti_p
      else
         vdwcut   = big
         dispcut  = big
         repcut   = big
         ctrncut  = big
         chgcut   = big
         mpolecut = big
      end if
      repcut        = 6.0_ti_p
      ctrncut       = 6.0_ti_p
      ewaldcut      = 7.0_ti_p
      ewaldshortcut = 5.0_ti_p
      repshortcut   = 4.0_ti_p
      dispshortcut  = 7.0_ti_p
      dewaldcut     = 7.0_ti_p
      dewaldshortcut= 5.0_ti_p
      ctrnshortcut  = 4.0_ti_p
      mpoleshortcut = 5.0_ti_p
      chgshortcut   = 5.0_ti_p
      vdwshortcut   = 7.0_ti_p
      ddcut         = 0.0_ti_p
      readBuf       = a_init
c
      shortheal = 0.5_ti_p
c
c     set defaults for tapering, neighbor buffers
c
      vdwtaper   = 0.90_ti_p
      reptaper   = 0.90_ti_p
      disptaper  = 0.90_ti_p
      ctrntaper  = 0.90_ti_p
      chgtaper   = 0.65_ti_p
      mpoletaper = 0.65_ti_p
c
c     set defaults for Ewald sum, tapering style and neighbor method
c
      use_ewald = .false.
      truncate  = .false.
      use_list  = .true.
      use_vlist = .true.
      use_dlist = .true.
      use_mlist = .true.
      use_clist = .true.
      use_shortvlist = .false.
      use_shortdlist = .false.
      use_shortmlist = .false.
      use_shortclist = .false.
      rebuild_lst    = .true.
c
c      set default for neighbor list update, reference = 20 for a 2fs time step
c
      ineigup = 20
c
c     set default buffer for nblist
c
#if TINKER_SINGLE_PREC+TINKER_MIXED_PREC
      defaultlbuffer  = merge( 0.4,0.7,(n.gt.95000) )
      defaultlbuffer1 = merge( 0.5,1.0,(n.gt.95000) )
#else
      defaultlbuffer  = 2.0
      defaultlbuffer1 = 2.0
#endif
c
c     search the keywords for various cutoff parameters
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:240)
c
c     get values related to use of Ewald summation
c
         if (keyword(1:6) .eq. 'EWALD ') then
            use_ewald = .true.
         else if (keyword(1:13) .eq. 'EWALD-CUTOFF ') then
            read (string,*,err=10,end=10)  ewaldcut
         else if (keyword(1:13) .eq. 'EWALDSHORT-CUTOFF ') then
            read (string,*,err=10,end=10)  ewaldshortcut
         else if (keyword(1:11) .eq. 'SHORT-HEAL ') then
            read (string,*,err=10,end=10)  shortheal
c
c     get values related to use of Ewald for dispersion
c
         else if (keyword(1:7) .eq. 'DEWALD ') then
            use_dewald = .true.
         else if (keyword(1:14) .eq. 'DEWALD-CUTOFF ') then
            read (string,*,err=10,end=10)  dewaldcut

c     get values for the tapering style and neighbor method
c
         else if (keyword(1:9) .eq. 'TRUNCATE ') then
            truncate = .true.
c
c     get the cutoff radii for potential energy functions
c
         else if (keyword(1:7) .eq. 'CUTOFF ') then
            read (string,*,err=10,end=10)  value
            vdwcut  = value
            mpolecut = value
            ewaldcut = value
            dispcut = value
            repcut  = value
            ctrncut = value
         else if (keyword(1:11) .eq. 'VDW-CUTOFF ') then
            read (string,*,err=10,end=10)  vdwcut
         else if (keyword(1:16) .eq. 'VDWSHORT-CUTOFF ') then
            read (string,*,err=10,end=10)  vdwshortcut
         else if (keyword(1:14) .eq. 'REPULS-CUTOFF ') then
            read (string,*,err=10,end=10)  repcut
         else if (keyword(1:19) .eq. 'REPULSSHORT-CUTOFF ') then
            read (string,*,err=10,end=10)  repshortcut
         else if (keyword(1:12) .eq. 'DISP-CUTOFF ') then
            read (string,*,err=10,end=10)  dispcut
         else if (keyword(1:12) .eq. 'DISPSHORT-CUTOFF ') then
            read (string,*,err=10,end=10)  dispshortcut
         else if (keyword(1:11) .eq. 'CHG-CUTOFF ') then
            read (string,*,err=10,end=10)  chgcut
         else if (keyword(1:11) .eq. 'CHGSHORT-CUTOFF ') then
            read (string,*,err=10,end=10)  chgshortcut
         else if (keyword(1:13) .eq. 'MPOLE-CUTOFF ') then
            read (string,*,err=10,end=10)  mpolecut
         else if (keyword(1:13) .eq. 'MPOLESHORT-CUTOFF ') then
            read (string,*,err=10,end=10)  mpoleshortcut
         else if (keyword(1:14) .eq. 'CHGTRN-CUTOFF ') then
            read (string,*,err=10,end=10)  ctrncut
         else if (keyword(1:19) .eq. 'CHGTRNSHORT-CUTOFF ') then
            read (string,*,err=10,end=10)  ctrnshortcut
c
c     get distance for initialization of energy switching
c
         else if (keyword(1:6) .eq. 'TAPER ') then
            read (string,*,err=10,end=10)  value
            vdwtaper  = value
            chgtaper  = value
            mpoletaper = value
            disptaper = value
            reptaper  = value
            ctrntaper = value
         else if (keyword(1:10) .eq. 'VDW-TAPER ') then
            read (string,*,err=10,end=10)  vdwtaper
         else if (keyword(1:13) .eq. 'REPULS-TAPER ') then
            read (string,*,err=10,end=10)  reptaper
         else if (keyword(1:11) .eq. 'DISP-TAPER ') then
            read (string,*,err=10,end=10)  disptaper
         else if (keyword(1:12) .eq. 'MPOLE-TAPER ') then
            read (string,*,err=10,end=10)  mpoletaper
         else if (keyword(1:13) .eq. 'CHGTRN-TAPER ') then
            read (string,*,err=10,end=10)  ctrntaper
c
c     get buffer width for use with pairwise neighbor lists
c
         else if (keyword(1:12) .eq. 'LIST-BUFFER ') then
            read (string,*,err=10,end=10) readBuf
c
c     fetch integrator
c
         else if (keyword(1:11) .eq. 'INTEGRATOR ') then
            call getword (record,integrate,next)
            call upcase (integrate)
         end if
   10    continue
      end do
c
c     keyword to include another (bigger than the non bonded ones) cutoff in the
c     spatial decomposition. This will add the processes within this cutoff in the
c     communication routines: positions, indexes, forces
c
      do i = 1, nkey
         next   = 1
         record = keyline(i)
         call upcase (record)
         call gettext (record,keyword,next)
         string = record(next:240)
         if (keyword(1:15) .eq. 'DD-CUTOFF ') then
            read (string,*,err=100,end=100)  ddcut
            if (rank.eq.0) write(iout,1000) ddcut
         end if
 1000    format (/,' Additional cutoff for spatial decompostion: ',
     $     F14.5)
  100    continue
      end do

      !Default buffer
      if (integrate.eq.'RESPA1'.or.integrate.eq.'BAOABRESPA1') then
         lbuffer = defaultlbuffer1
      else
         lbuffer = defaultlbuffer
      end if

      ! Change buffer if keyword
      if (readBuf .ne. a_init) lbuffer = readBuf

      ! Return an error if PME is not used
      if (.not.(use_ewald)) then
        if (rank.eq.0) 
     $   write (*,*) 'This program is only compatible with PME'
        call fatal
      end if

      !apply any Ewald cutoff to charge and multipole terms
      if (use_ewald) then
         mpolecut      = ewaldcut
         mpoleshortcut = ewaldshortcut
         chgcut        = ewaldcut
         chgshortcut   = ewaldshortcut
      end if
      if (use_dewald) then
         dispcut = dewaldcut
         dispshortcut = dewaldshortcut
      end if
c
c     convert any tapering percentages to absolute distances
c
      if (vdwtaper   .lt. 1.0_ti_p) vdwtaper   = vdwtaper  * vdwcut
      if (reptaper   .lt. 1.0_ti_p) reptaper   = reptaper  * repcut
      if (disptaper  .lt. 1.0_ti_p) disptaper  = disptaper * dispcut
      if (ctrntaper  .lt. 1.0_ti_p) ctrntaper  = ctrntaper * ctrncut
      if (mpoletaper .lt. 1.0_ti_p) mpoletaper = mpoletaper* mpolecut
      if (chgtaper   .lt. 1.0_ti_p) chgtaper   = chgtaper  * chgcut
c
c     apply truncation cutoffs if they were requested
c
      if (truncate) then
         vdwtaper   = big
         reptaper   = big
         disptaper  = big
         chgtaper   = big
         mpoletaper = big
         ctrntaper  = big
      end if
      call update_lbuffer(lbuffer)
      end

      subroutine update_lbuffer(list_buff)
      use argue
      use cutoff
      use neigh
      use tinheader
      implicit none
      real(t_p),intent(in):: list_buff
      real(r_p) dt
c
c     set buffer region limits for pairwise neighbor lists
c
      lbuffer    = list_buff
      lbuf2      = (0.5_ti_p*lbuffer)**2
      vbuf2      = (vdwcut+lbuffer)**2
      cbuf2      = (chgcut+lbuffer)**2
      dbuf2      = (dispcut+lbuffer)**2
      mbuf2      = (mpolecut+lbuffer)**2
      vshortbuf2 = (vdwshortcut+lbuffer)**2
      cshortbuf2 = (chgshortcut+lbuffer)**2
      mshortbuf2 = (mpoleshortcut+lbuffer)**2
      dshortbuf2 = (dispshortcut+lbuffer)**2

      ! Update ineigup
      read(arg(3),*,err=30,end=30) dt           ! Fetch timestep among arguments
      dt         = dt*0.001_re_p  ! Convert to picoseconds
      ineigup    = max( 1,int(((real(lbuffer,r_p)*0.02_re_p)/dt)+1d-3) )
 30   continue
      end subroutine
