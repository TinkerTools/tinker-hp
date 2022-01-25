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
      subroutine cutoffs
      use bound
      use cutoff
      use domdec
      use keys
      use iounit
      use neigh
      use potent
      implicit none
      integer i,next
      real*8 big,value
      logical truncate
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     set defaults for spherical energy cutoff distances
c
      big = 1.0d12
      if (use_bounds) then
         vdwcut = 9.0d0
         dispcut = 9.0d0
         repcut = 9.0d0
         ctrncut = 9.0d0
         chgcut = 9.0d0
         mpolecut = 9.0d0
      else
         vdwcut = big
         dispcut = big
         repcut = big
         ctrncut = big
         chgcut = big
         mpolecut = big
      end if
      repcut = 6.0d0
      ctrncut = 6.0d0
      ewaldcut = 7.0d0
      ewaldshortcut = 5.0d0
      repshortcut = 4.0d0
      dispshortcut = 7.0d0
      dewaldcut = 7.0d0
      dewaldshortcut = 5.0d0
      ctrnshortcut = 4.0d0
      mpoleshortcut = 5.0d0
      chgshortcut = 5.0d0
      vdwshortcut = 7.0d0
      ddcut = 0.0d0
c
      shortheal = 0.5d0
c
c     set defaults for tapering, neighbor buffers
c
      vdwtaper = 0.90d0
      reptaper = 0.90d0
      disptaper = 0.90d0
      ctrntaper = 0.90d0
      chgtaper = 0.65d0
      mpoletaper = 0.65d0
      lbuffer = 2.0d0
c
c     set defaults for Ewald sum, tapering style and neighbor method
c
      use_ewald = .false.
      use_dewald = .false.
      truncate = .false.
      use_list = .true.
      use_vlist = .true.
      use_dlist = .true.
      use_mlist = .true.
      use_clist = .true.
      use_shortvlist = .false.
      use_shortmlist = .false.
      use_shortclist = .false.
      use_shortdlist = .false.
c
c      set default for neighbor list update, reference = 20 for a 2fs time step
c
      ineigup = 20
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
            vdwcut = value
            mpolecut = value
            ewaldcut = value
            dispcut = value
            repcut = value
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
            vdwtaper = value
            chgtaper = value
            mpoletaper = value
            disptaper = value
            reptaper = value
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
            read (string,*,err=10,end=10)  lbuffer
         end if
   10    continue
      end do
c
c     keyword to include another (bigger than the non bonded ones) cutoff in the
c     spatial decomposition. This will add the processes within this cutoff in the
c     communication routines: positions, indexes, forces
c
      do i = 1, nkey
         next = 1
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
c
c     return an error if PME is not used
c
      if (.not.(use_ewald)) then
        if (rank.eq.0) 
     $   write(*,*) 'This program is only compatible with PME'
        call fatal
      end if
c
c     apply any Ewald cutoff to charge and multipole terms
c
      if (use_ewald) then
         mpolecut = ewaldcut
         mpoleshortcut = ewaldshortcut
         chgcut = ewaldcut
         chgshortcut = ewaldshortcut
      end if
      if (use_dewald) then
         dispcut = dewaldcut
         dispshortcut = dewaldshortcut
      end if
c
c     convert any tapering percentages to absolute distances
c
      if (vdwtaper .lt. 1.0d0)  vdwtaper = vdwtaper * vdwcut
      if (reptaper .lt. 1.0d0)  reptaper = reptaper * repcut
      if (disptaper .lt. 1.0d0)  disptaper = disptaper * dispcut
      if (ctrntaper .lt. 1.0d0)  ctrntaper = ctrntaper * ctrncut
      if (mpoletaper .lt. 1.0d0)  mpoletaper = mpoletaper * mpolecut
      if (chgtaper .lt. 1.0d0)  chgtaper = chgtaper * chgcut
c
c     apply truncation cutoffs if they were requested
c
      if (truncate) then
         vdwtaper = big
         reptaper = big
         disptaper = big
         chgtaper = big
         mpoletaper = big
         ctrntaper = big
      end if
c
c     set buffer region limits for pairwise neighbor lists
c
      lbuf2 = (0.5d0*lbuffer)**2
      vbuf2 = (vdwcut+lbuffer)**2
      cbuf2 = (chgcut+lbuffer)**2
      dbuf2 = (dispcut+lbuffer)**2
      mbuf2 = (mpolecut+lbuffer)**2
      vshortbuf2 = (vdwshortcut+lbuffer)**2
      cshortbuf2 = (chgshortcut+lbuffer)**2
      mshortbuf2 = (mpoleshortcut+lbuffer)**2
      dshortbuf2 = (dispshortcut+lbuffer)**2
      return
      end
