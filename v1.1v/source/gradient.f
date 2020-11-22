c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine gradient  --  find energy & gradient components  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "gradient" calls subroutines to calculate the potential energy
c     and first derivatives with respect to Cartesian coordinates
c
c
      subroutine gradient (energy,derivs)
      use atoms
      use sizes
      use deriv
      use domdec
      use energi
      use inter
      use iounit
      use potent
      use timestat
      use vdwpot
      use virial
      use mpi
      implicit none
      integer i,k
      real*8 energy
      real*8 derivs(3,nbloc)
      real*8 time0,time1,time2
      real*8 timevec0,timevec1
      logical isnan
c
c
c     zero out each of the potential energy components
c
      eb = 0.0d0
      ea = 0.0d0
      eba = 0.0d0
      eub = 0.0d0
      eaa = 0.0d0
      eopb = 0.0d0
      eopd = 0.0d0
      eid = 0.0d0
      eit = 0.0d0
      et = 0.0d0
      ept = 0.0d0
      ebt = 0.0d0
      ett = 0.0d0
      ev = 0.0d0
      ec = 0.0d0
      em = 0.0d0
      ep = 0.0d0
      eg = 0.0d0
      ex = 0.0d0
c
c     zero out each of the first derivative components
c
      deb   = 0.0d0
      dea   = 0.0d0
      deba  = 0.0d0
      deub  = 0.0d0
      deaa  = 0.0d0
      deopb = 0.0d0
      deopd = 0.0d0
      deid  = 0.0d0
      deit  = 0.0d0
      det   = 0.0d0
      dept  = 0.0d0
      debt  = 0.0d0
      dett  = 0.0d0
      dev   = 0.0d0
      dec   = 0.0d0
      dem   = 0.0d0
      dep   = 0.0d0
!!     dep1  = 0.0d0
!!     dep2  = 0.0d0
!!     dep3  = 0.0d0
      deg   = 0.0d0
      dex   = 0.0d0
c
c     zero out the virial and the intermolecular energy
c
      vir    = 0.0d0
      einter = 0.0d0
c
c     call the local geometry energy and gradient routines
c
      time0 = mpi_wtime()
      if (use_bond)  call ebond1
      if (use_strbnd)  call estrbnd1
      if (use_urey)  call eurey1
      if (use_angang)  call eangang1
      if (use_opbend)  call eopbend1
      if (use_opdist)  call eopdist1
      if (use_improp)  call eimprop1
      if (use_imptor)  call eimptor1
      if (use_tors)  call etors1
      if (use_pitors)  call epitors1
      if (use_strtor)  call estrtor1
      if (use_tortor)  call etortor1
      if (use_angle)  call eangle1
      time1 = mpi_wtime()
      timebonded = timebonded + time1-time0
c
c     call the van der Waals energy and gradient routines
c
      if (use_vdw) then
         timevec0 = mpi_wtime()
         if (vdwtyp .eq. 'LENNARD-JONES')  call elj1
         if (vdwtyp .eq. 'BUFFERED-14-7')  call ehal1
         timevec1 = mpi_wtime()
         if (rank.eq.0.and.tinkertime) 
     &        write(*,*) 'time vdw = ',timevec1-timevec0
      end if
c
c     call the electrostatic energy and gradient routines
c
      if (use_charge.or.use_mpole) then
         timevec0 = mpi_wtime()
         if (use_charge) call echarge1
         if (use_mpole)  call empole1
         timevec1 = mpi_wtime()
         if (rank.eq.0.and.tinkertime)
     &            write(*,*) 'time elec = ',timevec1-timevec0
      endif
      if (use_polar)  then
         timevec0 = mpi_wtime()
                      call epolar1
         timevec1 = mpi_wtime()
         if (rank.eq.0.and.tinkertime)
     &         write(*,*) 'time polar = ',timevec1-timevec0
      endif
c
c     call any miscellaneous energy and gradient routines
c
      if (use_geom)  call egeom1
      if (use_extra)  call extra1
      time2 = mpi_wtime()
      timenonbonded = timenonbonded + time2-time1
c
c     sum up to get the total energy and first derivatives
c
      esum = eit + eopd + eopb + eaa + eub + eba + ea + eb + em + ep
     $       + ec + ev + et + ept + ebt + ett + eg + ex + eid
      energy = esum
c
      desum = deaa + deba + dea + deb + dec + dem + deopb + deopd + deid
     $ + dep + dev + deub + deit + det + dept + debt + dett + deg + dex
!!     desum(1,:) =  deaa(1,:) + deba(1,:) + dea(1,:) + deb(1,:)
!!    &            + dec(1,:)  + dem(1,:) + deopb(1,:) + deopd(1,:)
!!    &            + deid(1,:) + dep1 + dev(1,:) + deub(1,:) + deit(1,:)
!!    &            + det(1,:) + dept(1,:) + debt(1,:) + dett(1,:)
!!    &            + deg(1,:) + dex(1,:)
!!
!!     desum(2,:) =  deaa(2,:) + deba(2,:) + dea(2,:) + deb(2,:)
!!    &            + dec(2,:)  + dem(2,:) + deopb(2,:) + deopd(2,:)
!!    &            + deid(2,:) + dep2 + dev(2,:) + deub(2,:) + deit(2,:)
!!    &            + det(2,:) + dept(2,:) + debt(2,:) + dett(2,:)
!!    &            + deg(2,:) + dex(2,:)
!!     desum(3,:) =  deaa(3,:) + deba(3,:) + dea(3,:) + deb(3,:)
!!    &            + dec(3,:)  + dem(3,:) + deopb(3,:) + deopd(3,:)
!!    &            + deid(3,:) + dep3 + dev(3,:) + deub(3,:) + deit(3,:)
!!    &            + det(3,:) + dept(3,:) + debt(3,:) + dett(3,:)
!!    &            + deg(3,:) + dex(3,:)
c
c
      debond = deaa + deba + dea + deb + deopb + deopd + deid 
     $ + deub + deit + det + dept + debt + dett + deg

c
      derivs(:,1:nloc) = derivs(:,1:nloc) + desum(:,1:nloc)
c
c     check for an illegal value for the total energy
c
      if (isnan(esum)) then
         write (iout,10)
   10    format (/,' GRADIENT  --  Illegal Value for the Total',
     &              ' Potential Energy')
         call fatal
      end if
      return
      end
