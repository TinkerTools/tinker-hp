c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine analysis  --  energy components and analysis  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "analysis" calls the series of routines needed to calculate
c     the potential energy and perform energy partitioning analysis
c     in terms of type of interaction or atom number
c
c
      subroutine analysis (energy)
      use action
      use analyz
      use domdec
      use energi
      use iounit
      use inter
      use potent
      use vdwpot
      use mpi
      implicit none
      integer ierr
      real*8 energy
      logical isnan
c
c     allocate arrays
c
      if (allocated(aec)) deallocate (aec)
      allocate (aec(nbloc))
      if (allocated(aem)) deallocate (aem)
      allocate (aem(nbloc))
      if (allocated(aep)) deallocate (aep)
      allocate (aep(nbloc))
      if (allocated(aev)) deallocate (aev)
      allocate (aev(nbloc))
      if (allocated(aeb)) deallocate (aeb)
      allocate (aeb(nbloc))
      if (allocated(aea)) deallocate (aea)
      allocate (aea(nbloc))
      if (allocated(aeba)) deallocate (aeba)
      allocate (aeba(nbloc))
      if (allocated(aub)) deallocate (aub)
      allocate (aub(nbloc))
      if (allocated(aeaa)) deallocate (aeaa)
      allocate (aeaa(nbloc))
      if (allocated(aeopb)) deallocate (aeopb)
      allocate (aeopb(nbloc))
      if (allocated(aeopd)) deallocate (aeopd)
      allocate (aeopd(nbloc))
      if (allocated(aeid)) deallocate (aeid)
      allocate (aeid(nbloc))
      if (allocated(aeit)) deallocate (aeit)
      allocate (aeit(nbloc))
      if (allocated(aet)) deallocate (aet)
      allocate (aet(nbloc))
      if (allocated(aept)) deallocate (aept)
      allocate (aept(nbloc))
      if (allocated(aebt)) deallocate (aebt)
      allocate (aebt(nbloc))
      if (allocated(aett)) deallocate (aett)
      allocate (aett(nbloc))
      if (allocated(aeg)) deallocate (aeg)
      allocate (aeg(nbloc))
      if (allocated(aex)) deallocate (aex)
      allocate (aex(nbloc))
      if (allocated(aerep)) deallocate (aerep)
      allocate (aerep(nbloc))
      if (allocated(aexdisp)) deallocate (aexdisp)
      allocate (aexdisp(nbloc))
      if (allocated(aect)) deallocate (aect)
      allocate (aect(nbloc))
      if (allocated(aesum)) deallocate (aesum)
      allocate (aesum(nbloc))
c
c     zero out each of the potential energy components
c
      ec = 0.0d0
      em = 0.0d0
      ep = 0.0d0
      eb = 0.0d0
      ev = 0.0d0
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
      eg = 0.0d0
      erep = 0.0d0
      exdisp = 0.0d0
      ect = 0.0d0
      ex = 0.0d0
c
c     zero out energy partitioning components for each atom
c
      aec = 0.0d0
      aem = 0.0d0
      aep = 0.0d0
      aev = 0.0d0
      aea = 0.0d0
      aeb = 0.0d0
      aub = 0.0d0
      aeba = 0.0d0
      aeaa = 0.0d0
      aeopb = 0.0d0
      aeopd = 0.0d0
      aeid = 0.0d0
      aeit = 0.0d0
      aet = 0.0d0
      aept = 0.0d0
      aebt = 0.0d0
      aett = 0.0d0
      aeg = 0.0d0
      aerep = 0.0d0
      aexdisp = 0.0d0
      aect = 0.0d0
      aex = 0.d0
c
c     zero out the total intermolecular energy
c
      einter = 0.0d0
c
c     call the local geometry energy component routines
c
      if (use_bond)  call ebond3
      if (use_angle)  call eangle3
      if (use_strbnd)  call estrbnd3
      if (use_urey)  call eurey3
      if (use_angang)  call eangang3
      if (use_opbend)  call eopbend3
      if (use_opdist)  call eopdist3
      if (use_improp)  call eimprop3
      if (use_imptor)  call eimptor3
      if (use_tors)  call etors3
      if (use_pitors)  call epitors3
      if (use_strtor)  call estrtor3
      if (use_tortor)  call etortor3
c
c     call the van der Waals energy component routines
c
      if (use_vdw) then
         if (vdwtyp .eq. 'LENNARD-JONES')  call elj3
         if (vdwtyp .eq. 'BUFFERED-14-7')  call ehal3
      end if
c
c     call the electrostatic energy component routines
c
      if (use_charge) call echarge3
      if (use_mpole)  call empole3
      if (use_polar)  call epolar3
c
      if (use_repulsion) call erepulsion3
      if (use_dispersion) call edispersion3
      if (use_ctransfer) call ectransfer3
c
c     call any miscellaneous energy component routines
c
      if (use_geom)  call egeom3
      if (use_extra)  call extra3
c
c     MPI : get total energy
c
      if (rank.eq.0) then
        call MPI_REDUCE(MPI_IN_PLACE,ec,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,nec,1,MPI_INT,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,em,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,nem,1,MPI_INT,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,ep,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,nep,1,MPI_INT,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,ev,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,nev,1,MPI_INT,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,eb,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,neb,1,MPI_INT,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,ea,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,nea,1,MPI_INT,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,eba,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,neba,1,MPI_INT,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,eub,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,neub,1,MPI_INT,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,eaa,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,neaa,1,MPI_INT,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,eopb,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,neopb,1,MPI_INT,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,eopd,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,neopd,1,MPI_INT,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,eid,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,neid,1,MPI_INT,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,eit,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,neit,1,MPI_INT,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,et,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,net,1,MPI_INT,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,ept,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,nept,1,MPI_INT,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,ebt,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,nebt,1,MPI_INT,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,ett,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,nett,1,MPI_INT,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,eg,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,neg,1,MPI_INT,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,erep,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,nerep,1,MPI_INT,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,exdisp,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,nexdisp,1,MPI_INT,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,ect,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,nect,1,MPI_INT,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,ex,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,nex,1,MPI_INT,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,einter,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_BEAD,ierr)
      else
        call MPI_REDUCE(ec,ec,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(nec,nec,1,MPI_INT,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(em,em,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(nem,nem,1,MPI_INT,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(ep,ep,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(nep,nep,1,MPI_INT,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(ev,ev,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(nev,nev,1,MPI_INT,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(eb,eb,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(neb,neb,1,MPI_INT,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(ea,ea,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(nea,nea,1,MPI_INT,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(eba,eba,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(neba,neba,1,MPI_INT,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(eub,eub,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(neub,neub,1,MPI_INT,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(eaa,eaa,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(neaa,neaa,1,MPI_INT,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(eopb,eopb,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(neopb,neopb,1,MPI_INT,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(eopd,eopd,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(neopd,neopd,1,MPI_INT,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(eid,eid,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(neid,neid,1,MPI_INT,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(eit,eit,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(neit,neit,1,MPI_INT,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(et,et,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(net,net,1,MPI_INT,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(ept,ept,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(nept,nept,1,MPI_INT,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(ebt,ebt,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(nebt,nebt,1,MPI_INT,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(ett,ett,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(nett,nett,1,MPI_INT,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(eg,eg,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(neg,neg,1,MPI_INT,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(erep,erep,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(nerep,nerep,1,MPI_INT,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(exdisp,exdisp,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(nexdisp,nexdisp,1,MPI_INT,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(ect,ect,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(nect,nect,1,MPI_INT,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(ex,ex,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(nex,nex,1,MPI_INT,MPI_SUM,0,
     $     COMM_BEAD,ierr)
        call MPI_REDUCE(einter,einter,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_BEAD,ierr)
      end if
c
c     sum up to give the total potential energy
c
      esum = eb + ea + eba + eub + eaa + eopb + eopd + eid + eit
     &          + et + ept + ebt + ett + ev + em
     &          + ec + ep +  eg + erep + exdisp + ect + ex
      energy = esum
c
c     sum up to give the total potential energy per atom
c
      aesum = aem + aec + aep + aev + aeb + aea + aeba + aub + aeaa +
     $ aeopb + aeopd + aeid + aeit + aet + aept + aebt + aett + aeg
     $ + aerep + aexdisp + aect + aex 
c
c     check for an illegal value for the total energy
c
      if (isnan(esum)) then
         write (iout,10)
   10    format (/,' ANALYSIS  --  Illegal Value for the Total',
     &              ' Potential Energy')
         call fatal
      end if
      return
      end
