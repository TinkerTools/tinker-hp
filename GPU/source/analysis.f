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
#include "tinker_precision.h"

      module analysis_inl
      contains
#include "convert.f.inc"
      end module

      subroutine analysis (energy)
      use action
      use analyz
      use analysis_inl
      use domdec
      use energi
      use iounit
      use inter
      use potent
      use vdwpot
      use timestat
      use tinheader
      use mpi
      use sizes
      implicit none
      integer ierr
      real(r_p) energy
      logical tinker_isnan_m
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
      if (allocated(aesum)) deallocate (aesum)
      allocate (aesum(nbloc))
!$acc enter data create(aec,aea,aeb,aub,aeopb,aeg,aeba,
!$acc&      aet,aept,aebt,aett) async
      call create_action_data_ondevice

 20   format(a,f12.9)
c
c     zero out each of the potential energy components
c
!$acc data present(eb,eba,eub,eopb,et,ept,ett,ebt,ea
!$acc&      ,eaa,eopd,eid,eit,ec,ev,em,ep,eg,ex,esum
!$acc&      ,ev_r,ec_r,em_r,ep_r,eb_r
!$acc&      ,emrec,eprec,nev,nec,nem,nep,nem_,nep_,nev_)

!$acc serial async
      ec    = 0.0_re_p
      ec_r  = 0
      em    = 0.0_re_p
      emrec = 0.0_re_p
      em_r  = 0
      ep    = 0.0_re_p
      eprec = 0.0_re_p
      ep_r  = 0
      eb    = 0.0_re_p
      ev    = 0.0_re_p
      ea    = 0.0_re_p
      eba   = 0.0_re_p
      eub   = 0.0_re_p
      eaa   = 0.0_re_p
      eopb  = 0.0_re_p
      eopd  = 0.0_re_p
      eid   = 0.0_re_p
      eit   = 0.0_re_p
      et    = 0.0_re_p
      ept   = 0.0_re_p
      ebt   = 0.0_re_p
      ett   = 0.0_re_p
      eg    = 0.0_re_p
      ev_r  = 0
      nec   = 0
      nev   = 0
      nem   = 0
      nep   = 0
      nev_  = 0.0
      nem_  = 0.0
      nep_  = 0.0
      nec_  = 0.0
!$acc end serial
c
c     zero out energy partitioning components for each atom
c

!$acc kernels async default(present)
      aec   = 0.0_ti_p
      aea   = 0.0_ti_p
      aeb   = 0.0_ti_p
      aub   = 0.0_ti_p
      aeba  = 0.0_ti_p
      aeopb = 0.0_ti_p
      aet   = 0.0_ti_p
      aept  = 0.0_ti_p
      aebt  = 0.0_ti_p
      aett  = 0.0_ti_p
      aeg   = 0.0_ti_p
!$acc end kernels
      aeaa  = 0.0_ti_p
      aeopd = 0.0_ti_p
      aeid  = 0.0_ti_p
      aeit  = 0.0_ti_p
      aem   = 0.0_ti_p
      aep   = 0.0_ti_p
      aev   = 0.0_ti_p
c
c     zero out the total intermolecular energy
c
      einter = 0.0_re_p
c
c     call the local geometry energy component routines
c
c     if (use_bond)    call ebond3
      if (use_bond)    call ebond3gpu
c     if (use_angle)   call eangle3
      if (use_bond)    call eangle3gpu
      if (use_strbnd)  call estrbnd3
c     if (use_urey)    call eurey3
      if (use_urey)    call eurey3gpu
      if (use_angang)  call eangang
c     if (use_opbend)  call eopbend3
      if (use_opbend)  call eopbend3gpu
      if (use_opdist)  call eopdist3
      if (use_improp)  call eimprop3
      if (use_imptor)  call eimptor3
      if (use_tors)    call etors3
      if (use_pitors)  call epitors3
      if (use_strtor)  call estrtor3
      if (use_tortor)  call etortor3
c
c     call the van der Waals energy component routines
c
      if (use_vdw) then
         if (vdwtyp .eq. 'LENNARD-JONES')  call elj3gpu
         call timer_enter( timer_ehal3 )
c        if (vdwtyp .eq. 'BUFFERED-14-7')  call ehal3
c        if (vdwtyp .eq. 'BUFFERED-14-7')  call ehal3vec
         if (vdwtyp .eq. 'BUFFERED-14-7')  call ehal3gpu
         call timer_exit( timer_ehal3 )
      end if
c
c     call the electrostatic energy component routines
c
      call MPI_BARRIER(hostcomm,ierr)
      if (use_charge) call echarge3gpu

      if (use_mpole) call timer_enter( timer_empole3 )
      if (use_mpole) call empole3gpu
      if (use_mpole) call timer_exit( timer_empole3 )

      if (use_polar) call timer_enter( timer_polar )
      if (use_polar) call epolar3gpu
      if (use_polar) call timer_exit( timer_polar )

c
c     call any miscellaneous energy component routines
c
      if (use_geom)   call egeom3gpu
      if (use_extra)  call extra3
c
c     Update data on host
c
!$acc wait
!$acc update host(eb,eba,eub,eopb,et,ept,ett,ebt,ea
!$acc&     ,eaa,eopd,eid,eit,ec,ev,em,ep,eg,ex,esum
!$acc&     ,eb_r,ev_r,ev_r,em_r,ep_r
!$acc&     ,nec,nev,nep,nem)

      ! get reducted contribution
      eb = eb + enr2en(eb_r)
      ev = ev + enr2en(ev_r)
c
c     MPI : get total energy
c
      if (rank.eq.0) then
        call MPI_REDUCE(MPI_IN_PLACE,ec,1,MPI_RPREC,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,nec,1,MPI_INT,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,em,1,MPI_RPREC,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,nem,1,MPI_INT,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,ep,1,MPI_RPREC,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,nep,1,MPI_INT,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,ev,1,MPI_RPREC,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,nev,1,MPI_INT,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,eb,1,MPI_RPREC,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,neb,1,MPI_INT,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,ea,1,MPI_RPREC,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,nea,1,MPI_INT,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,eba,1,MPI_RPREC,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,neba,1,MPI_INT,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,eub,1,MPI_RPREC,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,neub,1,MPI_INT,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,eaa,1,MPI_RPREC,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,neaa,1,MPI_INT,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,eopb,1,MPI_RPREC,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,neopb,1,MPI_INT,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,eopd,1,MPI_RPREC,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,neopd,1,MPI_INT,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,eid,1,MPI_RPREC,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,neid,1,MPI_INT,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,eit,1,MPI_RPREC,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,neit,1,MPI_INT,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,et,1,MPI_RPREC,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,net,1,MPI_INT,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,ept,1,MPI_RPREC,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,nept,1,MPI_INT,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,ebt,1,MPI_RPREC,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,nebt,1,MPI_INT,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,ett,1,MPI_RPREC,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,nett,1,MPI_INT,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,eg,1,MPI_RPREC,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,neg,1,MPI_INT,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,ex,1,MPI_RPREC,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,nex,1,MPI_INT,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,einter,1,MPI_RPREC,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
      else
        call MPI_REDUCE(ec,ec,1,MPI_RPREC,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(nec,nec,1,MPI_INT,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(em,em,1,MPI_RPREC,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(nem,nem,1,MPI_INT,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(ep,ep,1,MPI_RPREC,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(nep,nep,1,MPI_INT,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(ev,ev,1,MPI_RPREC,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(nev,nev,1,MPI_INT,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(eb,eb,1,MPI_RPREC,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(neb,neb,1,MPI_INT,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(ea,ea,1,MPI_RPREC,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(nea,nea,1,MPI_INT,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(eba,eba,1,MPI_RPREC,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(neba,neba,1,MPI_INT,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(eub,eub,1,MPI_RPREC,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(neub,neub,1,MPI_INT,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(eaa,eaa,1,MPI_RPREC,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(neaa,neaa,1,MPI_INT,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(eopb,eopb,1,MPI_RPREC,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(neopb,neopb,1,MPI_INT,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(eopd,eopd,1,MPI_RPREC,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(neopd,neopd,1,MPI_INT,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(eid,eid,1,MPI_RPREC,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(neid,neid,1,MPI_INT,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(eit,eit,1,MPI_RPREC,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(neit,neit,1,MPI_INT,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(et,et,1,MPI_RPREC,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(net,net,1,MPI_INT,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(ept,ept,1,MPI_RPREC,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(nept,nept,1,MPI_INT,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(ebt,ebt,1,MPI_RPREC,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(nebt,nebt,1,MPI_INT,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(ett,ett,1,MPI_RPREC,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(nett,nett,1,MPI_INT,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(eg,eg,1,MPI_RPREC,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(neg,neg,1,MPI_INT,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(ex,ex,1,MPI_RPREC,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(nex,nex,1,MPI_INT,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(einter,einter,1,MPI_RPREC,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
      end if
c
c     sum up to give the total potential energy
c
      esum = eb + ea + eba + eub + eaa + eopb + eopd + eid + eit
     &          + et + ept + ebt + ett + ev + em
     &          + ec + ep +  eg + ex
      energy = esum
c
c     sum up to give the total potential energy per atom
c
      aesum = aem + aec + aep + aev + aeb + aea + aeba + aub + aeaa +
     $ aeopb + aeopd + aeid + aeit + aet + aept + aebt + aett + aeg
     $ + aex

!$acc end data
c
      call delete_action_data_ondevice
!$acc exit data delete(aec,aea,aeb,aub,aeopb,aeg,aeba,
!$acc&      aet,aept,aebt,aett) async
c
c     check for an illegal value for the total energy
c
      if (tinker_isnan_m(esum)) then
         write (iout,10)
   10    format (/,' ANALYSIS  --  Illegal Value for the Total',
     &              ' Potential Energy')
         call fatal
      end if
      end
