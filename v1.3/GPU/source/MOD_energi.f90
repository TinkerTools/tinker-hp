!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #################################################################
!     ##                                                             ##
!     ##  module energi  --  individual potential energy components  ##
!     ##                                                             ##
!     #################################################################
!
!
!     calc_e switch to compute energy
!     esum   total potential energy of the system
!     eb     bond stretch potential energy of the system
!     ea     angle bend potential energy of the system
!     emlpot Machine learning potential energy of the system
!     eba    stretch-bend potential energy of the system
!     eub    Urey-Bradley potential energy of the system
!     eaa    angle-angle potential energy of the system
!     eopb   out-of-plane bend potential energy of the system
!     eopd   out-of-plane distance potential energy of the system
!     eid    improper dihedral potential energy of the system
!     eit    improper torsion potential energy of the system
!     et     torsional potential energy of the system
!     ept    pi-orbital torsion potential energy of the system
!     eat    angle-torsion potential energy of the system
!     ebt    stretch-torsion potential energy of the system
!     ett    torsion-torsion potential energy of the system
!     ev     van der Waals potential energy of the system
!     er     Pauli repulsion potential energy of the system
!     edsp   dispersion potential energy of the system
!     ec     charge-charge potential energy of the system
!     ecrec  reciproqual charge-charge potential energy of the system
!     em     atomic multipole potential energy of the system
!     emrec  reciproqual part of atomic multipole potential energy of the system
!     ep     polarization potential energy of the system
!     ect    charge transfer potential energy of the system
!     eprec  reciproqual part of polarization potential energy of the system
!     eg     geometric restraint potential energy of the system
!     ex     extra term potential energy of the system
!     esave  stored potential energy of the system
!     ensmd  extra term smd potential energy of the system
!     eDaMD  extra term aMD potential for dihedrals of the system
!     ePaMD  extra term aMD potential energy of the system
!     eW1aMD extra term aMD potential for bonds/angles of waters
!     eW2aMD extra term aMD potential energy of waters
!
!
#include "tinker_macro.h"
module energi
   implicit none
   logical   calc_e
   real(r_p) esum,eb,ea,eba
   real(r_p) emlpot
   real(r_p) eub,eaa,eopb,eopd
   real(r_p) eid,eit,et,ept
   real(r_p) eat,ebt,ett
   real(r_p) ev,ec,ecrec,em,emrec,ep,eprec,edsprec
   ener_rtyp er,edsp,ect
   real(r_p) eg,ex
   real(r_p) esave
   real(r_p) ensmd
   real(r_p) eDaMD,ePaMD,eW1aMD,eW2aMD
   ener_rtyp ev_r,ec_r,em_r,ep_r,eb_r
   real(r_p) epot_std10,epot_mean
   real(r_p) etot_std,etot_ave

   interface
      module subroutine create_energi_on_device()
      end subroutine
      module subroutine delete_energi_on_device()
      end subroutine
      module subroutine chk_energy_fluct(epot,ealt,abort)
         logical abort
         real(r_p) epot,ealt
      end subroutine
      module subroutine info_energy(rank)
         integer,intent(in)::rank
      end subroutine
      module subroutine allreduceen(epot)
         real(r_p), intent(inout) :: epot
      end subroutine allreduceen
      module subroutine reduceen(epot)
         real(r_p), intent(inout) :: epot
      end subroutine reduceen
   end interface

end module

submodule(energi) subEnergi
   logical ::chk_en=.true.
contains
#include "convert.inc.f90"

   module subroutine create_energi_on_device()
      implicit none
!$acc enter data create(esum,eb,ea,eba,eub,eaa,eopb,eopd, &
!$acc                  eid,eit,et,ept,ebt,ett,eg,ex,eat, &
!$acc                  ev,ec,ecrec,em,emrec,ep,eprec,esave,ensmd, &
!$acc                  er,edsp,edsprec,ect, &
!$acc                  eDaMD,ePaMD,eW1aMD,eW2aMD, &
!$acc                  emlpot, &
!$acc                  ev_r,ec_r,em_r,ep_r,eb_r)
   end subroutine
   module subroutine delete_energi_on_device()
      implicit none
!$acc exit data delete(esum,eb,ea,eba,eub,eaa,eopb,eopd, &
!$acc                  eid,eit,et,ept,ebt,ett,eg,ex,eat, &
!$acc                  ev,ec,ecrec,em,emrec,ep,eprec,esave,ensmd, &
!$acc                  er,edsp,edsprec,ect, &
!$acc                  eDaMD,ePaMD,eW1aMD,eW2aMD, &
!$acc                  emlpot, &
!$acc                  ev_r,ec_r,em_r,ep_r,eb_r)
   end subroutine

   module subroutine allreduceen(epot)
      use domdec, only: COMM_TINKER, nproc, ndir
      use inter
      use mpi
      use parallel, only: reduce_sum, reduction_all
      use potent, only: use_pmecore
      implicit none
      real(r_p), intent(inout) :: epot
      real(r_p) epot_buffer(1)
      real(r_p) buffer(22)

      if (nproc.eq.1.or.use_pmecore.and.ndir.eq.1) return

!$acc data create(buffer,epot_buffer)
!$acc serial
      buffer = [eba,ea,eb,ec,em,ep,ev,eub,eaa,eopb,eopd,eid,eit,et,ept,ebt,eat,ett,esum,eg,ex,einter]
      epot_buffer(1) = epot
!$acc end serial
!$acc host_data use_device(buffer,epot_buffer)
      call reduce_sum(buffer, reduction_all)
      call reduce_sum(epot_buffer, reduction_all)
!$acc end host_data
!$acc serial
      eba=buffer(1); ea=buffer(2); eb=buffer(3); ec=buffer(4); em=buffer(5); ep=buffer(6)
      ev=buffer(7); eub=buffer(8); eaa=buffer(9); eopb=buffer(10); eopd=buffer(11); eid=buffer(12)
      eit=buffer(13); et=buffer(14); ept=buffer(15); ebt=buffer(16); eat=buffer(17); ett=buffer(18)
      esum=buffer(19); eg=buffer(20); ex=buffer(21); einter=buffer(22); epot=epot_buffer(1)
!$acc end serial
!$acc end data
   end subroutine allreduceen

   module subroutine info_energy(rank)
      implicit none
      integer,intent(in):: rank
      real(r_p) ebonded

      if (rank.eq.0) then
!$acc wait
20       format ( 40('-'))
         print 20
30       format (1x,A,F18.6)
!$acc update host(esum,eb,ea,eba,eub,eaa,eopb,eopd, &
!$acc                  eid,eit,et,ept,ebt,ett,eg,ex,eat, &
!$acc                  ev,ec,ecrec,em,emrec,ep,eprec,esave,ensmd, &
!$acc                  er,edsp,edsprec,ect, &
!$acc                  eDaMD,ePaMD,eW1aMD,eW2aMD, &
!$acc                  emlpot, &
!$acc                  ev_r,ec_r,em_r,ep_r,eb_r)

         ebonded =&
            &eb+ea+eba+eub+eaa+eid+eit+et+ept+ebt+ett+eat+eopb+eopd+eg+ex
         if (eb   /=real(0,r_p)) print 30, 'eb     = ',eb
         !if (eb_r /=tp2enr(0.0)) print 30, 'eb_r = ',enr2en(eb_r)
         if (ea   /=real(0,r_p)) print 30, 'ea     = ',ea
         if (eba  /=real(0,r_p)) print 30, 'eba    = ',eba
         if (eub  /=real(0,r_p)) print 30, 'eub    = ',eub
         if (eaa  /=real(0,r_p)) print 30, 'eaa    = ',eaa
         if (eid  /=real(0,r_p)) print 30, 'eid    = ',eid
         if (eit  /=real(0,r_p)) print 30, 'eit    = ',eit
         if (et   /=real(0,r_p)) print 30, 'et     = ',et
         if (ept  /=real(0,r_p)) print 30, 'ept    = ',ept
         if (ebt  /=real(0,r_p)) print 30, 'ebt    = ',ebt
         if (ett  /=real(0,r_p)) print 30, 'ett    = ',ett
         if (eat  /=real(0,r_p)) print 30, 'eat    = ',eat
         if (eopb /=real(0,r_p)) print 30, 'eopb   = ',eopb
         if (eopd /=real(0,r_p)) print 30, 'eopd   = ',eopd
         if (eg   /=real(0,r_p)) print 30, 'eg     = ',eg
         if (ex   /=real(0,r_p)) print 30, 'ex     = ',ex
         if (emlpot /=real(0,r_p)) print 30, 'emlpot  =',emlpot
         if (ebonded/=real(0,r_p)) print 30, 'ebonded =',ebonded
         if (ec   /=real(0,r_p)) print 30, 'ec     = ',ec
         if (ev   /=real(0,r_p)) print 30, 'ev     = ',ev
         if (er   /=     0     ) print 30, 'er     = ',enr2en(er)
         if (edsp /=     0     ) print 30, 'edsp   = ',enr2en(edsp)
         if (em   /=real(0,r_p)) print 30, 'em     = ',em
         if (ep   /=real(0,r_p)) print 30, 'ep     = ',ep
         if (ect  /=     0     ) print 30, 'ect    = ',enr2en(ect)
         if (ensmd/=real(0,r_p)) print 30, 'ensmd  = ',ensmd
         if (eDaMD/=real(0,r_p)) print 30, 'eDaMD  = ',eDaMD
         if (ePaMD/=real(0,r_p)) print 30, 'ePaMD  = ',ePaMD
         if(eW1aMD/=real(0,r_p)) print 30,'eW1aMD  = ',eW1aMD
         !if (esum/=real(0,r_p)) print 30, 'esum    = ',esum
      end if
!$acc wait
   end subroutine

   module subroutine chk_energy_fluct(epot,ealt,abort)
      implicit none
      logical abort
      real(r_p) epot,ealt,edist

      if (emlpot /=real(0,r_p)) return

      if (chk_en) then
!$acc wait
!$acc update host(epot,ealt)
         epot_mean  = epot
         epot_std10 = abs(epot)
         chk_en     = .false.
      end if
!$acc serial async present(esum,eb,ea,eba,eub,eaa,eopb,eopd, &
!$acc             eid,eit,et,ept,ebt,ett,eg,ex,eat, &
!$acc             ev,ec,ecrec,em,emrec,ep,eprec,esave,ensmd, &
!$acc             er,edsp,edsprec,ect,epot,ealt &
!$acc            )   copy(abort)
      edist = abs(epot-epot_mean)
      if (edist.gt.epot_std10) then
         print*, ' --- Tinker-HP : '&
            &,'detected abnormal fluctuation on energy'
         print*, 'epot/ average/ 10*std',epot,epot_mean,epot_std10
         if (eb+ea  /=0.0) print*, 'eb/ea  = ',eb,ea
         if (eba+eub/=0.0) print*, 'eba/eub= ',eba,eub
         if (ealt   /=0.0) print*, 'ealt   = ',ealt
!        if (eaa  /=real(0,r_p)) print*, 'eaa    = ',eaa
!        if (eid+eit/=0.0) print*, 'eid/eit  = ',eid,eit
         if (et+ept /=0.0) print*, 'et/ept = ',et,ept
!        if (ebt  /=real(0,r_p)) print*, 'ebt    = ',ebt
!        if (eat  /=real(0,r_p)) print*, 'eat    = ',ett
!        if (ett  /=real(0,r_p)) print*, 'ett    = ',ett
!        if (eopb+eopd /=0.0) print*, 'eopb/eopd   = ',eopb,eopd
!        if (eg   /=real(0,r_p)) print*, 'eg     = ',eg
!        if (ex   /=real(0,r_p)) print*, 'ex     = ',eg
         if (ev+ec  /=0.0) print*, 'ev/ec  = ',ev,ec
         if (em+ep  /=0.0) print*, 'em+ep  = ',em+ep
         !if (ensmd/=real(0,r_p)) print*, 'ensmd  = ',ensmd
         if (edist.gt.5*epot_std10) abort = .true.
      end if
!$acc end serial
   end subroutine

!     subroutine reduceen : make the necessary reductions over the processes to get the
!     total values of the energies
!
module subroutine reduceen(epot)
   use domdec
   use inter
   use inform    ,only: deb_Energy
   use mdevice   ,only: reciprocal_work,sync_work
   use parallel  ,only: reduce_sum, reduction_all, reduction_root
   use potent
   use virial
   use mpi
   use timestat
   use tinheader
   use tinmemory ,only: mem_move
   implicit none
   integer ierr,commloc
   integer(8) nine
   real(r_p), intent(inout) :: epot
   real(r_p) buffer1(10),buffer2(20)
   real(r_p) vir_copy(3,3)
   parameter(nine=9)
!
   if (ndir.eq.1) return

   call timer_enter(timer_reduceen)
   if (use_pmecore) then
      commloc = comm_dir
   else
      commloc = COMM_TINKER
   end if
!$acc data create(buffer1,buffer2,vir_copy) &
!$acc     present_or_copy(epot) &
!$acc     present(eb,eba,eub,eopb,et,ept,ett,ebt,ea &
!$acc           ,eaa,eopd,eid,eit,ec,ev,er,edsp,em,ep,ect,eg,ex &
!$acc           ,esum,vir)
!
!$acc wait(reciprocal_work%queue)
!$acc serial
   buffer1(1) = ec
   buffer1(2) = em
   buffer1(3) = ep
   buffer1(4) = epot
   buffer1(5) = esum
   buffer1(6) = ev
   buffer1(7) = enr2en(edsp)
   buffer1(8) = enr2en(er )
   buffer1(9) = enr2en(ect)
!$acc end serial
!
!$acc host_data use_device(buffer1)
   call reduce_sum(buffer1(1:9), reduction_root)
!$acc end host_data
!
!$acc serial
   ec     = buffer1(1)
   em     = buffer1(2)
   ep     = buffer1(3)
   epot   = buffer1(4)
   esum   = buffer1(5)
   ev     = buffer1(6)
   edsp   = rp2enr(buffer1(7))
   er     = rp2enr(buffer1(8))
   ect    = rp2enr(buffer1(9))
!$acc end serial
!
!   MPI: get virial
!
   if (use_virial) then
      call mem_move(vir_copy(1:nine,1),vir(1:nine,1),nine,sync_work%stream)
!
!$acc host_data use_device(vir_copy)
      call reduce_sum(vir_copy, reduction_all)
!$acc end host_data
!
      call mem_move(vir(1:nine,1),vir_copy(1:nine,1),nine,sync_work%stream)
   end if
!
   if ((use_pmecore).and.(rank.gt.(ndir-1))) goto 10
   if (.not.deb_Energy) goto 10
!
!$acc serial
   buffer2(1)  = eba
   buffer2(2)  = ea
   buffer2(3)  = eb
   buffer2(4)  = einter
   buffer2(5)  = eub
   buffer2(6)  = eaa
   buffer2(7)  = eopb
   buffer2(8)  = eopd
   buffer2(9)  = eid
   buffer2(10) = eit
   buffer2(11) = et
   buffer2(12) = ept
   buffer2(13) = ebt
   buffer2(14) = ett
   buffer2(15) = eg
   buffer2(16) = ex
   buffer2(17) = eat
!$acc end serial
!
!$acc host_data use_device(buffer2)
   ! The debug energy subset follows the direct subgroup when PME cores are split.
   ! Keep this routing adapter until the domain privately owns that subgroup.
   if (rank.eq.0) then
      call MPI_REDUCE(MPI_IN_PLACE,buffer2,20,MPI_RPREC,MPI_SUM,0,commloc,ierr)
   else
      call MPI_REDUCE(buffer2,buffer2,20,MPI_RPREC,MPI_SUM,0,commloc,ierr)
   end if
!$acc end host_data
!
!$acc serial
   eba   = buffer2(1)
   ea    = buffer2(2)
   eb    = buffer2(3)
   einter= buffer2(4)
   eub   = buffer2(5)
   eaa   = buffer2(6)
   eopb  = buffer2(7)
   eopd  = buffer2(8)
   eid   = buffer2(9)
   eit   = buffer2(10)
   et    = buffer2(11)
   ept   = buffer2(12)
   ebt   = buffer2(13)
   ett   = buffer2(14)
   eg    = buffer2(15)
   ex    = buffer2(16)
   eat   = buffer2(17)
!$acc end serial

10 continue
!$acc end data
   call timer_exit(timer_reduceen,quiet_timers)
end subroutine reduceen

end submodule
