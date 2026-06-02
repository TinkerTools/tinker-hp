!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #############################################################
!     ##                                                         ##
!     ##  program analyze  --  energy partitioning and analysis  ##
!     ##                                                         ##
!     #############################################################
!
!
!     "analyze" computes and displays the total potential energy;
!
!
#include "tinker_macro.h"
module analyze_inl
contains
#include "convert.inc.f90"
end module

program analyze
   use mpi
   implicit none
   integer ierr!,nthreadsupport
   call MPI_INIT(ierr)
!      call MPI_INIT_THREAD(MPI_THREAD_MULTIPLE,nthreadsupport,ierr)
   call analyze_bis
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   call MPI_FINALIZE(ierr)
end
!
subroutine analyze_bis
   use dcdmod
   use domdec
   use files
   use inform
   use iounit
   use mpi
   implicit none
   integer i,ixyz,ierr
   integer frame
   integer freeunit
   integer trimtext
   real(r_p) energy
   logical doenergy,dodipoltot,dodipolmol
   logical exist
   character*1 letter
   character*240 string
   character*240 xyzfile
   character*240 dcdfile
   type(dcdinfo_t) :: dcdinfo
!
   ! Sign running program
   app_id = analyze_a
!
!     set up the structure and mechanics calculation
!
   call initial
   call initmpi
   call getxyz
   call unitcell
   call cutoffs
   call lattice
!
!     setup for MPI
!
   call drivermpi
   call reinitnl(0)
!
   call mechanic
!
!     call nblist(0)
!
!     get the desired types of analysis to be performed
!
   call nextarg (string,exist)
   if (.not. exist) then
      if (ranktot.eq.0) write (iout,10)
10    format (/,' The TINKER Analysis Facility can Provide :',&
              /,'   Total Potential Energy and its Components [E]',&
              /,'   Total Dipolar Moment [D]',&
              /,'   Molecular Dipolar Moments [M]')
      if (ranktot.eq.0) write (iout,30)
30    format (/,' You Need To Enter the Desired Analysis Types',&
         &' [E,D,M] :  ')
      call MPI_BARRIER(COMM_TINKER,ierr)
      __TINKER_FATAL__
   end if
!
!     set option control flags based desired analysis types
!
   doenergy   = .false.
   dodipoltot = .false.
   dodipolmol = .false.
   call upcase (string)
   do i = 1, trimtext(string)
      letter = string(i:i)
      if (letter .eq. 'E')  then
        doenergy   = .true.
      else if (letter .eq. 'D')  then
        dodipoltot = .true.
      else if (letter .eq. 'M')  then
        dodipolmol = .true.
      else
        if (ranktot.eq.0) write (iout,30)
        call MPI_BARRIER(COMM_TINKER,ierr)
        __TINKER_FATAL__
      end if
   end do
!
!     reopen the coordinates file and read the first structure
!
   frame = 0
   if (dcdio) then
      dcdfile = filename(1:leng)//'.dcd'
      call dcdfile_open(dcdinfo,dcdfile)
      call dcdfile_read_header(dcdinfo,.false.)
      call dcdfile_read_next(dcdinfo)
      call dcdfile_skip_next(dcdinfo,0)
   else
      ixyz = freeunit ()
      xyzfile = filename
      call suffix (xyzfile,'xyz','old')
      open (unit=ixyz,file=xyzfile,status ='old')
      rewind (unit=ixyz)
      call readxyz (ixyz)
   end if
!
!     perform analysis for each successive coordinate structure
!
   do while (.not. abort)
      frame = frame + 1
      if (frame .gt. 1) then
         if (rank.eq.0) write (iout,90)  frame
90       format (/,' Analysis for Archive Structure :',8x,i8)
      end if
!
!       setup for MPI
!
      call lattice
      call AllDirAssign
      call AllRecAssign
      call reinitnl(0)
      call mechanicstep(0)
      call nblist(0)
!
!     make the call to compute the potential energy
!
      if (doenergy)  call enrgyze
!
!     make the call to compute the induced dipoles and the total dipolar moments
!
      if (dodipoltot)  then
         call analysis(energy)
         call totaldipole
      end if
      if (dodipolmol)  then
         call analysis(energy)
         call dipolemol
      end if
!
!     energy partitioning by potential energy components
!
      if (doenergy) then
         if (rank.eq.0) call partyze
      end if
!
!     attempt to read next structure from the coordinate file
!
      if (.not.dcdio) then
         call readxyz (ixyz)
      else
         call dcdfile_read_next(dcdinfo)
         call dcdfile_skip_next(dcdinfo,0)
      end if
   end do
!
!     perform any final tasks before program exit
!
   if (dcdio) then
      call dcdfile_close(dcdinfo)
   else
      close (unit=ixyz)
   end if
   call final
end
!
!
!     ################################################################
!     ##                                                            ##
!     ##  subroutine enrgyze  --  compute & report energy analysis  ##
!     ##                                                            ##
!     ################################################################
!
!
!     "enrgyze" is an auxiliary routine for the analyze program
!     that performs the energy analysis and prints the total and
!     intermolecular energies
!
!
subroutine enrgyze
   use atoms
   use cutoff
   use domdec
   use inform
   use inter
   use iounit
   use molcul
   implicit none
   real(r_p) energy
   character*240 fstr
!
!
!     perform the energy analysis by atom and component
!
   call analysis (energy)
   if (rank.eq.0) then
!
!       print out the total potential energy of the system
!
      fstr='(/,'' Total Potential Energy :'',8x,f16.4,'' Kcal/mole'')'
      if (digits .ge. 6)  fstr(32:39) = '6x,f18.6'
      if (digits .ge. 8)  fstr(32:39) = '4x,f20.8'
      if (abs(energy) .ge. 1.0d10)  fstr(35:35) = 'd'
      write (iout,fstr)  energy
!
!       intermolecular energy for systems with multiple molecules
!
      fstr='(/,'' Intermolecular Energy :'',9x,f16.4,'' Kcal/mole'')'
      if (digits .ge. 6)  fstr(31:38) = '7x,f18.6'
      if (digits .ge. 8)  fstr(31:38) = '5x,f20.8'
      if (abs(einter) .ge. 1.0d10)  fstr(34:34) = 'd'
      if (nmol.gt.1 .and. nmol.lt.n .and. .not.use_ewald)&
         &write (iout,fstr)  einter
   end if
   return
end
!
!
!     ##############################################################
!     ##                                                          ##
!     ##  subroutine partyze  --  energy component decomposition  ##
!     ##                                                          ##
!     ##############################################################
!
!
!     "partyze" prints the energy component and number of
!     interactions for each of the potential energy terms
!
!
subroutine partyze
   use action
   use analyze_inl
   use energi
   use inform
   use iounit
   use mlpot
   implicit none
   character*12 form1
   character*12 form2
   character*240 fstr
   real(r_p) ect_,edsp_,er_
!
!
!     write out each energy component to the desired precision
!
   form1 = '5x,f18.4,i15'
   if (digits .ge. 6)  form1 = '3x,f18.6,i15'
   if (digits .ge. 8)  form1 = '1x,f20.8,i15'
   form2 = form1(1:3)//'d'//form1(5:12)
   fstr = '(/,'' Energy Component Breakdown :'',&
      &          11x,''Kcal/mole'',6x,''Interactions''/)'
   write (iout,fstr)
   if (neb.ne.0) then
      fstr = '('' Bond Stretching'',12x,'//form1//')'
      write (iout,fstr)  eb,neb
   end if
   if (nea.ne.0) then
      fstr = '('' Angle Bending'',14x,'//form1//')'
      write (iout,fstr)  ea,nea
   end if
   if (neba.ne.0) then
      fstr = '('' Stretch-Bend'',15x,'//form1//')'
      write (iout,fstr)  eba,neba
   end if
   if (neub.ne.0) then
      fstr = '('' Urey-Bradley'',15x,'//form1//')'
      write (iout,fstr)  eub,neub
   end if
   if (neaa.ne.0) then
      fstr = '('' Angle-Angle'',16x,'//form1//')'
      write (iout,fstr)  eaa,neaa
   end if
   if (neopb.ne.0) then
      fstr = '('' Out-of-Plane Bend'',10x,'//form1//')'
      write (iout,fstr)  eopb,neopb
   end if
   if (neopd.ne.0) then
      fstr = '('' Out-of-Plane Distance'',6x,'//form1//')'
      write (iout,fstr)  eopd,neopd
   end if
   if (neid.ne.0) then
      fstr = '('' Improper Dihedral'',10x,'//form1//')'
      write (iout,fstr)  eid,neid
   end if
   if (neit.ne.0) then
      fstr = '('' Improper Torsion'',11x,'//form1//')'
      write (iout,fstr)  eit,neit
   end if
   if (net.ne.0) then
      fstr = '('' Torsional Angle'',12x,'//form1//')'
      write (iout,fstr)  et,net
   end if
   if (nept.ne.0) then
      fstr = '('' Pi-Orbital Torsion'',9x,'//form1//')'
      write (iout,fstr)  ept,nept
   end if
   if (nebt.ne.0) then
      fstr = '('' Stretch-Torsion'',12x,'//form1//')'
      write (iout,fstr)  ebt,nebt
   end if
   if (neat.ne.0) then
      fstr = '('' Angle-Torsion'',14x,'//form1//')'
      write (iout,fstr)  eat,neat
   end if
   if (nett.ne.0) then
      fstr = '('' Torsion-Torsion'',12x,'//form1//')'
      write (iout,fstr)  ett,nett
   end if
   if (nev.ne.0) then
      if (abs(ev) .lt. 1.0d10) then
         fstr = '('' Van der Waals'',14x,'//form1//')'
      else
         fstr = '('' Van der Waals'',14x,'//form2//')'
      end if
      write (iout,fstr)  ev,nev
   end if
   if (ner.ne.0.or.er.ne.0) then
      er_ = enr2en(er)
      if (abs(er_) .lt. 1.0d10) then
         fstr = '('' Repulsion'',18x,'//form1//')'
      else
         fstr = '('' Repulsion'',18x,'//form2//')'
      end if
      write (iout,fstr)  er_,ner
   end if
   if (nedsp.ne.0.or.edsp.ne.0) then
      edsp_ = enr2en(edsp)
      fstr = '('' Dispersion'',17x,'//form1//')'
      write (iout,fstr)  edsp_,nedsp
   end if
   if (nec.ne.0) then
      if (abs(ec) .lt. 1.0d10) then
         fstr = '('' Charge-Charge'',14x,'//form1//')'
      else
         fstr = '('' Charge-Charge'',14x,'//form2//')'
      end if
      write (iout,fstr)  ec,nec
   end if
   if (nect.ne.0.or.ect.ne.0) then
      ect_ = enr2en(ect)
      if (abs(ect_) .lt. 1.0d10) then
         fstr = '('' Charge Transfer'',12x,'//form1//')'
      else
         fstr = '('' Charge Transfer'',12x,'//form2//')'
      end if
      write (iout,fstr)  ect_,nect
   end if
   if (nem.ne.0) then
      if (abs(em) .lt. 1.0d10) then
         fstr = '('' Atomic Multipoles'',10x,'//form1//')'
      else
         fstr = '('' Atomic Multipoles'',10x,'//form2//')'
      end if
      write (iout,fstr)  em,nem
   end if
   if (nep.ne.0) then
      if (abs(ep) .lt. 1.0d10) then
         fstr = '('' Polarization'',15x,'//form1//')'
      else
         fstr = '('' Polarization'',15x,'//form2//')'
      end if
      write (iout,fstr)  ep,nep
   end if
   if (neg.ne.0) then
      fstr = '('' Geometric Restraints'',7x,'//form1//')'
      write (iout,fstr)  eg,neg
   end if
   if (nex.ne.0) then
      fstr = '('' Extra Energy Terms'',9x,'//form1//')'
      write (iout,fstr)  ex,nex
   end if
   if (nemlpot.ne.0) then
      fstr = '('' ML Potential'',7x,'//form1//')'
      write (iout,fstr)  emlpot,nemlpot
   end if
   return
end
!
!     subroutine dipolemol: get the total dipole moment of the molecules system (permanent+induced)
!
subroutine dipolemol
   use atmlst
   use atoms
   use boxes
   use charge
   use domdec
   use iounit
   use molcul
   use mpole
   use polar
   use potent
   use units
   use mpi
   implicit none
   integer i
   integer iichg,iimol,kkpole,k
   real(t_p) q,xr,yr,zr
   real(r_p) dipx,dipy,dipz

   real(t_p) mux,muy,muz,mudx,mudy,mudz,mupx,mupy,mupz
1000 format(/'x dipolar moment of molecule',I5,' : ',F14.5)
1010 format(/'y dipolar moment of molecule',I5,' : ',F14.5)
1020 format(/'z dipolar moment of molecule',I5,' : ',F14.5)
1030 format(/'Norm of the dipolar moment of molecule',I5,' : ',F14.5)
!
!$acc wait
    if (nproc.gt.1) call molecule(.false.)

!$acc update host(molculeglob,rpole,uind,uinp)
   if (use_mpole) then
      do i = 1, nmoleloc
         iimol = molculeglob(i)
         dipx = 0
         dipy = 0
         dipz = 0
         do k = imol(1,iimol), imol(2,iimol)
            xr = x(k)
            yr = y(k)
            zr = z(k)
            kkpole = pollist(k)
            q = rpole(1,kkpole)
            mux = rpole(2,kkpole)
            muy = rpole(3,kkpole)
            muz = rpole(4,kkpole)
            mudx = uind(1,kkpole)
            mudy = uind(2,kkpole)
            mudz = uind(3,kkpole)
            mupx = uinp(1,kkpole)
            mupy = uinp(2,kkpole)
            mupz = uinp(3,kkpole)
            dipx = dipx + q*xr + mux + 0.5*(mudx+mupx)
            dipy = dipy + q*yr + muy + 0.5*(mudy+mupy)
            dipz = dipz + q*zr + muz + 0.5*(mudz+mupz)
         end do
         dipx = debye*dipx
         dipy = debye*dipy
         dipz = debye*dipz
         write(iout,1000) iimol,dipx
         write(iout,1010) iimol,dipy
         write(iout,1020) iimol,dipz
         write(iout,1030) iimol,sqrt(dipx**2+dipy**2+dipz**2)
      end do
   else if (use_charge) then
      do i = 1, nmoleloc
         dipx = 0
         dipy = 0
         dipz = 0
         iimol = molculeglob(i)
         do k = imol(1,iimol), imol(2,iimol)
            iichg = chglist(k)
            xr = x(k)
            yr = y(k)
            zr = z(k)
            q = pchg(iichg)
            dipx = dipx + q*xr
            dipy = dipy + q*yr
            dipz = dipz + q*zr
         end do
         dipx = debye*dipx
         dipy = debye*dipy
         dipz = debye*dipz
         write(iout,1000) iimol,dipx
         write(iout,1010) iimol,dipy
         write(iout,1020) iimol,dipz
         write(iout,1030) iimol,sqrt(dipx**2+dipy**2+dipz**2)
      end do
   end if
   return
end
