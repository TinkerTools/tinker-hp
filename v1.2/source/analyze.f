c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #############################################################
c     ##                                                         ##
c     ##  program analyze  --  energy partitioning and analysis  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "analyze" computes and displays the total potential energy;
c
c
      program analyze
      use mpi
      implicit none
      integer ierr
      call MPI_INIT(ierr)
      call analyze_bis
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_FINALIZE(ierr)
      end
c
      subroutine analyze_bis
      use dcdmod
      use domdec
      use files
      use inform
      use iounit
      use mpi
      use potent
      implicit none
      integer i,ixyz
      integer frame
      integer freeunit
      integer trimtext
      real*8 energy
      logical doenergy,dodipoltot,dodipolmol
      logical exist
      character*1 letter
      character*240 string
      character*240 xyzfile
      character*240 dcdfile
      type(dcdinfo_t) :: dcdinfo
c
c     initializing use_grad
c
      use_grad = .false.
c
c     set up the structure and mechanics calculation
c
      call initial
      call getxyz
      nproc = nproctot
      call initmpi
      call unitcell
      call cutoffs
      call lattice
c
c     get parameters
c
      call mechanic
c
c     setup for MPI
c
      call drivermpi
      call kewald_2
      call reinitnl(0)
      call mechanic_init_para
c
      call nblist(0)
c
c     get the desired types of analysis to be performed
c
      call nextarg (string,exist)
      if (.not. exist) then
         write (iout,10)
   10    format (/,' The TINKER Analysis Facility can Provide :',
     &           /,' Total Potential Energy and its Components [E]',
     &           /,' Total Dipolar Moment [D]',
     &           /,' Molecular Dipolar Moments [M]')
   20    continue
         write (iout,30)
   30    format (/,' Enter the Desired Analysis Types',
     &              ' [E,D,M] :  ',$)
         read (input,40,err=20)  string
   40    format (a240)
      end if
c
c     set option control flags based desired analysis types
c
      doenergy = .false.
      dodipoltot = .false.
      dodipolmol = .false.
      call upcase (string)
      do i = 1, trimtext(string)
         letter = string(i:i)
         if (letter .eq. 'E')  doenergy = .true.
         if (letter .eq. 'D')  dodipoltot = .true.
         if (letter .eq. 'M')  dodipolmol = .true.
      end do
c
c     reopen the coordinates file and read the first structure
c
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
c
c     perform analysis for each successive coordinate structure
c
      do while (.not. abort)
           frame = frame + 1
           if (frame .gt. 1) then
              if (rank.eq.0) write (iout,90)  frame
   90         format (/,' Analysis for Archive Structure :',8x,i8)
           end if
c
c        the box shape can change between frames
c
         call lattice
c
c       setup for MPI
c
         call ddpme3d
         call reassignpme(.true.)
         call reinitnl(0)
         call mechanic_up_para(0)
         call nblist(0)
c
c     make the call to compute the potential energy
c
         if (doenergy)  call enrgyze
c
c     make the call to compute the induced dipoles and the total dipolar moments
c
         if (dodipoltot)  then
           call analysis(energy)
           call totaldipole
         end if
         if (dodipolmol)  then
           call analysis(energy)
           call dipolemol
         end if
c
c     energy partitioning by potential energy components
c
         if (doenergy) then
            if (rank.eq.0) call partyze
         end if
c
c     attempt to read next structure from the coordinate file
c
         if (.not.dcdio) then
           call readxyz (ixyz)
         else 
           call dcdfile_read_next(dcdinfo)
           call dcdfile_skip_next(dcdinfo,0)
         end if
      end do
c
c     perform any final tasks before program exit
c
      if (dcdio) then
        call dcdfile_close(dcdinfo)
      else
        close (unit=ixyz)
      end if
      call final
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine enrgyze  --  compute & report energy analysis  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "enrgyze" is an auxiliary routine for the analyze program
c     that performs the energy analysis and prints the total and
c     intermolecular energies
c
c
      subroutine enrgyze
      use atoms
      use cutoff
      use domdec
      use inform
      use inter
      use iounit
      use molcul
      implicit none
      real*8 energy
      character*240 fstr
c
c
c     perform the energy analysis by atom and component
c
      call analysis (energy)
      if (rank.eq.0) then
c
c       print out the total potential energy of the system
c
        fstr='(/,'' Total Potential Energy :'',8x,f16.4,'' Kcal/mole'')'
        if (digits .ge. 6)  fstr(32:39) = '6x,f18.6'
        if (digits .ge. 8)  fstr(32:39) = '4x,f20.8'
        if (abs(energy) .ge. 1.0d10)  fstr(35:35) = 'd'
        write (iout,fstr)  energy
c
c       intermolecular energy for systems with multiple molecules
c
        fstr='(/,'' Intermolecular Energy :'',9x,f16.4,'' Kcal/mole'')'
        if (digits .ge. 6)  fstr(31:38) = '7x,f18.6'
        if (digits .ge. 8)  fstr(31:38) = '5x,f20.8'
        if (abs(einter) .ge. 1.0d10)  fstr(34:34) = 'd'
        if (nmol.gt.1 .and. nmol.lt.n .and. .not.use_ewald)
     &     write (iout,fstr)  einter
      end if
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine partyze  --  energy component decomposition  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "partyze" prints the energy component and number of
c     interactions for each of the potential energy terms
c
c
      subroutine partyze
      use action
      use energi
      use inform
      use iounit
      implicit none
      character*12 form1
      character*12 form2
      character*240 fstr
c
c
c     write out each energy component to the desired precision
c
      form1 = '5x,f16.4,i15'
      if (digits .ge. 6)  form1 = '3x,f18.6,i15'
      if (digits .ge. 8)  form1 = '1x,f20.8,i15'
      form2 = form1(1:3)//'d'//form1(5:12)
      fstr = '(/,'' Energy Component Breakdown :'',
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
      if (ner.ne.0.or.er.ne.0.0d0) then
         if (abs(er) .lt. 1.0d10) then
            fstr = '('' Repulsion'',18x,'//form1//')'
         else
            fstr = '('' Repulsion'',18x,'//form2//')'
         end if
         write (iout,fstr)  er,ner
      end if
      if (nedsp.ne.0.or.edsp.ne.0.0d0) then
         fstr = '('' Dispersion'',17x,'//form1//')'
         write (iout,fstr)  edsp,nedsp
      end if
      if (nec.ne.0) then
         if (abs(ec) .lt. 1.0d10) then
            fstr = '('' Charge-Charge'',14x,'//form1//')'
         else
            fstr = '('' Charge-Charge'',14x,'//form2//')'
         end if
         write (iout,fstr)  ec,nec
      end if
      if (nect.ne.0.or.ect.ne.0.0d0) then
         if (abs(ect) .lt. 1.0d10) then
            fstr = '('' Charge Transfer'',12x,'//form1//')'
         else
            fstr = '('' Charge Transfer'',12x,'//form2//')'
         end if
         write (iout,fstr)  ect,nect
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
      return
      end
c
c     subroutine dipolemol: get the total dipole moment of the molecules system (permanent+induced)
c
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
      real*8 q,xr,yr,zr
      real*8 dipx,dipy,dipz

      real*8 mux,muy,muz,mudx,mudy,mudz,mupx,mupy,mupz
 1000 format(/'x dipolar moment of molecule',I15,' : ',F14.5)
 1010 format(/'y dipolar moment of molecule',I15,' : ',F14.5)
 1020 format(/'z dipolar moment of molecule',I15,' : ',F14.5)
 1030 format(/'Norm of the dipolar moment of molecule',I15,' : ',F14.5)
c
      call molecule(.false.)
      if (use_mpole) then
        do i = 1, nmoleloc
          iimol = molculeglob(i)
          dipx = 0d0
          dipy = 0d0
          dipz = 0d0
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
          dipx = 0d0
          dipy = 0d0
          dipz = 0d0
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
