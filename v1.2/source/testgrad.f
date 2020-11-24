c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  program testgrad  --  derivative test; Cartesian version  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "testgrad" computes and compares the analytical and numerical
c     gradient vectors of the potential energy function with respect
c     to Cartesian coordinates
c
c
      program testgrad
      use atoms
      use deriv
      use domdec
      use energi
      use sizes
      use inform
      use iounit
      use mpole
      use usage
      use mpi
      implicit none
      integer i,j,next,iglob,ierr,nthreadsupport
      integer iloc
      real*8 etot,f,f0,eps,eps0,old,energy
      real*8 eb0,ea0,eba0,eub0,eaa0,eopb0
      real*8 eopd0,eid0,eit0,et0,ept0,ebt0
      real*8 ett0,ev0,ec0,em0,ep0
      real*8 eg0,ex0
      real*8 totnorm,ntotnorm,rms,nrms
      real*8, allocatable :: denorm(:)
      real*8, allocatable :: ndenorm(:)
      real*8, allocatable :: detot(:,:)
      real*8, allocatable :: ndetot(:,:)
      real*8, allocatable :: ndeb(:,:)
      real*8, allocatable :: ndea(:,:)
      real*8, allocatable :: ndeba(:,:)
      real*8, allocatable :: ndeub(:,:)
      real*8, allocatable :: ndeaa(:,:)
      real*8, allocatable :: ndeopb(:,:)
      real*8, allocatable :: ndeopd(:,:)
      real*8, allocatable :: ndeid(:,:)
      real*8, allocatable :: ndeit(:,:)
      real*8, allocatable :: ndet(:,:)
      real*8, allocatable :: ndept(:,:)
      real*8, allocatable :: ndebt(:,:)
      real*8, allocatable :: ndett(:,:)
      real*8, allocatable :: ndev(:,:)
      real*8, allocatable :: ndem(:,:)
      real*8, allocatable :: ndec(:,:)
      real*8, allocatable :: ndep(:,:)
      real*8, allocatable :: ndeg(:,:)
      real*8, allocatable :: ndex(:,:)
      logical exist,query
      logical doanalyt,donumer,dofull
      character*1 answer
      character*1 axis(3)
      character*120 record
      character*120 string
      data axis  / 'X','Y','Z' /
c
      call MPI_INIT_THREAD(MPI_THREAD_MULTIPLE,nthreadsupport,ierr)
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
      call drivermpi
      call reinitnl(0)
c
      dotstgrad = .true.
c
      call mechanic
      call nblist(0)
      call allocstep
c
c     decide whether to do an analytical gradient calculation
c
      doanalyt = .true.
      call nextarg (answer,exist)
      if (.not. exist) then
         write (iout,10)
   10    format (/,' Compute the Analytical Gradient Vector [Y] :  ',$)
         read (input,20)  record
   20    format (a120)
         next = 1
         call gettext (record,answer,next)
      end if
      call upcase (answer)
      if (answer .eq. 'N')  doanalyt = .false.
c
c     decide whether to do a numerical gradient calculation
c
      donumer = .true.
      call nextarg (answer,exist)
      if (.not. exist) then
         write (iout,30)
   30    format (/,' Compute the Numerical Gradient Vector [Y] :   ',$)
         read (input,40)  record
   40    format (a120)
         next = 1
         call gettext (record,answer,next)
      end if
      call upcase (answer)
      if (answer .eq. 'N')  donumer = .false.
c
c     get the stepsize for numerical gradient calculation
c
      if (donumer) then
         eps = -1.0d0
         eps0 = 1.0d-5
         query = .true.
         call nextarg (string,exist)
         if (exist) then
            read (string,*,err=50,end=50)  eps
            query = .false.
         end if
   50    continue
         if (query) then
            write (iout,60)  eps0
   60       format (/,' Enter a Numerical Stepsize [',d8.1,
     &                 ' Ang] :  ',$)
            read (input,70,err=50)  eps
   70       format (f20.0)
         end if
         if (eps .le. 0.0d0)  eps = eps0
      end if
c
c     decide whether to output results by gradient component
c
      dofull = .true.
      if (n .gt. 100) then
         dofull = .false.
         call nextarg (answer,exist)
         if (.not. exist) then
            write (iout,80)
   80       format (/,' Output Breakdown by Gradient Component',
     &                 ' [N] :  ',$)
            read (input,90)  record
   90       format (a120)
            next = 1
            call gettext (record,answer,next)
         end if
         call upcase (answer)
         if (answer .eq. 'Y')  dofull = .true.
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (detot(3,nbloc))
      detot = 0d0
      call ddpme3d
      call reassignpme(.false.)
      call reinitnl(0)
      call mechanicstep(0)
      call nblist(0)
c
c     compute the analytical gradient components
c
      if (doanalyt) then
         call gradient (etot,detot)
c
c    MPI : communicate analytical forces and torques
c
         call commforceststgrad(detot)
         call commforcesrec(detot)
c
c    MPI : get total energies from gradient routine
c
         call reduceen(etot)
      end if
c
c     print the total potential energy of the system
c
      if (doanalyt.and.rank.eq.0) then
         if (digits .ge. 8) then
            write (iout,100)  etot
  100       format (/,' Total Potential Energy :',8x,f20.8,' Kcal/mole')
         else if (digits .ge. 6) then
            write (iout,110)  etot
  110       format (/,' Total Potential Energy :',8x,f18.6,' Kcal/mole')
         else
            write (iout,120)  etot
  120       format (/,' Total Potential Energy :',8x,f16.4,' Kcal/mole')
         end if
c
c     print the energy breakdown over individual components
c
         write (iout,130)
  130    format (/,' Potential Energy Breakdown by Individual',
     &              ' Components :')
c         if (digits .ge. 8) then
c            write (iout,140)
c  140       format (/,'  Energy',7x,'EB',14x,'EA',14x,'EBA',13x,'EUB',
c     &              /,'  Terms',8x,'EAA',13x,'EOPB',12x,'EOPD',
c     &                 12x,'EID',
c     &              /,15x,'EIT',13x,'ET',14x,'EPT',13x,'EBT',
c     &              /,15x,'ETT',13x,'EV',14x,'EM',14x,'EP')
c            write (iout,150)  eb,ea,eba,eub,eaa,eopb,eopd,eid,
c     &                        eit,et,ept,ebt,ett,ev,em,ep
c  150       format (/,6x,4f16.8,/,6x,4f16.8,/,6x,4f16.8,/,6x,4f16.8,
c     &                 /,6x,4f16.8,/,6x,4f16.8)
c         else if (digits .ge. 6) then
c            write (iout,160)
c  160       format (/,'  Energy',7x,'EB',14x,'EA',14x,'EBA',13x,'EUB',
c     &              /,'  Terms',8x,'EAA',13x,'EOPB',12x,'EOPD',
c     &                 12x,'EID',
c     &              /,15x,'EIT',13x,'ET',14x,'EPT',13x,'EBT',
c     &              /,15x,'ETT',13x,'EV',14x,'EM',14x,'EP')
c            write (iout,170)  eb,ea,eba,eub,eaa,eopb,eopd,eid,
c     &                        eit,et,ept,ebt,ett,ev,em,ep
c  170       format (/,6x,5f14.6,/,6x,5f14.6,/,6x,5f14.6,/,6x,5f14.6,
c     &                 /,6x,4f14.6)
c         else
            write (iout,180)
  180       format (/,'  Energy',7x,'EB',14x,'EA',14x,'EBA',13x,'EUB',
     &              11x,'EX',
     &              /,'  Terms',8x,'EAA',13x,'EOPB',12x,'EOPD',
     &                 12x,'EID',11x,'EC',
     &              /,15x,'EIT',13x,'ET',14x,'EPT',13x,'EBT',
     &              11x,'EREP',/,15x,'ETT',13x,'EV',14x,'EM',14x,'EP',
     &               12x,'EXDISP')
            write (iout,190)  eb,ea,eba,eub,ex,eaa,eopb,eopd,eid,ec,
     &                        eit,et,ept,ebt,0d0,ett,ev,em,ep,0d0
  190       format (/,6x,5f15.4,/,6x,5f15.4,/,6x,5f15.4,/,6x,5f15.4)
c         end if
      end if
c
c     print a header for the gradients of individual potentials
c
      if (dofull.and.rank.eq.0) then
         write (iout,200)
  200    format (/,' Cartesian Gradient Breakdown by Individual',
     &                 ' Components :')
         if (digits .ge. 8) then
            write (iout,210)
  210       format (/,2x,'Atom',9x,'d EB',11x,'d EA',11x,'d EBA',
     &                 10x,'d EUB',8x,'d EX',
     &              /,2x,'Axis',9x,'d EAA',9x,'d EOPB',10x,'d EOPD',
     &                 9x,'d EID',8x,'d EC',
     &              /,2x,'Type',9x,'d EIT',10x,'d ET',11x,'d EPT',
     &                 10x,'d EBT',8x,'d EREP',
     &              /,15x,'d ETT',10x,'d EV',11x,'d EM',11x,'d EP',9x,
     &               'd EXDISP')
         else if (digits .ge. 6) then
            write (iout,220)
  220       format (/,2x,'Atom',9x,'d EB',11x,'d EA',11x,'d EBA',
     &                 10x,'d EUB',8x,' d EX',
     &              /,2x,'Axis',9x,'d EAA',9x,'d EOPB',10x,'d EOPD',
     &                 9x,'d EID',8x,'d EC',
     &              /,2x,'Type',9x,'d EIT',10x,'d ET',11x,'d EPT',
     &                 10x,'d EBT',8x,'d EREP'
     &              /,15x,'d ETT',10x,'d EV',11x,'d EM',11x,'d EP',9x,
     &                'd EXDISP')
         else
            write (iout,230)
  230       format (/,2x,'Atom',9x,'d EB',11x,'d EA',11x,'d EBA',
     &                 10x,'d EUB',8x,'d EX',
     &              /,2x,'Axis',9x,'d EAA',10x,'d EOPB',9x,'d EOPD',
     &                 9x,'d EID',8x,'d EC',
     &              /,2x,'Type',9x,'d EIT',10x,'d ET',11x,'d EPT',
     &                 10x,'d EBT',8x,'d EREP',
     &              /,15x,'d ETT',10x,'d EV',11x,'d EM',11x,'d EP',9x,
     &               'd EXDISP')
         end if
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (ndetot(3,nloc))
      allocate (ndeb(3,nloc))
      allocate (ndea(3,nloc))
      allocate (ndeba(3,nloc))
      allocate (ndeub(3,nloc))
      allocate (ndeaa(3,nloc))
      allocate (ndeopb(3,nloc))
      allocate (ndeopd(3,nloc))
      allocate (ndeid(3,nloc))
      allocate (ndeit(3,nloc))
      allocate (ndet(3,nloc))
      allocate (ndept(3,nloc))
      allocate (ndebt(3,nloc))
      allocate (ndett(3,nloc))
      allocate (ndev(3,nloc))
      allocate (ndem(3,nloc))
      allocate (ndec(3,nloc))
      allocate (ndep(3,nloc))
      allocate (ndeg(3,nloc))
      allocate (ndex(3,nloc))
      ndetot = 0d0
      ndeb = 0d0
      ndea = 0d0
      ndeba = 0d0
      ndeub = 0d0
      ndeaa = 0d0
      ndeopb = 0d0
      ndeopd = 0d0
      ndeid = 0d0
      ndet = 0d0
      ndept = 0d0
      ndebt = 0d0
      ndett = 0d0
      ndev = 0d0
      ndem = 0d0
      ndec = 0d0
      ndep = 0d0
      ndeg = 0d0
      ndex = 0d0
c
c     get the Cartesian component two-sided numerical gradients
c
      do i = 1, n
         if (repart(i).eq.rank) iloc = loc(i)
         call MPI_BARRIER(COMM_TINKER,ierr)
         if (donumer .and. use(i)) then
            do j = 1, 3
               if (j .eq. 1) then
                  old = x(i)
                  x(i) = x(i) - 0.5d0*eps
               else if (j .eq. 2) then
                  old = y(i)
                  y(i) = y(i) - 0.5d0*eps
               else if (j .eq. 3) then
                  old = z(i)
                  z(i) = z(i) - 0.5d0*eps
               end if
               f0 = energy ()
               call allreduceen(f0)
               eb0 = eb
               ea0 = ea
               eba0 = eba
               eub0 = eub
               eaa0 = eaa
               eopb0 = eopb
               eopd0 = eopd
               eid0 = eid
               eit0 = eit
               et0 = et
               ept0 = ept
               ebt0 = ebt
               ett0 = ett
               ev0 = ev
               ec0 = ec
               em0 = em
               ep0 = ep
               eg0 = eg
               ex0 = ex
               eg0 = eg
               if (j .eq. 1) then
                  x(i) = x(i) + eps
               else if (j .eq. 2) then
                  y(i) = y(i) + eps
               else if (j .eq. 3) then
                  z(i) = z(i) + eps
               end if
               f = energy ()
               call allreduceen(f)
               if (j .eq. 1) then
                  x(i) = old
               else if (j .eq. 2) then
                  y(i) = old
               else if (j .eq. 3) then
                  z(i) = old
               end if
c
c   si i dans domenloc alors nde et affichage
c
               if (repart(i).eq.rank) then
                 ndetot(j,iloc) = (f - f0) / eps
                 ndeb(j,iloc) = (eb - eb0) / eps
                 ndea(j,iloc) = (ea - ea0) / eps
                 ndeba(j,iloc) = (eba - eba0) / eps
                 ndeub(j,iloc) = (eub - eub0) / eps
                 ndeaa(j,iloc) = (eaa - eaa0) / eps
                 ndeopb(j,iloc) = (eopb - eopb0) / eps
                 ndeopd(j,iloc) = (eopd - eopd0) / eps
                 ndeid(j,iloc) = (eid - eid0) / eps
                 ndeit(j,iloc) = (eit - eit0) / eps
                 ndet(j,iloc) = (et - et0) / eps
                 ndept(j,iloc) = (ept - ept0) / eps
                 ndebt(j,iloc) = (ebt - ebt0) / eps
                 ndett(j,iloc) = (ett - ett0) / eps
                 ndev(j,iloc) = (ev - ev0) / eps
                 ndem(j,iloc) = (em - em0) / eps
                 ndec(j,iloc) = (ec - ec0) / eps
                 ndep(j,iloc) = (ep - ep0) / eps
                 ndeg(j,iloc) = (eg - eg0) / eps
                 ndex(j,iloc) = (ex - ex0) / eps
               end if
            end do
         end if
c
c     print analytical gradients of each energy term for each atom
c
         if (dofull .and. use(i)) then
            do j = 1, 3
               if (doanalyt.and.repart(i).eq.rank) then
c                  if (digits .ge. 8) then
c                     write (iout,240) iglob,deb(j,i),dea(j,i),deba(j,i),
c     &                                deub(j,i),axis(j),deaa(j,i),
c     &                                deopb(j,i),deopd(j,i),deid(j,i),
c     &                                deit(j,i),det(j,i),dept(j,i),
c     &                                debt(j,i),dett(j,i),dev(j,i),
c     &                                dem(j,i),dep(j,i)
c  240                format (/,i6,4f16.8,/,5x,a1,4f16.8,
c     &                          /,' Anlyt',4f16.8,/,6x,4f16.8,
c     &                          /,6x,4f16.8,/,6x,4f16.8)
c                  else if (digits .ge. 6) then
c                     write (iout,250) iglob,deb(j,i),dea(j,i),deba(j,i),
c     &                                deub(j,i),deaa(j,i),axis(j),
c     &                                deopb(j,i),deopd(j,i),deid(j,i),
c     &                                deit(j,i),det(j,i),dept(j,i),
c     &                                debt(j,i),dett(j,i),dev(j,i),
c     &                                dem(j,i),dep(j,i)
c  250                format (/,i6,5f14.6,/,5x,a1,5f14.6,/,' Anlyt',
c     &                          5f14.6,/,6x,5f14.6,/,6x,4f14.6)
c                  else
                     write (iout,260) i,deb(j,iloc),dea(j,iloc),
     &                        deba(j,iloc),deub(j,iloc),dex(j,iloc),
     &                        axis(j),
     &                        deaa(j,iloc),deopb(j,iloc),deopd(j,iloc),
     &                        deid(j,iloc),dec(j,iloc),deit(j,iloc),
     &                        det(j,iloc),
     &                        dept(j,iloc),debt(j,iloc),0d0,
     &                        dett(j,iloc),dev(j,iloc),dem(j,iloc),
     &                        dep(j,iloc),0d0
  260                format (/,i6,5f15.4,/,5x,a1,5f15.4,/,' Anlyt',
     &                          5f15.4,/,6x,5f15.4)
c                  end if
               end if
c
c     print numerical gradients of each energy term for each atom
c
               if (donumer.and.repart(i).eq.rank) then
c                  if (digits .ge. 8) then
c                     write (iout,270) iglob,ndeb(j,i),ndea(j,i),
c     &                                ndeba(j,i),ndeub(j,i),ndeaa(j,i),
c     &                                ndeopb(j,i),axis(j),
c     &                                ndeaa(j,i),ndeopb(j,i),
c     &                                ndeopd(j,i),ndeid(j,i),ndeit(j,i), c     &                                ndet(j,i),ndept(j,i),ndebt(j,i),
c     &                                ndett(j,i),ndev(j,i),
c     &                                ndem(j,i),ndep(j,i)
c  270                format (/,i6,4f16.8,/,5x,a1,4f16.8,
c     &                          /,' Numer',4f16.8,/,6x,4f16.8,
c     &                          /,6x,4f16.8,/,6x,4f16.8)
c                  else if (digits .ge. 6) then
c                     write (iout,280) iglob,ndeb(j,i),ndea(j,i),
c     &                                ndeba(j,i),ndeub(j,i),ndeaa(j,i),
c     &                                ndeopb(j,i),axis(j),
c     &                                ndeaa(j,i),ndeopb(j,i),
c     &                                ndeopd(j,i),ndeid(j,i),ndeit(j,i),
c     &                                ndet(j,i),ndept(j,i),ndebt(j,i),
c     &                                ndett(j,i),ndev(j,i),
c     &                                ndem(j,i),ndep(j,i)
c  280                format (/,i6,5f14.6,/,5x,a1,5f14.6,/,' Numer',
c     &                          5f14.6,/,6x,5f14.6,/,6x,4f14.6)
c                  else
                     write (iout,290) i,ndeb(j,iloc),ndea(j,iloc),
     &                           ndeba(j,iloc),ndeub(j,iloc),
     &                           ndex(j,iloc),axis(j),
     &                           ndeaa(j,iloc),ndeopb(j,iloc),
     &                           ndeopd(j,iloc),ndeid(j,iloc),
     &                           ndec(j,iloc),ndeit(j,iloc),
     &                           ndet(j,iloc),
     &                           ndept(j,iloc),ndebt(j,iloc),
     &                           0d0,ndett(j,iloc),
     &                           ndev(j,iloc),ndem(j,iloc),ndep(j,iloc),
     &                           0d0

  290                format (/,i6,5f15.4,/,5x,a1,5f15.4,/,' Numer',
     &                          5f15.4,/,6x,5f15.4)
c                  end if
               end if
            end do
         end if
      end do
      call MPI_BARRIER(COMM_TINKER,ierr)
c
c     perform deallocation of some local arrays
c
      deallocate (ndeb)
      deallocate (ndea)
      deallocate (ndeba)
      deallocate (ndeub)
      deallocate (ndeaa)
      deallocate (ndeopb)
      deallocate (ndeopd)
      deallocate (ndeid)
      deallocate (ndeit)
      deallocate (ndet)
      deallocate (ndept)
      deallocate (ndebt)
      deallocate (ndett)
      deallocate (ndev)
      deallocate (ndem)
      deallocate (ndec)
      deallocate (ndep)
      deallocate (ndeg)
      deallocate (ndex)
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (denorm(nloc))
      allocate (ndenorm(nloc))
      denorm = 0d0
      ndenorm = 0d0
c
      call MPI_BARRIER(COMM_TINKER,ierr)
c
c     print the total gradient components for each atom
c
      if (rank.eq.0) then
        write (iout,300)
  300 format (/,' Cartesian Gradient Breakdown over Individual Atoms :')
        if (digits .ge. 8) then
           write (iout,310)
  310      format (/,2x,'Type',4x,'Atom',10x,'dE/dX',11x,'dE/dY',
     &                11x,'dE/dZ',11x,'Norm',/)
        else if (digits .ge. 6) then
           write (iout,320)
  320      format (/,2x,'Type',6x,'Atom',11x,'dE/dX',9x,'dE/dY',
     &                9x,'dE/dZ',11x,'Norm',/)
        else
           write (iout,330)
  330      format (/,2x,'Type',6x,'Atom',14x,'dE/dX',7x,'dE/dY',
     &                7x,'dE/dZ',10x,'Norm',/)
        end if
      end if
      totnorm = 0.0d0
      ntotnorm = 0.0d0
      do i = 1, nloc
         iglob = glob(i)
         if (doanalyt .and. use(iglob)) then
            denorm(i) = detot(1,i)**2 + detot(2,i)**2 +
     &                       detot(3,i)**2
            totnorm = totnorm + denorm(i)
            denorm(i) = sqrt(denorm(i))
c            if (digits .ge. 8) then
c               write (iout,340)  iglob,(detot(j,i),j=1,3),denorm(i)
c  340          format (' Anlyt',i8,1x,3f16.8,f16.8)
c            else if (digits .ge. 6) then
c               write (iout,350)  iglob,(detot(j,i),j=1,3),denorm(i)
c  350          format (' Anlyt',2x,i8,3x,3f14.6,2x,f14.6)
c            else
               write (iout,360)  iglob,(detot(j,i),j=1,3),denorm(i)
  360          format (' Anlyt',2x,i8,7x,3f12.4,2x,f12.4)
c            end if
         end if
         if (donumer .and. use(iglob)) then
            ndenorm(i) = ndetot(1,i)**2 + ndetot(2,i)**2 +
     &                        ndetot(3,i)**2
            ntotnorm = ntotnorm + ndenorm(i)
            ndenorm(i) = sqrt(ndenorm(i))
c            if (digits .ge. 8) then
c               write (iout,370)  iglob,(ndetot(j,i),j=1,3),ndenorm(i)
c  370          format (' Numer',i8,1x,3f16.8,f16.8)
c            else if (digits .ge. 6) then
c               write (iout,380)  iglob,(ndetot(j,i),j=1,3),ndenorm(i)
c  380          format (' Numer',2x,i8,3x,3f14.6,2x,f14.6)
c            else
               write (iout,390)  iglob,(ndetot(j,i),j=1,3),ndenorm(i)
  390          format (' Numer',2x,i8,7x,3f12.4,2x,f12.4)
c            end if
         end if
      end do
      call MPI_BARRIER(COMM_TINKER,ierr)
      if (rank.eq.0) then
       call MPI_REDUCE(MPI_IN_PLACE,totnorm,1,MPI_REAL8,MPI_SUM,0,
     $    COMM_TINKER,ierr)
      else
       call MPI_REDUCE(totnorm,totnorm,1,MPI_REAL8,MPI_SUM,0,
     $    COMM_TINKER,ierr)
      end if
      if (rank.eq.0) then
       call MPI_REDUCE(MPI_IN_PLACE,ntotnorm,1,MPI_REAL8,MPI_SUM,
     $    0,COMM_TINKER,ierr)
      else
       call MPI_REDUCE(ntotnorm,ntotnorm,1,MPI_REAL8,MPI_SUM,0,
     $    COMM_TINKER,ierr)
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (detot)
      deallocate (ndetot)
      deallocate (denorm)
      deallocate (ndenorm)
c
c
c     print the total norm for the analytical gradient
c
      if (rank.eq.0) then
        write (iout,400)
  400   format (/,' Total Gradient Norm and RMS Gradient per Atom :')
        if (doanalyt) then
           totnorm = sqrt(totnorm)
c           if (digits .ge. 8) then
c              write (iout,410)  totnorm
c  410         format (/,' Anlyt',6x,'Total Gradient Norm Value',6x,f20.8)
c           else if (digits .ge. 6) then
c              write (iout,420)  totnorm
c  420         format (/,' Anlyt',6x,'Total Gradient Norm Value',6x,f18.6)
c           else
             write (iout,430)  totnorm
  430        format (/,' Anlyt',6x,'Total Gradient Norm Value',6x,f16.4)
c           end if
        end if
      end if
c
c     print the total norm for the numerical gradient
c
      if (donumer) then
         ntotnorm = sqrt(ntotnorm)
c         if (digits .ge. 8) then
c            write (iout,440)  ntotnorm
c  440       format (' Numer',6x,'Total Gradient Norm Value',6x,f20.8)
c         else if (digits .ge. 6) then
c            write (iout,450)  ntotnorm
c  450       format (' Numer',6x,'Total Gradient Norm Value',6x,f18.6)
c         else
            if (rank.eq.0) write (iout,460)  ntotnorm
  460       format (' Numer',6x,'Total Gradient Norm Value',6x,f16.4)
c         end if
      end if
c
c     print the rms per atom norm for the analytical gradient
c
      if (doanalyt) then
         rms = totnorm / sqrt(dble(nuse))
c         if (digits .ge. 8) then
c            write (iout,470)  rms
c  470       format (/,' Anlyt',6x,'RMS Gradient over All Atoms',
c     &                 4x,f20.8)
c         else if (digits .ge. 6) then
c            write (iout,480)  rms
c  480       format (/,' Anlyt',6x,'RMS Gradient over All Atoms',
c     &                 4x,f18.6)
c         else
            if (rank.eq.0) write (iout,490)  rms
  490       format (/,' Anlyt',6x,'RMS Gradient over All Atoms',
     &                 4x,f16.4)
c         end if
      end if
c
c     print the rms per atom norm for the numerical gradient
c
      if (donumer) then
         nrms = ntotnorm / sqrt(dble(nuse))
c         if (digits .ge. 8) then
c            write (iout,500)  nrms
c  500       format (' Numer',6x,'RMS Gradient over All Atoms',4x,f20.8)
c         else if (digits .ge. 6) then
c            write (iout,510)  nrms
c  510       format (' Numer',6x,'RMS Gradient over All Atoms',4x,f18.6)
c         else
            if (rank.eq.0) write (iout,520)  nrms
  520       format (' Numer',6x,'RMS Gradient over All Atoms',4x,f16.4)
c         end if
      end if
c
c     perform any final tasks before program exit
c
      call final
      call MPI_FINALIZE(ierr)
      end
