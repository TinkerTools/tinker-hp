c
c          ###    ##    ##    ###
c         #       # #  # #    #  ##
c        #        #  ##  #    #    #
c         ###     #      #    #    #
c           #     #      #    #    #
c          #      #      #    #  ##
c       ###       #      #    ###
c
c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine ksmd  --  initialize smd procedure            ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "ksmd" checks the smd activation and assign the necessary
c     parameters enters in the key file
c
c
#include "tinker_precision.h"
      subroutine ksmd(init)
      use angle
      use atomsMirror
      use atmtyp
      use bitor
      use bond
      use deriv
      use domdec
      use energi
      use keys
      use iounit
      use improp
      use inform ,only:deb_path
      use mdstuf
      use msmd
      use mpi
      use potent
      use sizes
      use tinheader
      use tors
      implicit none


      character(len=100)::file_name1
      character(len=1)::num1
      character(len=2)::num2
      character(len=3)::num3
      character(len=4)::num4
      character(len=5)::num5
      character(len=6)::num6
      character(len=7)::num7
      integer i,j,k,cpt1
      integer ii, jj
      integer next
      integer nsmdcap
      integer ismd !temporary SMD atom number during the processors distrubution
      integer, allocatable :: reqrec(:)
      logical, allocatable :: smdproc(:)
      integer :: reqsend,ierr,status(MPI_STATUS_SIZE),tagmpi
      real r
      logical init, quit
      logical calccom, manualcom
      logical warning

      real(t_p) eps,dalt,dalt2
      integer nalt, nalt2
      character*20 keyword
      character*240 record
      character*240 string      
c
c     Check if initialization or not (LOOP1)
c
      if (init) then
c
c     Initialisation parameters
c
      use_smd_velconst = .false.
      use_smd_forconst = .false.
      use_smdk2 = .false.
      use_atfol = .false.
      calccom   = .true.
      manualcom = .false.
      warning   = .false.
      smdprocprint = -1
      quit      = .true.
      SMDk      = 0.0_ti_p
      SMDk2     = 0.0_ti_p
      SMDVel    = 0.0_ti_p
      SMDFor    = 0.0_ti_p
      SMDoutputFreq = 1
      stock_dedx    = 0.0_ti_p
      stock_dedy    = 0.0_ti_p
      stock_dedz    = 0.0_ti_p
      stock_ensmd   = 0.0_ti_p
      nalt    = 1
      nalt2   = 1
      atfol   = 0
         xcom = 0.0_ti_p
         ycom = 0.0_ti_p
         zcom = 0.0_ti_p
         xdir = 0.0_ti_p
         ydir = 0.0_ti_p
         zdir = 0.0_ti_p
      mtotcom = 0.0_ti_p
      ensmd   = 0.0_ti_p
      tsmd    = 0.0_ti_p
      SMDdt   = 0.001
      tpass   = 0

      if (deb_path) print*,'ksmd'
c
c######################################################################
c
c     Check the smd key activation for constant velocity SMD
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:240)
         if (keyword(1:6) .eq. 'CVSMD ') then
            use_smd_velconst = .true.
            ismdout = 3
            if (rank.eq.0) then
                open (ismdout,file='SMD_output.dat',action="write")
                call promosmd(1)
            end if
         end if
      end do
c
c     Check the smd key activation for constant force SMD
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:240)
         if (keyword(1:6) .eq. 'CFSMD ') then
            use_smd_forconst = .true.
            ismdout = 3
            if (rank.eq.0) then
                open (ismdout,file='SMD_output.dat',action="write")
                call promosmd(1)
            end if
         end if
      end do
c
c######################################################################
c
c     Assign the desired parameters for the SMD (LOOP2)
c
      if (use_smd_velconst .or. use_smd_forconst) then
c
c     Elastic constant asignment SMDk (own to the both SMD)
c
         if (use_smd_velconst) then
         do i = 1, nkey
             next = 1
             record = keyline(i)
             call gettext (record,keyword,next)
             call upcase (keyword)
             string = record(next:240)
             if (keyword(1:5) .eq. 'SMDK ') then
                read (string,*,err=10) SMDk
             endif
         end do
         end if
c
c    Transverse elastic constant SMDk2 (own to the both SMD)
c
         if (use_smd_velconst) then
         do i = 1, nkey
             next = 1
             record = keyline(i)
             call gettext (record,keyword,next)
             call upcase (keyword)
             string = record(next:240)
             if (keyword(1:6) .eq. 'SMDK2 ') then
                use_smdk2 = .true.
                read (string,*,err=10) SMDk2
             endif
         end do
         end if
c
c    Constant velocity (only own to the constant velocity SMD)
c
         if (use_smd_velconst) then
         do i = 1, nkey
             next = 1
             record = keyline(i)
             call gettext (record,keyword,next)
             call upcase (keyword)
             string = record(next:240)
             if (keyword(1:7) .eq. 'SMDVEL ') then
                read (string,*,err=10) SMDVel
             endif
         end do
         end if
c
c    Constant force (only own to the constant force SMD)
c
         if (use_smd_forconst) then
         do i = 1, nkey
             next = 1
             record = keyline(i)
             call gettext (record,keyword,next)
             call upcase (keyword)
             string = record(next:240)
             if (keyword(1:7) .eq. 'SMDFOR ') then
                read (string,*,err=10) SMDFor
             endif
         end do
         end if
c
c    Time step of the SMD force
c
         do i = 1, nkey
             next = 1
             record = keyline(i)
             call gettext (record,keyword,next)
             call upcase (keyword)
             string = record(next:240)
             if (keyword(1:6) .eq. 'SMDDT ') then
                read (string,*,err=10) SMDdt
c
c    Rescaling of the SMD velocity if integrate /= VERLET
c
             do j = 1, nkey
             next = 1
             record = keyline(j)
             call gettext (record,keyword,next)
             call upcase (keyword)
             string = record(next:240)
             if (keyword(1:11) .eq. 'INTEGRATOR ') then
                call getword (record,integrate,next)
                call upcase (integrate)
             end if
             end do
             if ((integrate .eq. 'RESPA') .or. (integrate .eq. 
     $       'BAOABRESPA')) then
                dalt = 0.001_ti_p
                do j = 1, nkey
                   next = 1
                   record = keyline(j)
                   call gettext (record,keyword,next)
                   call upcase (keyword)
                   string = record(next:240)
                   if (keyword(1:7) .eq. 'DSHORT ') then
                      read (string,*,err=10) dalt
                   end if
                end do
                nalt = int(abs(SMDdt)/(dalt+2*dalt*prec_eps)) + 1
                naltsmd = nalt
                SMDdt = SMDdt/nalt
             end if
             if ((integrate .eq. 'RESPA1') .or. (integrate .eq. 
     $       'BAOABRESPA1')) then
                dalt  = 0.002
                dalt2 = 0.00025
                do j = 1, nkey
                   next = 1
                   record = keyline(j)
                   call gettext (record,keyword,next)
                   call upcase (keyword)
                   string = record(next:240)
                   if (keyword(1:7) .eq. 'DSHORT ') then
                      read (string,*,err=10) dalt
                   end if
                end do
                do j = 1, nkey
                   next = 1
                   record = keyline(j)
                   call gettext (record,keyword,next)
                   call upcase (keyword)
                   string = record(next:240)
                   if (keyword(1:7) .eq. 'DINTER ') then
                      read (string,*,err=10) dalt2
                   end if
                end do
                nalt = int(abs(SMDdt)/(dalt+2*dalt*prec_eps)) + 1
                nalt2 = int(dalt/(dalt2+2*dalt*prec_eps)) + 1
                naltsmd = nalt*nalt2
                SMDdt = SMDdt/(nalt*nalt2)
             end if
             endif
         end do
c
c    Number of atoms assigned to the center of mass (only to the 
c    constant velocity SMD) or number of atoms assigned to the 
c    constant force (only to the constant force SMD)
c
c    set defaults for the numbers and lists of SMD atoms
c
         ncsmd = -1000
         cpt1 = 1
         allocate (tcsmd(n))
!$acc enter data create(tcsmd)
         do i = 1, n
              tcsmd(i) = 0
         end do
         do i = 1, nkey
             next = 1
             record = keyline(i)
             call gettext (record,keyword,next)
             call upcase (keyword)
             string = record(next:240)
             if (keyword(1:9) .eq. 'SMDATOMS ') then
                read (string,*,err=10) ncsmd, j
             end if
         end do
c
c      1st case: SMD atoms are depicted one by one
c
         if (j .gt. 0) then
         do i = 1, nkey
             next = 1
             record = keyline(i)
             call gettext (record,keyword,next)
             call upcase (keyword)
             string = record(next:240)
             if (keyword(1:9) .eq. 'SMDATOMS ') then
                 read (string,*,err=10) ncsmd, (tcsmd(j), j=1,ncsmd)
             end if
         end do
         end if
c
c      2nd case: SMD atoms are depicted like a group: -1 100
c      => Useful to use it for a large number of SMD atoms
c
         if (j .lt. 0) then
         do i = 1, nkey
             next = 1
             record = keyline(i)
             call gettext (record,keyword,next)
             call upcase (keyword)
             string = record(next:240)
             if (keyword(1:9) .eq. 'SMDATOMS ') then
                 read (string,*,err=10) ncsmd, j, k
             end if
         end do
         jj = -j
         ii = 1
         do i = jj, k
             tcsmd(ii) = i
             ii = ii + 1
         end do
         end if
c
c    Manual parametrization of the center of mass (only for the constant
c    velocity SMD
c
         if (use_smd_velconst) then
         do i = 1, nkey
             next = 1
             record = keyline(i)
             call gettext (record,keyword,next)
             call upcase (keyword)
             string = record(next:240)
             if (keyword(1:10) .eq. 'MANUALCOM ') then
                read (string,*,err=10) xcom, ycom, zcom
                calccom = .false.
                manualcom = .true.
             end if
         end do
         end if
c
c    Direction of the steered motion (own to the both SMDs)
c
         do i = 1, nkey
             next = 1
             record = keyline(i)
             call gettext (record,keyword,next)
             call upcase (keyword)
             string = record(next:240)
             if (keyword(1:7) .eq. 'SMDDIR ') then
                read (string,*,err=10) xdir, ydir, zdir
             endif
         end do
c
c    Atom to print in a specific file during the simulation (own to the both SMDs)
c
         do i = 1, nkey
             next = 1
             record = keyline(i)
             call gettext (record,keyword,next)
             call upcase (keyword)
             string = record(next:240)
             if (keyword(1:7) .eq. 'SMDFOL ') then
                if (associated(tabatfol)) deallocate(tabatfol)
                allocate(tabatfol(n))
                use_atfol = .true.
             endif
         end do
         if (use_atfol) then
             do i = 1, ncsmd
                tabatfol(i) = tcsmd(i)
             end do
             atfol = ncsmd
         else
             use_atfol = .false.
         end if
c
c    Assign the frequency output of the data (own to the both SMDs)
c
         do i = 1, nkey
             next = 1
             record = keyline(i)
             call gettext (record,keyword,next)
             call upcase (keyword)
             string = record(next:240)
             if (keyword(1:14) .eq. 'SMDOUTPUTFREQ ') then
                read (string,*,err=10) SMDoutputFreq
             endif
         end do
c
c     Need to adjust it if RESPA, RESPA1, BAOABRESPA or
c     BAOABRESPA1 is used
c
         if ((integrate .eq. 'RESPA') .or. (integrate .eq. 
     $   'BAOABRESPA')) SMDoutputFreq = SMDoutputFreq*nalt
         if ((integrate .eq. 'RESPA1') .or. (integrate .eq. 
     $   'BAOABRESPA1')) SMDoutputFreq = SMDoutputFreq*nalt*nalt2
c
c    Restart process of the SMD procedure
c
         do i = 1, nkey
             next = 1
             record = keyline(i)
             call gettext (record,keyword,next)
             call upcase (keyword)
             string = record(next:240)
             if (keyword(1:12) .eq. 'RESTART-SMD ') then
                read (string,*,err=10) tsmd
             endif
         end do
c
c######################################################################
c
c    Check the validity of the SMD parameters
c
         quit = .false.
   10    continue
c #### Check here if any problem of lectures in key file ####
         if (quit) then
             if (rank.eq.0) then
             write(ismdout,11)
   11        format(/,'SMDinfo: Error in the lecture of the SMD values!'
     &       ,/)
             call promosmd(3)
             call fatal
             end if
         end if
c
c #### Check the validity of mandatory variables ####
c
c Titantron
c
         if (rank.eq.0) write (ismdout,110)
  110    format(/,3x,'*************** WARNING ISSUES ***************')
c
c CVSMD case
c
         if (use_smd_velconst) then
             if (SMDk .lt. 0.0_ti_p .and. rank.eq.0) then
                write (ismdout,21)
   21           format(3x,'Problem with SMDk; Negative ',
     &          3x,'Value => has to be a positive real or nul !')
                call promosmd(3)
                call fatal
             end if
             if (SMDk == 0.0_ti_p .and. rank.eq.0) then
                write (ismdout,22) 
   22           format(3x,'Problem with SMDk: equal to 0',
     &  /,3x,'=> Check the key file !')
                call promosmd(3)
                call fatal
             end if
             if (SMDVel == 0.0_ti_p .and. rank.eq.0) then
                write (ismdout,23) 
   23           format(3x,'WARNING about SMDVel:',
     & 'equal to 0',/,3x,'=> COM will not move !',
     & /,3x,'Note that : ',/,5x,'1/ If only k1, so it has no sens ...',
     & /,5x,'2/ If k1 AND k2, it corresponds to',/,5x, 
     & 'a restrain-position on the COM')
                warning = .true.
c => We find here a single position restrain !
             end if
             if (SMDVel .lt. 0.0_ti_p .and. rank.eq.0) then
                write (ismdout,24) 
   24           format(3x,'Problem with SMDVel: Negative ',
     & /,3x,'Value => has to be a positive real or nul !')
                call promosmd(3)
                call fatal
             end if
         end if
c CFSMD case
         if (use_smd_forconst) then
             if (SMDFor == 0.0_ti_p .and. rank.eq.0) then
                write (ismdout,25) 
   25           format(3x,'Problem with SMDFor : ',
     & 'equal to 0',/,3x,'=> Check the key file !')
                call promosmd(3)
                call fatal
             end if
             if (SMDFor .lt. 0.0_ti_p .and. rank.eq.0) then
                write (ismdout,20)
   20           format(3x,'Problem with SMDFor : ',
     & 'has to be upper than 0',/,3x,'=> Check the key file !')
                call promosmd(3)
                call fatal
             end if
         end if
c CVSMD and CFSMD both cases

         if (ncsmd .eq. -1000 .and. rank.eq.0) then
             write (ismdout,26)
   26        format(3x,'Problem with SMDatoms :',/,3x,
     &  'Missing keyword in the key file !',/,3x,
     &  '=> Check the key file !')
             call promosmd(3)
             call fatal
         end if
         if (ncsmd .le. 0 .and. rank.eq.0) then
             write (ismdout,27) 
   27        format(3x,'Problem with SMDatoms :',/,3x, 
     &  'the first number has to be strictely positive !',/,3x, 
     &  '=> Check the key file !')
             call promosmd(3)
             call fatal
         end if
         do i = 1, ncsmd
             if (tcsmd(i) .le. 0 .and. rank.eq.0) then
                write (ismdout,28) 
   28           format(3x,'Problem with the SMDatoms : ',/,3x,
     & 'bad number of atoms',/,3x,'=> Check the key file !')
                call promosmd(3)
                call fatal
             end if
         end do
c
c #### Check the validity of optional variables ####
c CVSMD case
         if (use_smdk2) then
             if (SMDk2 .lt. 0.0_ti_p .and. rank.eq.0) then
                write (ismdout,200) 
  200           format(3x,'Problem with SMDk2 : has to be',
     & 'positive',/,3x,'=> Check the key file !')
                call promosmd(3)
                call fatal
             end if
         end if
         if (tsmd .lt. 0 .and. rank.eq.0) then
             write (ismdout,201)
  201        format(3x,'Problem with Restart-SMD: '/,3x,'has to be 
     & strictely positive',/,3x,'=> Check the key file !')
             call promosmd(3)
             call fatal
         end if
c CVSMD and CFSMD both cases
         if (SMDdt .le. 0 .and. rank.eq.0) then
             write (ismdout,202)
  202        format(3x,'Problem with SMDdt: '/,3x,'has to be 
     & strictely positive',/,3x,'=> Check the key file !')
             call promosmd(3)
             call fatal
         end if
         if (SMDoutputFreq .le. 0 .and. rank.eq.0) then
             write (ismdout,203) 
  203        format(3x,'Problem with SMDoutputFreq : ',/,3x,'has to be
     & strictely positive',/,3x,'=> Check the key file !')
             call promosmd(3)
             call fatal
         end if
c
c  End of the check point
c
         if ((warning .eqv. .false.) .and. (rank .eq. 0)) then
             write (ismdout,207)
  207        format(3x,'No WARNING was observed',/,3x, 
     &       '=> all parameters seem to be good :D')
             call promosmd(2)
         end if
         if (rank.eq.0) write (ismdout,208)
  208        format(3x,46('*'),/,/)
c
c#############################################################
c    Sum of the used parameters within the smd output file
c#############################################################
c
       if (rank.eq.0) then
         write (ismdout,500) 
  500    format(3x,44('*'),/,
     &   3x,10('*'),' Initial SMD parameters ',10('*'),/,
     &   3x,44('*'),/)
c
c    Type of SMD used (CVSMD or CFSMD)
c
         if (use_smd_velconst) write (ismdout,501) 
  501    format(/,3x,'Type of SMD used : Constant velocity SMD')
         if (use_smd_forconst) write (ismdout,502) 
  502    format(/,3x,'Type of SMD used : Constant force SMD')
c
c    Type of integrator used (and rescaling if RESPA is used)
c
         write (ismdout,503) integrate 
  503    format(/,3x,'Integrator used : ', a)
         if (integrate .eq. 'RESPA' .or. integrate .eq. 
     $   'BAOABRESPA') then
             write (ismdout,504) 
  504        format(6x,'Rescaling of the SMDdt',
     $ ' according to the time step used for the RESPA and BAOABRESPA', 
     $ ' simulation ! ')
             write (ismdout,505) SMDdt*nalt
  505        format(6x,'SMDdt chosen : ',f12.6) 
             write (ismdout,506) nalt
  506        format(6x,'RESPA/BAOABRESPA factor : ',i2) 
             write (ismdout,507) SMDdt
  507        format(6x,'New SMDdt : ', f12.6)
         end if
         if (integrate .eq. 'RESPA1' .or. integrate .eq. 
     $   'BAOABRESPA1') then
             write (ismdout,5072)
 5072        format(6x,'Rescaling of the SMDdt',
     $ ' according to the time step used for the RESPA1 and', 
     $ ' BAOABRESPA1 simulation ! ')
             write (ismdout,5073) SMDdt*nalt*nalt2
 5073        format(6x,'SMDdt chosen : ',f12.6)
             write (ismdout,5074) nalt, nalt2
 5074        format(6x,'RESPA1/BAOABRESPA1 factor : ',/,9x,'nalt = '
     $       ,i2,/,9x,'nalt2 = ',i2)
c             write (ismdout,5075) SMDdt
c 5075        format(6x,'New SMDdt : ', f12.6)
         end if
c
c    Restart process
c
         if(tsmd .ne. 0.0_ti_p) write (ismdout,5070) tsmd
 5070    format(/,3x,10('*'),'SMD Restart process activated !',10('*'),
     & /,3x,'Time step considered : ',f15.5,' ps')
c
c    Periodicity of the printage in the output file
c
         if (integrate .eq. 'VERLET') write (ismdout,5071) SMDoutputFreq
         if (integrate .eq. 'RESPA' .or. integrate .eq.
     $   'BAOABRESPA') then
         write (ismdout,5071) SMDoutputFreq/nalt
         end if
         if (integrate .eq. 'RESPA1' .or. integrate .eq.
     $   'BAOABRESPA1') then 
         write (ismdout,5071) SMDoutputFreq/(nalt*nalt2)
         end if
 5071    format(/,3x,'Frequency for the printing : ', i10)
c
c    Parameters for the CVSMD
c
         if (use_smd_velconst) then 
             write (ismdout,508) SMDk
  508        format(/,3x,'k = ', f6.3, ' kcal/mol.A^2')
             if (use_smdk2) write (ismdout,509) SMDk2
  509        format(3x,'Additional transverse spring added k2 = ', f6.3,
     & ' kcal/mol.A^2')
             write (ismdout,510) SMDVel
  510        format(3x,'Velocity of the pulling SMDVel = ', f12.6, 
     & ' A/ps')
             write (ismdout,511) SMDdt*nalt*nalt2
  511        format(3x,'Time step SMDdt of the SMD = ', f12.6, ' ps')
             write (ismdout,512) SMDVel*SMDdt*nalt*nalt2*10**6
  512        format(3x,'Total predicted stretching distance: ', 
     &       f12.6, ' Angstrom/ns') 
         end if
c
c    Parameters for the CFSMD
c
         if (use_smd_forconst) write (ismdout,513) SMDFor 
  513    format(3x,'SMDFor = ', f6.3, ' kcal/mol.A^2')
c
c    Parameters for the COM and fixed atom
c
         if (use_smd_velconst) then
             write (ismdout,516) 
  516        format(/,3x,'************* COM AND DIRECTION VECTOR ',
     &       'DEFINITION *************')
             if(.not. calccom) write(ismdout,517) 
  517        format(3x,'COM has been defined manually !')
             if(.not. calccom) write(ismdout,518) xcom, ycom, zcom
  518        format(3x,'Manual coordinates read :',/,3x,3f8.3)
             write (ismdout,519) ncsmd
  519        format(3x,'Number of SMD atoms considered : ', i6)
             cpt1 = 1
             do i = 1, ncsmd
                write (ismdout,520) tcsmd(cpt1)
  520           format(3x,'Atoms ID : ', i6)
                cpt1 = cpt1 + 1
             end do
         end if
         if (use_smd_forconst) then
             write (ismdout,516)
             write (ismdout,521)
  521        format(3x,'Number of SMD atoms considered : ', i6)
             cpt1 = 1
             do i = 1, ncsmd
                write (ismdout,522) tcsmd(cpt1)
  522           format(3x,'Atoms ID : ', i6)
                cpt1 = cpt1 + 1
             end do
         end if
             write (ismdout,523) xdir, ydir, zdir 
  523        format(3x,'Initial direction vector read for the SMD : ',
     &       /,3x,3f8.3)
         if (xdir.eq.0 .and. ydir.eq.0 .and. zdir.eq.0) then
             write (ismdout,524)
  524        format(3x,'WARNING about the direction vector !',/,3x, 
     &       '=> Is equal to 0, so your COM will be restrain at its',
     &       /,3x,'initial position !')
         end if
             write (ismdout,525)
  525        format(3x,63('*'))
c
c    ATFOL procedure sum up
c
         if (use_smd_velconst) then
             write (ismdout,526)
  526        format(/,3x,'*************** SMDfol ',
     & 'procedure ****************')
             if (.not. use_atfol) then
                write(ismdout,530) 
  530           format(3x,'No SMDatfol keyword has been recognized !')
             end if
             if (use_atfol) then
                write(ismdout,531) ncsmd
  531           format(3x,'Number of selected SMDatfol: ', i6)
                do i = 1, ncsmd
                    write(ismdout,532) tabatfol(i)
  532               format(3x,'Atoms ID: ', i6)
                end do
             end if
             write (ismdout,536)
  536        format(3x,49('*'),/)
         end if
       end if ! End rank 0 for writting outputfile
c
c    End of the init condition (LOOP2)
c
!$acc update device(tcsmd)
      end if
c
c    End of the init condition (LOOP1)
c
      end if
c
c######################################################################
c
c    Assignation of the corresponded atoms to the right processor
c
c######################################################################
c
c    Allocation and deallocation of tables
c
      if (init) then ! LOOP1
        if (use_smd_velconst .or. use_smd_forconst) then ! LOOP2
            if (associated(nsmdglob)) deallocate(nsmdglob)
            if (associated(dedx)) deallocate(dedx)
            if (associated(dedy)) deallocate(dedy)
            if (associated(dedz)) deallocate(dedz)
            allocate(nsmdglob(ncsmd))
            allocate(dedx(n),dedy(n),dedz(n))
!$acc enter data create(nsmdglob)
c
c      Initialization of the tables
c
            do i = 1, n
                  dedx(i) = 0.0_re_p
                  dedy(i) = 0.0_re_p
                  dedz(i) = 0.0_re_p
            end do
        end if ! LOOP2
      end if ! LOOP1
c
c    => FROM HERE, PART WHICH IS UPLOADED EVERY STEP
c       OF THE DYNAMIC PROCEDURE
c
c###############################################
c################ MPI PROCEDURE ################
c###############################################
c
c    SMD atoms
c
      if (use_smd_velconst .or. use_smd_forconst) then ! LOOP1
          smdprocprint = -1
          nsmdloc = 0
!$acc parallel loop async default(present)
!$acc&         copy(smdprocprint,nsmdloc)
          do i = 1, ncsmd
             ismd = tcsmd(i)
             if (repart(ismd).eq.rank) then
!$acc atomic capture
                nsmdloc = nsmdloc + 1
                nsmdcap = nsmdloc
!$acc end atomic
                nsmdglob(nsmdcap) = ismd
                if (i==1) then
!$acc atomic write
                smdprocprint = rank
                end if
             end if
          end do
c
c     the rank of the process that will compute the com data is sent to the master
c
         if (rank.eq.0) then
c
           allocate(reqrec(nproc))
           allocate(smdproc(nproc))
           smdproc = .false.
           do i = 1, nproc-1
             tagmpi = i + 1
             call MPI_IRECV(smdproc(i),1,MPI_LOGICAL,i,tagmpi,
     $        COMM_TINKER,reqrec(i),ierr)
           end do
           do i = 1, nproc-1
             tagmpi = i + 1
             call MPI_WAIT(reqrec(i),status,ierr)
           end do
c
           do i = 1, nproc-1
             if (smdproc(i)) then
               smdprocprint = i
             end if
           end do
c
           deallocate(reqrec,smdproc)
c
         else
           tagmpi = rank+1
           call MPI_ISEND((smdprocprint.eq.rank),1,MPI_LOGICAL,0,tagmpi,
     $     COMM_TINKER,reqsend,ierr)
           call MPI_WAIT(reqsend,status,ierr)
         end if
      end if ! LOOP1
c
c#####################################################
c### Current COM calculation for the SMD procedure ###
c#####################################################
c
      if (use_smd_velconst .or. use_smd_forconst) then ! LOOP1
          cur_xcom = 0.0_ti_p
          cur_ycom = 0.0_ti_p
          cur_zcom = 0.0_ti_p
          curmtotcom = 0.0_ti_p
!$acc parallel loop default(present) async
!$acc&    reduction(+:cur_xcom,cur_ycom,cur_zcom,curmtotcom)
          do i = 1, ncsmd
             j = tcsmd(i)
             cur_xcom   = cur_xcom + x(j)*mass(j)
             cur_ycom   = cur_ycom + y(j)*mass(j)
             cur_zcom   = cur_zcom + z(j)*mass(j)
             curmtotcom = curmtotcom + mass(j)
          end do
!$acc wait
          cur_xcom = cur_xcom/curmtotcom
          cur_ycom = cur_ycom/curmtotcom
          cur_zcom = cur_zcom/curmtotcom
      end if ! LOOP1
c
c    => HERE IS THE END OF THE UPLOADED PROCEDURE !    
c
c######################################################################
c
c    Calculation of the initial center of mass coordinates
c
c######################################################################
c
c    Initial center of mass coordinates
c
      if (init) then ! LOOP1
      if (use_smd_velconst .and. calccom) then ! LOOP2
         xcom = 0._ti_p
         ycom = 0._ti_p
         zcom = 0._ti_p
         mtotcom = 0._ti_p
         do i = 1, ncsmd
             j = tcsmd(i)
             xcom = xcom + x(j)*mass(j)
             ycom = ycom + y(j)*mass(j)
             zcom = zcom + z(j)*mass(j)
             mtotcom = mtotcom + mass(j)
         end do 
         xcom = xcom/mtotcom
         ycom = ycom/mtotcom
         zcom = zcom/mtotcom
c
c    Update the coordinates within the smd outputfile
c
         if (rank.eq.0) then
             write (ismdout,525)
             write (ismdout,1000) xcom, ycom, zcom 
 1000        format(3x,'No manual COM has been specified.',
     & /,3x,'An automatic COM is generated according to the COM of',
     & ' each'/,3x,'SMD atoms considered.',
     & /,3x,'Automatic coordinates of the SMD COM : ',/,5x,
     &       3f12.3)
             write (ismdout,525) 
         end if
      end if ! LOOP2
      if (use_smd_velconst .and. .not. calccom) then ! LOOP2
         do i = 1, ncsmd
             j = tcsmd(i)
             mtotcom = mtotcom + mass(j)
         end do
      end if ! LOOP2
      end if ! LOOP1
c
c    Normalisation of the steered motion vector
c
      if (init) then ! LOOP1
        if (    use_smd_velconst
     &    .or.use_smd_forconst .and. calccom) then ! LOOP2
          r = sqrt(xdir*xdir + ydir*ydir + zdir*zdir)
          if (r .ne. 1 .and. r .ne. 0) then ! LOOP3
             xdir = xdir/r 
             ydir = ydir/r
             zdir = zdir/r
             if (rank.eq.0) then 
                 write (ismdout,1001)
 1001            format(/,3x,'********** NORMALISATION OF THE SMD ',
     &           'VECTOR **********')
                 write (ismdout,1002) 1.0/sqrt(r)
 1002            format(3x,'Normalization of the motion directory ',
     &           'vector done, factor = ',f6.3)
                 write (ismdout,1003) xdir, ydir, zdir
 1003            format(3x,'Normalized motion directory vector: ',/,5x,
     &           3g16.9)
                 write (ismdout,1004)
 1004            format(3x,'***************************************',
     &           '**************')
             end if
          else
             if (rank.eq.0) then
                 write (ismdout,1001)
                 write (ismdout,1005)
                 write (ismdout,1004)
 1005            format(3x,'No normalization of the direction vector',
     & ' is needed !')
             end if
          end if ! LOOP3
c
c    End of the conditional loop and the initializationn
c
        end if ! LOOP2
      end if ! LOOP1
      return
      end

      subroutine ksmd_recom
      use atoms    ,only:x,y,z
      use atmtyp   ,only:mass
      use tinheader,only:ti_p
      use potent   ,only:use_smd_velconst,use_smd_forconst
      use msmd
      implicit none
      integer i,j
c
c#####################################################
c### Current COM calculation for the SMD procedure ###
c#####################################################
c
      if (use_smd_velconst .or. use_smd_forconst) then ! LOOP1
          cur_xcom   = 0.0_re_p
          cur_ycom   = 0.0_re_p
          cur_zcom   = 0.0_re_p
          curmtotcom = 0.0_re_p
!$acc parallel loop default(present) async
!$acc&    reduction(+:cur_xcom,cur_ycom,cur_zcom,curmtotcom)
          do i = 1, ncsmd
             j = tcsmd(i)
             cur_xcom   = cur_xcom + x(j)*mass(j)
             cur_ycom   = cur_ycom + y(j)*mass(j)
             cur_zcom   = cur_zcom + z(j)*mass(j)
             curmtotcom = curmtotcom + mass(j)
          end do
!$acc wait
          cur_xcom = cur_xcom/curmtotcom
          cur_ycom = cur_ycom/curmtotcom
          cur_zcom = cur_zcom/curmtotcom
      end if ! LOOP1
      end
