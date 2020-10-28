c
c           ##          ##    ##    ###
c          #  #         # #  # #    #  ##
c         #    #        #  ##  #    #    #
c        ########       #      #    #    #
c       #        #      #      #    #    #
c      #          #     #      #    #  ##
c     #            #    #      #    ###
c
c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  kamd.f  --  parameters                                     ## 
c     ##                  for aMD procedure                          ##
c     ##                                                             ##
c     #################################################################
c
c
c     "kamd" checks the aMD activation and assign the necessary
c     parameters enters in the key file
c
c
      subroutine kamd(init)
      use atoms
      use deriv
      use domdec
      use energi
      use iounit
      use keys
      use mamd
      use mdstuf
      use mpi
      use potent
      use sizes
      use tors
      implicit none
      character(len=100)::ligne
      character(len=100)::file_name1
      character(len=3)::cpt
      integer i,j,k,cpt1
      integer ii, jj
      real*8 eps,dalt,dalt2
      integer nalt,nalt2
      integer next
      logical init, quit, warning
      character*20 keyword
      character*240 record
      character*240 string
c
      if (init) then    ! LOOP INIT
c
c     Initialisation parameters
c
      use_amd_dih = .false.
      use_amd_ene = .false.
      use_amd_wat1 = .false.
      use_gamd = .false.
      amd_dih_ene = 0.0d0
      amd_dih_alpha = 1.0
      amd_ep_ene = 0.0d0
      amd_ep_alpha = 1.0
      amd_factor_dih = 1.0
      amd_factor_tot = 1.0
      amdoutputfreq = 1
      amdoutputfreq_dih = 1
      amdboostavg = 0.0d0
      amdboostavg_dih = 0.0d0
      aMDdt = 0.001
      aMDdt_dih = 0.001
      tamd = 0
      amdtpass = 0
      amdtpass_dih = 0
      nalt = 1
      nalt2 = 1
      gamd_ie = 1
      gamd_cmdprepsteps = 200000
      gamd_cmdsteps = 1000000
      gamd_eqprepsteps = 200000
      gamd_eqsteps = 1000000
      gamd_sigma0P = 6.0
      gamd_sigma0D = 6.0
      use_gamd_restart = .false.
      gamd_restartfile = "GaMD_restart.dat"
      VavgD = 0.0d0
      VavgP = 0.0d0
      VavgW1 = 0.0d0
      M2D = 0.0d0
      M2P = 0.0d0
      M2W1 = 0.0d0
      sigmaVD = 0.0d0
      sigmaVP = 0.0d0
      sigmaVW1 = 0.0d0
      cptgamdD = 1
      cptgamdP = 1
      cptgamdW1 = 1
      gamdED = 0.0d0
      gamdEP = 0.0d0
      gamdEW1 = 0.0d0
      gamdkD = 0.0d0
      gamdkP = 0.0d0
      gamdkW1 = 0.0d0
      gamd_factor_dih = 1.0
      gamd_factor_tot = 1.0
      gamd_factor_wat1 = 1.0
      aMDwattype(1) = 349
      aMDwattype(2) = 350
      quit = .true.
      warning = .false.
!$acc update device(aMDwattype)
c
c######################################################################
c
c     Check the aMD key activation for:
c               - a dihedral potential boost
c               - a potential energy boost
c               - both dihedral and potential boost
c               - a bond/angle potential boost on water
c               - a potential energy boost on water
c               - use or not of the GaMD boost potential
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:240)
         if (keyword(1:12) .eq. 'AMDDIHEDRAL ') then
            use_amd_dih = .true.
            iamdout = 500
            if (rank.eq.0) then
                open (iamdout,file='aMD_output.dat',action="write")
                call promoamd(1)
            end if
         end if
         if (keyword(1:10) .eq. 'AMDENERGY ') then
            use_amd_ene = .true.
            iamdout = 500
            if (.not. use_amd_dih .and. rank.eq.0) then
                open (iamdout,file='aMD_output.dat',action="write")
                call promoamd(1)
            end if
         end if
         if (keyword(1:8) .eq. 'AMDDUAL ') then
            use_amd_dih = .true.
            use_amd_ene = .true.
            iamdout = 500
            if (rank.eq.0) then
                open (iamdout,file='aMD_output.dat',action="write")
                call promoamd(1)
            end if
         end if
         if (keyword(1:10) .eq. 'AMDWATER1 ') then
            use_amd_wat1 = .true.
            iamdout = 500
            if (.not. use_amd_dih .and. rank.eq.0) then
                open (iamdout,file='aMD_output.dat',action="write")
                call promoamd(1)
            end if
         end if
         if (keyword(1:5) .eq. 'GAMD ') then
            use_gamd = .true.
            iamdout = 500
            igamdrestart = 501
         end if
      end do
c
c######################################################################
c
c     Assign the desired parameters for the aMD (LOOP2)
c
      if (use_amd_dih .or. use_amd_ene .or. use_amd_wat1) then
c
c     Energetic threshold asignment for aMDdihE (aMD on dihedrals)
c
         do i = 1, nkey
             next = 1
             record = keyline(i)
             call gettext (record,keyword,next)
             call upcase (keyword)
             string = record(next:240)
             if (keyword(1:8) .eq. 'AMDDIHE ') then
                read (string,*,err=10) amd_dih_ene
             endif
         end do
c
c     Energetic threshold asignment for aMDePE (aMD on potential energy)
c
         do i = 1, nkey
             next = 1
             record = keyline(i)
             call gettext (record,keyword,next)
             call upcase (keyword)
             string = record(next:240)
             if (keyword(1:7) .eq. 'AMDEPE ') then
                read (string,*,err=10) amd_ep_ene
             endif
         end do
c
c     alpha acceleration factor aMDdihalpha (aMD on dihedrals)
c
         do i = 1, nkey
             next = 1
             record = keyline(i)
             call gettext (record,keyword,next)
             call upcase (keyword)
             string = record(next:240)
             if (keyword(1:12) .eq. 'AMDDIHALPHA ') then
                read (string,*,err=10) amd_dih_alpha
             endif
         end do
c     
c     alpha acceleration factor aMDePalpha (aMD on potential energy)
c 
         do i = 1, nkey
             next = 1
             record = keyline(i)
             call gettext (record,keyword,next)
             call upcase (keyword)
             string = record(next:240)
             if (keyword(1:11) .eq. 'AMDEPALPHA ') then
                read (string,*,err=10) amd_ep_alpha
             endif
         end do
c     
c     aMD outputFreq
c     
         do i = 1, nkey
             next = 1
             record = keyline(i)
             call gettext (record,keyword,next)
             call upcase (keyword)
             string = record(next:240)
             if (keyword(1:14) .eq. 'AMDOUTPUTFREQ ') then
                read (string,*,err=10) amdoutputFreq
             end if
         end do
         if (use_amd_dih) amdoutputFreq_dih = amdoutputFreq
         if (use_amd_wat1) amdoutputFreq_dih = amdoutputFreq
c
c     Reading of the timestep
c
         do i = 1, nkey
             next = 1
             record = keyline(i)
             call gettext (record,keyword,next)
             call upcase (keyword)
             string = record(next:240)
             if (keyword(1:6) .eq. 'AMDDT ') then
                read (string,*,err=10) aMDdt
             end if
         end do
         if (use_amd_dih) aMDdt_dih = aMDdt
c
c     Rescaling of the aMDoutputFreq_dih for dihedrals
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
     $   'BAOABRESPA')) then
            eps =  0.00000001d0
            dalt = 0.00025
            nalt = int(abs(aMDdt)/(dalt+eps)) + 1
            aMDdt = aMDdt/nalt
         end if
         if ((integrate .eq. 'RESPA1') .or. (integrate .eq.
     $   'BAOABRESPA1')) then
            eps =  0.00000001d0
            dalt = 0.002
            dalt2 = 0.00025
            do j = 1, nkey
               next = 1
               record = keyline(j)
               call gettext (record,keyword,next)
               call upcase (keyword)
               string = record(next:240)
               if (keyword(1:7) .eq. 'DINTER ') then
                  read (string,*,err=10) dalt
               end if
            end do
            do j = 1, nkey
               next = 1
               record = keyline(j)
               call gettext (record,keyword,next)
               call upcase (keyword)
               string = record(next:240)
               if (keyword(1:7) .eq. 'DSHORT ') then
                  read (string,*,err=10) dalt2
               end if
            end do
            nalt = int(abs(aMDdt)/(dalt+eps)) + 1
            nalt2 = int(dalt/(dalt2+eps)) + 1
            aMDdt = aMDdt/(nalt*nalt2)
         end if
c
c     Reading of the number of water molecules
c
         if (use_amd_wat1) then
         ncamd = -1000
         cpt1 = 1
         allocate (tcamd(n))
         do i = 1, n
              tcamd(i) = 0
         end do
         do i = 1, nkey
             next = 1
             record = keyline(i)
             call gettext (record,keyword,next)
             call upcase (keyword)
             string = record(next:240)
             if (keyword(1:12) .eq. 'AMDWATATOMS ') then
                read (string,*,err=10) ncamd, j
             end if
         end do
         if (j .lt. 0) then
         do i = 1, nkey
             next = 1
             record = keyline(i)
             call gettext (record,keyword,next)
             call upcase (keyword)
             string = record(next:240)
             if (keyword(1:12) .eq. 'AMDWATATOMS ') then
                 read (string,*,err=10) ncamd, j, k
             end if
         end do
         jj = -j
         ii = 1
         do i = jj, k
             tcamd(ii) = i
             ii = ii + 1
         end do
         end if
         end if
c
c     Reading of the type of water molecules
c
         if (use_amd_wat1) then
         do i = 1, nkey
             next = 1
             record = keyline(i)
             call gettext (record,keyword,next)
             call upcase (keyword)
             string = record(next:240)
             if (keyword(1:12) .eq. 'AMDWATTYPES ') then
                read (string,*,err=10) amdwattype(1), amdwattype(2)
             end if
         end do
!$acc update device(aMDwattype)
         end if
c
         if ((integrate .eq. 'RESPA') .or. (integrate .eq.
     $   'BAOABRESPA')) aMDoutputFreq_dih = aMDoutputFreq_dih*nalt
         if ((integrate .eq. 'RESPA1') .or. (integrate .eq.
     $   'BAOABRESPA1')) aMDoutputFreq_dih =aMDoutputFreq_dih*nalt*nalt2
c
c########################################
c#### Specific of the GaMD procedure ####
c########################################
c
c     Threshold energy for adding boost potential
c
         do i = 1, nkey
             next = 1
             record = keyline(i)
             call gettext (record,keyword,next)
             call upcase (keyword)
             string = record(next:240)
             if (keyword(1:7) .eq. 'GAMDIE ') then
                read (string,*,err=10) gamd_ie
             endif
         end do
c
c     1/ Preparation steps for GaMD procedure
c
         do i = 1, nkey
             next = 1
             record = keyline(i)
             call gettext (record,keyword,next)
             call upcase (keyword)
             string = record(next:240)
             if (keyword(1:18) .eq. 'GAMD-CMDPREPSTEPS ') then
                read (string,*,err=10) gamd_cmdprepsteps
             endif
         end do
c
c     2/ cMD steps
c
         do i = 1, nkey
             next = 1
             record = keyline(i)
             call gettext (record,keyword,next)
             call upcase (keyword)
             string = record(next:240)
             if (keyword(1:14) .eq. 'GAMD-CMDSTEPS ') then
                read (string,*,err=10) gamd_cmdsteps
             endif
         end do
c
c     3/ Preparatory equilibration steps
c
         do i = 1, nkey
             next = 1
             record = keyline(i)
             call gettext (record,keyword,next)
             call upcase (keyword)
             string = record(next:240)
             if (keyword(1:17) .eq. 'GAMD-EQPREPSTEPS ') then
                read (string,*,err=10) gamd_eqprepsteps
             endif
         end do
c
c     4/ Total equilibration steps
c
         do i = 1, nkey
             next = 1
             record = keyline(i)
             call gettext (record,keyword,next)
             call upcase (keyword)
             string = record(next:240)
             if (keyword(1:13) .eq. 'GAMD-EQSTEPS ') then
                read (string,*,err=10) gamd_eqsteps
             endif
         end do
c
c     Upper limit of SD of the potential energy
c
         do i = 1, nkey
             next = 1
             record = keyline(i)
             call gettext (record,keyword,next)
             call upcase (keyword)
             string = record(next:240)
             if (keyword(1:13) .eq. 'GAMD-SIGMA0P ') then
                read (string,*,err=10) gamd_sigma0P
             endif
         end do
c
c     Upper limit of SD of the dihedral potential boost
c
         do i = 1, nkey
             next = 1
             record = keyline(i)
             call gettext (record,keyword,next)
             call upcase (keyword)
             string = record(next:240)
             if (keyword(1:13) .eq. 'GAMD-SIGMA0D ') then
                read (string,*,err=10) gamd_sigma0D
             endif
         end do
c
c     Upper limit of SD of the bond/angle water potential boost
c
         do i = 1, nkey
             next = 1
             record = keyline(i)
             call gettext (record,keyword,next)
             call upcase (keyword)
             string = record(next:240)
             if (keyword(1:14) .eq. 'GAMD-SIGMA0W1 ') then
                read (string,*,err=10) gamd_sigma0W1
             endif
         end do
c
c     Restart procedure
c
         do i = 1, nkey
             next = 1
             record = keyline(i)
             call gettext (record,keyword,next)
             call upcase (keyword)
             string = record(next:240)
             if (keyword(1:13) .eq. 'GAMD-RESTART ') then
                use_gamd_restart = .true.
             endif
         end do
c
c     Restart file name
c
         if (use_gamd_restart) then
         do i = 1, nkey
             next = 1
             record = keyline(i)
             call gettext (record,keyword,next)
             call upcase (keyword)
             string = record(next:240)
             if (keyword(1:18) .eq. 'GAMD-RESTARTFILE ') then
                read (string,*,err=10) gamd_restartfile
                igamdrestart = 5
                open (igamdrestart,file=gamd_restartfile,action="read",
     &          status='old')
                read (igamdrestart,*)
                if (use_amd_dih) then 
                   read (igamdrestart,*) tamd 
                   read (igamdrestart,*) VminD
                   read (igamdrestart,*) VmaxD 
                   read (igamdrestart,*) VavgD
                   read (igamdrestart,*) sigmaVD
                   read (igamdrestart,*) M2D
                   read (igamdrestart,*) gamdED
                   read (igamdrestart,*) gamdk0D
                   read (igamdrestart,*) gamdkD
                   read (igamdrestart,*) amdtpass
                   read (igamdrestart,*) cptgamdD
                   amdtpass_dih = amdtpass
                end if
                if (.not. use_amd_dih .and. use_amd_ene) then 
                   read (igamdrestart,*) tamd 
                   read (igamdrestart,*) VminP 
                   read (igamdrestart,*) VmaxP
                   read (igamdrestart,*) VavgP
                   read (igamdrestart,*) sigmaVP
                   read (igamdrestart,*) M2P
                   read (igamdrestart,*) gamdEP
                   read (igamdrestart,*) gamdk0P
                   read (igamdrestart,*) gamdkP
                   read (igamdrestart,*) amdtpass
                   read (igamdrestart,*) cptgamdP
                end if
                if (use_amd_wat1 .and. use_amd_dih) then
                   read (igamdrestart,*)
                   read (igamdrestart,*) tamd
                   read (igamdrestart,*) VminW1
                   read (igamdrestart,*) VmaxW1
                   read (igamdrestart,*) VavgW1
                   read (igamdrestart,*) sigmaVW1
                   read (igamdrestart,*) M2W1
                   read (igamdrestart,*) gamdEW1
                   read (igamdrestart,*) gamdk0W1
                   read (igamdrestart,*) gamdkW1
                   read (igamdrestart,*) amdtpass
                   read (igamdrestart,*) cptgamdW1
                end if
                if (use_amd_wat1 .and. .not. use_amd_dih) then
                   read (igamdrestart,*) tamd
                   read (igamdrestart,*) VminW1
                   read (igamdrestart,*) VmaxW1
                   read (igamdrestart,*) VavgW1
                   read (igamdrestart,*) sigmaVW1
                   read (igamdrestart,*) M2W1
                   read (igamdrestart,*) gamdEW1
                   read (igamdrestart,*) gamdk0W1
                   read (igamdrestart,*) gamdkW1
                   read (igamdrestart,*) amdtpass
                   read (igamdrestart,*) cptgamdW1
                end if
                if (use_amd_dih .and. use_amd_ene) then
                   read (igamdrestart,*)
                   read (igamdrestart,*) tamd
                   read (igamdrestart,*) VminP
                   read (igamdrestart,*) VmaxP
                   read (igamdrestart,*) VavgP
                   read (igamdrestart,*) sigmaVP
                   read (igamdrestart,*) M2P
                   read (igamdrestart,*) gamdEP
                   read (igamdrestart,*) gamdk0P
                   read (igamdrestart,*) gamdkP
                   read (igamdrestart,*) amdtpass
                   read (igamdrestart,*) cptgamdP
                end if
                close(igamdrestart)
             endif
         end do
         end if
c
c######################################################################
c
c    Check the validity of the aMD parameters
c
         quit = .false.
   10    continue
c #### Check here if any problem of lectures in key file ####
         if (quit) then
             if (rank.eq.0) then
             write(iamdout,11)
   11        format(/,'aMDinfo: Error in the lecture of the aMD values!'
     &       ,/)
             call promoamd(3)
             call fatal
             end if
         end if
c
c #### Check the validity of mandatory variables ####
c
         if (rank.eq.0) then            ! LOOP1
         write (iamdout,110)
  110    format(/,3x,'*************** WARNING ISSUES ***************')
c
         if (use_amd_dih) then
         if (amd_dih_alpha .le. 0.0d0) then
             write (iamdout,21)
   21        format(3x,'Problem with aMDdihalpha; Negative ',
     &       3x,'Value => has to be a positive real !')
             call promoamd(3)
             call fatal
         end if
         end if
         if (use_amd_ene) then
         if (amd_ep_alpha .le. 0.0d0) then
             write (iamdout,22)
   22        format(3x,'Problem with aMDePalpha; Negative ',
     &       3x,'Value => has to be a positive real !')
             call promoamd(3)
             call fatal
         end if
         if (integrate .eq. "RESPA" .or. integrate .eq. "RESPA1"
     &   .or. integrate .eq. "BAOABRESPA" .or. integrate .eq. 
     &   "BAOABRESPA1") then
             write (iamdout,23)
   23        format(3x,'aMD and GaMD are not available ', 
     &       /,3x,'with multitimsteping integrators !')
             call promoamd(3)
             call fatal
         end if
         end if
         if (amdoutputfreq .le. 0) then
             write (iamdout,25)
   25        format(3x,'Problem with aMDoutputFreq; Negative ',
     &       3x,'Value => has to be a positive integer !')
             call promoamd(3)
             call fatal
         end if
c
c  End of the check point
c
         if (.not. warning) then
             write (iamdout,26)
   26        format(3x,'No WARNING was observed',/,3x,
     &       '=> all parameters seem to be good :D')
             call promoamd(2)
         end if
         write (iamdout,27)
   27    format(3x,46('*'),/,/)
c
         end if         ! END LOOP1
c
c#############################################################
c    Sum of the used parameters within the smd output file
c#############################################################
c
         if (rank.eq.0) then            ! LOOP1
         write (iamdout,499)
  499    format(3x,44('*'),/,
     &   3x,10('*'),' Initial aMD parameters ',10('*'),/,
     &   3x,44('*'),/)
c
c    Type of aMD used (aMDdihE / aMDePE / dual / water / GaMD)
c
         write (iamdout,500)
  500    format(3x,"#####################")
         if (use_amd_dih) cpt = "YES"
         if (.not. use_amd_dih) cpt = "NO"
         write (iamdout,501) cpt
  501    format(3x,'Use of aMDdihE : ',a)
         if (use_amd_ene) cpt = "YES"
         if (.not. use_amd_ene) cpt = "NO"
         write (iamdout,502) cpt
  502    format(3x,'Use of aMDePE : ',a)
         if (use_amd_dih .and. use_amd_ene) then
             write (iamdout,503)
  503        format(3x,'aMD is used in a dual mode !')
         end if
         if (use_amd_wat1) then
             write (iamdout,5003)
 5003        format(3x,'aMD is used on bonds/angles ',/,
     $       3x,'of water molecules')
         end if
         if (use_gamd) then
             write (iamdout,504)
  504        format(3x,'Use of GaMD : YES',/,
     &       3x,'GaMDinfo: GaMD potential will be applied instead ',
     &       /,3x,'of the aMD classical boost !')
         else
             write (iamdout,505)
  505        format(3x,'Use of GaMD : NO')
         end if
         write (iamdout,500)
c
c    Restart procedure for GaMD
c
         if (use_gamd_restart) then
             write (iamdout,506) gamd_restartfile
  506        format(/,3x,'#################### GaMD Restart ',
     &       'procedure activated ! ####################',
     &       /,3x,'GaMDinfo: Name of the restart file = ',a)
         end if
c
c    Parameters for each aMD
c
         if (use_amd_dih .and. .not. use_gamd) 
     &   write (iamdout,507) amd_dih_ene
         if (use_amd_ene .and. .not. use_gamd) 
     &   write (iamdout,508) amd_ep_ene
  507    format(/,3x,'aMDinfo: Threshold energy for dihedrals   = '
     &   ,f12.3)
  508    format(3x,'aMDinfo: Threshold energy for potential energy = '
     &   ,f12.3)
         if (use_amd_dih .and. .not. use_gamd) 
     &   write (iamdout,509) amd_dih_alpha
         if (use_amd_ene .and. .not. use_gamd) 
     &   write (iamdout,510) amd_ep_alpha
  509    format(3x,'aMDinfo: alpha acceleration for dihedrals = '
     &   ,f12.3)
  510    format(3x,'aMDinfo: alpha acceleration for potential energy = '
     &   ,f12.3)
         write (iamdout,511) amdoutputFreq
  511    format(/,3x,'aMDinfo: Frequency for the printing = ', i10)
         if (use_amd_dih) write (iamdout,512) ntors
  512    format(/,3x,'aMDinfo: Number of torsional angles = ', i10)
         if (use_amd_wat1) write(iamdout,5112) ncamd/3
 5112    format(/,3x,'aMDinfo: Number of water molecules  = ', i10)
         if (.not. use_gamd) write (iamdout,513)
  513    format(/)
         if (use_gamd) then
             write (iamdout,514)
  514        format(/,3x,"#############################################"
     &       ,"#####",/,3x,'GaMDinfo: Sum of the the GaMD parameters')
             write (iamdout,515) gamd_ie
  515        format(6x,"GaMD flag for energy threshold E   = ",i10)
             write (iamdout,516) gamd_cmdprepsteps, gamd_cmdsteps,
     &       gamd_eqprepsteps, gamd_eqsteps
  516        format(6x,"GaMD preparation steps for cMD     = ",i10,
     &              /,6x,"GaMD number of steps for cMD       = ",i10,
     &              /,6x,"GaMD preparation steps for eq GaMD = ",i10,
     &              /,6x,"GaMD number of steps for eq GaMD   = ",i10)
             if (use_amd_dih) write(iamdout,517) gamd_sigma0D
             if (use_amd_ene) write(iamdout,518) gamd_sigma0P
             if (use_amd_wat1) write(iamdout,5188) gamd_sigma0W1
  517        format(6x,"GaMD max SD allowed for dihedral   = ",f10.3)
  518        format(6x,"GaMD max SD allowed for potential  = ",f10.3)
 5188        format(6x,"GaMD max SD allowed for water1     = ",f10.3)
             write (iamdout,519)
  519        format(3x,"############################################",
     &       "######")
             write (iamdout,520)
  520        format(3x,"GaMDinfo: All the number of steps have been ",
     &       /,3x,"adapted according to the read restart file if it ",
     &       /,3x,"exists ! ",/,3x,
     &       "##################################################")
             if (use_gamd_restart .and. use_amd_dih) write (iamdout,521)
     &       tamd, VminD, VmaxD, VavgD, sigmaVD, M2D, gamdED, gamdk0D, 
     &       gamdkD
  521        format(/,3x,"#############################################"
     &       ,"#####",/,3x,'GaMDinfo: Restart procedure is activated !',
     &       /,6x,'Read parameters for dihedrals: ',/,6x,'Timestep = ',
     &       f12.3,/,6x,'VminD    = ',f12.3,/,6x,'VmaxD    = ',f12.3,/,
     &       6x,
     &       'VavgD    = ',f12.3,/,6x,'SigmaVD  = ',f12.3,/,6x,
     &       'M2D      = ',f12.3,/,6x,'ED       = ',f12.3,/,6x,
     &       'k0D      = ',f12.9,/,6x,'kD       = ',f12.9,/,
     &       3x,"############################################",
     &       "######")
             if (use_gamd_restart .and. use_amd_ene) write (iamdout,522)
     &       tamd, VminP, VmaxP, VavgP, sigmaVP, M2P, gamdEP, gamdk0P,
     &       gamdkP
  522        format(/,3x,"#############################################"
     &       ,"#####",/,3x,'GaMDinfo: Restart procedure is activated !',
     &       /,6x,'Read parameters for total potential energy: ',/,6x,
     &       'Timestep = ',f12.3
     &       ,/,6x,'VminP    = ',f12.3,/,6x,'VmaxP    = ',f12.3,/,6x,
     &       'VavgP    = ',f12.3,/,6x,'SigmaVP  = ',f12.3,/,6x,
     &       'M2P      = ',f12.3,/,6x,'EP       = ',f12.3,/,6x,
     &       'k0P      = ',f12.3,/,6x,'kP       = ',f12.3,/,
     &       3x,"############################################",
     &       "######")
             if (use_gamd_restart .and. use_amd_wat1) write(iamdout,523)
     &       tamd, VminW1, VmaxW1, VavgW1, sigmaVW1, M2W1, gamdEW1, 
     &       gamdk0W1, gamdkW1
  523        format(/,3x,"#############################################"
     &       ,"#####",/,3x,'GaMDinfo: Restart procedure is activated !',
     &       /,6x,'Read parameters for water molecules: ',/,6x,
     &       'Timestep = ',f12.3
     &       ,/,6x,'VminW    = ',f12.3,/,6x,'VmaxW    = ',f12.3,/,6x,
     &       'VavgW    = ',f12.3,/,6x,'SigmaVW  = ',f12.3,/,6x,
     &       'M2W      = ',f12.3,/,6x,'EW       = ',f12.3,/,6x,
     &       'k0W      = ',f12.3,/,6x,'kW       = ',f12.3,/,
     &       3x,"############################################",
     &       "######")
             write (iamdout,524)
  524        format(/,3x,"GaMD output file for Tinker-HP: ",/,/,
     &       13x,"STEP",4x,"iG",10x,"dV",7x,"dVavg",8x,"Vmax",8x,"Vmin",
     &       8x,"Vavg",6x,"SigmaV",11x,"E",10x,"k0",11x,"k",/)
         end if
         if (.not. use_gamd) write (iamdout,525)
  525    format(7x,"STEP",9x,"dV",8x,"dVAVG")
c
         end if         ! END LOOP1
c
      end if
c
c     Adjusting of the GaMD steps
c
      gamd_cmdprepsteps = gamd_cmdprepsteps*nalt*nalt2
      gamd_cmdsteps = gamd_cmdsteps*nalt*nalt2 + gamd_cmdprepsteps
      gamd_eqprepsteps = gamd_eqprepsteps*nalt*nalt2 + gamd_cmdsteps
      gamd_eqsteps = gamd_eqsteps*nalt*nalt2 + gamd_eqprepsteps
c
      end if            !END LOOP INIT
c
      end


