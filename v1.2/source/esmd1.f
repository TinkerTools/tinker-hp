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
c     ##  subroutine esmd1  --  make smd calculations              ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "esmd1" assigns the biased forced on each concerned atoms and
c     stores the smd and derivatives energy components
c
c     literature reference:
c     Massively Parallel Implementation of Steered Molecular Dynamics in Tinker-HP:
c     Comparisons of Polarizable and Non-Polarizable Simulations of Realistic Systems
c     F Celerse, L Lagardere, E Derat, JP Piquemal
c     Journal of Chemical Theory and Computation 15 (6), 3694-3709
c
c
      subroutine esmd1
      use msmd
      use potent
      use deriv
      use energi
      use atmtyp
      use atoms
      use domdec
      use inform
      use iounit
      use virial
      use mpi
      implicit none
      character(len=100)::file_name1
      character(len=1)::num1
      character(len=2)::num2
      character(len=3)::num3
      character(len=4)::num4
      character(len=5)::num5
      character(len=6)::num6
      character(len=7)::num7      

      real*8 diffx, diffy, diffz
      real*8 diff, diff2, vel
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8 buffers(7),buffer(7)
      integer ierr,reqrec,reqsend,STATUS(MPI_STATUS_SIZE),tagmpi
      integer a,b,c,j



      integer ib
c
      if (deb_Path) write(iout,*), 'esmd1 '
c
c
c     SMD storage: initialization
c
      stock_dedx = 0.0d0
      stock_dedy = 0.0d0
      stock_dedz = 0.0d0
      stock_ensmd = 0.0d0
c
c###########################################
c###########################################
c#### Case of the constant velocity SMD ####
c###########################################
c###########################################
c
c####################################
c#### Case where only k1 is used ####
c####################################
c
      if (use_smd_velconst .and. .not. use_smdk2) then  ! LOOP1 CVSMDK1
c
c     Initial SMD force calculation according to the current COM
c
          diffx = (cur_xcom - xcom)*xdir
          diffy = (cur_ycom - ycom)*ydir
          diffz = (cur_zcom - zcom)*zdir
          diff = diffx + diffy + diffz
c
          com_dedx = SMDk*xdir*(SMDVel*tsmd - diff)
          com_dedy = SMDk*ydir*(SMDVel*tsmd - diff)
          com_dedz = SMDk*zdir*(SMDVel*tsmd - diff)
c
c     Initial SMD energy calculation according to the current COM
c
          vel = xdir*SMDVel + ydir*SMDVel + zdir*SMDVel
          ensmd = 0.5*SMDk*(vel*tsmd - diff)**2          
c
c     Storage for printing
c     => Mandatory if integrator != VERLET (so multitimestep)
c
          stock_dedx = stock_dedx + com_dedx
          stock_dedy = stock_dedy + com_dedy
          stock_dedz = stock_dedz + com_dedz
          stock_ensmd = stock_ensmd + ensmd
c
c     the proc of rank 'smdprocprint' sends the info to the rank 0
c
          if (rank.eq.smdprocprint) then
            buffers(1) = cur_xcom
            buffers(2) = cur_ycom
            buffers(3) = cur_zcom
            buffers(4) = stock_dedx
            buffers(5) = stock_dedy
            buffers(6) = stock_dedz
            buffers(7) = curmtotcom
            tagmpi = rank+1
            call MPI_ISEND(buffers,7,MPI_REAL8,0,tagmpi,COMM_TINKER,
     $       reqsend,ierr)
            call MPI_WAIT(reqsend,status,ierr)
          else if (rank.eq.0) then
            tagmpi = smdprocprint+1
            call MPI_IRECV(buffer,7,MPI_REAL8,smdprocprint,tagmpi,
     $        COMM_TINKER,reqrec,ierr)
            call MPI_WAIT(reqrec,status,ierr)
            cur_xcom = buffer(1)
            cur_ycom = buffer(2)
            cur_zcom = buffer(3)
            stock_dedx = buffer(4)
            stock_dedy = buffer(5)
            stock_dedz = buffer(6)
            curmtotcom = buffer(7)
          end if
c
c     SMD output data
c
          if (rank==0) then
              c = modulo(tpass,SMDoutputFreq)
              if (c == 0) then          ! LOOPA
              open (ismdout,file='SMD_output.dat',position='append',
     &        action="write",status='old')
              write (ismdout,45)  tsmd,cur_xcom,cur_ycom,cur_zcom,
     &        stock_dedx,stock_dedy,stock_dedz
   45         format ('SMD',1x,f15.5,3f10.4,3f10.4)
              close(ismdout)
              end if                    ! END LOOPB
          end if
c
c     Redistribution on each SMD atoms forming the COM
c
          do a = 1, nsmdloc     ! LOOP2 CVSMDK2
              b = nsmdglob(a)
              ib = loc(b)
c
              dedx(b) = (mass(b)/curmtotcom)*com_dedx
              dedy(b) = (mass(b)/curmtotcom)*com_dedy
              dedz(b) = (mass(b)/curmtotcom)*com_dedz
c
c     Incrementation of the derivative terms
c
              desmd(1,ib) = desmd(1,ib) - dedx(b)
              desmd(2,ib) = desmd(2,ib) - dedy(b)
              desmd(3,ib) = desmd(3,ib) - dedz(b)
c
c     Increment the internal virial tensor components
c
              vxx = (cur_xcom - xcom)*dedx(b)
              vyx = (cur_ycom - ycom)*dedx(b)
              vzx = (cur_zcom - zcom)*dedx(b)
              vyy = (cur_ycom - ycom)*dedy(b)
              vzy = (cur_zcom - zcom)*dedy(b)
              vzz = (cur_zcom - zcom)*dedz(b)
              vir(1,1) = vir(1,1) + vxx
              vir(2,1) = vir(2,1) + vyx
              vir(3,1) = vir(3,1) + vzx
              vir(1,2) = vir(1,2) + vyx
              vir(2,2) = vir(2,2) + vyy
              vir(3,2) = vir(3,2) + vzy
              vir(1,3) = vir(1,3) + vzx
              vir(2,3) = vir(2,3) + vzy
              vir(3,3) = vir(3,3) + vzz
          end do
c
c##########################################
c#### Write results in the outputfile ####
c##########################################
c
          do a = 1, ncsmd
            b = tcsmd(a)
c
c     Check of the atom
c
              do j = 1, atfol   ! LOOP3 CVSMDK1
              if (b.eq.tabatfol(j)) then        ! LOOP4 CVSMDK1
c
c     Checking with respect of the SMDOUTPUTFREQ
c
                  c = modulo(tpass,SMDoutputFreq)
                  if (c == 0) then      ! LOOP5 CVSMDK1
c
c     Write in the SMD_output file
c
                  if (use_atfol) then ! LOOP 6 
                  if (b .lt. 10) then
                     write(num1,'(i1)') b
                     file_name1 = "SMD"//num1//".dat"
                     open (1000+b,file=file_name1,position='append',
     &               action="write",status='old')
                     write (1000+b,fmt='(i6,f15.5,f12.5,f12.5,
     &               f12.5,f12.5,f12.5,f12.5)')
     &               type(b), tsmd, x(b), y(b), z(b),
     &               (mass(b)/curmtotcom)*stock_dedx,
     &               (mass(b)/curmtotcom)*stock_dedy,
     &               (mass(b)/curmtotcom)*stock_dedz
                     close(1000+b)
                  end if
                  if (b .ge. 10 .and. b .lt. 100) then
                     write(num2,'(i2)') b
                     file_name1 = "SMD"//num2//".dat"
                     open (1000+b,file=file_name1,position='append',
     &               action="write",status='old')
                     write (1000+b,fmt='(i6,f15.5,f12.5,f12.5,
     &               f12.5,f12.5,f12.5,f12.5)')
     &               type(b), tsmd, x(b), y(b), z(b),
     &               (mass(b)/curmtotcom)*stock_dedx,
     &               (mass(b)/curmtotcom)*stock_dedy,
     &               (mass(b)/curmtotcom)*stock_dedz                  
                     close(1000+b)
                  end if
                  if (b .ge. 100 .and. b .lt. 1000) then
                     write(num3,'(i3)') b
                     file_name1 = "SMD"//num3//".dat"
                     open (1000+b,file=file_name1,position='append',
     &               action="write",status='old')
                     write (1000+b,fmt='(i6,f15.5,f12.5,f12.5,
     &               f12.5,f12.5,f12.5,f12.5)')
     &               type(b), tsmd, x(b), y(b), z(b),
     &               (mass(b)/curmtotcom)*stock_dedx,
     &               (mass(b)/curmtotcom)*stock_dedy,
     &               (mass(b)/curmtotcom)*stock_dedz
                     close(1000+b)
                  end if
                  if (b .ge. 1000 .and. b .lt. 10000) then
                     write(num4,'(i4)') b
                     file_name1 = "SMD"//num4//".dat"
                     open (1000+b,file=file_name1,position='append',
     &               action="write",status='old')
                     write (1000+b,fmt='(i6,f15.5,f12.5,f12.5,
     &               f12.5,f12.5,f12.5,f12.5)')
     &               type(b), tsmd, x(b), y(b), z(b),
     &               (mass(b)/curmtotcom)*stock_dedx,
     &               (mass(b)/curmtotcom)*stock_dedy,
     &               (mass(b)/curmtotcom)*stock_dedz
                     close(1000+b)
                  end if
                  if (b .ge. 10000 .and. b .lt. 100000) then
                     write(num5,'(i5)') b
                     file_name1 = "SMD"//num5//".dat"
                     open (1000+b,file=file_name1,position='append',
     &               action="write",status='old')
                     write (1000+b,fmt='(i6,f15.5,f12.5,f12.5,
     &               f12.5,f12.5,f12.5,f12.5)')
     &               type(b), tsmd, x(b), y(b), z(b),
     &               (mass(b)/curmtotcom)*stock_dedx,
     &               (mass(b)/curmtotcom)*stock_dedy,
     &               (mass(b)/curmtotcom)*stock_dedz
                     close(1000+b)
                  end if
                  if (b .ge. 100000 .and. b .lt. 1000000) then
                     write(num6,'(i6)') b
                     file_name1 = "SMD"//num6//".dat"
                     open (1000+b,file=file_name1,position='append',
     &               action="write",status='old')
                     write (1000+b,fmt='(i6,f15.5,f12.5,f12.5,
     &               f12.5,f12.5,f12.5,f12.5)')
     &               type(b), tsmd, x(b), y(b), z(b),
     &               (mass(b)/curmtotcom)*stock_dedx,
     &               (mass(b)/curmtotcom)*stock_dedy,
     &               (mass(b)/curmtotcom)*stock_dedz
                     close(1000+b)
                  end if
                  if (b .ge. 1000000 .and. b .lt. 10000000) then
                     write(num7,'(i7)') b
                     file_name1 = "SMD"//num7//".dat"
                     open (1000+b,file=file_name1,position='append',
     &               action="write",status='old')
                     write (1000+b,fmt='(i6,f15.5,f12.5,f12.5,
     &               f12.5,f12.5,f12.5,f12.5)')
     &               type(b), tsmd, x(b), y(b), z(b),
     &               (mass(b)/curmtotcom)*stock_dedx,
     &               (mass(b)/curmtotcom)*stock_dedy,
     &               (mass(b)/curmtotcom)*stock_dedz
                     close(1000+b)
                  end if
                  end if ! END LOOP 6
                  end if ! END LOOP5 CVSMDK1
              end if    ! END LOOP4 CVSMDK1
          end do        ! END LOOP3 CVSMDK1
          end do        ! END LOOP2 CVSMDK1
c
      end if    ! END LOOP1 CVSMDK1
c
c#######################################
c#### Case where k2 is used with k1 ####
c#######################################
c
      if (use_smd_velconst .and. use_smdk2) then
c
c     Initial SMD force calculation according to the current COM
c
          diffx = (cur_xcom - xcom)*xdir
          diffy = (cur_ycom - ycom)*ydir
          diffz = (cur_zcom - zcom)*zdir
          diff = diffx + diffy + diffz
c
          com_dedx = SMDk*xdir*(SMDVel*tsmd - diff)
          com_dedy = SMDk*ydir*(SMDVel*tsmd - diff)
          com_dedz = SMDk*zdir*(SMDVel*tsmd - diff)
c
c     Initial SMD energy calculation according to the current COM
c
          vel = xdir*SMDVel + ydir*SMDVel + zdir*SMDVel
          ensmd = 0.5*SMDk*(vel*tsmd - diff)**2
c
c     Calculation of the derivative components of k2
c
c     #Part 1: Contribution without the motion#
          com_dedx = com_dedx - SMDk2*(cur_xcom - xcom)
          com_dedy = com_dedy - SMDk2*(cur_ycom - ycom)
          com_dedz = com_dedz - SMDk2*(cur_zcom - zcom)
c
c     #Part 2: Contribution with the motion#
c
          com_dedx = com_dedx + (SMDk2*diff*xdir)
          com_dedy = com_dedy + (SMDk2*diff*ydir)
          com_dedz = com_dedz + (SMDk2*diff*zdir)
c
c     SMD energy calculation according to the current COM and k2
c
          diff2 = sqrt((cur_xcom - xcom)**2 + (cur_ycom - ycom)**2 + 
     &    (cur_zcom - zcom)**2)
          ensmd = ensmd + 0.5*SMDk2*(diff2 - (diff)**2)
c
c     Storage for printing
c     => Mandatory if integrator != VERLET (so multitimestep)
c
          stock_dedx = stock_dedx + com_dedx
          stock_dedy = stock_dedy + com_dedy
          stock_dedz = stock_dedz + com_dedz
          stock_ensmd = stock_ensmd + ensmd        
c
c     the proc of rank 'smdprocprint' sends the info to the rank 0
c
          if (rank.eq.smdprocprint) then
            buffers(1) = cur_xcom
            buffers(2) = cur_ycom
            buffers(3) = cur_zcom
            buffers(4) = stock_dedx
            buffers(5) = stock_dedy
            buffers(6) = stock_dedz
            buffers(7) = curmtotcom
            tagmpi = rank+1
            call MPI_ISEND(buffers,7,MPI_REAL8,0,tagmpi,COMM_TINKER,
     $       reqsend,ierr)
            call MPI_WAIT(reqsend,status,ierr)
          else if (rank.eq.0) then
            tagmpi = smdprocprint+1
            call MPI_IRECV(buffer,7,MPI_REAL8,smdprocprint,tagmpi,
     $        COMM_TINKER,reqrec,ierr)
            call MPI_WAIT(reqrec,status,ierr)
            cur_xcom = buffer(1)
            cur_ycom = buffer(2)
            cur_zcom = buffer(3)
            stock_dedx = buffer(4)
            stock_dedy = buffer(5)
            stock_dedz = buffer(6)
            curmtotcom = buffer(7)
          end if
c
c     SMD output data
c
          if (rank==0) then
              c = modulo(tpass,SMDoutputFreq)
              if (c == 0) then          ! LOOPA
              open (ismdout,file='SMD_output.dat',position='append',
     &        action="write",status='old')
              write (ismdout,48)  tsmd,cur_xcom,cur_ycom,cur_zcom,
     &        stock_dedx,stock_dedy,stock_dedz
   48         format ('SMD',1x,f15.5,3f10.4,3f10.4)
              close(ismdout)
              end if                    ! END LOOPB
          end if
c
c     Redistribution on each SMD atoms forming the COM
c
          do a = 1, nsmdloc     ! LOOP2 CVSMDK2
              b = nsmdglob(a)
              ib = loc(b)
c
              dedx(b) = (mass(b)/curmtotcom)*com_dedx
              dedy(b) = (mass(b)/curmtotcom)*com_dedy
              dedz(b) = (mass(b)/curmtotcom)*com_dedz
c
c     Incrementation of the derivative terms
c
              desmd(1,ib) = desmd(1,ib) - dedx(b)
              desmd(2,ib) = desmd(2,ib) - dedy(b)
              desmd(3,ib) = desmd(3,ib) - dedz(b)
c
c     Increment the internal virial tensor components
c
              vxx = (cur_xcom - xcom)*dedx(b)
              vyx = (cur_ycom - ycom)*dedx(b)
              vzx = (cur_zcom - zcom)*dedx(b)
              vyy = (cur_ycom - ycom)*dedy(b)
              vzy = (cur_zcom - zcom)*dedy(b)
              vzz = (cur_zcom - zcom)*dedz(b)
              vir(1,1) = vir(1,1) + vxx
              vir(2,1) = vir(2,1) + vyx
              vir(3,1) = vir(3,1) + vzx
              vir(1,2) = vir(1,2) + vyx
              vir(2,2) = vir(2,2) + vyy
              vir(3,2) = vir(3,2) + vzy
              vir(1,3) = vir(1,3) + vzx
              vir(2,3) = vir(2,3) + vzy
              vir(3,3) = vir(3,3) + vzz
          end do
c
c##########################################
c#### Writte results in the outputfile ####
c##########################################
c
          do a = 1, ncsmd
            b = tcsmd(a)
c
c     Check of the use of CVSMD
c
              do j = 1, atfol
              if (b.eq.tabatfol(j)) then
c
c     Checking with respect of the SMDOUTPUTFREQ
c
                  c = modulo(tpass,SMDoutputFreq)
                  if (c == 0) then
c
c     Write in the SMD_output file
c
                  if (use_atfol) then
                  if (b .lt. 10) then
                     write(num1,'(i1)') b
                     file_name1 = "SMD"//num1//".dat"
                     open (1000+b,file=file_name1,position='append',
     &               action="write",status='old')
                     write (1000+b,fmt='(i6,f15.5,f12.5,f12.5,
     &               f12.5,f12.5,f12.5,f12.5)')
     &               type(b), tsmd, x(b), y(b), z(b),
     &               (mass(b)/curmtotcom)*stock_dedx,
     &               (mass(b)/curmtotcom)*stock_dedy,
     &               (mass(b)/curmtotcom)*stock_dedz
                     close(1000+b)
                  end if
                  if (b .ge. 10 .and. b .lt. 100) then
                     write(num2,'(i2)') b
                     file_name1 = "SMD"//num2//".dat"
                     open (1000+b,file=file_name1,position='append',
     &               action="write",status='old')
                     write (1000+b,fmt='(i6,f15.5,f12.5,f12.5,
     &               f12.5,f12.5,f12.5,f12.5)')
     &               type(b), tsmd, x(b), y(b), z(b),
     &               (mass(b)/curmtotcom)*stock_dedx,
     &               (mass(b)/curmtotcom)*stock_dedy,
     &               (mass(b)/curmtotcom)*stock_dedz
                     close(1000+b)
                  end if
                  if (b .ge. 100 .and. b .lt. 1000) then
                     write(num3,'(i3)') b
                     file_name1 = "SMD"//num3//".dat"
                     open (1000+b,file=file_name1,position='append',
     &               action="write",status='old')
                     write (1000+b,fmt='(i6,f15.5,f12.5,f12.5,
     &               f12.5,f12.5,f12.5,f12.5)')
     &               type(b), tsmd, x(b), y(b), z(b),
     &               (mass(b)/curmtotcom)*stock_dedx,
     &               (mass(b)/curmtotcom)*stock_dedy,
     &               (mass(b)/curmtotcom)*stock_dedz
                     close(1000+b)
                  end if
                  if (b .ge. 1000 .and. b .lt. 10000) then
                     write(num4,'(i4)') b
                     file_name1 = "SMD"//num4//".dat"
                     open (1000+b,file=file_name1,position='append',
     &               action="write",status='old')
                     write (1000+b,fmt='(i6,f15.5,f12.5,f12.5,
     &               f12.5,f12.5,f12.5,f12.5)')
     &               type(b), tsmd, x(b), y(b), z(b),
     &               (mass(b)/curmtotcom)*stock_dedx,
     &               (mass(b)/curmtotcom)*stock_dedy,
     &               (mass(b)/curmtotcom)*stock_dedz
                     close(1000+b)
                  end if
                  if (b .ge. 10000 .and. b .lt. 100000) then
                     write(num5,'(i5)') b
                     file_name1 = "SMD"//num5//".dat"
                     open (1000+b,file=file_name1,position='append',
     &               action="write",status='old')
                     write (1000+b,fmt='(i6,f15.5,f12.5,f12.5,
     &               f12.5,f12.5,f12.5,f12.5)')
     &               type(b), tsmd, x(b), y(b), z(b),
     &               (mass(b)/curmtotcom)*stock_dedx,
     &               (mass(b)/curmtotcom)*stock_dedy,
     &               (mass(b)/curmtotcom)*stock_dedz
                     close(1000+b)
                  end if
                  if (b .ge. 100000 .and. b .lt. 1000000) then
                     write(num6,'(i6)') b
                     file_name1 = "SMD"//num6//".dat"
                     open (1000+b,file=file_name1,position='append',
     &               action="write",status='old')
                     write (1000+b,fmt='(i6,f15.5,f12.5,f12.5,
     &               f12.5,f12.5,f12.5,f12.5)')
     &               type(b), tsmd, x(b), y(b), z(b),
     &               (mass(b)/curmtotcom)*stock_dedx,
     &               (mass(b)/curmtotcom)*stock_dedy,
     &               (mass(b)/curmtotcom)*stock_dedz
                     close(1000+b)
                  end if
                  if (b .ge. 1000000 .and. b .lt. 10000000) then
                     write(num7,'(i7)') b
                     file_name1 = "SMD"//num7//".dat"
                     open (1000+b,file=file_name1,position='append',
     &               action="write",status='old')
                     write (1000+b,fmt='(i6,f15.5,f12.5,f12.5,
     &               f12.5,f12.5,f12.5,f12.5)')
     &               type(b), tsmd, x(b), y(b), z(b),
     &               (mass(b)/curmtotcom)*stock_dedx,
     &               (mass(b)/curmtotcom)*stock_dedy,
     &               (mass(b)/curmtotcom)*stock_dedz
                     close(1000+b)
                  end if
                  end if
                  end if
              end if
              end do
          end do
      end if
c
c##########################################
c#### End of the Constant Velocity SMD ####
c##########################################
c
c########################################
c########################################
c#### Case of the constant force SMD ####
c########################################
c########################################
c
c
      if (use_smd_forconst) then
c
c     Initial SMD force calculation according to the current COM
c
          com_dedx = xdir*SMDFor
          com_dedy = ydir*SMDFor
          com_dedz = zdir*SMDFor
c
c     Initial SMD energy calculation according to the current COM
c
          ensmd = xdir*SMDFor*cur_xcom + ydir*SMDFor*cur_ycom + 
     &    zdir*SMDFor*cur_zcom
c
c     SMD output data
c
          if (rank==0) then
              c = modulo(tpass,SMDoutputFreq)
              if (c == 0) then          ! LOOPA
              open (ismdout,file='SMD_output.dat',position='append',
     &        action="write",status='old')
              write (ismdout,51)  tsmd,cur_xcom,cur_ycom,cur_zcom
   51         format ('SMD',1x,f15.5,3f10.4)
              close (ismdout)
              end if                    ! END LOOPB
          end if
c
c     Redistribution on each SMD atoms forming the COM
c
          do a = 1, nsmdloc
              b = nsmdglob(a)
              ib = loc(b)
              dedx(b) = (mass(b)/curmtotcom)*com_dedx 
              dedy(b) = (mass(b)/curmtotcom)*com_dedy
              dedz(b) = (mass(b)/curmtotcom)*com_dedz
c
c     Incrementation of the derivative terms
c
              desmd(1,ib) = desmd(1,ib) - dedx(b)
              desmd(2,ib) = desmd(2,ib) - dedy(b)
              desmd(3,ib) = desmd(3,ib) - dedz(b)
c
c     Incrementation of the Virial
c
              vxx = x(b) * dedx(b)
              vyx = x(b) * dedx(b)
              vzx = x(b) * dedx(b)
              vyy = y(b) * dedy(b)
              vzy = y(b) * dedy(b)
              vzz = z(b) * dedz(b)
              vir(1,1) = vir(1,1) + vxx
              vir(2,1) = vir(2,1) + vyx
              vir(3,1) = vir(3,1) + vzx
              vir(1,2) = vir(1,2) + vyx
              vir(2,2) = vir(2,2) + vyy
              vir(3,2) = vir(3,2) + vzy
              vir(1,3) = vir(1,3) + vzx
              vir(2,3) = vir(2,3) + vzy
              vir(3,3) = vir(3,3) + vzz
          end do
c
c#######################################
c#### End of the Constant Force SMD ####
c#######################################
      end if
c
c##############################
c#### END OF THE PROCEDURE ####
c##############################
c
c
c     Incrementation of the tpass counter
c     => +1 each times one enters in !
c
      tpass = tpass + 1
c
      tsmd = tsmd + SMDdt
      return
      end



