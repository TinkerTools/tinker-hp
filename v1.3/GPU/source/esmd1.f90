!
!          ###    ##    ##    ###
!         #       # #  # #    #  ##
!        #        #  ##  #    #    #
!         ###     #      #    #    #
!           #     #      #    #    #
!          #      #      #    #  ##
!       ###       #      #    ###
!
!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###############################################################
!     ##                                                           ##
!     ##  subroutine esmd1  --  make smd calculations              ##
!     ##                                                           ##
!     ###############################################################
!
!
!     "esmd1" assigns the biased forced on each concerned atoms and
!     stores the smd and derivatives energy components
!
!     literature reference:
!     Massively Parallel Implementation of Steered Molecular Dynamics in Tinker-HP:
!     Comparisons of Polarizable and Non-Polarizable Simulations of Realistic Systems
!     F Célerse, L Lagardère, E Derat, JP Piquemal
!     Journal of Chemical Theory and Computation 15 (6), 3694-3709
!
!
#include "tinker_macro.h"
module esmd1_inl
contains
#include "convert.inc.f90"
end module

subroutine esmd1
   use atmtyp
   use atoms
   use deriv
   use domdec
   use energi
   use esmd1_inl
   use inform    ,only: deb_Path
   use mpi
   use msmd
   use potent
   use sizes     ,only: tinkerdebug
   use tinheader ,only:ti_p,re_p
   use virial
   implicit none
   character(len=100)::file_name1
   character(len=1)::num1
   character(len=2)::num2
   character(len=3)::num3
   character(len=4)::num4
   character(len=5)::num5
   character(len=6)::num6
   character(len=7)::num7

   real(t_p) diffx, diffy, diffz
   real(t_p) diff, diff2
   real(t_p) vxx,vyy,vzz
   real(t_p) vyx,vzx,vzy
   real(t_p) buffers(7),buffer(7)
   real(r_p) ded_x,ded_y,ded_z,x_b,y_b,z_b
   integer ierr,reqrec,reqsend,STATUS(MPI_STATUS_SIZE),tagmpi
   integer a,b,c,j
   integer ib
   logical IsParallel

   if (deb_Path) print*,"esmd1"
!
!     SMD storage: initialization
!
   stock_dedx  = 0.0_re_p
   stock_dedy  = 0.0_re_p
   stock_dedz  = 0.0_re_p
   stock_ensmd = 0.0_re_p
   IsParallel  = merge(.true.,.false.,(ndir.ge.1))
!
!###########################################
!###########################################
!#### Case of the constant velocity SMD ####
!###########################################
!###########################################
!
!####################################
!#### Case where only k1 is used ####
!####################################
!
   if (use_smd_velconst .and. .not. use_smdk2) then  ! LOOP1 CVSMDK1
!
!     Initial SMD force calculation according to the current COM
!
      diffx = (cur_xcom - xcom)*xdir
      diffy = (cur_ycom - ycom)*ydir
      diffz = (cur_zcom - zcom)*zdir
      diff = diffx + diffy + diffz
!
      com_dedx = SMDk*xdir*(SMDVel*tsmd - diff)
      com_dedy = SMDk*ydir*(SMDVel*tsmd - diff)
      com_dedz = SMDk*zdir*(SMDVel*tsmd - diff)
!
!     Initial SMD energy calculation according to the current COM
!
      diff = sqrt(xcom**2 + ycom**2 + zcom**2)
      ensmd = 0.5*SMDk*(SMDVel*tsmd - diff)**2
!$acc update device(ensmd) async
!
!     Storage for printing
!     => Mandatory if integrator != VERLET (so multitimestep)
!
      stock_dedx = stock_dedx + com_dedx
      stock_dedy = stock_dedy + com_dedy
      stock_dedz = stock_dedz + com_dedz
      stock_ensmd = stock_ensmd + ensmd
!
!     the proc of rank 'smdprocprint' sends the info to the rank 0
!
      if (IsParallel) then
         if (rank.eq.smdprocprint) then
            buffers(1) = cur_xcom
            buffers(2) = cur_ycom
            buffers(3) = cur_zcom
            buffers(4) = stock_dedx
            buffers(5) = stock_dedy
            buffers(6) = stock_dedz
            buffers(7) = curmtotcom
            tagmpi = rank+1
            call MPI_ISEND(buffers,7,MPI_TPREC,0,tagmpi,COMM_TINKER,&
               &reqsend,ierr)
            call MPI_WAIT(reqsend,status,ierr)
         else if (rank.eq.0) then
            tagmpi = smdprocprint+1
            call MPI_IRECV(buffer,7,MPI_TPREC,smdprocprint,tagmpi,&
               &COMM_TINKER,reqrec,ierr)
            call MPI_WAIT(reqrec,status,ierr)
            cur_xcom = buffer(1)
            cur_ycom = buffer(2)
            cur_zcom = buffer(3)
            stock_dedx = buffer(4)
            stock_dedy = buffer(5)
            stock_dedz = buffer(6)
            curmtotcom = buffer(7)
         end if
      end if
!
!     SMD output data
!
      if (rank==0) then
         c = modulo(tpass,SMDoutputFreq)
         if (c == 0) then          ! LOOPA
            open (ismdout,file='SMD_output.dat',position='append',&
               &action="write",status='old')
            write (ismdout,45)  tsmd,cur_xcom,cur_ycom,cur_zcom,&
               &stock_dedx,stock_dedy,stock_dedz
45          format ('SMD',1x,f15.5,3f10.4,3f10.4)
            close(ismdout)
         end if                    ! END LOOPB
      end if
!
!     Redistribution on each SMD atoms forming the COM
!
!$acc parallel loop async default(present) &
!$acc         present(g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
      do a = 1, nsmdloc     ! LOOP2 CVSMDK2
         b  = nsmdglob(a)
         ib = loc(b)
!
         ded_x = (mass(b)/curmtotcom)*com_dedx
         ded_y = (mass(b)/curmtotcom)*com_dedy
         ded_z = (mass(b)/curmtotcom)*com_dedz
!
!     Incrementation of the derivative terms
!
         desmd(1,ib) = desmd(1,ib) - rp2mdr(ded_x)
         desmd(2,ib) = desmd(2,ib) - rp2mdr(ded_y)
         desmd(3,ib) = desmd(3,ib) - rp2mdr(ded_z)
!
!     Increment the internal virial tensor components
!
         g_vxx = g_vxx + (cur_xcom - xcom)*ded_x
         g_vxy = g_vxy + (cur_ycom - ycom)*ded_x
         g_vxz = g_vxz + (cur_zcom - zcom)*ded_x
         g_vyy = g_vyy + (cur_ycom - ycom)*ded_y
         g_vyz = g_vyz + (cur_zcom - zcom)*ded_y
         g_vzz = g_vzz + (cur_zcom - zcom)*ded_z
      end do
!
   end if    ! END LOOP1 CVSMDK1
!
!#######################################
!#### Case where k2 is used with k1 ####
!#######################################
!
   if (use_smd_velconst .and. use_smdk2) then
!
!     Initial SMD force calculation according to the current COM
!
      diffx = (cur_xcom - xcom)*xdir
      diffy = (cur_ycom - ycom)*ydir
      diffz = (cur_zcom - zcom)*zdir
      diff = diffx + diffy + diffz
!
      com_dedx = SMDk*xdir*(SMDVel*tsmd - diff)
      com_dedy = SMDk*ydir*(SMDVel*tsmd - diff)
      com_dedz = SMDk*zdir*(SMDVel*tsmd - diff)
!
!     Initial SMD energy calculation according to the current COM
!
      diff = sqrt(xcom**2 + ycom**2 + zcom**2)
      ensmd = 0.5*SMDk*(SMDVel*tsmd - diff)**2
!
!     Calculation of the derivative components of k2
!
!     #Part 1: Contribution without the motion#
      com_dedx = com_dedx - SMDk2*(cur_xcom - xcom)
      com_dedy = com_dedy - SMDk2*(cur_ycom - ycom)
      com_dedz = com_dedz - SMDk2*(cur_zcom - zcom)
!
!     #Part 2: Contribution with the motion#
!
      com_dedx = com_dedx + (SMDk2*diff*xdir)
      com_dedy = com_dedy + (SMDk2*diff*ydir)
      com_dedz = com_dedz + (SMDk2*diff*zdir)
!
!     SMD energy calculation according to the current COM and k2
!
      diff2 = sqrt((cur_xcom - xcom)**2 + (cur_ycom - ycom)**2 +&
         &(cur_zcom - zcom)**2)
      ensmd = ensmd + 0.5*SMDk2*(diff2 - (diff)**2)
!$acc update device(ensmd) async
!
!     Storage for printing
!     => Mandatory if integrator != VERLET (so multitimestep)
!
      stock_dedx = stock_dedx + com_dedx
      stock_dedy = stock_dedy + com_dedy
      stock_dedz = stock_dedz + com_dedz
      stock_ensmd = stock_ensmd + ensmd
!
!     the proc of rank 'smdprocprint' sends the info to the rank 0
!
      if (IsParallel) then
         if (rank.eq.smdprocprint) then
            buffers(1) = cur_xcom
            buffers(2) = cur_ycom
            buffers(3) = cur_zcom
            buffers(4) = stock_dedx
            buffers(5) = stock_dedy
            buffers(6) = stock_dedz
            buffers(7) = curmtotcom
            tagmpi = rank+1
            call MPI_ISEND(buffers,7,MPI_TPREC,0,tagmpi,COMM_TINKER,&
               &reqsend,ierr)
            call MPI_WAIT(reqsend,status,ierr)
         else if (rank.eq.0) then
            tagmpi = smdprocprint+1
            call MPI_IRECV(buffer,7,MPI_TPREC,smdprocprint,tagmpi,&
               &COMM_TINKER,reqrec,ierr)
            call MPI_WAIT(reqrec,status,ierr)
            cur_xcom = buffer(1)
            cur_ycom = buffer(2)
            cur_zcom = buffer(3)
            stock_dedx = buffer(4)
            stock_dedy = buffer(5)
            stock_dedz = buffer(6)
            curmtotcom = buffer(7)
         end if
      end if
!
!     SMD output data
!
      if (rank==0) then
         c = modulo(tpass,SMDoutputFreq)
         if (c == 0) then          ! LOOPA
            open (ismdout,file='SMD_output.dat',position='append',&
               &action="write",status='old')
            write (ismdout,48)  tsmd,cur_xcom,cur_ycom,cur_zcom,&
               &stock_dedx,stock_dedy,stock_dedz
48          format ('SMD',1x,f15.5,3f10.4,3f10.4)
            close(ismdout)
         end if                    ! END LOOPB
      end if
!
!     Redistribution on each SMD atoms forming the COM
!
!$acc parallel loop async default(present) &
!$acc         present(g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
      do a = 1, nsmdloc     ! LOOP2 CVSMDK2
         b  = nsmdglob(a)
         ib = loc(b)
!
         ded_x = (mass(b)/curmtotcom)*com_dedx
         ded_y = (mass(b)/curmtotcom)*com_dedy
         ded_z = (mass(b)/curmtotcom)*com_dedz
!
!     Incrementation of the derivative terms
!
         desmd(1,ib) = desmd(1,ib) - rp2mdr(ded_x)
         desmd(2,ib) = desmd(2,ib) - rp2mdr(ded_y)
         desmd(3,ib) = desmd(3,ib) - rp2mdr(ded_z)
!
!     Increment the internal virial tensor components
!
         g_vxx = g_vxx + (cur_xcom - xcom)*ded_x
         g_vxy = g_vxy + (cur_ycom - ycom)*ded_x
         g_vxz = g_vxz + (cur_zcom - zcom)*ded_x
         g_vyy = g_vyy + (cur_ycom - ycom)*ded_y
         g_vyz = g_vyz + (cur_zcom - zcom)*ded_y
         g_vzz = g_vzz + (cur_zcom - zcom)*ded_z
      end do
!
   end if
!
!##########################################
!#### End of the Constant Velocity SMD ####
!##########################################
!
!########################################
!########################################
!#### Case of the constant force SMD ####
!########################################
!########################################
!
!
   if (use_smd_forconst) then
!
!     Initial SMD force calculation according to the current COM
!
      com_dedx = xdir*SMDFor
      com_dedy = ydir*SMDFor
      com_dedz = zdir*SMDFor
!
!     Initial SMD energy calculation according to the current COM
!
      ensmd = xdir*SMDFor*cur_xcom + ydir*SMDFor*cur_ycom +&
         &zdir*SMDFor*cur_zcom
!$acc update device(ensmd) async
!
!     SMD output data
!
      if (rank==0) then
         c = modulo(tpass,SMDoutputFreq)
         if (c == 0) then          ! LOOPA
            open (ismdout,file='SMD_output.dat',position='append',&
               &action="write",status='old')
            write (ismdout,51)  tsmd,cur_xcom,cur_ycom,cur_zcom
51          format ('SMD',1x,f15.5,3f10.4)
            close (ismdout)
         end if                    ! END LOOPB
      end if
!
!     Redistribution on each SMD atoms forming the COM
!
!$acc parallel loop async default(present) &
!$acc         present(g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
      do a = 1, nsmdloc
         b   = nsmdglob(a)
         ib   = loc(b)
         x_b   = x(b)
         y_b   = y(b)
         z_b   = z(b)
         ded_x = (mass(b)/curmtotcom)*com_dedx
         ded_y = (mass(b)/curmtotcom)*com_dedy
         ded_z = (mass(b)/curmtotcom)*com_dedz
!
!     Incrementation of the derivative terms
!
         desmd(1,ib) = desmd(1,ib) - rp2mdr(ded_x)
         desmd(2,ib) = desmd(2,ib) - rp2mdr(ded_y)
         desmd(3,ib) = desmd(3,ib) - rp2mdr(ded_z)
!
!     Incrementation of the Virial
!
         g_vxx = g_vxx + x_b * ded_x
         g_vxy = g_vxy + x_b * ded_x
         g_vxz = g_vxz + x_b * ded_x
         g_vyy = g_vyy + y_b * ded_y
         g_vyz = g_vyz + y_b * ded_y
         g_vzz = g_vzz + z_b * ded_z
      end do
!
!#######################################
!#### End of the Constant Force SMD ####
!#######################################
   end if
!
!##############################
!#### END OF THE PROCEDURE ####
!##############################
!
!
!     Incrementation of the tpass counter
!     => +1 each times one enters in !
!
   tpass = tpass + 1
!
   tsmd = tsmd + SMDdt
end
