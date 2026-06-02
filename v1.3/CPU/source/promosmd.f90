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
!     ##  subroutine promosmd  --  SMD input file design           ##
!     ##                                                           ##
!     ###############################################################
!
!
!     "promosmd" writes a short message containing information
!     about the nature of the calculation performed and some designs
!
!
!> @brief 
!> "promosmd" writes a short message containing information
!> about the nature of the calculation performed and some designs
!> @param no params
subroutine promosmd(i)
   use iounit
   use msmd
   implicit none
   integer i
!
!     print out the informational header message
!
   if (i==1) write (ismdout,10)

10 format (/,6x,'Sorbonne University',&
   &/,6x,'Washington University in Saint Louis',&
   &/,6x,'University of Texas at Austin',&
   &/,/,/,&
   &11x,3('#'),4x,2('#'),4x,2('#'),4x,3('#'),&
   &/,10x,1('#'),7x,1('#'),1x,1('#'),2x,1('#'),1x,&
   &1('#'),4x,1('#'),2x,2('#'),&
   &/,9x,1('#'),8x,1('#'),2x,2('#'),2x,1('#'),4x,&
   &1('#'),4x,1('#'),&
   &/,10x,3('#'),5x,1('#'),6x,1('#'),4x,1('#'),4x,1('#'),&
   &/,12x,1('#'),5x,1('#'),6x,1('#'),4x,1('#'),4x,1('#'),&
   &/,11x,1('#'),6x,1('#'),6x,1('#'),4x,1('#'),2x,2('#'),&
   &/,8x,3('#'),7x,1('#'),6x,1('#'),4x,3('#'),/&
   &/,/,6x,38('#'),&
   &/,6x,1('#'),36x,1('#'),&
   &/,6x,1('#'),1x,'TINKER-HP SMD MODULE : VERSION 1.2',&
   &1x,1('#'),&
   &/,6x,1('#'),36x,1('#'),&
   &/,6x,38('#'),&
   &/,/)
!
!     print out no warning
!
   if (i==2) write (ismdout,20)

20 format (/,5x,'^',2x,'^',/,&
   &6x,'--',/,&
   &5x,'(^^)\________',/&
   &5x,'(__)\',8x,')\/\...SMD',/&
   &9x,'||-----w |',/&
   &9x,'||      ||',&
   &/)
!
!     print out dead message
!
   if (i==3) write (ismdout,30)

30 format (/,&
   &'                            /-----------\',/&
   &'                           /#############\',/&
   &'                          /###############\',/&
   &'                         /#################\',/&
   &'                         |##/~~~\###/~~~\##|',/&
   &'                         |##| o |###| o |##|',/&
   &'                         |##\---/XXX\___|##|',/&
   &'                         |########X########|',/&
   &'                         \#######XXX#######/',/&
   &'                          \#####XXXXX#####/',/&
   &'                           |#|%%%%%%%%|##|',/&
   &'                           |#|OHOHOHOH|##|',/&
   &'            XX             |#|OHOHOHOH|##|            XX',&
   &/,&
   &'          XX..X            \#############/',11x,'XX..X',/,&
   &'        XX.....X            --$$$$$$$$$--',11x,'XX.....X',&
   &/,&
   &'   XXXXX.....XX               \-------/',13x,&
   &'XX.....XXXXX',/&
   &'  X |.......XX%,.@',32x,'@#%XX.......| X',/&
   &' X |......X  @#%,.@',30x,'@#%,.@X......| X',/&
   &'  X  \..X      @#%,.@',26x,'@#%,.@     X../  X',/&
   &'   X# \.X       @#%,.@',24x,'@#%,.@      X./ #X',/&
   &'    ##  X        @#%,.@',22x,'@#%,.@       X  ##',/&
   &'    X# #X         @#%,.@',20x,'@#%,.@        X# #X',/&
   &'     X###X         @#%,.@',18x,'@#%,.@        X###X',/&
   &'      X ###         @#%,.@',16x,'@#%,.@        ### X',/&
   &'        X;X          @#%,.@',14x,'@#%,.@         X;X',/&
   &'        X             @#%,.@',12x,'@#%,.@            X',/&
   &'                       @#%,.@          @#%,.@',/&
   &'                        @#%,.@        @#%,.@',/&
   &'                         @#%,.@      @#%,.@',/&
   &'                          @#%,.@    @#%,.@',/&
   &'                           @#%,.@  @#%,.@',/&
   &'                            @#%,.@@#%,.@',/&
   &'            #-------------------------------------------#',/&
   &'            #                                           #',/&
   &'            #              The SMD caclulation          #',/&
   &'            #                     is                    #',/&
   &'            #                    DEAD                   #',/&
   &'            #___________________________________________#',/&
   &)
!
!     print out final message
!
   if (i==4) then
      open (ismdout,file='SMD_output.dat',position='append',&
      &action="write",status='old')
      write (ismdout,40)
      close (ismdout)
   end if

40 format (/,/,'SMDinfo: Calculation is ended normally.',/,&
   &'SMDinfo: No bugs were detected during the procedure.',/,&
   &'SMDinfo: End of the procedure !',/,/,&
   &'                ###',/&
   &'              #######',/&
   &'             #########',/&
   &'            ###########',/&
   &'   ^  ^',/&
   &'    --',/&
   &'   (^^)\________',/&
   &'   (__)\        )\/\...SMD',/&
   &'       ||-----w |',/&
   &'       ||      ||',/&
   &'#########################################',/&
   &'# ______      ____      _      _____    #',/&
   &'# |  ___|    |  _ \    | |    |  _  \   #',/&
   &'# | |        | | \ \   | |    | | \  \  #',/&
   &'# | |___     | |  \ \  | |    | |  \  \ #',/&
   &'# | ____|    | |   \ \ | |    | |   | | #',/&
   &'# | |        | |    \ \| |    | |   / / #',/&
   &'# | |___     | |     \   |    | |_ / /  #',/&
   &'# |_____|    |_|      \__|    |_____/   #',/&
   &'#                                       #',/&
   &'#########################################',/&
   &)
!
   flush (ismdout)
   return
end
