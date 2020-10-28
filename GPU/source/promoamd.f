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
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine promoamd  --  aMD input file design           ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "promoamd" writes a short message containing information
c     about the nature of the calculation performed and some designs
c
c
      subroutine promoamd(i)
      use iounit
      use mamd
      implicit none
      integer i
c
c     print out the informational header message
c
      if (i==1) write (iamdout,10)

   10 format (/,6x,'Sorbonne University',
     &        /,6x,'Washington University in Saint Louis',
     &        /,6x,'University of Texas at Austin',
     &        /,/,/,
     &        '         ##          ##    ##    ###',/,
     &        '        #  #         # #  # #    #  ##',/,
     &        '       #    #        #  ##  #    #    #',/,
     &        '      ########       #      #    #    #',/,
     &        '     #        #      #      #    #    #',/,
     &        '    #          #     #      #    #  ##',/,
     &        '   #            #    #      #    ###',/,
     &        /,/,6x,38('#'),
     &        /,6x,1('#'),36x,1('#'),
     &        /,6x,1('#'),1x,'TINKER-HP aMD MODULE : VERSION 1.3',
     &        1x,1('#'),
     &        /,6x,1('#'),36x,1('#'),
     &        /,6x,38('#'),
     &        /,/)
c
c     print out no warning 
c
      if (i==2) write (iamdout,20)

   20 format (/,5x,'^',2x,'^',/,
     &        6x,'--',/,
     &        5x,'(^^)\________',/
     &        5x,'(__)\',8x,')\/\...aMD',/
     &        9x,'||-----w |',/
     &        9x,'||      ||',
     &        /)
c
c     print out dead message
c
      if (i==3) write (iamdout,30)

   30 format (/,
     &       '                            /-----------\',/
     &       '                           /#############\',/
     &       '                          /###############\',/
     &       '                         /#################\',/
     &       '                         |##/~~~\###/~~~\##|',/
     &       '                         |##| o |###| o |##|',/
     &       '                         |##\---/XXX\___|##|',/
     &       '                         |########X########|',/
     &       '                         \#######XXX#######/',/
     &       '                          \#####XXXXX#####/',/
     &       '                           |#|%%%%%%%%|##|',/
     &       '                           |#|OHOHOHOH|##|',/
     &       '            XX             |#|OHOHOHOH|##|            XX',
     &       /,
     &       '          XX..X            \#############/',11x,'XX..X',/,
     &       '        XX.....X            --$$$$$$$$$--',11x,'XX.....X',
     &       /,
     &       '   XXXXX.....XX               \-------/',13x,
     &                                                  'XX.....XXXXX',/
     &       '  X |.......XX%,.@',32x,'@#%XX.......| X',/
     &       ' X |......X  @#%,.@',28x,'@#%,.@X......| X',/
     &       '  X  \..X      @#%,.@',26x,'@#%,.@     X../  X',/
     &       '   X# \.X       @#%,.@',24x,'@#%,.@      X./ #X',/
     &       '    ##  X        @#%,.@',22x,'@#%,.@       X  ##',/
     &       '    X# #X         @#%,.@',20x,'@#%,.@        X# #X',/
     &       '     X###X         @#%,.@',18x,'@#%,.@        X###X',/
     &       '      X ###         @#%,.@',16x,'@#%,.@        ### X',/
     &       '        X;X          @#%,.@',14x,'@#%,.@         X;X',/
     &       '        X             @#%,.@',12x,'@#%,.@            X',/
     &       '                       @#%,.@          @#%,.@',/
     &       '                        @#%,.@        @#%,.@',/
     &       '                         @#%,.@      @#%,.@',/
     &       '                          @#%,.@    @#%,.@',/
     &       '                           @#%,.@  @#%,.@',/
     &       '                            @#%,.@@#%,.@',/
     &    '            #-------------------------------------------#',/
     &    '            #                                           #',/
     &    '            #              The aMD caclulation          #',/
     &    '            #                     is                    #',/
     &    '            #                    DEAD                   #',/
     &    '            #___________________________________________#',/
     &       )
c     
c     print out final message
c
      if (i==4) then
          open (iamdout,file='aMD_output.dat',position='append',
     &    action="write",status='old')
          write (iamdout,40)
          close (iamdout)
      end if

   40 format (/,/,'aMDinfo: Calculation is ended normally.',/,
     &       'aMDinfo: No bugs were detected during the procedure.',/,
     &       'aMDinfo: End of the procedure !',/,/,
     &       '                ###',/
     &       '              #######',/
     &       '             #########',/
     &       '            ###########',/
     &       '   ^  ^',/
     &       '    --',/
     &       '   (^^)\________',/
     &       '   (__)\        )\/\...aMD',/
     &       '       ||-----w |',/
     &       '       ||      ||',/
     &       '#########################################',/
     &       '# ______      ____      _      _____    #',/
     &       '# |  ___|    |  _ \    | |    |  _  \   #',/
     &       '# | |        | | \ \   | |    | | \  \  #',/
     &       '# | |___     | |  \ \  | |    | |  \  \ #',/
     &       '# | ____|    | |   \ \ | |    | |   | | #',/
     &       '# | |        | |    \ \| |    | |   / / #',/
     &       '# | |___     | |     \   |    | |_ / /  #',/
     &       '# |_____|    |_|      \__|    |_____/   #',/
     &       '#                                       #',/
     &       '#########################################',/
     &       )
c
c      flush (ismdout)
      return
      end
