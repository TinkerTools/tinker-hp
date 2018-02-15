c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine promo  --  version info and copywrite notice  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "promo" writes a short message containing information
c     about the Tinker-HP version number and the copyright notice
c
c
      subroutine promo
      implicit none
      include 'iounit.i'
      include 'openmp.i'
c
c
c     print out the informational header message
c
      write (iout,10)
   10 format (/,5x,70('#'),
     &        /,3x,74('#'),
     &        /,2x,'###',70x,'###',
     &        /,1x,'###',12x,'Tinker-HP  ---  Software Tools for',
     &           ' Molecular Design',9x,'###',
     &        /,1x,'##',74x,'##',
     &        /,1x,'##',24x,'Version 1.0  October 2017',25x,'##',
     &        /,1x,'##',74x,'##',
     &        /,1x,'##',5x,'Copyright (c) Washington University',
     &          ' in Saint Louis (WU)',14x,'##',
     &        /,1x,'##',19x,'The University of Texas at Austin',
     &        22x,'##',
     &        /,1x,'##',19x,'Sorbonne Universites, UPMC (Sorbonne)',
     &        18x,'##',
     &        /,1x,'##',28x,'  1990-2017',35x,'##',
     &        /,1x,'###',23x,'All Rights Reserved',30x,'###',
     &        /,2x,'###',70x,'###',
     &        /,3x,74('#'),
     &        /,5x,70('#'),/,/,'License Number : LICENSENUMBER'
     &        /)

      return
      end
