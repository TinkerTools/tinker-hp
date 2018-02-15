c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################
c     ##                                                           ##
c     ##  files.i  --  name and number of current structure files  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     nprior     number of previously existing cycle files
c     ldir       length in characters of the directory name
c     leng       length in characters of the base filename
c     filename   base filename used by default for all files
c     outfile    output filename used for intermediate results
c
c
      integer nprior,ldir,leng
      character*120 filename,outfile
      common /files/ nprior,ldir,leng,filename,outfile
