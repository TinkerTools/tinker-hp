c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###################################################################
c     ##                                                               ##
c     ##  module  sizes  --  parameter values to set array dimensions  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     "sizes.f" sets values for critical array dimensions used
c     throughout the software
c
c     parameter:      maximum allowed number of:
c
c     maxvalue          atoms directly bonded to an atom
c     maxgrp          user-defined groups of atoms
c     maxtyp          force field atom type definitions
c     maxclass        force field atom class definitions
c     maxprm          lines in the parameter file
c     maxkey          lines in the keyword file
c     maxvlst         neighbors in van der Waals pair list
c     maxelst         neighbors in electrostatics pair list
c     maxfft          grid points in each FFT dimension
c     maxring         3-, 4-, or 5-membered rings
c     maxbio          biopolymer atom definitions
c     maxres          residues in the macromolecule
c     maxele          elements in periodic table
c     maxamino        amino acid residue types
c     maxnuc          nucleic acid residue types
c     maxtors         torsional angles in molecular system
c     maxbitor        bitorsions in molecular system
c
c
      module sizes
      implicit none
      integer maxvalue,maxgrp
      integer maxtyp,maxclass
      integer maxprm,maxkey
      integer maxopt
      integer maxvlst,maxelst
      integer maxfft
      integer maxcell,maxref
      integer maxring,maxbio,maxres
      integer maxele,maxamino,maxnuc
      integer maxbnd,maxang,maxtors
      integer maxbitor,maxlp
      parameter (maxvalue=8)
      parameter (maxgrp=1000)
      parameter (maxref=10)
      parameter (maxtyp=5000)
      parameter (maxclass=1000)
      parameter (maxprm=25000)
      parameter (maxkey=5000)
      parameter (maxvlst=2500)
      parameter (maxelst=1200)
      parameter (maxfft=864)
      parameter (maxring=10000)
      parameter (maxbio=10000)
      parameter (maxres=10000)
      parameter (maxele=112)
      parameter (maxamino=38)
      parameter (maxnuc=12)
      save
      end
