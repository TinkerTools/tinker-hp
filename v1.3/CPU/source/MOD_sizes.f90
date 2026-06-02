!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###################################################################
!     ##                                                               ##
!     ##  module  sizes  --  parameter values to set array dimensions  ##
!     ##                                                               ##
!     ###################################################################
!
!
!     "sizes.f" sets values for critical array dimensions used
!     throughout the software
!
!     parameter:      maximum allowed number of:
!
!     maxvalue        !<max allowed number of atoms directly bonded to an atom
!     maxgrp          !<max allowed number of user-defined groups of atoms
!     maxtyp          !<max allowed number of force field atom type definitions
!     maxclass        !<max allowed number of force field atom class definitions
!     maxprm          !<max allowed number of lines in the parameter file
!     maxkey          !<max allowed number of lines in the keyword file
!     maxvlst         !<max allowed number of neighbors in van der Waals pair list
!     maxelst         !<max allowed number of neighbors in electrostatics pair list
!     maxfft          !<max allowed number of grid points in each FFT dimension
!     maxring         !<max allowed number of 3-, 4-, or 5-membered rings
!     maxbio          !<max allowed number of biopolymer atom definitions
!     maxres          !<max allowed number of residues in the macromolecule
!     maxele          !<max allowed number of elements in periodic table
!     maxamino        !<max allowed number of amino acid residue types
!     maxnuc          !<max allowed number of nucleic acid residue types
!     maxtors         !<max allowed number of torsional angles in molecular system
!     maxbitor        !<max allowed number of bitorsions in molecular system
!
!
module sizes
   implicit none
   integer :: maxvalue !<max allowed number of atoms directly bonded to an atom
   integer :: maxgrp !<max allowed number of user-defined groups of atoms
   integer :: maxtyp !<max allowed number of force field atom type definitions
   integer :: maxclass !<max allowed number of force field atom class definitions
   integer :: maxprm !<max allowed number of lines in the parameter file
   integer :: maxkey !<max allowed number of lines in the keyword file
   integer :: maxvlst !<max allowed number of neighbors in van der Waals pair listt
   integer :: maxelst !<max allowed number of neighbors in electrostatics pair list
   integer :: maxfft !<max allowed number of grid points in each FFT dimension
   integer :: maxring !<max allowed number of 3-, 4-, or 5-membered rings
   integer :: maxbio !<max allowed number of biopolymer atom definitions
   integer :: maxele !<max allowed number of elements in periodic table
   integer :: maxamino !<max allowed number of amino acid residue types
   integer :: maxnuc !<max allowed number of nucleic acid residue types
   integer :: maxbnd !<max allowed number of bonds
   integer :: maxang !<max allowed number of angles
   integer :: maxtors !<max allowed number of torsions
   integer :: maxbitor !<max allowed number of bitorsions
   parameter (maxvalue=8)
   parameter (maxgrp=1000)
   parameter (maxtyp=5000)
   parameter (maxclass=2000)
   parameter (maxprm=25000)
   parameter (maxkey=5000)
   parameter (maxvlst=2500)
   parameter (maxelst=1200)
   parameter (maxfft=864)
   parameter (maxring=10000)
   parameter (maxbio=10000)
   parameter (maxele=112)
   parameter (maxamino=38)
   parameter (maxnuc=12)
   integer :: tinkerdebug !<integer defining level of output for debugging purposes
   save
end
