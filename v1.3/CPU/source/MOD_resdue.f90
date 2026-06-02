!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #####################################################################
!     ##                                                                 ##
!     ##  module resdue  --  biopolymer residue names and biotype codes  ##
!     ##                                                                 ##
!     #####################################################################
!
!
!
module resdue
   use sizes
   implicit none
   integer :: ntyp(maxamino) !<biotypes for mid-chain peptide backbone N atoms
   integer :: catyp(maxamino) !<biotypes for mid-chain peptide backbone CA atoms
   integer :: ctyp(maxamino) !<biotypes for mid-chain peptide backbone C atoms
   integer :: hntyp(maxamino) !<biotypes for mid-chain peptide backbone HN atoms
   integer :: otyp(maxamino) !<biotypes for mid-chain peptide backbone O atoms
   integer :: hatyp(maxamino) !<biotypes for mid-chain peptide backbone HA atoms
   integer :: cbtyp(maxamino) !<biotypes for mid-chain peptide backbone CB atoms
   integer :: nntyp(maxamino) !<biotypes for N-terminal peptide backbone N atoms
   integer :: cantyp(maxamino) !<biotypes for N-terminal peptide backbone CA atoms
   integer :: cntyp(maxamino) !<biotypes for N-terminal peptide backbone C atoms
   integer :: hnntyp(maxamino) !<biotypes for N-terminal peptide backbone HN atoms
   integer :: ontyp(maxamino) !<biotypes for N-terminal peptide backbone O atoms
   integer :: hantyp(maxamino) !<biotypes for N-terminal peptide backbone HA atoms
   integer :: nctyp(maxamino) !<biotypes for C-terminal peptide backbone N atoms
   integer :: cactyp(maxamino) !<biotypes for C-terminal peptide backbone CA atoms
   integer :: cctyp(maxamino) !<biotypes for C-terminal peptide backbone C atoms
   integer :: hnctyp(maxamino) !<biotypes for C-terminal peptide backbone HN atoms
   integer :: octyp(maxamino) !<biotypes for C-terminal peptide backbone O atoms
   integer :: hactyp(maxamino) !<biotypes for C-terminal peptide backbone HA atoms
   integer :: o5typ(maxnuc) !<biotypes for nucleotide backbone and sugar O5' atoms
   integer :: c5typ(maxnuc) !<biotypes for nucleotide backbone and sugar C5' atoms
   integer :: h51typ(maxnuc) !<biotypes for nucleotide backbone and sugar H5' atoms
   integer :: h52typ(maxnuc) !<biotypes for nucleotide backbone and sugar H5'' atoms
   integer :: c4typ(maxnuc) !<biotypes for nucleotide backbone and sugar C4' atoms
   integer :: h4typ(maxnuc) !<biotypes for nucleotide backbone and sugar H4' atoms
   integer :: o4typ(maxnuc) !<biotypes for nucleotide backbone and sugar O4' atoms
   integer :: c1typ(maxnuc) !<biotypes for nucleotide backbone and sugar C1' atoms
   integer :: h1typ(maxnuc) !<biotypes for nucleotide backbone and sugar H1' atoms
   integer :: c3typ(maxnuc) !<biotypes for nucleotide backbone and sugar C3' atoms
   integer :: h3typ(maxnuc) !<biotypes for nucleotide backbone and sugar H3' atoms
   integer :: c2typ(maxnuc) !<biotypes for nucleotide backbone and sugar C2' atoms
   integer :: h21typ(maxnuc) !<biotypes for nucleotide backbone and sugar H2' atoms
   integer :: o2typ(maxnuc) !<biotypes for nucleotide backbone and sugar O2' atoms
   integer :: h22typ(maxnuc) !<biotypes for nucleotide backbone and sugar H2'' atoms
   integer :: o3typ(maxnuc) !<biotypes for nucleotide backbone and sugar O3' atoms
   integer :: ptyp(maxnuc) !<biotypes for nucleotide backbone and sugar P atoms
   integer :: optyp(maxnuc) !<biotypes for nucleotide backbone and sugar OP atoms
   integer :: h5ttyp(maxnuc) !<biotypes for nucleotide backbone and sugar H5T atoms
   integer :: h3ttyp(maxnuc) !<biotypes for nucleotide backbone and sugar H3T atoms
   character*1 :: amino1(maxamino) !<three-letter abbreviations for amino acids types
   character*1 :: nuclz1(maxnuc) !<three-letter abbreviations for nucleic acids types
   character*3 :: amino(maxamino) !<one-letter abbreviations for amino acids types
   character*3 :: nuclz(maxnuc) !<one-letter abbreviations for nucleic acids types
   save
end
