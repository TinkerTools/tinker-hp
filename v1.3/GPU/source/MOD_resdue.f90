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
!     ntyp     biotypes for mid-chain peptide backbone N atoms
!     catyp    biotypes for mid-chain peptide backbone CA atoms
!     ctyp     biotypes for mid-chain peptide backbone C atoms
!     hntyp    biotypes for mid-chain peptide backbone HN atoms
!     otyp     biotypes for mid-chain peptide backbone O atoms
!     hatyp    biotypes for mid-chain peptide backbone HA atoms
!     cbtyp    biotypes for mid-chain peptide backbone CB atoms
!     nntyp    biotypes for N-terminal peptide backbone N atoms
!     cantyp   biotypes for N-terminal peptide backbone CA atoms
!     cntyp    biotypes for N-terminal peptide backbone C atoms
!     hnntyp   biotypes for N-terminal peptide backbone HN atoms
!     ontyp    biotypes for N-terminal peptide backbone O atoms
!     hantyp   biotypes for N-terminal peptide backbone HA atoms
!     nctyp    biotypes for C-terminal peptide backbone N atoms
!     cactyp   biotypes for C-terminal peptide backbone CA atoms
!     cctyp    biotypes for C-terminal peptide backbone C atoms
!     hnctyp   biotypes for C-terminal peptide backbone HN atoms
!     octyp    biotypes for C-terminal peptide backbone O atoms
!     hactyp   biotypes for C-terminal peptide backbone HA atoms
!     o5typ    biotypes for nucleotide backbone and sugar O5' atoms
!     c5typ    biotypes for nucleotide backbone and sugar C5' atoms
!     h51typ   biotypes for nucleotide backbone and sugar H5' atoms
!     h52typ   biotypes for nucleotide backbone and sugar H5'' atoms
!     c4typ    biotypes for nucleotide backbone and sugar C4' atoms
!     h4typ    biotypes for nucleotide backbone and sugar H4' atoms
!     o4typ    biotypes for nucleotide backbone and sugar O4' atoms
!     c1typ    biotypes for nucleotide backbone and sugar C1' atoms
!     h1typ    biotypes for nucleotide backbone and sugar H1' atoms
!     c3typ    biotypes for nucleotide backbone and sugar C3' atoms
!     h3typ    biotypes for nucleotide backbone and sugar H3' atoms
!     c2typ    biotypes for nucleotide backbone and sugar C2' atoms
!     h21typ   biotypes for nucleotide backbone and sugar H2' atoms
!     o2typ    biotypes for nucleotide backbone and sugar O2' atoms
!     h22typ   biotypes for nucleotide backbone and sugar H2'' atoms
!     o3typ    biotypes for nucleotide backbone and sugar O3' atoms
!     ptyp     biotypes for nucleotide backbone and sugar P atoms
!     optyp    biotypes for nucleotide backbone and sugar OP atoms
!     h5ttyp   biotypes for nucleotide backbone and sugar H5T atoms
!     h3ttyp   biotypes for nucleotide backbone and sugar H3T atoms
!     amino    three-letter abbreviations for amino acids types
!     nuclz    three-letter abbreviations for nucleic acids types
!     amino1   one-letter abbreviations for amino acids types
!     nuclz1   one-letter abbreviations for nucleic acids types
!
!
module resdue
   use sizes
   implicit none
   integer ntyp(maxamino),catyp(maxamino)
   integer ctyp(maxamino),hntyp(maxamino)
   integer otyp(maxamino),hatyp(maxamino),cbtyp(maxamino)
   integer nntyp(maxamino),cantyp(maxamino),cntyp(maxamino)
   integer hnntyp(maxamino),ontyp(maxamino),hantyp(maxamino)
   integer nctyp(maxamino),cactyp(maxamino),cctyp(maxamino)
   integer hnctyp(maxamino),octyp(maxamino),hactyp(maxamino)
   integer o5typ(maxnuc),c5typ(maxnuc),h51typ(maxnuc)
   integer h52typ(maxnuc),c4typ(maxnuc),h4typ(maxnuc)
   integer o4typ(maxnuc),c1typ(maxnuc),h1typ(maxnuc)
   integer c3typ(maxnuc),h3typ(maxnuc),c2typ(maxnuc)
   integer h21typ(maxnuc),o2typ(maxnuc),h22typ(maxnuc)
   integer o3typ(maxnuc),ptyp(maxnuc),optyp(maxnuc)
   integer h5ttyp(maxnuc),h3ttyp(maxnuc)
   character*1 amino1(maxamino),nuclz1(maxnuc)
   character*3 amino(maxamino),nuclz(maxnuc)
   save
end
