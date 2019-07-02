# LiveCOMS Paper
Raising the Performance of the Tinker-HP Molecular Modeling Package on Intel's HPC Architectures: a Living Review [Article v1.0]

Versions:

Current Github version: Preprint version 1.0 (07/02/2019)


This living paper reviews the present High Performance Computing (HPC)
capabilities of the Tinker-HP molecular modeling package. We focus here on the
reference, double precision, massively parallel molecular dynamics engine
present in Tinker-HP and dedicated to perform large scale simulations. We show
how it can be adapted to recent Intel®Central Processing Unit (CPU) petascale
architectures. First, we discuss the new set of Intel®Advanced Vector
Extensions 512 (Intel AVX-512) instructions present in recent Intel processors
(e.g., the Intel®Xeon®Scalable and Intel®Xeon Phi#2nd generation processors)
allowing for larger vectorization enhancements. These instructions constitute
the central source of potential computational gains when using the latest
processors, justifying important vectorization efforts for developers. We then
brie2y review the organization of the Tinker-HP code and identify the
computational hotspots which require Intel AVX-512 optimization and we propose
a general and optimal strategy to vectorize those particular parts of the code.
We intended to present our optimization strategy in a pedagogical way so it
could bene1t to other researchers and students interested in gaining
performances in their own software.  Finally we present the performance
enhancements obtained compared to the unoptimized code both sequentially and at
the scaling limit in parallel for classical non-polarizable (CHARMM) and
polarizable force 1elds (AMOEBA). This paper never ceases to be updated as we
accumulate new data on the associated Github repository between new versions of
this living paper.
