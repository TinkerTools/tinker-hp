#
#     Sorbonne University
#     Washington University in Saint Louis
#     University of Texas at Austin
#
 
fft_fftw3.o: decomp_2d.o
fft_fftw3_f03.o: decomp_2d.o
fft_generic.o: decomp_2d.o glassman.o
fft_mkl.o: decomp_2d.o mkl_dfti.o
glassman.o: decomp_2d.o 
glassman.mod: decomp_2d.o
io.o: decomp_2d.o
decomp_2d_io.mod:decomp_2d.o
