#
#  Sorbonne University
#  Washington University in Saint Louis
#  University of Texas at Austin
#

decomp_2d.o: transpose_x_to_y.f90 transpose_y_to_x.f90 transpose_y_to_z.f90 transpose_z_to_y.f90 factor.f90 halo.f90 alloc.f90
fft_acml.o:  decomp_2d.o fft_common.f90 acml_plan.f90
fft_essl.o:  decomp_2d.o glassman.o fft_common.f90 fft_common_3d.f90
fft_fftw3.o: decomp_2d.o fft_common.f90
fft_fftw3_f03.o:  decomp_2d.o
fft_fftpack5.f90: decomp_2d.o fft_common.f90
fft_ffte.o:  decomp_2d.o glassman.o factor.f90 fft_common.f90 fft_common_3d.f90
fft_mkl.o:   decomp_2d.o fft_common.f90
fft_cufft.o: decomp_2d.o fft_$(FFT).o
fft_generic.o: decomp_2d.o glassman.o fft_common.f90 fft_common_3d.f90
glassman.o:  decomp_2d.o
io.o: decomp_2d.o io_read_one.f90 io_read_var.f90 io_write_one.f90 io_write_var.f90 io_write_every.f90 io_write_plane.f90
