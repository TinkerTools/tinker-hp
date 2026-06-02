INC_DIR=include_s
LIB_NAME=lib2decomp_fft_s.a

echo "-----------------------------------------------"
echo " Re-installing 2decomp_fft single precision in "
echo "  include directory : $INC_DIR  "
echo "  library's name    : $LIB_NAME "
echo "-----------------------------------------------"

cd ../source
make 2decomp_fft_rebuild_single

cd ../2decomp_fft/include
mkdir -p ../$INC_DIR
rm -f ../$INC_DIR/*
mv *.mod ../$INC_DIR/
cd ../lib
rm -f $LIB_NAME
mv lib2decomp_fft.a $LIB_NAME
cd ../../
