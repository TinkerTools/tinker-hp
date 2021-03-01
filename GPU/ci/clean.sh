#!/bin/bash
cd 2decomp_fft && make clean
cd ../plumed2  && make clean
cd ../wrappers && make clean
cd ../source   && make clean
cd ../bin      && rm *
