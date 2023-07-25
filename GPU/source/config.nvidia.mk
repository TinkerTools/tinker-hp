# -- Include MKL ?
add_mkl__        := no
ifeq ($(REL_BUILD),0)
   add_mkl__     := yes
endif
ifeq ($(arch),$(filter $(arch),host cpu))
   add_mkl__     := yes
endif

# -- arithmetic precision
ifeq ($(prec),$(filter $(prec),single s))
   main_prec        := s
   secn_prec        := s
else ifeq ($(prec),$(filter $(prec),mixed m))
   main_prec        := s
   secn_prec        := m
else
   main_prec        := d
   secn_prec        := d
endif

# -- Verbosity control
ifeq ($(VERBOSE),0)
   AT := @
endif

# -- Program suffix configuration
bin_suffix_       :=
ifeq ($(NVSHMEM_SUPPORT),1)
   bin_suffix_    := $(bin_suffix_)_sh
endif
ifeq ($(PLUMED_SUPPORT),1)
   bin_suffix_    := $(bin_suffix_)_plm
endif
ifeq ($(COLVARS_SUPPORT),1)
   bin_suffix_    := $(bin_suffix_)_clv
endif
ifeq ($(NN_SUPPORT),1)
   bin_suffix_    := $(bin_suffix_)_nn
endif
ifeq ($(ORTHO_BOX_ONLY_SUPPORT),1)
   bin_suffix_    := $(bin_suffix_)_ob
endif
ifeq ($(prec),$(filter $(prec),single s))
   bin_suffix_    := $(bin_suffix_).single
else ifeq ($(prec),$(filter $(prec),mixed m))
ifeq ($(FPA_SUPPORT),1)
   bin_suffix_    := $(bin_suffix_).fpa
else
   bin_suffix_    := $(bin_suffix_).mixed
endif
else
   ifeq ($(arch),$(filter $(arch),device gpu))
      bin_suffix_ := $(bin_suffix_).gpu
   else
      bin_suffix_ := $(bin_suffix_).cpu
   endif
endif
ifeq ($(prog_suffix),$(empty__))
   bin_sufx       := $(bin_suffix_)
else
   bin_sufx       := $(prog_suffix)
endif


#
# ------------------
# Preprocessor Flags
# ------------------
pp_flags_common_    :=
pp_flags_f_         := -cpp
pp_flags_cxx_       :=
pp_flags_cuda_c_    :=
pp_flags_cuda_f_    := -cpp
ifeq ($(opt),debug+)
   pp_flags_common_ += -DTINKER_DEBUG
endif
ifeq ($(prec),$(filter $(prec),single s))
   pp_flags_common_ += -DSINGLE
   pp_flags_cuda_c_ += -DUSE_ERFC_HASTINGS
else ifeq ($(prec),$(filter $(prec),mixed m))
   pp_flags_common_ += -DMIXED
   pp_flags_cuda_c_ += -DUSE_ERFC_HASTINGS
endif
ifeq ($(NVSHMEM_SUPPORT),1)
ifeq ($(arch),$(filter $(arch),device gpu))
   pp_flags_common_ += -DUSE_NVSHMEM
endif
endif
ifeq ($(NVTX_SUPPORT),1)
   pp_flags_common_ += -DUSE_NVTX
endif
ifeq ($(ORTHO_BOX_ONLY_SUPPORT),1)
   pp_flags_common_ += -DORTHOGONAL_BOX_SHAPE_ONLY
endif
ifeq ($(NO_MUTATION),1)
   pp_flags_common_ += -DTINKER_NO_MUTATE
endif
ifeq ($(FPA_SUPPORT),1)
ifeq ($(prec),$(filter $(prec),mixed m))
   pp_flags_common_ += -DUSE_DETERMINISTIC_REDUCTION
endif
endif
ifeq ($(PLUMED_SUPPORT),1)
   pp_flags_common_ += -DPLUMED
endif
ifeq ($(COLVARS_SUPPORT),1)
   pp_flags_common_ += -DCOLVARS
   ifeq ($(TCL_SUPPORT),1)
      pp_flags_common_ += -DCOLVARS_TCL
   endif
endif
ifeq ($(NN_SUPPORT),1)
   pp_flags_common_ += -DNN_SUPPORT
endif
ifneq ($(add_options),$(empty__))
   pp_flags_common_ +=$(add_options)
endif
ifneq ($(add_options_f),$(empty__))
   pp_flags_f_      += $(add_options_f)
endif
ifneq ($(add_options_cuf),$(empty__))
   pp_flags_cuda_f_ += $(add_options_cuf)
endif
ifneq ($(add_options_cu),$(empty__))
   pp_flags_cuda_c_ += $(add_options_cu)
endif
pp_flags_f_         := $(pp_flags_f_) $(pp_flags_common_)
pp_flags_cxx_       += $(pp_flags_common_)
pp_flags_cuda_c_    += $(pp_flags_common_)
pp_flags_cuda_f_    += $(pp_flags_common_)


#
# -----------------
# Compilation Flags
# -----------------
legacy_flags        := off
cxx_std             := 14
compute_capability  := 60,70
cuda_version        := 11.7
device_c_comp       := -ccbin $(RunCXX)

# -- Nvidia Device's compute capability
cc_list_            := $(subst $(comma__), ,$(compute_capability))
# 60 70 => cc60 cc70
nvidia_cc_f__       := $(foreach var,$(cc_list_),cc$(var))
# cc60 cc70 => cc60,cc70
nvidia_cc_f__       := $(subst $(space__),$(comma__),$(nvidia_cc_f__))
# 60 70 => -gencode arch=compute_60,code=sm_60 -gencode arch=compute_70,code=sm_70
nvidia_cc_cxx__     := $(foreach var,$(cc_list_),-gencode arch=compute_$(var)$(comma__)code=sm_$(var))


# -- GPU Device property
ifeq ($(prec),$(filter $(prec),mixed m single s))
   nvidia_prop_sufx__  := ,fastmath
endif
ifeq ($(legacy_flags),off)
   nvidia_prop_prfx__  := -gpu=
endif
nvidia_prop_f_      := $(nvidia_cc_f__),cuda$(cuda_version),unroll$(nvidia_prop_sufx__)
nvidia_prop_cxx_    := $(nvidia_cc_f__),cuda$(cuda_version),unroll$(nvidia_prop_sufx__)
amd_prop_f_         :=

ifeq ($(legacy_flags),on)
   oacc_flags_f_      := -ta=tesla:$(nvidia_prop_f_)
   oacc_flags_cxx_    := -ta=tesla:$(nvidia_prop_cxx_)
   omp_flags_f_       :=
   omp_flags_cxx_     :=
   cuda_flags_f_      := -Mcuda=$(nvidia_prop_f_)
   cuda_flags_c_      := $(device_c_comp) $(nvidia_cc_cxx__)
else
   device_flags_f_    := $(nvidia_prop_prfx__)$(nvidia_prop_f_)
   device_flags_cxx_  := $(nvidia_prop_prfx__)$(nvidia_prop_cxx_)
   oacc_flags_f_      := -acc
   oacc_flags_cxx_    := -acc
   omp_flags_f_       := -mp=gpu
   omp_flags_cxx_     := -mp=gpu
   cuda_flags_f_      := -cuda
   cuda_flags_c_      := $(device_c_comp) $(nvidia_cc_cxx__)
   hip_flags_f_       :=
   hip_flags_cxx_     :=
endif
ifeq ($(NVSHMEM_SUPPORT), 1)
   cuda_flags_c_      += $(cuda_c_compiler) -rdc=true -I$(INC_NVSHMEM)
endif
ifeq ($(prec),$(filter $(prec),mixed m single s))
   cuda_flags_c_      += --use_fast_math
endif

# -- CPU Optimisation flags
inline_f__            := -Minline=maxsize:340
ifeq ($(opt), release)
   host_flags_f_      := -traceback -g -fast -Mdalign
   host_flags_cxx_    := -traceback -g -fast -Mdalign -std=c++$(cxx_std)
   host_flags_cuda_c_ := -std=c++$(cxx_std) -g -O3
   host_flags_cuda_f_ := -traceback -g -fast -Mdalign $(inline_f__)
else ifeq ($(opt), debug)
   host_flags_f_      := -traceback -g -O0 -Mbounds
   host_flags_cxx_    := -std=c++$(cxx_std) -traceback -g -O0 -Mbounds
   host_flags_cuda_c_ := -std=c++$(cxx_std) -g -O0
   host_flags_cuda_f_ := -traceback -g -O0 -Mbounds
else ifeq ($(opt), debug+)
   host_flags_f_      := -traceback -g -O0 -C -Mbounds
   host_flags_cxx_    := -std=c++$(cxx_std) -traceback -g -O0 -C -Mbounds
   host_flags_cuda_c_ := -std=c++$(cxx_std) -g -O0
   host_flags_cuda_f_ := -traceback -g -O0 -C -Mbounds
endif
# -- Static/Dynamic binaries
ifneq ($(library_link__),static)
   host_flags_f_      += -fPIC
   host_flags_cxx_    += -fPIC
   host_flags_cuda_c_ += --compiler-options -fPIC
endif
# -- Fortran Arightmetic precision
ifeq ($(prec),$(filter $(prec),mixed m single s))
   host_flags_f_      += -r4
   host_flags_cuda_f_ += -r4
else
   host_flags_f_      += -r8
   host_flags_cuda_f_ += -r8
endif
# -- Additional flags
ifneq ($(add_host),$(empty__))
   host_flags_f_      += $(add_host)
   host_flags_cxx_    += $(add_host)
   host_flags_cuda_c_ += $(add_host)
   host_flags_cuda_f_ += $(add_host)
endif
ifneq ($(add_host_f),$(empty__))
   host_flags_f_      += $(add_host_f)
endif
ifneq ($(add_host_cu),$(empty__))
   host_flags_cuda_c_ += $(add_host_cu)
endif
ifneq ($(add_host_cuf),$(empty__))
   host_flags_cuda_f_ += $(add_host_cuf)
endif

# -- Build final flags
comp_flags_f_         := $(pp_flags_f_) $(host_flags_f_)
decp_flags_f_         := $(host_flags_f_)
ifeq ($(arch),$(filter $(arch),device gpu))
ifeq ($(OPENMP_SUPPORT),1)
   comp_flags_f_      += $(omp_flags_f_)
   decp_flags_f_      += $(omp_flags_f_)
else
   comp_flags_f_      += $(oacc_flags_f_)
   decp_flags_f_      += $(oacc_flags_f_)
endif
ifneq ($(device_flags_f_),$(empty__))
   comp_flags_f_      += $(device_flags_f_)
   decp_flags_f_      += $(device_flags_f_)
endif
endif
comp_flags_mod_f_     := $(comp_flags_f_) $(inline_f__)
ifeq ($(arch),$(filter $(arch),device gpu))
   ifeq ($(HIP_SUPPORT),0)
   comp_flags_mod_f_     += $(cuda_flags_f_)
   endif
endif
comp_flags_mod_f_     += -I$(INC_WRAP)
comp_flags_gpu_f_     := $(comp_flags_mod_f_)
comp_flags_cxx_       := $(pp_flags_cxx_) $(host_flags_cxx_)
ifeq ($(arch),$(filter $(arch),device gpu))
ifeq ($(OPENMP_SUPPORT),1)
   comp_flags_cxx_    += $(omp_flags_cxx_)
else
   comp_flags_cxx_    += $(oacc_flags_cxx_)
endif
ifneq ($(device_flags_cxx_),$(empty__))
   comp_flags_cxx_    += $(device_flags_cxx_)
endif
endif
ifeq ($(COLVARS_SUPPORT),1)
   comp_flags_cxx_    += -I$(INC_COLVARS)
   ifeq ($(TCL_SUPPORT),1)
      comp_flags_cxx_ += -I$(INC_TCL)
   endif
endif
comp_flags_cuda_f_    := $(pp_flags_cuda_f_) $(host_flags_cuda_f_) $(cuda_flags_f_)
ifneq ($(device_flags_f_),$(empty__))
comp_flags_cuda_f_    += $(device_flags_f_)
endif
comp_flags_cuda_c_    := $(pp_flags_cuda_c_) $(host_flags_cuda_c_) $(cuda_flags_c_)
comp_flags_hip_c_     := $(pp_flags_cuda_c_) $(host_flags_hip_c_)

# --- Finalize compilation flags
FFLAGS                 = $(comp_flags_f_)
DECOMP_FFLAGS          = $(decp_flags_f_)
GPUFLAGS               = $(comp_flags_gpu_f_)
CUFFLAGS               = $(comp_flags_cuda_f_)
CXXFLAGS               = $(comp_flags_cxx_)
COLVARS_CXX_F          = $(host_flags_cxx_)
ifeq ($(HIP_SUPPORT),0)
   CUCFLAGS            = $(comp_flags_cuda_c_)
else
   CUCFLAGS            = $(comp_flags_hip_c_)
endif
FFLAGS2                = $(comp_flags_gpu_f_)

#
#  ------------------------------------
#  Configure libraries and dependencies
#  ------------------------------------
ifeq ($(library_link__),static)
   LIBS                = libtinker.a
   fobj_lext__        := o
else
   LIBS                = -L. -ltinker
   fobj_lext__        := f
endif
ifeq ($(NN_SUPPORT),1)
   LIBS               += lib$(mlplugin_name).dylib
endif
LDLIBS                 = -lm -L$(LIB_FFTDECOMP) -L$(LIB_CPP)
ifeq ($(add_mkl__),yes)
   LDLIBS             += -L$(LIB_MKL)
else
   LDLIBS             += $(LIB_LAPACK)
endif

dev_ldlibs            :=
ifeq ($(HIP_SUPPORT),0)
   ifeq ($(legacy_flags),on)
      dev_ldlibs      += -Mcudalib=curand,cufft,cublas -Mcuda=$(acc_cc_flag__) -lcusolver
   else
      dev_ldlibs      += -cudalib=curand,cufft,cublas -lcusolver
   endif
   ifeq ($(findstring NVTX, $(OPTIONS)),NVTX)
      dev_ldlibs      += -lnvToolsExt
   endif
endif
dev_ldlibs            += -L$(LIB_WRAP)
depend_targets         = 2decomp_fft
ifeq ($(PLUMED_SUPPORT),1)
   LDLIBS             += -L$(LIB_PLUMED)
   depend_targets     += plumed
endif
ifeq ($(COLVARS_SUPPORT),1)
   LDLIBS             += -L$(LIB_COLVARS)
   ifeq ($(TCL_SUPPORT),1)
      ifeq ($(TCL_HOME),$(TCL_HOME__))
         depend_targets += tcl
      endif
      LDLIBS          += -L$(LIB_TCL)
   endif
   depend_targets     += colvars
   CULIBS += -lmpi_cxx
endif
ifeq ($(NVSHMEM_SUPPORT),1)
   dev_ldlibs         += -L$(LIB_NVSHMEM) -L/$(LIB_CUDA_DRIVER) -lmpi_cxx
endif
ifneq ($(FFTW_SUPPORT),0)
   LDLIBS             += -L$(LIB_FFT)
endif
ifeq ($(arch),$(filter $(arch),device gpu))
   LDLIBS             += $(dev_ldlibs)
   depend_targets     += thrust
   ifeq ($(NN_SUPPORT),1)
      depend_targets  += $(mlplugin_name)
   endif
endif
