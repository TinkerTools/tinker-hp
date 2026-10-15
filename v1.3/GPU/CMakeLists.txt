cmake_minimum_required(VERSION 3.31)

option(TINKERHP_INTERNAL_VALIDATE_ONLY "Internal syntax-only configure path." OFF)
mark_as_advanced(TINKERHP_INTERNAL_VALIDATE_ONLY)

set(_tinkerhp_languages C CXX Fortran)
if(NOT TINKERHP_INTERNAL_VALIDATE_ONLY)
  list(APPEND _tinkerhp_languages CUDA)
endif()

project(TinkerHP VERSION 0.1.0 LANGUAGES ${_tinkerhp_languages})

foreach(_lang CXX Fortran CUDA)
  set(_release_flags "${CMAKE_${_lang}_FLAGS_RELEASE}")
  if(_release_flags)
    string(REGEX REPLACE "(^| )-O[0-9s]([ ]|$)" " " _release_flags "${_release_flags}")
    if(_lang STREQUAL "Fortran")
      string(REGEX REPLACE "(^| )-fast([ ]|$)" " " _release_flags "${_release_flags}")
    endif()
    string(STRIP "${_release_flags}" _release_flags)
    set(CMAKE_${_lang}_FLAGS_RELEASE "${_release_flags}" CACHE STRING "Release flags for ${_lang}" FORCE)
  endif()
endforeach()

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}" CACHE PATH "Install prefix" FORCE)
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE OFF)

if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

include(GNUInstallDirs)
include(cmake/TinkerHPSourceLists.cmake)
include(cmake/TinkerHPPreflight.cmake)
include(cmake/TinkerHPBuildReport.cmake)
include(cmake/TinkerHPDeepHP.cmake)
include(cmake/TinkerHPColvars.cmake)
include(cmake/TinkerHPPlumed.cmake)

set(TINKERHP_BACKEND "GPU" CACHE STRING "Requested execution backend.")
set_property(CACHE TINKERHP_BACKEND PROPERTY STRINGS GPU CPU)

set(TINKERHP_PRECISION "DOUBLE" CACHE STRING "Floating-point precision mode.")
set_property(CACHE TINKERHP_PRECISION PROPERTY STRINGS DOUBLE MIXED SINGLE)

set(TINKERHP_LINK_MODE "STATIC" CACHE STRING "Internal Tinker-HP library link mode.")
set_property(CACHE TINKERHP_LINK_MODE PROPERTY STRINGS STATIC SHARED)

set(TINKERHP_FFT_BACKEND "GENERIC" CACHE STRING "2decomp FFT backend.")
set_property(CACHE TINKERHP_FFT_BACKEND PROPERTY STRINGS GENERIC FFTW3 MKL)

set(TINKERHP_GPU_ARCHITECTURES "native" CACHE STRING "NVIDIA GPU architectures: 'native', 'all', or a semicolon-separated SM list such as '70;80'.")
set_property(CACHE TINKERHP_GPU_ARCHITECTURES PROPERTY STRINGS native all)
set(TINKERHP_BINARY_SUFFIX "" CACHE STRING "Override executable suffix, including the leading dot.")
set(TINKERHP_INSTALL_BINDIR "${CMAKE_INSTALL_BINDIR}" CACHE STRING "Install destination for executables.")

option(TINKERHP_ENABLE_PLUMED "Enable Plumed integration." OFF)
option(TINKERHP_ENABLE_COLVARS "Enable Colvars integration." OFF)
option(TINKERHP_ENABLE_ORTHO_BOX_ONLY "Restrict the GPU build to orthogonal box shapes." OFF)
option(TINKERHP_ENABLE_FPA "Enable Fixed Precision Arithmetic in mixed precision." ON)
option(TINKERHP_ENABLE_NN "Enable neural-network plugin integration." OFF)
option(TINKERHP_ENABLE_NVTX "Enable NVTX markers." OFF)
option(TINKERHP_ENABLE_PARALLEL_SCIENTIFIC_TESTS
  "Register opt-in MPI scientific contracts." OFF)
option(TINKERHP_ENABLE_PERFORMANCE_TESTS
  "Register opt-in performance observations for maintained scientific workloads." OFF)
option(TINKERHP_ENABLE_OPENMP "Enable OpenMP GPU offload mode." OFF)
option(TINKERHP_ENABLE_HIP "Enable HIP compilation path." OFF)
option(TINKERHP_ENABLE_NO_MUTATION "Enable no-mutation mode." OFF)
option(TINKERHP_ENABLE_NVSHMEM "Enable NVSHMEM support." OFF)

option(TINKERHP_ENABLE_NVIDIA_21X_WORKAROUND "Apply the NVHPC 21.3-21.9 workaround automatically." ON)
set(TINKERHP_NVIDIA_21X_WORKAROUND_FLAG "-Mnollvm" CACHE STRING "Compiler flag used for NVHPC 21.3-21.9.")
set(TINKERHP_NVFORTRAN_OPT_LEVEL "FAST" CACHE STRING
  "NVFORTRAN Fortran optimization level (FAST, O1, or O0).")
set_property(CACHE TINKERHP_NVFORTRAN_OPT_LEVEL PROPERTY STRINGS FAST O1 O0)
string(TOUPPER "${TINKERHP_NVFORTRAN_OPT_LEVEL}" TINKERHP_NVFORTRAN_OPT_LEVEL_NORMALIZED)
if(TINKERHP_NVFORTRAN_OPT_LEVEL_NORMALIZED STREQUAL "FAST")
  set(TINKERHP_NVFORTRAN_OPT_FLAG -fast)
elseif(TINKERHP_NVFORTRAN_OPT_LEVEL_NORMALIZED STREQUAL "O1")
  set(TINKERHP_NVFORTRAN_OPT_FLAG -O1)
elseif(TINKERHP_NVFORTRAN_OPT_LEVEL_NORMALIZED STREQUAL "O0")
  set(TINKERHP_NVFORTRAN_OPT_FLAG -O0)
else()
  message(FATAL_ERROR
    "TINKERHP_NVFORTRAN_OPT_LEVEL must be FAST, O1, or O0; got "
    "${TINKERHP_NVFORTRAN_OPT_LEVEL}")
endif()
option(TINKERHP_ENABLE_NVFORTRAN_NOINLINE
  "Disable NVFORTRAN inlining for all Fortran targets." OFF)

tinkerhp_validate_requested_configuration()
if(TINKERHP_LINK_MODE STREQUAL "SHARED")
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
else()
  set(CMAKE_POSITION_INDEPENDENT_CODE OFF)
endif()
tinkerhp_preflight_toolchain(TINKERHP_DETECTED_CUDA_VERSION TINKERHP_WORKAROUND_FLAGS)
tinkerhp_configure_deep_hp()

set(TINKERHP_GPU_ARCHITECTURES_REQUESTED "${TINKERHP_GPU_ARCHITECTURES}")
if(TINKERHP_GPU_ARCHITECTURES_REQUESTED STREQUAL "native")
  set(TINKERHP_GPU_ARCHITECTURES_REPORT "native")
  set(TINKERHP_CUDA_ARCHITECTURES_VALUE native)
  set(TINKERHP_NVIDIA_FORTRAN_ARCH_LIST ccnative)
elseif(TINKERHP_GPU_ARCHITECTURES_REQUESTED STREQUAL "all")
  set(TINKERHP_GPU_ARCHITECTURES_REPORT "all")
  set(TINKERHP_CUDA_ARCHITECTURES_VALUE all)
  set(TINKERHP_NVIDIA_FORTRAN_ARCH_LIST ccall)
else()
  set(TINKERHP_CUDA_ARCHITECTURES_VALUE "${TINKERHP_GPU_ARCHITECTURES}")
  set(TINKERHP_NVIDIA_FORTRAN_ARCHES "")
  foreach(_arch IN LISTS TINKERHP_GPU_ARCHITECTURES)
    if(NOT _arch MATCHES "^[0-9]+$")
      message(FATAL_ERROR
        "TINKERHP_GPU_ARCHITECTURES must be 'native', 'all', or a semicolon-separated list of numeric SM values; got '${_arch}'."
      )
    endif()
    list(APPEND TINKERHP_NVIDIA_FORTRAN_ARCHES "cc${_arch}")
  endforeach()
  if(NOT TINKERHP_NVIDIA_FORTRAN_ARCHES)
    message(FATAL_ERROR
      "TINKERHP_GPU_ARCHITECTURES must not be empty; use 'native', 'all', or a semicolon-separated SM list."
    )
  endif()
  list(JOIN TINKERHP_GPU_ARCHITECTURES "," TINKERHP_GPU_ARCHITECTURES_REPORT)
  list(JOIN TINKERHP_NVIDIA_FORTRAN_ARCHES "," TINKERHP_NVIDIA_FORTRAN_ARCH_LIST)
endif()

if(TINKERHP_ENABLE_PLUMED AND TINKERHP_ENABLE_COLVARS)
  message(WARNING
    "TINKERHP_ENABLE_PLUMED and TINKERHP_ENABLE_COLVARS are both ON. "
    "This uncommon combined configuration is supported but has no dedicated preset and is not validated."
  )
endif()

if(TINKERHP_BINARY_SUFFIX STREQUAL "")
  set(TINKERHP_EFFECTIVE_BINARY_SUFFIX "")
  if(TINKERHP_ENABLE_NN)
    string(APPEND TINKERHP_EFFECTIVE_BINARY_SUFFIX "_nn")
  endif()
  if(TINKERHP_ENABLE_PLUMED)
    string(APPEND TINKERHP_EFFECTIVE_BINARY_SUFFIX "_plm")
  endif()
  if(TINKERHP_ENABLE_COLVARS)
    string(APPEND TINKERHP_EFFECTIVE_BINARY_SUFFIX "_clv")
  endif()
  if(TINKERHP_PRECISION STREQUAL "MIXED")
    if(TINKERHP_ENABLE_FPA)
      string(APPEND TINKERHP_EFFECTIVE_BINARY_SUFFIX ".fpa")
    else()
      string(APPEND TINKERHP_EFFECTIVE_BINARY_SUFFIX ".mixed")
    endif()
  else()
    string(APPEND TINKERHP_EFFECTIVE_BINARY_SUFFIX ".gpu")
  endif()
else()
  set(TINKERHP_EFFECTIVE_BINARY_SUFFIX "${TINKERHP_BINARY_SUFFIX}")
endif()

set(TINKERHP_MODULE_DIR "${CMAKE_BINARY_DIR}/mod")
set(TINKERHP_MAIN_MODULE_DIR "${TINKERHP_MODULE_DIR}/tinkerhp")
set(TINKERHP_WRAPPER_MODULE_DIR "${TINKERHP_MODULE_DIR}/wrappers")
set(TINKERHP_2DECOMP_MODULE_DIR "${TINKERHP_MODULE_DIR}/2decomp_fft")

set(TINKERHP_NVIDIA_GPU_FLAG
  "-gpu=${TINKERHP_NVIDIA_FORTRAN_ARCH_LIST},cuda${TINKERHP_DETECTED_CUDA_VERSION},unroll"
)
if(TINKERHP_NVFORTRAN_OPT_LEVEL_NORMALIZED STREQUAL "O0")
  # NVHPC 24.7 emits invalid NVVM debug metadata for O0 OpenACC code.
  string(APPEND TINKERHP_NVIDIA_GPU_FLAG ",nodebug")
endif()

if(TINKERHP_PRECISION STREQUAL "MIXED")
  string(APPEND TINKERHP_NVIDIA_GPU_FLAG ",fastmath")
  set(TINKERHP_FORTRAN_REAL_FLAG -r4)
  set(TINKERHP_MIXED_CUDA_FLAGS --use_fast_math)
else()
  set(TINKERHP_FORTRAN_REAL_FLAG -r8)
  set(TINKERHP_MIXED_CUDA_FLAGS)
endif()

set(TINKERHP_COMMON_FORTRAN_FLAGS
  -traceback
  -g
  ${TINKERHP_NVFORTRAN_OPT_FLAG}
  -Mdalign
  ${TINKERHP_FORTRAN_REAL_FLAG}
  -acc
  "${TINKERHP_NVIDIA_GPU_FLAG}"
  ${TINKERHP_WORKAROUND_FLAGS}
)

set(TINKERHP_GPU_INLINING_FLAGS)
if(NOT TINKERHP_ENABLE_NVFORTRAN_NOINLINE)
  list(APPEND TINKERHP_GPU_INLINING_FLAGS -Minline=maxsize:340)
endif()
set(TINKERHP_GPU_FORTRAN_EXTRA_FLAGS
  -cuda
  ${TINKERHP_GPU_INLINING_FLAGS}
)

set(TINKERHP_CUDA_FORTRAN_FLAGS
  -traceback
  -g
  ${TINKERHP_NVFORTRAN_OPT_FLAG}
  -Mdalign
  ${TINKERHP_FORTRAN_REAL_FLAG}
  -cuda
  "${TINKERHP_NVIDIA_GPU_FLAG}"
  ${TINKERHP_GPU_INLINING_FLAGS}
  ${TINKERHP_WORKAROUND_FLAGS}
)

set(TINKERHP_HOST_CXX_FLAGS
  -traceback
  -g
  -fast
  -Mdalign
  -std=c++14
  ${TINKERHP_WORKAROUND_FLAGS}
)

set(TINKERHP_CXX_FLAGS
  ${TINKERHP_HOST_CXX_FLAGS}
  -acc
  "${TINKERHP_NVIDIA_GPU_FLAG}"
)

set(TINKERHP_CUDA_FLAGS
  -std=c++14
  -g
  -O3
  -ccbin
  "${CMAKE_CXX_COMPILER}"
  ${TINKERHP_MIXED_CUDA_FLAGS}
)

set(TINKERHP_GPU_LINK_OPTIONS
  -acc
  "${TINKERHP_NVIDIA_GPU_FLAG}"
  -cuda
  -cudalib=curand,cufft,cublas
  -lcusolver
)

# NVHPC 24.7 can ICE when CMake's Unix Makefiles generator passes an absolute
# Fortran source path.  The compiler launcher keeps the generated graph intact
# while making project Fortran source paths relative to the build directory.
# It also removes the redundant _OPENACC definition from NVHPC invocations:
# -acc already provides it, whereas CMake still retains the definition for its
# Makefile conditional-module dependency scan.
if(CMAKE_GENERATOR STREQUAL "Unix Makefiles"
    AND CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
  set(TINKERHP_MAKEFILE_LAUNCHER_SOURCE_DIR "${CMAKE_SOURCE_DIR}")
  set(TINKERHP_NVHPC_MAKEFILE_FORTRAN_LAUNCHER
    "${CMAKE_BINARY_DIR}/cmake/tinkerhp-nvhpc-make-fortran-launcher.sh")
  file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/cmake")
  configure_file(
    "${CMAKE_SOURCE_DIR}/cmake/TinkerHPNVHPCMakeFortranLauncher.sh.in"
    "${TINKERHP_NVHPC_MAKEFILE_FORTRAN_LAUNCHER}"
    @ONLY
  )
  file(CHMOD "${TINKERHP_NVHPC_MAKEFILE_FORTRAN_LAUNCHER}"
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
  )
  set(CMAKE_Fortran_COMPILER_LAUNCHER "${TINKERHP_NVHPC_MAKEFILE_FORTRAN_LAUNCHER}")
endif()

# Fortran_PREPROCESS generates an intermediate used by CMake's module-dependency
# scanner.  Keep macro definitions in COMPILE_DEFINITIONS: CMake supplies them
# to that preprocessing step, but does not pass them to compilation of the
# already-preprocessed intermediate.
function(tinkerhp_add_preprocessor_definitions target_name)
  get_target_property(_target_type ${target_name} TYPE)
  if(_target_type STREQUAL "INTERFACE_LIBRARY")
    set(_definition_scope INTERFACE)
  else()
    set(_definition_scope PRIVATE)
  endif()
  if(NOT "${ARGN}" STREQUAL "")
    target_compile_definitions(${target_name} ${_definition_scope} ${ARGN})
  endif()
endfunction()

set(TINKERHP_MAIN_PREPROCESS_DEFINITIONS)
if(TINKERHP_PRECISION STREQUAL "MIXED")
  list(APPEND TINKERHP_MAIN_PREPROCESS_DEFINITIONS MIXED)
  if(TINKERHP_ENABLE_FPA)
    list(APPEND TINKERHP_MAIN_PREPROCESS_DEFINITIONS USE_DETERMINISTIC_REDUCTION)
  endif()
endif()
if(TINKERHP_ENABLE_ORTHO_BOX_ONLY)
  list(APPEND TINKERHP_MAIN_PREPROCESS_DEFINITIONS ORTHOGONAL_BOX_SHAPE_ONLY)
endif()
if(TINKERHP_ENABLE_NVTX)
  list(APPEND TINKERHP_MAIN_PREPROCESS_DEFINITIONS USE_NVTX)
endif()
if(TINKERHP_ENABLE_NO_MUTATION)
  list(APPEND TINKERHP_MAIN_PREPROCESS_DEFINITIONS TINKER_NO_MUTATE)
endif()
if(TINKERHP_ENABLE_NN)
  list(APPEND TINKERHP_MAIN_PREPROCESS_DEFINITIONS NN_SUPPORT)
endif()

set(TINKERHP_CUDA_PREPROCESS_DEFINITIONS)
if(TINKERHP_PRECISION STREQUAL "MIXED")
  list(APPEND TINKERHP_CUDA_PREPROCESS_DEFINITIONS USE_ERFC_HASTINGS)
endif()

# NVHPC defines these macros while compiling with -acc and -cuda.  The Unix
# Makefiles dependency scanner does not infer them from compile options, so
# pass them explicitly there to retain the conditional module graph.  Apply
# _CUDA to every Fortran target configured with -cuda. Ninja receives NVHPC's
# preprocessed scan input and must not receive redefinitions.
set(TINKERHP_OPENACC_PREPROCESS_DEFINITIONS)
set(TINKERHP_CUDA_FORTRAN_PREPROCESS_DEFINITIONS)
if(CMAKE_GENERATOR STREQUAL "Unix Makefiles")
  list(APPEND TINKERHP_OPENACC_PREPROCESS_DEFINITIONS _OPENACC)
  list(APPEND TINKERHP_CUDA_FORTRAN_PREPROCESS_DEFINITIONS _CUDA)
endif()

if(TINKERHP_PRECISION STREQUAL "DOUBLE")
  set(TINKERHP_2DECOMP_PREPROCESS_DEFINITIONS DOUBLE_PREC OVERWRITE)
else()
  set(TINKERHP_2DECOMP_PREPROCESS_DEFINITIONS OVERWRITE)
endif()

set(TINKERHP_PLUMED_PREPROCESS_DEFINITIONS)
if(TINKERHP_ENABLE_PLUMED)
  list(APPEND TINKERHP_PLUMED_PREPROCESS_DEFINITIONS PLUMED)
endif()

set(TINKERHP_COLVARS_PREPROCESS_DEFINITIONS)
if(TINKERHP_ENABLE_COLVARS)
  list(APPEND TINKERHP_COLVARS_PREPROCESS_DEFINITIONS COLVARS COLVARS_MPI COLVARS_TCL)
endif()

if(NOT TINKERHP_INTERNAL_VALIDATE_ONLY)
  set(CMAKE_CUDA_ARCHITECTURES "${TINKERHP_CUDA_ARCHITECTURES_VALUE}")
  set(CMAKE_CUDA_RUNTIME_LIBRARY None)
endif()

tinkerhp_configure_plumed()
tinkerhp_configure_colvars(${TINKERHP_HOST_CXX_FLAGS})
tinkerhp_add_deep_hp_plugin()

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

set(TINKERHP_BLAS_LINK_LIBRARIES
  -lblas_lp64
)
set(TINKERHP_LAPACK_LINK_LIBRARIES
  -llapack_lp64
)
set(TINKERHP_NVTX_LINK_LIBRARIES)
if(TINKERHP_ENABLE_NVTX)
  list(APPEND TINKERHP_NVTX_LINK_LIBRARIES -lnvToolsExt)
endif()

if(TINKERHP_INTERNAL_VALIDATE_ONLY)
  add_library(tinkerhp_wrapper INTERFACE)
  target_include_directories(tinkerhp_wrapper INTERFACE wrappers)
else()
  add_subdirectory(wrappers)
endif()

add_subdirectory(2decomp_fft)

tinkerhp_partition_fortran_sources(TINKERHP_FORTRAN_COMMON_SRCS TINKERHP_FORTRAN_GPU_SRCS)

# CMake preprocesses each source to discover Fortran module dependencies.  The
# compiler must consume that generated *.f90-pp.f90 file as Fortran, without a
# second -cpp pass (NVFORTRAN otherwise crashes on its CPP line markers).
# Compile submodules after their base and CUDA Fortran module providers.
set(TINKERHP_SUBMODULE_SRCS)
foreach(_source IN LISTS MOD_SRCS)
  if(_source MATCHES "/sm_.*\\.f90$")
    list(APPEND TINKERHP_SUBMODULE_SRCS "${_source}")
  endif()
endforeach()
list(REMOVE_ITEM MOD_SRCS ${TINKERHP_SUBMODULE_SRCS})
list(APPEND TINKERHP_FORTRAN_GPU_SRCS ${TINKERHP_SUBMODULE_SRCS})

if(TINKERHP_ENABLE_NVFORTRAN_NOINLINE)
  set(_tinkerhp_nvhpc_noinline_sources
    ${MOD_SRCS}
    ${TINKERHP_FORTRAN_COMMON_SRCS}
    ${TINKERHP_FORTRAN_GPU_SRCS}
    ${DEV_F_SRCS}
  )
  list(REMOVE_DUPLICATES _tinkerhp_nvhpc_noinline_sources)
  set_property(SOURCE ${_tinkerhp_nvhpc_noinline_sources}
    APPEND PROPERTY COMPILE_OPTIONS "-Mnoinline"
  )
endif()

set(TINKERHP_DIRECT_FORTRAN_PREPROCESS_SRCS
  ${MOD_SRCS}
  ${TINKERHP_FORTRAN_COMMON_SRCS}
  ${TINKERHP_FORTRAN_GPU_SRCS}
  ${DEV_F_SRCS}
)
foreach(_program IN LISTS TINKERHP_PROGRAM_TARGETS)
  list(APPEND TINKERHP_DIRECT_FORTRAN_PREPROCESS_SRCS "source/${_program}.f90")
endforeach()
list(REMOVE_DUPLICATES TINKERHP_DIRECT_FORTRAN_PREPROCESS_SRCS)
set_source_files_properties(${TINKERHP_DIRECT_FORTRAN_PREPROCESS_SRCS} PROPERTIES
  Fortran_PREPROCESS ON
)

function(tinkerhp_configure_fortran_target target_name module_dir extra_flags)
  set_target_properties(${target_name} PROPERTIES
    Fortran_MODULE_DIRECTORY "${module_dir}"
    LINKER_LANGUAGE Fortran
  )
  target_compile_options(${target_name} PRIVATE
    $<$<COMPILE_LANGUAGE:Fortran>:${TINKERHP_COMMON_FORTRAN_FLAGS};${extra_flags}>
  )
  target_include_directories(${target_name} PRIVATE
    "${TINKERHP_WRAPPER_MODULE_DIR}"
    "${CMAKE_SOURCE_DIR}/source"
    "${TINKERHP_MAIN_MODULE_DIR}"
    "${TINKERHP_2DECOMP_MODULE_DIR}"
  )
  target_link_libraries(${target_name} PRIVATE
    tinkerhp_2decomp_fft
    tinkerhp_wrapper
  )
  target_link_options(${target_name} PRIVATE ${TINKERHP_GPU_LINK_OPTIONS})
endfunction()

add_library(tinkerhp_modules OBJECT ${MOD_SRCS})
tinkerhp_configure_fortran_target(tinkerhp_modules "${TINKERHP_MAIN_MODULE_DIR}" "${TINKERHP_GPU_FORTRAN_EXTRA_FLAGS}")
add_dependencies(tinkerhp_modules tinkerhp_wrapper tinkerhp_2decomp_fft)

add_library(tinkerhp_fortran_common OBJECT ${TINKERHP_FORTRAN_COMMON_SRCS})
tinkerhp_configure_fortran_target(tinkerhp_fortran_common "${TINKERHP_MAIN_MODULE_DIR}" "")
add_dependencies(tinkerhp_fortran_common tinkerhp_wrapper)

add_library(tinkerhp_fortran_gpu OBJECT ${TINKERHP_FORTRAN_GPU_SRCS})
tinkerhp_configure_fortran_target(tinkerhp_fortran_gpu "${TINKERHP_MAIN_MODULE_DIR}" "${TINKERHP_GPU_FORTRAN_EXTRA_FLAGS}")
add_dependencies(tinkerhp_fortran_gpu tinkerhp_wrapper)

add_library(tinkerhp_cuda_fortran OBJECT ${DEV_F_SRCS})
set_target_properties(tinkerhp_cuda_fortran PROPERTIES
  Fortran_MODULE_DIRECTORY "${TINKERHP_MAIN_MODULE_DIR}"
)
target_compile_options(tinkerhp_cuda_fortran PRIVATE
  $<$<COMPILE_LANGUAGE:Fortran>:${TINKERHP_CUDA_FORTRAN_FLAGS}>
)
target_include_directories(tinkerhp_cuda_fortran PRIVATE
  "${TINKERHP_WRAPPER_MODULE_DIR}"
  "${CMAKE_SOURCE_DIR}/source"
  "${TINKERHP_MAIN_MODULE_DIR}"
)
add_dependencies(tinkerhp_cuda_fortran tinkerhp_wrapper)
target_link_libraries(tinkerhp_cuda_fortran PRIVATE
  tinkerhp_modules
  tinkerhp_2decomp_fft
  tinkerhp_wrapper
)
set_source_files_properties(${DEV_F_SRCS} PROPERTIES LANGUAGE Fortran)
set_source_files_properties(
  source/promoamd.f90
  source/promosmd.f90
  PROPERTIES COMPILE_OPTIONS "-Mbackslash"
)

target_link_libraries(tinkerhp_fortran_common PRIVATE
  tinkerhp_modules
  tinkerhp_cuda_fortran
)
target_link_libraries(tinkerhp_fortran_gpu PRIVATE
  tinkerhp_modules
  tinkerhp_cuda_fortran
)

if(TINKERHP_ENABLE_COLVARS)
  list(APPEND CC_SRCS source/colvarproxy_tinkerhp.cpp)
endif()

add_library(tinkerhp_cxx OBJECT ${CC_SRCS})
target_compile_options(tinkerhp_cxx PRIVATE
  $<$<COMPILE_LANGUAGE:CXX>:${TINKERHP_CXX_FLAGS}>
)
target_include_directories(tinkerhp_cxx PRIVATE
  "${CMAKE_SOURCE_DIR}/source"
)
if(TINKERHP_ENABLE_COLVARS)
  target_include_directories(tinkerhp_cxx PRIVATE
    "${COLVARS_HOME}"
    "${TINKERHP_COLVARS_TCL_INCLUDE_DIR}"
  )
endif()

if(NOT TINKERHP_INTERNAL_VALIDATE_ONLY)
  add_library(tinkerhp_cuda OBJECT ${DEV_CC_SRCS})
  set_target_properties(tinkerhp_cuda PROPERTIES
    CUDA_ARCHITECTURES "${TINKERHP_CUDA_ARCHITECTURES_VALUE}"
  )
  target_compile_options(tinkerhp_cuda PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:${TINKERHP_CUDA_FLAGS}>
  )
  target_include_directories(tinkerhp_cuda PRIVATE
    "${CMAKE_SOURCE_DIR}/source"
  )
endif()

set(TINKERHP_COMMON_COMPILE_TARGETS
  tinkerhp_modules
  tinkerhp_fortran_common
  tinkerhp_fortran_gpu
  tinkerhp_cuda_fortran
  tinkerhp_cxx
)
if(TARGET tinkerhp_cuda)
  list(APPEND TINKERHP_COMMON_COMPILE_TARGETS tinkerhp_cuda)
endif()

function(tinkerhp_apply_main_preprocessor_definitions target_name)
  tinkerhp_add_preprocessor_definitions(${target_name}
    ${TINKERHP_MAIN_PREPROCESS_DEFINITIONS}
  )
endfunction()

foreach(_target IN LISTS TINKERHP_COMMON_COMPILE_TARGETS)
  tinkerhp_apply_main_preprocessor_definitions(${_target})
endforeach()
foreach(_target
    tinkerhp_2decomp_fft
    tinkerhp_modules
    tinkerhp_fortran_common
    tinkerhp_fortran_gpu)
  tinkerhp_add_preprocessor_definitions(${_target}
    ${TINKERHP_OPENACC_PREPROCESS_DEFINITIONS}
  )
endforeach()
foreach(_target
    tinkerhp_fortran_gpu
    tinkerhp_cuda_fortran)
  tinkerhp_add_preprocessor_definitions(${_target}
    ${TINKERHP_CUDA_FORTRAN_PREPROCESS_DEFINITIONS}
  )
endforeach()
if(TINKERHP_PRECISION STREQUAL "MIXED" AND TARGET tinkerhp_cuda)
  tinkerhp_add_preprocessor_definitions(tinkerhp_cuda
    ${TINKERHP_CUDA_PREPROCESS_DEFINITIONS}
  )
endif()

if(TINKERHP_ENABLE_PLUMED)
  set(_tinkerhp_plumed_compile_targets
    tinkerhp_modules
    tinkerhp_fortran_common
    tinkerhp_fortran_gpu
    tinkerhp_cuda_fortran
  )
  if(TARGET tinkerhp_cuda)
    list(APPEND _tinkerhp_plumed_compile_targets tinkerhp_cuda)
  endif()
  foreach(_target IN LISTS _tinkerhp_plumed_compile_targets)
    tinkerhp_add_preprocessor_definitions(${_target}
      ${TINKERHP_PLUMED_PREPROCESS_DEFINITIONS}
    )
  endforeach()
endif()

if(TINKERHP_ENABLE_COLVARS)
  set(_tinkerhp_colvars_compile_targets
    tinkerhp_modules
    tinkerhp_fortran_common
    tinkerhp_fortran_gpu
    tinkerhp_cuda_fortran
    tinkerhp_cxx
  )
  if(NOT TINKERHP_INTERNAL_VALIDATE_ONLY)
    list(APPEND _tinkerhp_colvars_compile_targets tinkerhp_cuda)
  endif()
  foreach(_target IN LISTS _tinkerhp_colvars_compile_targets)
    tinkerhp_add_preprocessor_definitions(${_target}
      ${TINKERHP_COLVARS_PREPROCESS_DEFINITIONS}
    )
  endforeach()
endif()

add_library(tinkerhp_core ${TINKERHP_LINK_MODE}
  $<TARGET_OBJECTS:tinkerhp_modules>
  $<TARGET_OBJECTS:tinkerhp_fortran_common>
  $<TARGET_OBJECTS:tinkerhp_fortran_gpu>
  $<TARGET_OBJECTS:tinkerhp_cuda_fortran>
  $<TARGET_OBJECTS:tinkerhp_cxx>
)

if(NOT TINKERHP_INTERNAL_VALIDATE_ONLY)
  target_sources(tinkerhp_core PRIVATE $<TARGET_OBJECTS:tinkerhp_cuda>)
endif()

target_link_libraries(tinkerhp_core PUBLIC
  tinkerhp_2decomp_fft
  tinkerhp_wrapper
  ${TINKERHP_BLAS_LINK_LIBRARIES}
  ${TINKERHP_LAPACK_LINK_LIBRARIES}
  ${TINKERHP_NVTX_LINK_LIBRARIES}
)

if(TINKERHP_LINK_MODE STREQUAL "SHARED")
  set_target_properties(tinkerhp_core PROPERTIES
    LINKER_LANGUAGE Fortran
    BUILD_WITH_INSTALL_RPATH TRUE
  )
  target_link_options(tinkerhp_core PRIVATE ${TINKERHP_GPU_LINK_OPTIONS})
endif()

if(TINKERHP_ENABLE_NN)
  target_link_libraries(tinkerhp_core PUBLIC tinkerhp_mlinterface)
endif()

if(TINKERHP_ENABLE_COLVARS)
  add_dependencies(tinkerhp_core tinkerhp_colvars)
  target_link_libraries(tinkerhp_core PUBLIC tinkerhp_colvars_link)
endif()

if(TINKERHP_ENABLE_PLUMED)
  add_dependencies(tinkerhp_core tinkerhp_plumed)
  target_link_libraries(tinkerhp_core PUBLIC tinkerhp_plumed_link)
endif()

target_include_directories(tinkerhp_core PUBLIC
  "${CMAKE_SOURCE_DIR}/source"
  "${TINKERHP_MAIN_MODULE_DIR}"
)

set_target_properties(tinkerhp_core PROPERTIES
  ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
  LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
)

set(_tinkerhp_runtime_rpaths)
if(TINKERHP_LINK_MODE STREQUAL "SHARED" OR TINKERHP_ENABLE_NN)
  if(APPLE)
    set(_tinkerhp_runtime_library_rpath "@loader_path/../${CMAKE_INSTALL_LIBDIR}")
  else()
    set(_tinkerhp_runtime_library_rpath "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
  endif()
  list(APPEND _tinkerhp_runtime_rpaths "${_tinkerhp_runtime_library_rpath}")
  if(TINKERHP_ENABLE_NN)
    list(APPEND _tinkerhp_runtime_rpaths "${TINKERHP_DEEPHP_PYTHON_LIBRARY_DIR}")
  endif()
endif()

foreach(_program IN LISTS TINKERHP_PROGRAM_TARGETS)
  add_executable(${_program} "source/${_program}.f90")
  tinkerhp_apply_main_preprocessor_definitions(${_program})
  set_target_properties(${_program} PROPERTIES
    Fortran_MODULE_DIRECTORY "${TINKERHP_MAIN_MODULE_DIR}"
    LINKER_LANGUAGE Fortran
    OUTPUT_NAME "${_program}${TINKERHP_EFFECTIVE_BINARY_SUFFIX}"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
  )
  set(_tinkerhp_program_fortran_flags ${TINKERHP_COMMON_FORTRAN_FLAGS})
  if(TINKERHP_LINK_MODE STREQUAL "SHARED" AND _program STREQUAL "radial")
    # radial is a serial analysis utility.  Its OpenACC compilation creates
    # device references to module data that NVHPC cannot resolve from a shared
    # core, so keep this utility host-only in shared-core builds.
    list(REMOVE_ITEM _tinkerhp_program_fortran_flags
      -acc
      "${TINKERHP_NVIDIA_GPU_FLAG}"
    )
  else()
    tinkerhp_add_preprocessor_definitions(${_program}
      ${TINKERHP_OPENACC_PREPROCESS_DEFINITIONS}
    )
  endif()
  target_compile_options(${_program} PRIVATE
    $<$<COMPILE_LANGUAGE:Fortran>:${_tinkerhp_program_fortran_flags}>
  )
  target_link_options(${_program} PRIVATE ${TINKERHP_GPU_LINK_OPTIONS})
  target_include_directories(${_program} PRIVATE
    "${CMAKE_SOURCE_DIR}/source"
  )
  add_dependencies(${_program}
    tinkerhp_modules
    tinkerhp_cuda_fortran
  )
  target_link_libraries(${_program} PRIVATE
    tinkerhp_core
    ${TINKERHP_BLAS_LINK_LIBRARIES}
    ${TINKERHP_LAPACK_LINK_LIBRARIES}
  )
  if(TINKERHP_ENABLE_NN)
    add_dependencies(${_program} tinkerhp_mlinterface_build)
  endif()
  if(TINKERHP_LINK_MODE STREQUAL "SHARED" OR TINKERHP_ENABLE_NN)
    set_target_properties(${_program} PROPERTIES
      BUILD_RPATH "${_tinkerhp_runtime_rpaths}"
      INSTALL_RPATH "${_tinkerhp_runtime_rpaths}"
      BUILD_WITH_INSTALL_RPATH TRUE
    )
  endif()
  install(TARGETS ${_program} RUNTIME DESTINATION "${TINKERHP_INSTALL_BINDIR}")
endforeach()

if(TINKERHP_INTERNAL_VALIDATE_ONLY AND TINKERHP_ENABLE_NN)
  set(_tinkerhp_runtime_rpath_probe_source
    "${CMAKE_CURRENT_BINARY_DIR}/tinkerhp_deephp_runtime_rpath_probe.c")
  file(WRITE "${_tinkerhp_runtime_rpath_probe_source}" "int main(void) { return 0; }\n")
  add_executable(tinkerhp-deephp-runtime-rpath-probe EXCLUDE_FROM_ALL
    "${_tinkerhp_runtime_rpath_probe_source}")
  set_target_properties(tinkerhp-deephp-runtime-rpath-probe PROPERTIES
    BUILD_RPATH "${_tinkerhp_runtime_rpaths}"
    INSTALL_RPATH "${_tinkerhp_runtime_rpaths}"
    BUILD_WITH_INSTALL_RPATH TRUE
  )
  install(TARGETS tinkerhp-deephp-runtime-rpath-probe
    RUNTIME DESTINATION "${TINKERHP_INSTALL_BINDIR}"
    COMPONENT tinkerhp-deephp-runtime-rpath-contract
  )
endif()

install(TARGETS tinkerhp_core tinkerhp_2decomp_fft
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
)

if(TINKERHP_ENABLE_NN)
  install(FILES "${TINKERHP_DEEPHP_PLUGIN}" DESTINATION "${CMAKE_INSTALL_LIBDIR}")
endif()

if(NOT TINKERHP_INTERNAL_VALIDATE_ONLY)
  install(TARGETS tinkerhp_wrapper
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  )
endif()

tinkerhp_generate_build_report(
  OUTPUT "${CMAKE_BINARY_DIR}/tinkerhp-build-report.txt"
  CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}"
  CMAKE_C_COMPILER "${CMAKE_C_COMPILER}"
  CMAKE_CXX_COMPILER "${CMAKE_CXX_COMPILER}"
  CMAKE_Fortran_COMPILER "${CMAKE_Fortran_COMPILER}"
  CMAKE_CUDA_COMPILER "${CMAKE_CUDA_COMPILER}"
  CMAKE_CXX_COMPILER_ID "${CMAKE_CXX_COMPILER_ID}"
  CMAKE_CXX_COMPILER_VERSION "${CMAKE_CXX_COMPILER_VERSION}"
  CMAKE_Fortran_COMPILER_ID "${CMAKE_Fortran_COMPILER_ID}"
  CMAKE_Fortran_COMPILER_VERSION "${CMAKE_Fortran_COMPILER_VERSION}"
  CMAKE_CUDA_COMPILER_ID "${CMAKE_CUDA_COMPILER_ID}"
  CMAKE_CUDA_COMPILER_VERSION "${CMAKE_CUDA_COMPILER_VERSION}"
  TINKERHP_BACKEND "${TINKERHP_BACKEND}"
  TINKERHP_PRECISION "${TINKERHP_PRECISION}"
  TINKERHP_LINK_MODE "${TINKERHP_LINK_MODE}"
  TINKERHP_NVFORTRAN_OPT_LEVEL "${TINKERHP_NVFORTRAN_OPT_LEVEL_NORMALIZED}"
  TINKERHP_ENABLE_ORTHO_BOX_ONLY "${TINKERHP_ENABLE_ORTHO_BOX_ONLY}"
  TINKERHP_ENABLE_FPA "${TINKERHP_ENABLE_FPA}"
  TINKERHP_ENABLE_NVFORTRAN_NOINLINE "${TINKERHP_ENABLE_NVFORTRAN_NOINLINE}"
  TINKERHP_ENABLE_NVTX "${TINKERHP_ENABLE_NVTX}"
  TINKERHP_ENABLE_NO_MUTATION "${TINKERHP_ENABLE_NO_MUTATION}"
  TINKERHP_ENABLE_NN "${TINKERHP_ENABLE_NN}"
  TINKERHP_DEEPHP_ENVIRONMENT_ORIGIN "${TINKERHP_DEEPHP_ENVIRONMENT_ORIGIN}"
  TINKERHP_DEEPHP_PYTHON "${TINKERHP_DEEPHP_PYTHON}"
  TINKERHP_DEEPHP_PYTHON_VERSION "${TINKERHP_DEEPHP_PYTHON_VERSION}"
  TINKERHP_DEEPHP_PYTHON_LIBRARY "${TINKERHP_DEEPHP_PYTHON_LIBRARY}"
  TINKERHP_DEEPHP_CUDA_PROFILE "${TINKERHP_DEEPHP_CUDA_PROFILE}"
  TINKERHP_EFFECTIVE_BINARY_SUFFIX "${TINKERHP_EFFECTIVE_BINARY_SUFFIX}"
  TINKERHP_GPU_ARCHITECTURES "${TINKERHP_GPU_ARCHITECTURES_REPORT}"
  TINKERHP_DETECTED_CUDA_VERSION "${TINKERHP_DETECTED_CUDA_VERSION}"
  TINKERHP_ENABLE_PLUMED "${TINKERHP_ENABLE_PLUMED}"
  PLUMED_SOURCE_HOME "${PLUMED_SOURCE_HOME}"
  PLUMED_HOME "${PLUMED_HOME}"
  TINKERHP_PLUMED_LIBRARY "${TINKERHP_PLUMED_LIBRARY}"
  TINKERHP_PLUMED_KERNEL_LIBRARY "${TINKERHP_PLUMED_KERNEL_LIBRARY}"
  TINKERHP_ENABLE_COLVARS "${TINKERHP_ENABLE_COLVARS}"
  COLVARS_HOME "${COLVARS_HOME}"
  TCL_HOME "${TCL_HOME}"
  TINKERHP_COLVARS_TCL_INCLUDE_DIR "${TINKERHP_COLVARS_TCL_INCLUDE_DIR}"
  TINKERHP_COLVARS_TCL_LIBRARY "${TINKERHP_COLVARS_TCL_LIBRARY}"
  TINKERHP_COLVARS_MPI_VERSION "${TINKERHP_COLVARS_MPI_VERSION}"
  TINKERHP_COLVARS_MPI_CXX_LIBRARY "${TINKERHP_COLVARS_MPI_CXX_LIBRARY}"
  TINKERHP_COMMON_FORTRAN_FLAGS "${TINKERHP_COMMON_FORTRAN_FLAGS}"
  TINKERHP_CUDA_FORTRAN_FLAGS "${TINKERHP_CUDA_FORTRAN_FLAGS}"
  TINKERHP_CXX_FLAGS "${TINKERHP_CXX_FLAGS}"
  TINKERHP_CUDA_FLAGS "${TINKERHP_CUDA_FLAGS}"
  TINKERHP_MAIN_PREPROCESS_DEFINITIONS "${TINKERHP_MAIN_PREPROCESS_DEFINITIONS}"
  TINKERHP_CUDA_PREPROCESS_DEFINITIONS "${TINKERHP_CUDA_PREPROCESS_DEFINITIONS}"
  TINKERHP_2DECOMP_PREPROCESS_DEFINITIONS "${TINKERHP_2DECOMP_PREPROCESS_DEFINITIONS}"
  TINKERHP_PLUMED_PREPROCESS_DEFINITIONS "${TINKERHP_PLUMED_PREPROCESS_DEFINITIONS}"
  TINKERHP_COLVARS_PREPROCESS_DEFINITIONS "${TINKERHP_COLVARS_PREPROCESS_DEFINITIONS}"
  CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}"
  TINKERHP_INSTALL_BINDIR "${TINKERHP_INSTALL_BINDIR}"
)

option(BUILD_TESTING "Build Tinker-HP tests." OFF)
include(CTest)
if(BUILD_TESTING)
  find_package(Python3 REQUIRED COMPONENTS Interpreter)
  # Dimension-prefixed labels (scope, claim, method) follow the test vocabulary
  # in GLOSSARY.md. APPEND keeps labels already set on a test.
  function(tinkerhp_label_test name)
    set_property(TEST "${name}" APPEND PROPERTY LABELS ${ARGN})
  endfunction()
  add_test(
    NAME tinkerhp-runtime-refactor-oracle-contract
    COMMAND "${Python3_EXECUTABLE}" -m unittest discover
      -s "${CMAKE_SOURCE_DIR}/tests"
      -p test_runtime_refactor_oracle.py
  )
  tinkerhp_label_test(tinkerhp-runtime-refactor-oracle-contract scope:unit claim:functional method:conformance)
  add_test(
    NAME tinkerhp-generated-dependency-closure
    COMMAND "${CMAKE_COMMAND}"
      -DTINKERHP_BUILD_DIR=${CMAKE_BINARY_DIR}
      -P "${CMAKE_SOURCE_DIR}/cmake/TinkerHPGeneratedDependencyValidation.cmake"
  )
  tinkerhp_label_test(tinkerhp-generated-dependency-closure scope:static claim:functional method:conformance)
  add_test(
    NAME tinkerhp-cmake-configuration-contract
    COMMAND "${CMAKE_COMMAND}"
      -DTINKERHP_SOURCE_DIR=${CMAKE_SOURCE_DIR}
      -DTINKERHP_VALIDATION_DIR=${CMAKE_BINARY_DIR}/cmake-configuration-validation
      -P "${CMAKE_SOURCE_DIR}/cmake/TinkerHPConfigurationValidation.cmake"
  )
  tinkerhp_label_test(tinkerhp-cmake-configuration-contract scope:static claim:functional method:conformance)
  add_executable(tinkerhp-dynamic-cli-contract tests/test_dynamic_cli.f90)
  set_source_files_properties(tests/test_dynamic_cli.f90 PROPERTIES
    Fortran_PREPROCESS ON)
  tinkerhp_configure_fortran_target(tinkerhp-dynamic-cli-contract
    "${TINKERHP_MAIN_MODULE_DIR}" "")
  target_link_libraries(tinkerhp-dynamic-cli-contract PRIVATE tinkerhp_core)
  add_test(
    NAME tinkerhp-dynamic-cli-contract
    COMMAND "$<TARGET_FILE:tinkerhp-dynamic-cli-contract>"
  )
  tinkerhp_label_test(tinkerhp-dynamic-cli-contract scope:integration claim:functional method:conformance)
  add_test(
    NAME tinkerhp-dynamic-energy-format-fixed-contract
    COMMAND "$<TARGET_FILE:tinkerhp-dynamic-cli-contract>" fixed
  )
  tinkerhp_label_test(tinkerhp-dynamic-energy-format-fixed-contract scope:integration claim:functional method:conformance)
  add_test(
    NAME tinkerhp-dynamic-energy-format-scientific-contract
    COMMAND "$<TARGET_FILE:tinkerhp-dynamic-cli-contract>" scientific
  )
  tinkerhp_label_test(tinkerhp-dynamic-energy-format-scientific-contract scope:integration claim:functional method:conformance)
  set_tests_properties(
    tinkerhp-dynamic-cli-contract
    tinkerhp-dynamic-energy-format-fixed-contract
    tinkerhp-dynamic-energy-format-scientific-contract
    PROPERTIES
    ENVIRONMENT "TINKERHP_SOURCE_DIR=${CMAKE_SOURCE_DIR};TINKER_PARAMS=."
    RUN_SERIAL TRUE
  )
  if(TINKERHP_ENABLE_COLVARS)
    add_executable(tinkerhp-colvars-lifecycle-contract
      tests/test_colvars_lifecycle.f90)
    tinkerhp_configure_fortran_target(tinkerhp-colvars-lifecycle-contract
      "${TINKERHP_MAIN_MODULE_DIR}" "")
    target_link_libraries(tinkerhp-colvars-lifecycle-contract PRIVATE
      tinkerhp_core)
    foreach(rank_count IN ITEMS 1 2)
      add_test(
        NAME tinkerhp-colvars-lifecycle-contract-${rank_count}
        COMMAND mpirun -n ${rank_count}
          "$<TARGET_FILE:tinkerhp-colvars-lifecycle-contract>"
      )
      set_tests_properties(tinkerhp-colvars-lifecycle-contract-${rank_count}
        PROPERTIES
          ENVIRONMENT
            "OMPI_MCA_oob_tcp_if_include=lo;OMPI_MCA_btl=self,vader"
          PROCESSORS ${rank_count}
          RUN_SERIAL TRUE
      )
      tinkerhp_label_test(tinkerhp-colvars-lifecycle-contract-${rank_count}
        scope:integration claim:functional method:conformance)
    endforeach()
    function(tinkerhp_add_colvars_runtime_test name ranks)
      add_test(
        NAME "${name}"
        COMMAND "${Python3_EXECUTABLE}"
          "${CMAKE_SOURCE_DIR}/tests/test_colvars_runtime.py"
          --source-dir "${CMAKE_SOURCE_DIR}"
          --executable "$<TARGET_FILE:dynamic>"
          --ranks "${ranks}"
          ${ARGN}
      )
      set_tests_properties("${name}" PROPERTIES
        ENVIRONMENT "OMPI_MCA_oob_tcp_if_include=lo;OMPI_MCA_btl=self,vader"
        PROCESSORS "${ranks}"
        RUN_SERIAL TRUE
        TIMEOUT 180
      )
      if("--failure" IN_LIST ARGN)
        tinkerhp_label_test("${name}" scope:integration claim:functional method:failure-robustness)
      else()
        tinkerhp_label_test("${name}" scope:integration claim:functional method:conformance)
      endif()
    endfunction()
    tinkerhp_add_colvars_runtime_test(tinkerhp-colvars-runtime-activation-1 1)
    tinkerhp_add_colvars_runtime_test(tinkerhp-colvars-runtime-activation-2 2)
    tinkerhp_add_colvars_runtime_test(tinkerhp-colvars-runtime-repeated-activation 1
      --repeat-finalize)
    tinkerhp_add_colvars_runtime_test(tinkerhp-colvars-runtime-failure-1 1 --failure)
    tinkerhp_add_colvars_runtime_test(tinkerhp-colvars-runtime-failure-2 2 --failure)
  endif()
  if(TINKERHP_ENABLE_PLUMED)
    function(tinkerhp_add_plumed_runtime_test name)
      add_test(
        NAME "${name}"
        COMMAND "${Python3_EXECUTABLE}"
          "${CMAKE_SOURCE_DIR}/tests/test_plumed_runtime.py"
          --source-dir "${CMAKE_SOURCE_DIR}"
          --executable "$<TARGET_FILE:dynamic>"
      )
      set_tests_properties("${name}" PROPERTIES
        ENVIRONMENT
          "LD_LIBRARY_PATH=${PLUMED_HOME}/lib:$ENV{LD_LIBRARY_PATH};OMPI_MCA_oob_tcp_if_include=lo;OMPI_MCA_btl=self,vader"
        PROCESSORS 1
        RUN_SERIAL TRUE
        TIMEOUT 180
      )
      tinkerhp_label_test("${name}" scope:integration claim:functional method:conformance)
    endfunction()
    tinkerhp_add_plumed_runtime_test(tinkerhp-plumed-runtime-activation-1)
  endif()
  add_test(
    NAME tinkerhp-minimize-xyz-contract
    COMMAND "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/tests/test_minimize_xyz.py"
      --source-dir "${CMAKE_SOURCE_DIR}"
      --executable "$<TARGET_FILE:minimize>"
      --build-report "${CMAKE_BINARY_DIR}/tinkerhp-build-report.txt"
      --report "${CMAKE_BINARY_DIR}/minimize-xyz-contract/report.json"
      --precision "${TINKERHP_PRECISION}"
      --backend "${TINKERHP_BACKEND}"
      --mpi-layout 1x1x1
  )
  tinkerhp_label_test(tinkerhp-minimize-xyz-contract scope:integration claim:functional method:conformance)
  set_tests_properties(tinkerhp-minimize-xyz-contract PROPERTIES
    ENVIRONMENT "OMPI_MCA_oob_tcp_if_include=lo;OMPI_MCA_btl=self,vader"
    RUN_SERIAL TRUE
    SKIP_RETURN_CODE 77
    TIMEOUT 180
  )
  if(TINKERHP_PRECISION STREQUAL "MIXED")
    add_test(
      NAME tinkerhp-mdstat-skipped-energy-boundary-mixed
      COMMAND "${Python3_EXECUTABLE}"
        "${CMAKE_SOURCE_DIR}/tests/test_mdstat_energy_boundary.py"
        --source-dir "${CMAKE_SOURCE_DIR}"
        --executable "$<TARGET_FILE:dynamic>"
        --mpi-launcher mpirun
    )
    tinkerhp_label_test(tinkerhp-mdstat-skipped-energy-boundary-mixed scope:integration claim:functional method:conformance)
    set_tests_properties(tinkerhp-mdstat-skipped-energy-boundary-mixed PROPERTIES
      ENVIRONMENT "OMPI_MCA_oob_tcp_if_include=lo;OMPI_MCA_btl=self,vader"
      PROCESSORS 1
      RUN_SERIAL TRUE
      SKIP_RETURN_CODE 77
      TIMEOUT 300
    )
  endif()
  add_executable(tinkerhp-restart-contract tests/test_restart.f90)
  set_source_files_properties(tests/test_restart.f90 PROPERTIES
    Fortran_PREPROCESS ON)
  tinkerhp_configure_fortran_target(tinkerhp-restart-contract
    "${TINKERHP_MAIN_MODULE_DIR}" "")
  target_link_libraries(tinkerhp-restart-contract PRIVATE tinkerhp_core)
  add_test(
    NAME tinkerhp-restart-contract
    COMMAND "$<TARGET_FILE:tinkerhp-restart-contract>"
  )
  tinkerhp_label_test(tinkerhp-restart-contract scope:integration claim:functional method:conformance)
  add_executable(tinkerhp-mdtraj-contract tests/test_mdtraj.f90)
  set_source_files_properties(tests/test_mdtraj.f90 PROPERTIES
    Fortran_PREPROCESS ON)
  tinkerhp_configure_fortran_target(tinkerhp-mdtraj-contract
    "${TINKERHP_MAIN_MODULE_DIR}" "")
  target_link_libraries(tinkerhp-mdtraj-contract PRIVATE tinkerhp_core)
  add_test(
    NAME tinkerhp-mdtraj-contract
    COMMAND "$<TARGET_FILE:tinkerhp-mdtraj-contract>"
  )
  tinkerhp_label_test(tinkerhp-mdtraj-contract scope:integration claim:functional method:conformance method:failure-robustness)
  set_tests_properties(tinkerhp-mdtraj-contract PROPERTIES
    ENVIRONMENT "TINKERHP_SOURCE_DIR=${CMAKE_SOURCE_DIR};TINKER_PARAMS=."
  )
  add_test(
    NAME tinkerhp-replica-restart-contract
    COMMAND "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/tests/test_replica_restart_contract.py"
      --executable "$<TARGET_FILE:dynamic>"
      --source-dir "${CMAKE_SOURCE_DIR}"
      --mpi-launcher mpirun
  )
  tinkerhp_label_test(tinkerhp-replica-restart-contract scope:integration claim:functional method:conformance method:failure-robustness)
  set_tests_properties(tinkerhp-replica-restart-contract PROPERTIES
    ENVIRONMENT "OMPI_MCA_oob_tcp_if_include=lo;OMPI_MCA_btl=self,vader"
    RUN_SERIAL TRUE
    TIMEOUT 180
  )
  add_executable(tinkerhp-initialization-contract tests/test_initialization.f90)
  set_source_files_properties(tests/test_initialization.f90 PROPERTIES
    Fortran_PREPROCESS ON)
  tinkerhp_configure_fortran_target(tinkerhp-initialization-contract
    "${TINKERHP_MAIN_MODULE_DIR}" "")
  target_link_libraries(tinkerhp-initialization-contract PRIVATE tinkerhp_core)
  foreach(phase RANGE 1 12)
    add_test(
      NAME tinkerhp-initialization-contract-phase-${phase}
      COMMAND mpirun -n 1 "$<TARGET_FILE:tinkerhp-initialization-contract>" ${phase}
    )
    set_tests_properties(tinkerhp-initialization-contract-phase-${phase} PROPERTIES
      WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/examples"
      ENVIRONMENT "TINKERHP_SOURCE_DIR=${CMAKE_SOURCE_DIR};TINKER_PARAMS=.;OMPI_MCA_oob_tcp_if_include=lo;OMPI_MCA_btl=self,vader"
    )
    tinkerhp_label_test(tinkerhp-initialization-contract-phase-${phase}
      scope:integration claim:functional method:conformance
      method:failure-robustness)
  endforeach()
  add_executable(tinkerhp-parallel-contract tests/test_parallel.f90)
  set_source_files_properties(tests/test_parallel.f90 PROPERTIES
    Fortran_PREPROCESS ON)
  tinkerhp_configure_fortran_target(tinkerhp-parallel-contract
    "${TINKERHP_MAIN_MODULE_DIR}" "")
  tinkerhp_add_preprocessor_definitions(tinkerhp-parallel-contract
    ${TINKERHP_MAIN_PREPROCESS_DEFINITIONS}
    ${TINKERHP_OPENACC_PREPROCESS_DEFINITIONS}
  )
  target_link_libraries(tinkerhp-parallel-contract PRIVATE tinkerhp_core)
  add_executable(tinkerhp-parallel-topology-contract tests/test_parallel_topology.f90)
  set_source_files_properties(tests/test_parallel_topology.f90 PROPERTIES
    Fortran_PREPROCESS ON)
  tinkerhp_configure_fortran_target(tinkerhp-parallel-topology-contract
    "${TINKERHP_MAIN_MODULE_DIR}" "")
  tinkerhp_add_preprocessor_definitions(tinkerhp-parallel-topology-contract
    ${TINKERHP_MAIN_PREPROCESS_DEFINITIONS}
    ${TINKERHP_OPENACC_PREPROCESS_DEFINITIONS}
  )
  target_link_libraries(tinkerhp-parallel-topology-contract PRIVATE tinkerhp_core)
  add_test(
    NAME tinkerhp-parallel-contract
    COMMAND mpirun -n 1 "$<TARGET_FILE:tinkerhp-parallel-contract>"
  )
  tinkerhp_label_test(tinkerhp-parallel-contract scope:integration claim:functional method:conformance method:failure-robustness)
  set_tests_properties(tinkerhp-parallel-contract PROPERTIES
    ENVIRONMENT "TINKERHP_SOURCE_DIR=${CMAKE_SOURCE_DIR};TINKER_PARAMS=.;TINKER_ASYNC_COVER=1;OMPI_MCA_oob_tcp_if_include=lo;OMPI_MCA_btl=self,vader"
  )
  add_test(
    NAME tinkerhp-parallel-halo-contract
    COMMAND mpirun -n 2 "$<TARGET_FILE:tinkerhp-parallel-contract>"
  )
  tinkerhp_label_test(tinkerhp-parallel-halo-contract scope:integration claim:functional method:conformance method:failure-robustness)
  set_tests_properties(tinkerhp-parallel-halo-contract PROPERTIES
    ENVIRONMENT "TINKERHP_SOURCE_DIR=${CMAKE_SOURCE_DIR};TINKER_PARAMS=.;TINKER_ASYNC_COVER=1;OMPI_MCA_oob_tcp_if_include=lo;OMPI_MCA_btl=self,vader"
    PROCESSORS 2
    RUN_SERIAL TRUE
  )
  add_test(
    NAME tinkerhp-parallel-root-exchange-contract
    COMMAND mpirun -n 4 "$<TARGET_FILE:tinkerhp-parallel-topology-contract>"
  )
  tinkerhp_label_test(tinkerhp-parallel-root-exchange-contract scope:integration claim:functional method:conformance method:failure-robustness)
  set_tests_properties(tinkerhp-parallel-root-exchange-contract PROPERTIES
    ENVIRONMENT "TINKERHP_SOURCE_DIR=${CMAKE_SOURCE_DIR};TINKER_PARAMS=.;OMPI_MCA_oob_tcp_if_include=lo;OMPI_MCA_btl=self,vader"
    PROCESSORS 4
    RUN_SERIAL TRUE
  )
  add_executable(tinkerhp-parallel-release-contract tests/test_parallel_release.f90)
  set_source_files_properties(tests/test_parallel_release.f90 PROPERTIES
    Fortran_PREPROCESS ON)
  tinkerhp_configure_fortran_target(tinkerhp-parallel-release-contract
    "${TINKERHP_MAIN_MODULE_DIR}" "")
  tinkerhp_add_preprocessor_definitions(tinkerhp-parallel-release-contract
    ${TINKERHP_MAIN_PREPROCESS_DEFINITIONS}
    ${TINKERHP_OPENACC_PREPROCESS_DEFINITIONS}
  )
  target_link_libraries(tinkerhp-parallel-release-contract PRIVATE tinkerhp_core)
  add_test(
    NAME tinkerhp-parallel-release-contract
    COMMAND mpirun -n 2 "$<TARGET_FILE:tinkerhp-parallel-release-contract>"
  )
  tinkerhp_label_test(tinkerhp-parallel-release-contract scope:integration claim:functional method:conformance method:failure-robustness)
  set_tests_properties(tinkerhp-parallel-release-contract PROPERTIES
    ENVIRONMENT "TINKERHP_SOURCE_DIR=${CMAKE_SOURCE_DIR};TINKER_PARAMS=.;OMPI_MCA_oob_tcp_if_include=lo;OMPI_MCA_btl=self,vader"
    PROCESSORS 2
    RUN_SERIAL TRUE
  )
  add_executable(tinkerhp-mdevice-contract tests/test_mdevice.f90)
  set_source_files_properties(tests/test_mdevice.f90 PROPERTIES
    Fortran_PREPROCESS ON)
  tinkerhp_configure_fortran_target(tinkerhp-mdevice-contract
    "${TINKERHP_MAIN_MODULE_DIR}" "")
  target_link_libraries(tinkerhp-mdevice-contract PRIVATE tinkerhp_core)
  add_test(
    NAME tinkerhp-mdevice-contract
    COMMAND mpirun -n 1 "$<TARGET_FILE:tinkerhp-mdevice-contract>"
  )
  tinkerhp_label_test(tinkerhp-mdevice-contract scope:integration claim:functional method:conformance method:failure-robustness)
  set_tests_properties(tinkerhp-mdevice-contract PROPERTIES
    ENVIRONMENT "TINKERHP_SOURCE_DIR=${CMAKE_SOURCE_DIR};TINKER_PARAMS=.;OMPI_MCA_oob_tcp_if_include=lo;OMPI_MCA_btl=self,vader"
  )
  add_executable(tinkerhp-workspace-lifecycle-contract tests/test_workspace_lifecycle.f90)
  set_source_files_properties(tests/test_workspace_lifecycle.f90 PROPERTIES
    Fortran_PREPROCESS ON)
  tinkerhp_configure_fortran_target(tinkerhp-workspace-lifecycle-contract
    "${TINKERHP_MAIN_MODULE_DIR}" "")
  target_link_libraries(tinkerhp-workspace-lifecycle-contract PRIVATE tinkerhp_core)
  add_test(
    NAME tinkerhp-workspace-lifecycle-contract
    COMMAND mpirun -n 1 "$<TARGET_FILE:tinkerhp-workspace-lifecycle-contract>"
  )
  tinkerhp_label_test(tinkerhp-workspace-lifecycle-contract scope:integration claim:functional method:conformance)
  set_tests_properties(tinkerhp-workspace-lifecycle-contract PROPERTIES
    ENVIRONMENT "TINKERHP_SOURCE_DIR=${CMAKE_SOURCE_DIR};TINKER_PARAMS=.;OMPI_MCA_oob_tcp_if_include=lo;OMPI_MCA_btl=self,vader"
  )
  add_executable(tinkerhp-pme-lifecycle-contract tests/test_pme_lifecycle.f90)
  set_source_files_properties(tests/test_pme_lifecycle.f90 PROPERTIES
    Fortran_PREPROCESS ON)
  tinkerhp_configure_fortran_target(tinkerhp-pme-lifecycle-contract
    "${TINKERHP_MAIN_MODULE_DIR}" "")
  target_link_libraries(tinkerhp-pme-lifecycle-contract PRIVATE tinkerhp_core)
  add_test(
    NAME tinkerhp-pme-lifecycle-contract
    COMMAND mpirun -n 1 "$<TARGET_FILE:tinkerhp-pme-lifecycle-contract>"
  )
  tinkerhp_label_test(tinkerhp-pme-lifecycle-contract scope:integration claim:functional method:conformance)
  set_tests_properties(tinkerhp-pme-lifecycle-contract PROPERTIES
    ENVIRONMENT "TINKERHP_SOURCE_DIR=${CMAKE_SOURCE_DIR};TINKER_PARAMS=.;OMPI_MCA_oob_tcp_if_include=lo;OMPI_MCA_btl=self,vader"
  )
  add_executable(tinkerhp-structure-contract tests/test_structure.f90)
  set_source_files_properties(tests/test_structure.f90 PROPERTIES
    Fortran_PREPROCESS ON)
  tinkerhp_configure_fortran_target(tinkerhp-structure-contract
    "${TINKERHP_MAIN_MODULE_DIR}" "")
  target_link_libraries(tinkerhp-structure-contract PRIVATE tinkerhp_core)
  add_test(
    NAME tinkerhp-structure-contract
    COMMAND mpirun -n 1 "$<TARGET_FILE:tinkerhp-structure-contract>"
  )
  tinkerhp_label_test(tinkerhp-structure-contract scope:integration claim:functional method:conformance method:failure-robustness)
  set_tests_properties(tinkerhp-structure-contract PROPERTIES
    ENVIRONMENT "TINKERHP_SOURCE_DIR=${CMAKE_SOURCE_DIR};TINKER_PARAMS=.;OMPI_MCA_oob_tcp_if_include=lo;OMPI_MCA_btl=self,vader"
  )
  add_executable(tinkerhp-vcouple-contract
    tests/test_vcouple_contract.f90 source/MOD_mutant.f90)
  tinkerhp_configure_fortran_target(tinkerhp-vcouple-contract
    "${TINKERHP_MAIN_MODULE_DIR}" "")
  set_target_properties(tinkerhp-vcouple-contract PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod/vcouple-contract")
  add_test(
    NAME tinkerhp-vcouple-contract
    COMMAND "${CMAKE_COMMAND}"
      "-Dcontract=$<TARGET_FILE:tinkerhp-vcouple-contract>"
      -P "${CMAKE_SOURCE_DIR}/tests/run_vcouple_contract.cmake"
  )
  tinkerhp_label_test(tinkerhp-vcouple-contract scope:unit claim:functional method:conformance method:failure-robustness)
  function(tinkerhp_add_dhfr_analyze_baseline suffix precision)
    set(_tinkerhp_analyze_features)
    if(precision STREQUAL "MIXED")
      list(APPEND _tinkerhp_analyze_features --feature FPA)
    endif()
    add_test(
      NAME "tinkerhp-scientific-analyze-dhfr-${suffix}-sequential"
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/tests/runtime_refactor_oracle.py" run
        --profile analyze-baseline
        --executable "$<TARGET_FILE:analyze>"
        --build-report "${CMAKE_BINARY_DIR}/tinkerhp-build-report.txt"
        --precision "${precision}"
        ${_tinkerhp_analyze_features}
        --backend GPU
        --mpi-launcher mpirun
        --mpi-ranks 1
        --mpi-layout 1x1x1
        --fixture-dir "${CMAKE_SOURCE_DIR}/examples"
        --key dhfr2.key
        --primary-xyz dhfr2.xyz
        --replace-output-dir
        --output-dir "${CMAKE_BINARY_DIR}/scientific-test-contract/dhfr-analyze-${suffix}-sequential"
        -- dhfr2 -k dhfr2 E
    )
    set_tests_properties("tinkerhp-scientific-analyze-dhfr-${suffix}-sequential" PROPERTIES
      ENVIRONMENT "OMPI_MCA_oob_tcp_if_include=lo;OMPI_MCA_btl=self,vader"
      PROCESSORS 1
      RUN_SERIAL TRUE
      TIMEOUT 180
      FIXTURES_SETUP "dhfr_analyze_${suffix}_sequential"
    )
    tinkerhp_label_test("tinkerhp-scientific-analyze-dhfr-${suffix}-sequential"
      scope:integration claim:functional method:conformance)
  endfunction()
  if(TINKERHP_PRECISION STREQUAL "DOUBLE")
    tinkerhp_add_dhfr_analyze_baseline(double DOUBLE)
    add_test(
      NAME tinkerhp-scientific-dynamic-bench6-double
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/tests/runtime_refactor_oracle.py" run
        --profile dynamic
        --executable "$<TARGET_FILE:dynamic>"
        --build-report "${CMAKE_BINARY_DIR}/tinkerhp-build-report.txt"
        --precision DOUBLE
        --backend GPU
        --mpi-launcher mpirun
        --mpi-ranks 1
        --mpi-layout 1x1x1
        --fixture-dir "${CMAKE_SOURCE_DIR}/examples"
        --key bench6_r.key
        --primary-xyz bench6.xyz
        --restart bench6.dyn
        --restart-reference "${CMAKE_SOURCE_DIR}/examples/bench6.dyn"
        --format-reference "${CMAKE_SOURCE_DIR}/examples/bench6.dyn"
        --reference "${CMAKE_SOURCE_DIR}/examples/bench6_r.double.dyn.ref"
        --output-root "${CMAKE_BINARY_DIR}/scientific-test-contract"
        -- bench6 50 2 20 2 300 -k bench6_r
    )
    tinkerhp_label_test(tinkerhp-scientific-dynamic-bench6-double scope:integration claim:functional claim:numerical-regression method:conformance)
    set_tests_properties(tinkerhp-scientific-dynamic-bench6-double PROPERTIES
      ENVIRONMENT "OMPI_MCA_oob_tcp_if_include=lo;OMPI_MCA_btl=self,vader"
      PROCESSORS 1
      RUN_SERIAL TRUE
      TIMEOUT 180
    )
    add_test(
      NAME tinkerhp-scientific-dynamic-dhfr-b-double
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/tests/runtime_refactor_oracle.py" run
        --profile dynamic
        --executable "$<TARGET_FILE:dynamic>"
        --build-report "${CMAKE_BINARY_DIR}/tinkerhp-build-report.txt"
        --precision DOUBLE
        --backend GPU
        --mpi-launcher mpirun
        --mpi-ranks 1
        --mpi-layout 1x1x1
        --fixture-dir "${CMAKE_SOURCE_DIR}/examples"
        --key dhfr2dyn_B.key
        --primary-xyz dhfr2dyn.xyz
        --restart dhfr2dyn.dyn
        --restart-reference "${CMAKE_SOURCE_DIR}/examples/dhfr2dyn.dyn"
        --format-reference "${CMAKE_SOURCE_DIR}/examples/dhfr2dyn.dyn"
        --reference "${CMAKE_SOURCE_DIR}/examples/dhfr2dyn_B.double.dyn.ref"
        --energy-horizon 25
        --temperature-horizon 20
        --output-root "${CMAKE_BINARY_DIR}/scientific-test-contract"
        -- dhfr2dyn 25 10 200 4 300 1.0 -k dhfr2dyn_B
    )
    tinkerhp_label_test(tinkerhp-scientific-dynamic-dhfr-b-double scope:integration claim:functional claim:numerical-regression method:conformance)
    set_tests_properties(tinkerhp-scientific-dynamic-dhfr-b-double PROPERTIES
      ENVIRONMENT "OMPI_MCA_oob_tcp_if_include=lo;OMPI_MCA_btl=self,vader"
      PROCESSORS 1
      RUN_SERIAL TRUE
      TIMEOUT 180
    )
    add_test(
      NAME tinkerhp-scientific-dynamic-watersmall-amoebap-double
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/tests/runtime_refactor_oracle.py" run
        --profile dynamic
        --executable "$<TARGET_FILE:dynamic>"
        --build-report "${CMAKE_BINARY_DIR}/tinkerhp-build-report.txt"
        --precision DOUBLE
        --backend GPU
        --mpi-launcher mpirun
        --mpi-ranks 1
        --mpi-layout 1x1x1
        --fixture-dir "${CMAKE_SOURCE_DIR}/examples"
        --key watersmall_amoebap.key
        --primary-xyz watersmall.xyz
        --restart watersmall.dyn
        --restart-reference "${CMAKE_SOURCE_DIR}/examples/watersmall.dyn"
        --format-reference "${CMAKE_SOURCE_DIR}/examples/watersmall.dyn"
        --reference "${CMAKE_SOURCE_DIR}/examples/watersmall_amoebap.double.dyn.ref"
        --energy-horizon 50
        --temperature-horizon 20
        --output-root "${CMAKE_BINARY_DIR}/scientific-test-contract"
        -- watersmall 50 2 200 2 300 -k watersmall_amoebap.key
    )
    tinkerhp_label_test(tinkerhp-scientific-dynamic-watersmall-amoebap-double scope:integration claim:functional claim:numerical-regression method:conformance)
    set_tests_properties(tinkerhp-scientific-dynamic-watersmall-amoebap-double PROPERTIES
      ENVIRONMENT "OMPI_MCA_oob_tcp_if_include=lo;OMPI_MCA_btl=self,vader"
      PROCESSORS 1
      RUN_SERIAL TRUE
      TIMEOUT 180
    )
    add_test(
      NAME tinkerhp-scientific-minimize-iceii-double
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/tests/runtime_refactor_oracle.py" run
        --profile minimize
        --executable "$<TARGET_FILE:minimize>"
        --build-report "${CMAKE_BINARY_DIR}/tinkerhp-build-report.txt"
        --precision DOUBLE
        --backend GPU
        --mpi-launcher mpirun
        --mpi-ranks 1
        --mpi-layout 1x1x1
        --fixture-dir "${CMAKE_SOURCE_DIR}/examples"
        --key iceII-large-amoeba_b.key
        --primary-xyz iceII-large-amoeba.xyz
        --reference "${CMAKE_SOURCE_DIR}/examples/iceII-large-amoeba_b.double.min.ref"
        --output-root "${CMAKE_BINARY_DIR}/scientific-test-contract"
        -- iceII-large-amoeba.xyz 1.0 -k iceII-large-amoeba_b
    )
    tinkerhp_label_test(tinkerhp-scientific-minimize-iceii-double scope:integration claim:functional claim:numerical-regression method:conformance)
    set_tests_properties(tinkerhp-scientific-minimize-iceii-double PROPERTIES
      ENVIRONMENT "OMPI_MCA_oob_tcp_if_include=lo;OMPI_MCA_btl=self,vader"
      PROCESSORS 1
      RUN_SERIAL TRUE
      TIMEOUT 180
    )
    add_test(
      NAME tinkerhp-scientific-minimize-watersmall-double
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/tests/runtime_refactor_oracle.py" run
        --profile minimize-convergence
        --executable "$<TARGET_FILE:minimize>"
        --build-report "${CMAKE_BINARY_DIR}/tinkerhp-build-report.txt"
        --precision DOUBLE
        --backend GPU
        --mpi-launcher mpirun
        --mpi-ranks 1
        --mpi-layout 1x1x1
        --fixture-dir "${CMAKE_SOURCE_DIR}/examples"
        --key watersmall_amoebap.key
        --primary-xyz watersmall.xyz
        --replace-output-dir
        --output-dir "${CMAKE_BINARY_DIR}/scientific-test-contract/watersmall-minimize-double-sequential"
        -- watersmall 1.0 -k watersmall_amoebap.key
    )
    tinkerhp_label_test(tinkerhp-scientific-minimize-watersmall-double scope:integration claim:functional method:conformance)
    set_tests_properties(tinkerhp-scientific-minimize-watersmall-double PROPERTIES
      ENVIRONMENT "OMPI_MCA_oob_tcp_if_include=lo;OMPI_MCA_btl=self,vader"
      PROCESSORS 1
      RUN_SERIAL TRUE
      TIMEOUT 300
      FIXTURES_SETUP watersmall_minimize_double_sequential
    )
    set(_tinkerhp_bar_output_root
      "${CMAKE_BINARY_DIR}/bar-scientific-test-contract")
    set(_tinkerhp_bar_fixtures
      --fixture watersmall1.xyz
      --fixture watersmall1.key
      --fixture watersmall1.dcd
      --fixture watersmall2.xyz
      --fixture watersmall2.key
      --fixture watersmall2.dcd
      --fixture water.prm)
    add_test(
      NAME tinkerhp-scientific-bar-mode1-double
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/tests/runtime_refactor_oracle.py" run
        --profile bar
        --bar-mode mode1
        --executable "$<TARGET_FILE:bar>"
        --build-report "${CMAKE_BINARY_DIR}/tinkerhp-build-report.txt"
        --precision DOUBLE
        --backend GPU
        --mpi-launcher mpirun
        --mpi-ranks 1
        --mpi-layout 1x1x1
        --fixture-dir "${CMAKE_SOURCE_DIR}/examples/bar"
        --key watersmall1.key
        --primary-xyz watersmall1.xyz
        ${_tinkerhp_bar_fixtures}
        --replace-output-dir
        --output-dir "${_tinkerhp_bar_output_root}/serial-mode1"
    )
    tinkerhp_label_test(tinkerhp-scientific-bar-mode1-double scope:integration claim:functional method:conformance)
    set_tests_properties(tinkerhp-scientific-bar-mode1-double PROPERTIES
      ENVIRONMENT "OMPI_MCA_oob_tcp_if_include=lo;OMPI_MCA_btl=self,vader"
      PROCESSORS 1
      RUN_SERIAL TRUE
      TIMEOUT 900
    )
    add_test(
      NAME tinkerhp-scientific-bar-mode2-double
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/tests/runtime_refactor_oracle.py" run
        --profile bar
        --bar-mode mode2
        --executable "$<TARGET_FILE:bar>"
        --build-report "${CMAKE_BINARY_DIR}/tinkerhp-build-report.txt"
        --precision DOUBLE
        --backend GPU
        --mpi-launcher mpirun
        --mpi-ranks 1
        --mpi-layout 1x1x1
        --fixture-dir "${CMAKE_SOURCE_DIR}/examples/bar"
        --key watersmall1.key
        --primary-xyz watersmall1.xyz
        --fixture watersmall1.key
        --bar-input "${_tinkerhp_bar_output_root}/serial-mode1/workspace/watersmall1.bar"
        --reference "${CMAKE_SOURCE_DIR}/examples/bar/watersmall.outbar2"
        --replace-output-dir
        --output-dir "${_tinkerhp_bar_output_root}/serial-mode2"
    )
    tinkerhp_label_test(tinkerhp-scientific-bar-mode2-double scope:integration claim:functional claim:numerical-regression method:conformance)
    set_tests_properties(tinkerhp-scientific-bar-mode2-double PROPERTIES
      DEPENDS tinkerhp-scientific-bar-mode1-double
      ENVIRONMENT "OMPI_MCA_oob_tcp_if_include=lo;OMPI_MCA_btl=self,vader"
      PROCESSORS 1
      RUN_SERIAL TRUE
      TIMEOUT 300
    )
    add_test(
      NAME tinkerhp-scientific-bar-parallel-consume-double
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/tests/runtime_refactor_oracle.py" run
        --profile bar
        --bar-mode parallel-consume
        --executable "$<TARGET_FILE:bar>"
        --build-report "${CMAKE_BINARY_DIR}/tinkerhp-build-report.txt"
        --precision DOUBLE
        --backend GPU
        --mpi-launcher mpirun
        --mpi-ranks 2
        --mpi-layout 1x1x2
        --fixture-dir "${CMAKE_SOURCE_DIR}/examples/bar"
        --key watersmall1.key
        --primary-xyz watersmall1.xyz
        ${_tinkerhp_bar_fixtures}
        --reference "${CMAKE_SOURCE_DIR}/examples/bar/watersmall.outbar2"
        --bar-reference-file "${_tinkerhp_bar_output_root}/serial-mode1/workspace/watersmall1.bar"
        --replace-output-dir
        --output-dir "${_tinkerhp_bar_output_root}/parallel-consume"
    )
    tinkerhp_label_test(tinkerhp-scientific-bar-parallel-consume-double scope:integration claim:functional claim:numerical-regression method:conformance method:differential)
    set_tests_properties(tinkerhp-scientific-bar-parallel-consume-double PROPERTIES
      DEPENDS tinkerhp-scientific-bar-mode1-double
      ENVIRONMENT "OMPI_MCA_oob_tcp_if_include=lo;OMPI_MCA_btl=self,vader"
      PROCESSORS 2
      RUN_SERIAL TRUE
      TIMEOUT 900
    )
    set(_tinkerhp_pibar_output_root
      "${CMAKE_BINARY_DIR}/pibar-scientific-test-contract")
    set(_tinkerhp_pibar_fixtures
      --fixture watersmall1.xyz
      --fixture watersmall1.key
      --fixture watersmall1.dcd
      --fixture watersmall1_beads001.dcd
      --fixture watersmall2.xyz
      --fixture watersmall2.key
      --fixture watersmall2.dcd
      --fixture watersmall2_beads001.dcd
      --fixture water.prm)
    add_test(
      NAME tinkerhp-scientific-pibar-mode1-double
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/tests/runtime_refactor_oracle.py" run
        --profile pibar
        --bar-mode mode1
        --executable "$<TARGET_FILE:pibar>"
        --build-report "${CMAKE_BINARY_DIR}/tinkerhp-build-report.txt"
        --precision DOUBLE
        --backend GPU
        --mpi-launcher mpirun
        --mpi-ranks 1
        --mpi-layout 1x1x1
        --fixture-dir "${CMAKE_SOURCE_DIR}/examples/bar"
        --key watersmall1.key
        --primary-xyz watersmall1.xyz
        ${_tinkerhp_pibar_fixtures}
        --replace-output-dir
        --output-dir "${_tinkerhp_pibar_output_root}/mode1"
    )
    tinkerhp_label_test(tinkerhp-scientific-pibar-mode1-double scope:integration claim:functional method:conformance)
    set_tests_properties(tinkerhp-scientific-pibar-mode1-double PROPERTIES
      ENVIRONMENT "OMPI_MCA_oob_tcp_if_include=lo;OMPI_MCA_btl=self,vader"
      PROCESSORS 1
      RUN_SERIAL TRUE
      TIMEOUT 900
    )
    add_test(
      NAME tinkerhp-scientific-pibar-mode2-double
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/tests/runtime_refactor_oracle.py" run
        --profile pibar
        --bar-mode mode2
        --executable "$<TARGET_FILE:pibar>"
        --build-report "${CMAKE_BINARY_DIR}/tinkerhp-build-report.txt"
        --precision DOUBLE
        --backend GPU
        --mpi-launcher mpirun
        --mpi-ranks 1
        --mpi-layout 1x1x1
        --fixture-dir "${CMAKE_SOURCE_DIR}/examples/bar"
        --key watersmall1.key
        --primary-xyz watersmall1.xyz
        --fixture watersmall1.xyz
        --fixture watersmall1.key
        --bar-input "${_tinkerhp_pibar_output_root}/mode1/workspace/watersmall1.bar"
        --replace-output-dir
        --output-dir "${_tinkerhp_pibar_output_root}/mode2"
    )
    tinkerhp_label_test(tinkerhp-scientific-pibar-mode2-double scope:integration claim:functional method:conformance)
    set_tests_properties(tinkerhp-scientific-pibar-mode2-double PROPERTIES
      DEPENDS tinkerhp-scientific-pibar-mode1-double
      ENVIRONMENT "OMPI_MCA_oob_tcp_if_include=lo;OMPI_MCA_btl=self,vader"
      PROCESSORS 1
      RUN_SERIAL TRUE
      TIMEOUT 300
    )
    add_test(
      NAME tinkerhp-functional-pibar-invalid-mode-double
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/tests/runtime_refactor_oracle.py" run
        --profile pibar
        --bar-mode mode1
        --executable "$<TARGET_FILE:pibar>"
        --build-report "${CMAKE_BINARY_DIR}/tinkerhp-build-report.txt"
        --precision DOUBLE
        --backend GPU
        --mpi-launcher mpirun
        --mpi-ranks 1
        --mpi-layout 1x1x1
        --fixture-dir "${CMAKE_SOURCE_DIR}/examples/bar"
        --key watersmall1.key
        --primary-xyz watersmall1.xyz
        --fixture watersmall1.xyz
        --fixture watersmall1.key
        --expected-returncode -1
        --expected-diagnostic "PIBAR mode must be 1 or 2"
        --replace-output-dir
        --output-dir "${_tinkerhp_pibar_output_root}/invalid-mode"
        -- 3
    )
    tinkerhp_label_test(tinkerhp-functional-pibar-invalid-mode-double scope:integration claim:functional method:failure-robustness)
    set_tests_properties(tinkerhp-functional-pibar-invalid-mode-double PROPERTIES
      ENVIRONMENT "OMPI_MCA_oob_tcp_if_include=lo;OMPI_MCA_btl=self,vader"
      PROCESSORS 1
      RUN_SERIAL TRUE
      TIMEOUT 300
    )
    add_test(
      NAME tinkerhp-functional-pibar-missing-bar-double
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/tests/runtime_refactor_oracle.py" run
        --profile pibar
        --bar-mode mode2
        --executable "$<TARGET_FILE:pibar>"
        --build-report "${CMAKE_BINARY_DIR}/tinkerhp-build-report.txt"
        --precision DOUBLE
        --backend GPU
        --mpi-launcher mpirun
        --mpi-ranks 1
        --mpi-layout 1x1x1
        --fixture-dir "${CMAKE_SOURCE_DIR}/examples/bar"
        --key watersmall1.key
        --primary-xyz watersmall1.xyz
        --fixture watersmall1.xyz
        --fixture watersmall1.key
        --bar-input "${_tinkerhp_pibar_output_root}/mode1/workspace/watersmall1.bar"
        --expected-returncode -1
        --expected-diagnostic "BAR file not found"
        --replace-output-dir
        --output-dir "${_tinkerhp_pibar_output_root}/missing-bar"
        -- 2 missing false 1 100 1 1 100 1
    )
    tinkerhp_label_test(tinkerhp-functional-pibar-missing-bar-double scope:integration claim:functional method:failure-robustness)
    set_tests_properties(tinkerhp-functional-pibar-missing-bar-double PROPERTIES
      DEPENDS tinkerhp-scientific-pibar-mode1-double
      ENVIRONMENT "OMPI_MCA_oob_tcp_if_include=lo;OMPI_MCA_btl=self,vader"
      PROCESSORS 1
      RUN_SERIAL TRUE
      TIMEOUT 300
    )
    add_test(
      NAME tinkerhp-functional-pibar-malformed-bar-double
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/tests/runtime_refactor_oracle.py" run
        --profile pibar
        --bar-mode mode2
        --executable "$<TARGET_FILE:pibar>"
        --build-report "${CMAKE_BINARY_DIR}/tinkerhp-build-report.txt"
        --precision DOUBLE
        --backend GPU
        --mpi-launcher mpirun
        --mpi-ranks 1
        --mpi-layout 1x1x1
        --fixture-dir "${CMAKE_SOURCE_DIR}/examples/bar"
        --key watersmall1.key
        --primary-xyz watersmall1.xyz
        --fixture watersmall1.xyz
        --fixture watersmall1.key
        --bar-input "${CMAKE_SOURCE_DIR}/examples/bar/watersmall1_malformed.bar"
        --expected-returncode -1
        --expected-diagnostic "Invalid BAR file"
        --replace-output-dir
        --output-dir "${_tinkerhp_pibar_output_root}/malformed-bar"
        -- 2 watersmall1_malformed false 1 100 1 1 100 1
    )
    tinkerhp_label_test(tinkerhp-functional-pibar-malformed-bar-double scope:integration claim:functional method:failure-robustness)
    set_tests_properties(tinkerhp-functional-pibar-malformed-bar-double PROPERTIES
      ENVIRONMENT "OMPI_MCA_oob_tcp_if_include=lo;OMPI_MCA_btl=self,vader"
      PROCESSORS 1
      RUN_SERIAL TRUE
      TIMEOUT 300
    )
  endif()
  if(TINKERHP_PRECISION STREQUAL "MIXED" AND TINKERHP_ENABLE_FPA)
    tinkerhp_add_dhfr_analyze_baseline(fpa MIXED)
    add_test(
      NAME tinkerhp-scientific-minimize-watersmall-fpa
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/tests/runtime_refactor_oracle.py" run
        --profile minimize-convergence
        --executable "$<TARGET_FILE:minimize>"
        --build-report "${CMAKE_BINARY_DIR}/tinkerhp-build-report.txt"
        --precision MIXED
        --backend GPU
        --mpi-launcher mpirun
        --mpi-ranks 1
        --mpi-layout 1x1x1
        --feature FPA
        --fixture-dir "${CMAKE_SOURCE_DIR}/examples"
        --key watersmall_amoebap.key
        --primary-xyz watersmall.xyz
        --replace-output-dir
        --output-dir "${CMAKE_BINARY_DIR}/scientific-test-contract/watersmall-minimize-fpa-sequential"
        -- watersmall 1.0 -k watersmall_amoebap.key
    )
    tinkerhp_label_test(tinkerhp-scientific-minimize-watersmall-fpa scope:integration claim:functional method:conformance)
    set_tests_properties(tinkerhp-scientific-minimize-watersmall-fpa PROPERTIES
      ENVIRONMENT "OMPI_MCA_oob_tcp_if_include=lo;OMPI_MCA_btl=self,vader"
      PROCESSORS 1
      RUN_SERIAL TRUE
      TIMEOUT 300
      FIXTURES_SETUP watersmall_minimize_fpa_sequential
    )
    add_test(
      NAME tinkerhp-scientific-dynamic-bench6-fpa
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/tests/runtime_refactor_oracle.py" run
        --profile dynamic
        --executable "$<TARGET_FILE:dynamic>"
        --build-report "${CMAKE_BINARY_DIR}/tinkerhp-build-report.txt"
        --precision MIXED
        --backend GPU
        --mpi-launcher mpirun
        --mpi-ranks 1
        --mpi-layout 1x1x1
        --feature FPA
        --dynamic-energy-threshold 1e-4
        --fixture-dir "${CMAKE_SOURCE_DIR}/examples"
        --key bench6_r.key
        --primary-xyz bench6.xyz
        --restart bench6.dyn
        --restart-reference "${CMAKE_SOURCE_DIR}/examples/bench6.dyn"
        --format-reference "${CMAKE_SOURCE_DIR}/examples/bench6.dyn"
        --reference "${CMAKE_SOURCE_DIR}/examples/bench6_r.fpa.dyn.ref"
        --output-root "${CMAKE_BINARY_DIR}/scientific-test-contract"
        -- bench6 50 2 20 2 300 -k bench6_r
    )
    tinkerhp_label_test(tinkerhp-scientific-dynamic-bench6-fpa scope:integration claim:functional claim:numerical-regression method:conformance)
    set_tests_properties(tinkerhp-scientific-dynamic-bench6-fpa PROPERTIES
      ENVIRONMENT "OMPI_MCA_oob_tcp_if_include=lo;OMPI_MCA_btl=self,vader"
      PROCESSORS 1
      RUN_SERIAL TRUE
      TIMEOUT 180
    )
  endif()
  function(tinkerhp_add_pimd_watersmall_scientific_test
      name ranks precision output_root target_temperature)
    set(extra_arguments)
    if(precision STREQUAL "MIXED" AND TINKERHP_ENABLE_FPA)
      list(APPEND extra_arguments --feature FPA)
    endif()
    add_test(
      NAME "${name}"
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/tests/runtime_refactor_oracle.py" run
        --profile pimd
        --executable "$<TARGET_FILE:pimd>"
        --build-report "${CMAKE_BINARY_DIR}/tinkerhp-build-report.txt"
        --precision "${precision}"
        --backend GPU
        --mpi-launcher mpirun
        --mpi-ranks "${ranks}"
        --mpi-layout auto
        ${extra_arguments}
        --fixture-dir "${CMAKE_SOURCE_DIR}/examples/NQE/liquid_water/RPMD"
        --key watersmall.key
        --primary-xyz watersmall.xyz
        --reference "${CMAKE_SOURCE_DIR}/examples/NQE/liquid_water/RPMD/watersmall.double.pimd.ref"
        --centroid-temperature-target "${target_temperature}"
        --output-root "${output_root}"
        -- watersmall 1000 2 20 2 300
    )
    set_tests_properties("${name}" PROPERTIES
      ENVIRONMENT "OMPI_MCA_oob_tcp_if_include=lo;OMPI_MCA_btl=self,vader"
      PROCESSORS "${ranks}"
      RUN_SERIAL TRUE
      TIMEOUT 120
    )
    tinkerhp_label_test("${name}"
      scope:integration claim:functional claim:scientific-result method:conformance)
  endfunction()
  function(tinkerhp_add_pimd_layout_test name precision reference output_root)
    set(feature_arguments ${ARGN})
    add_test(
      NAME "${name}"
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/tests/runtime_refactor_oracle.py" run
        --profile pimd-layout
        --executable "$<TARGET_FILE:pimd>"
        --build-report "${CMAKE_BINARY_DIR}/tinkerhp-build-report.txt"
        --precision "${precision}"
        --backend GPU
        --mpi-launcher mpirun
        --mpi-ranks 4
        --mpi-layout 2-bead-groups/2-spatial-ranks
        --fixture-dir "${CMAKE_SOURCE_DIR}/examples/NQE/liquid_water/RPMD"
        --key watersmall_b2.key
        --primary-xyz watersmall.xyz
        --reference "${CMAKE_SOURCE_DIR}/examples/NQE/liquid_water/RPMD/${reference}"
        ${feature_arguments}
        --energy-horizon 50
        --output-root "${output_root}"
        -- watersmall 50 2 20 2 300 -k watersmall_b2
    )
    set_tests_properties("${name}" PROPERTIES
      ENVIRONMENT "OMPI_MCA_oob_tcp_if_include=lo;OMPI_MCA_btl=self,vader"
      PROCESSORS 4
      RUN_SERIAL TRUE
      TIMEOUT 180
    )
    tinkerhp_label_test("${name}"
      scope:integration claim:functional claim:numerical-regression
      method:conformance method:differential)
  endfunction()
  function(tinkerhp_add_analyze_beads_watersmall_functional_test name precision output_root)
    set(extra_arguments)
    if(precision STREQUAL "MIXED")
      list(APPEND extra_arguments --feature FPA)
    endif()
    add_test(
      NAME "${name}"
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/tests/runtime_refactor_oracle.py" run
        --profile analyze-beads
        --executable "$<TARGET_FILE:analyze_beads>"
        --build-report "${CMAKE_BINARY_DIR}/tinkerhp-build-report.txt"
        --precision "${precision}"
        --backend GPU
        --mpi-launcher mpirun
        --mpi-ranks 1
        --mpi-layout auto
        ${extra_arguments}
        --fixture-dir "${CMAKE_SOURCE_DIR}/examples/NQE/liquid_water/RPMD"
        --key watersmall.key
        --primary-xyz watersmall.xyz
        --bead-archive-source watersmall.xyz
        --bead-archive-count 32
        --output-root "${output_root}"
        -- watersmall E 300
    )
    set_tests_properties("${name}" PROPERTIES
      ENVIRONMENT "OMPI_MCA_oob_tcp_if_include=lo;OMPI_MCA_btl=self,vader"
      PROCESSORS 1
      RUN_SERIAL TRUE
      TIMEOUT 120
    )
    tinkerhp_label_test("${name}" scope:integration claim:functional method:conformance)
  endfunction()
  function(tinkerhp_add_mlpot_only_analyze_integration_test name precision ranks output_root fixture)
    set(feature_arguments)
    set(comparison_arguments)
    if(precision STREQUAL "MIXED" AND TINKERHP_ENABLE_FPA)
      list(APPEND feature_arguments --feature FPA)
    endif()
    if(TINKERHP_ENABLE_COLVARS)
      list(APPEND feature_arguments --feature COLVARS)
    endif()

    if(ranks EQUAL 1)
      set(profile analyze-baseline)
      set(mpi_layout 1x1x1)
      set(output_dir "${output_root}/one-rank")
      set(labels "scope:integration;claim:functional;method:conformance")
      set(fixture_property FIXTURES_SETUP "${fixture}")
    else()
      set(profile analyze)
      set(mpi_layout 1x1x2)
      set(output_dir "${output_root}/two-ranks")
      set(labels
        "scope:integration;claim:functional;claim:numerical-regression;method:differential")
      set(fixture_property FIXTURES_REQUIRED "${fixture}")
      set(comparison_arguments
        --sequential-report "${output_root}/one-rank/report.json"
        --expected-decomposition auto)
    endif()

    add_test(
      NAME "${name}"
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/tests/runtime_refactor_oracle.py" run
        --profile "${profile}"
        --executable "$<TARGET_FILE:analyze>"
        --build-report "${CMAKE_BINARY_DIR}/tinkerhp-build-report.txt"
        --precision "${precision}"
        ${feature_arguments}
        --backend GPU
        --mpi-launcher mpirun
        --mpi-ranks "${ranks}"
        --mpi-layout "${mpi_layout}"
        ${comparison_arguments}
        --fixture-dir "${CMAKE_SOURCE_DIR}/examples"
        --key watersmall_ml.key
        --primary-xyz watersmall.xyz
        --fixture watersmall.xyz
        --fixture watersmall_ml.key
        --model-file "${CMAKE_SOURCE_DIR}/ml_models/ani2x.fnx"
        --replace-output-dir
        --output-dir "${output_dir}"
        -- watersmall E -k watersmall_ml.key
    )
    set_tests_properties("${name}" PROPERTIES
      ENVIRONMENT "OMPI_MCA_oob_tcp_if_include=lo;OMPI_MCA_btl=self,vader"
      PROCESSORS "${ranks}"
      RUN_SERIAL TRUE
      TIMEOUT 180
      ${fixture_property}
    )
    tinkerhp_label_test("${name}" ${labels})
  endfunction()
  function(tinkerhp_add_pimd_piqtb_functional_test name output_root)
    add_test(
      NAME "${name}"
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/tests/runtime_refactor_oracle.py" run
        --profile pimd-piqtb
        --executable "$<TARGET_FILE:pimd>"
        --build-report "${CMAKE_BINARY_DIR}/tinkerhp-build-report.txt"
        --precision MIXED
        --backend GPU
        --mpi-launcher mpirun
        --mpi-ranks 4
        --mpi-layout auto
        --fixture-dir "${CMAKE_SOURCE_DIR}/examples/NQE/liquid_water/QTB"
        --key watersmall_pimd.key
        --primary-xyz watersmall.xyz
        --output-root "${output_root}"
        -- watersmall 1 2 2 2 300 -k watersmall_pimd
    )
    set_tests_properties("${name}" PROPERTIES
      ENVIRONMENT "OMPI_MCA_oob_tcp_if_include=lo;OMPI_MCA_btl=self,vader"
      PROCESSORS 4
      RUN_SERIAL TRUE
      TIMEOUT 120
    )
    tinkerhp_label_test("${name}" scope:integration claim:functional method:conformance)
  endfunction()
  function(tinkerhp_add_pimd_restart_test name ranks precision case output_root)
    set(extra_arguments)
    set(command_tail watersmall 1 2 2 2 300)
    if(TINKERHP_ENABLE_COLVARS)
      list(APPEND extra_arguments --feature COLVARS)
    endif()
    if(precision STREQUAL "MIXED" AND TINKERHP_ENABLE_FPA)
      list(APPEND extra_arguments --feature FPA)
    endif()
    if(case STREQUAL "resumed")
      set(command_tail watersmall 2 2 0.002 2 300)
      list(APPEND extra_arguments
        --expected-bead-loads 32
        --expected-bead-rewrites 32)
    elseif(case STREQUAL "missing")
      list(APPEND extra_arguments
        --bead-restart-missing 17
        --expected-bead-loads 0
        --expected-returncode -1
        --expected-diagnostic "not all beads can be restarted")
    elseif(case STREQUAL "malformed")
      list(APPEND extra_arguments
        --bead-restart-malformed 17
        --expected-bead-loads 17
        --expected-returncode -1
        --expected-diagnostic "Malformed DYN restart")
    else()
      list(APPEND extra_arguments --expected-bead-loads 32)
    endif()
    add_test(
      NAME "${name}"
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/tests/runtime_refactor_oracle.py" run
        --profile pimd-restart
        --executable "$<TARGET_FILE:pimd>"
        --build-report "${CMAKE_BINARY_DIR}/tinkerhp-build-report.txt"
        --precision "${precision}"
        --backend GPU
        --mpi-launcher mpirun
        --mpi-ranks "${ranks}"
        --mpi-layout auto
        ${extra_arguments}
        --fixture-dir "${CMAKE_SOURCE_DIR}/examples/NQE/liquid_water/RPMD"
        --key watersmall.key
        --primary-xyz watersmall.xyz
        --bead-restart-count 32
        --output-root "${output_root}"
        -- ${command_tail}
    )
    set_tests_properties("${name}" PROPERTIES
      ENVIRONMENT "OMPI_MCA_oob_tcp_if_include=lo;OMPI_MCA_btl=self,vader"
      PROCESSORS "${ranks}"
      RUN_SERIAL TRUE
      TIMEOUT 120
    )
    if(case STREQUAL "missing" OR case STREQUAL "malformed")
      tinkerhp_label_test("${name}" scope:integration claim:functional method:failure-robustness)
    else()
      tinkerhp_label_test("${name}" scope:integration claim:functional method:conformance)
    endif()
  endfunction()
  function(tinkerhp_add_pimd_dcd_contract name)
    add_test(
      NAME "${name}"
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/tests/test_pimd_dcd_contract.py"
        --source-dir "${CMAKE_SOURCE_DIR}"
        --executable "$<TARGET_FILE:pimd>"
        --mpi-launcher mpirun
        --mpi-ranks 2
    )
    set_tests_properties("${name}" PROPERTIES
      ENVIRONMENT "OMPI_MCA_oob_tcp_if_include=lo;OMPI_MCA_btl=self,vader"
      PROCESSORS 2
      RUN_SERIAL TRUE
      TIMEOUT 180
    )
    tinkerhp_label_test("${name}" scope:integration claim:functional method:conformance)
  endfunction()
  if(TINKERHP_PRECISION STREQUAL "DOUBLE")
    tinkerhp_add_pimd_dcd_contract(tinkerhp-pimd-dcd-contract)
    tinkerhp_add_pimd_watersmall_scientific_test(
      tinkerhp-scientific-pimd-watersmall-double 1 DOUBLE
      "${CMAKE_BINARY_DIR}/scientific-test-contract" 300)
    tinkerhp_add_analyze_beads_watersmall_functional_test(
      tinkerhp-functional-analyze-beads-watersmall-double DOUBLE
      "${CMAKE_BINARY_DIR}/functional-test-contract")
    tinkerhp_add_pimd_restart_test(
      tinkerhp-functional-pimd-restart-complete-double 1 DOUBLE complete
      "${CMAKE_BINARY_DIR}/functional-test-contract")
    tinkerhp_add_pimd_restart_test(
      tinkerhp-functional-pimd-restart-resumed-double 1 DOUBLE resumed
      "${CMAKE_BINARY_DIR}/functional-test-contract")
    tinkerhp_add_pimd_restart_test(
      tinkerhp-functional-pimd-restart-missing-double 1 DOUBLE missing
      "${CMAKE_BINARY_DIR}/functional-test-contract")
    tinkerhp_add_pimd_restart_test(
      tinkerhp-functional-pimd-restart-malformed-double 1 DOUBLE malformed
      "${CMAKE_BINARY_DIR}/functional-test-contract")
  elseif(TINKERHP_PRECISION STREQUAL "MIXED")
    tinkerhp_add_pimd_dcd_contract(tinkerhp-pimd-dcd-contract-mixed)
    tinkerhp_add_pimd_watersmall_scientific_test(
      tinkerhp-scientific-pimd-watersmall-mixed 1 MIXED
      "${CMAKE_BINARY_DIR}/scientific-test-contract" 300)
    if(TINKERHP_ENABLE_FPA)
      tinkerhp_add_analyze_beads_watersmall_functional_test(
        tinkerhp-functional-analyze-beads-watersmall-fpa MIXED
        "${CMAKE_BINARY_DIR}/functional-test-contract")
    endif()
    tinkerhp_add_pimd_restart_test(
      tinkerhp-functional-pimd-restart-complete-mixed 1 MIXED complete
      "${CMAKE_BINARY_DIR}/functional-test-contract")
    tinkerhp_add_pimd_restart_test(
      tinkerhp-functional-pimd-restart-resumed-mixed 1 MIXED resumed
      "${CMAKE_BINARY_DIR}/functional-test-contract")
    tinkerhp_add_pimd_restart_test(
      tinkerhp-functional-pimd-restart-missing-mixed 1 MIXED missing
      "${CMAKE_BINARY_DIR}/functional-test-contract")
    tinkerhp_add_pimd_restart_test(
      tinkerhp-functional-pimd-restart-malformed-mixed 1 MIXED malformed
      "${CMAKE_BINARY_DIR}/functional-test-contract")
  endif()
  if(TINKERHP_ENABLE_NN)
    if(TINKERHP_PRECISION STREQUAL "DOUBLE")
      set(_tinkerhp_mlpot_analyze_precision_tag double)
    elseif(TINKERHP_ENABLE_FPA)
      set(_tinkerhp_mlpot_analyze_precision_tag fpa)
    else()
      set(_tinkerhp_mlpot_analyze_precision_tag mixed)
    endif()
    set(_tinkerhp_mlpot_analyze_variant "nn-${_tinkerhp_mlpot_analyze_precision_tag}")
    if(TINKERHP_ENABLE_COLVARS)
      set(_tinkerhp_mlpot_analyze_variant
        "nn-colvars-${_tinkerhp_mlpot_analyze_precision_tag}")
    endif()
    set(_tinkerhp_mlpot_analyze_test_prefix
      "tinkerhp-functional-analyze-mlpot-only-watersmall-${_tinkerhp_mlpot_analyze_variant}")
    set(_tinkerhp_mlpot_analyze_output_root
      "${CMAKE_BINARY_DIR}/mlpot-only-analyze-${_tinkerhp_mlpot_analyze_variant}")
    set(_tinkerhp_mlpot_analyze_fixture
      "mlpot_only_analyze_${_tinkerhp_mlpot_analyze_variant}_one_rank")
    tinkerhp_add_mlpot_only_analyze_integration_test(
      "${_tinkerhp_mlpot_analyze_test_prefix}-one-rank"
      "${TINKERHP_PRECISION}" 1 "${_tinkerhp_mlpot_analyze_output_root}"
      "${_tinkerhp_mlpot_analyze_fixture}")
    tinkerhp_add_mlpot_only_analyze_integration_test(
      "${_tinkerhp_mlpot_analyze_test_prefix}-two-ranks"
      "${TINKERHP_PRECISION}" 2 "${_tinkerhp_mlpot_analyze_output_root}"
      "${_tinkerhp_mlpot_analyze_fixture}")
  endif()
  if(TINKERHP_ENABLE_PERFORMANCE_TESTS
      AND TINKERHP_PRECISION STREQUAL "MIXED"
      AND TINKERHP_ENABLE_FPA)
    add_test(
      NAME tinkerhp-performance-dynamic-dhfr-b-fpa
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/tests/runtime_refactor_oracle.py" run
        --profile dynamic-performance
        --executable "$<TARGET_FILE:dynamic>"
        --build-report "${CMAKE_BINARY_DIR}/tinkerhp-build-report.txt"
        --precision MIXED
        --backend GPU
        --mpi-launcher mpirun
        --mpi-ranks 1
        --mpi-layout 1x1x1
        --feature FPA
        --fixture-dir "${CMAKE_SOURCE_DIR}/examples"
        --key dhfr2dyn_B_perf.key
        --primary-xyz dhfr2dyn.xyz
        --restart dhfr2dyn.dyn
        --reference "${CMAKE_SOURCE_DIR}/examples/dhfr2dyn_B_perf.fpa.dyn.ref"
        --output-root "${CMAKE_BINARY_DIR}/performance-test-contract"
        -- dhfr2dyn 5000 10 200 2 300 -k dhfr2dyn_B_perf
    )
    set_tests_properties(tinkerhp-performance-dynamic-dhfr-b-fpa PROPERTIES
      ENVIRONMENT "OMPI_MCA_oob_tcp_if_include=lo;OMPI_MCA_btl=self,vader"
      PROCESSORS 1
      RUN_SERIAL TRUE
      TIMEOUT 900
    )
    tinkerhp_label_test(tinkerhp-performance-dynamic-dhfr-b-fpa
      scope:integration claim:functional claim:performance method:conformance
      performance)
  endif()
  if(TINKERHP_ENABLE_PARALLEL_SCIENTIFIC_TESTS)
    set(_tinkerhp_parallel_scientific_output_root
      "${CMAKE_BINARY_DIR}/parallel-scientific-test-contract")
    function(tinkerhp_add_parallel_dhfr_analyze_test suffix precision sequential_fixture)
      set(_tinkerhp_analyze_features)
      if(precision STREQUAL "MIXED")
        list(APPEND _tinkerhp_analyze_features --feature FPA)
      endif()
      add_test(
        NAME "tinkerhp-parallel-scientific-analyze-dhfr-${suffix}"
        COMMAND "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/tests/runtime_refactor_oracle.py" run
          --profile analyze
          --executable "$<TARGET_FILE:analyze>"
          --build-report "${CMAKE_BINARY_DIR}/tinkerhp-build-report.txt"
          --precision "${precision}"
          ${_tinkerhp_analyze_features}
          --backend GPU
          --mpi-launcher mpirun
          --mpi-ranks 2
          --mpi-layout 1x1x2
          --sequential-report "${CMAKE_BINARY_DIR}/scientific-test-contract/dhfr-analyze-${suffix}-sequential/report.json"
          --fixture-dir "${CMAKE_SOURCE_DIR}/examples"
          --key dhfr2.key
          --primary-xyz dhfr2.xyz
          --expected-decomposition auto
          --output-root "${_tinkerhp_parallel_scientific_output_root}"
          -- dhfr2 -k dhfr2 E
      )
      set_tests_properties("tinkerhp-parallel-scientific-analyze-dhfr-${suffix}" PROPERTIES
        ENVIRONMENT "OMPI_MCA_oob_tcp_if_include=lo;OMPI_MCA_btl=self,vader"
        PROCESSORS 2
        RUN_SERIAL TRUE
        TIMEOUT 180
        FIXTURES_REQUIRED "${sequential_fixture}"
      )
      tinkerhp_label_test("tinkerhp-parallel-scientific-analyze-dhfr-${suffix}"
        scope:integration claim:functional claim:numerical-regression
        method:differential)
    endfunction()
    function(tinkerhp_add_parallel_minimize_watersmall_test name precision sequential_dir fixture_name)
      set(_tinkerhp_minimize_features)
      if(precision STREQUAL "MIXED")
        list(APPEND _tinkerhp_minimize_features --feature FPA)
      endif()
      add_test(
        NAME "${name}"
        COMMAND "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/tests/runtime_refactor_oracle.py" run
          --profile minimize-convergence
          --executable "$<TARGET_FILE:minimize>"
          --build-report "${CMAKE_BINARY_DIR}/tinkerhp-build-report.txt"
          --precision "${precision}"
          ${_tinkerhp_minimize_features}
          --backend GPU
          --mpi-launcher mpirun
          --mpi-ranks 2
          --mpi-layout 1x1x2
          --sequential-report "${sequential_dir}/report.json"
          --fixture-dir "${CMAKE_SOURCE_DIR}/examples"
          --key watersmall_amoebap.key
          --primary-xyz watersmall.xyz
          --output-root "${_tinkerhp_parallel_scientific_output_root}"
          -- watersmall 1.0 -k watersmall_amoebap.key
      )
      set_tests_properties("${name}" PROPERTIES
        ENVIRONMENT "OMPI_MCA_oob_tcp_if_include=lo;OMPI_MCA_btl=self,vader"
        PROCESSORS 2
        RUN_SERIAL TRUE
        TIMEOUT 300
        FIXTURES_REQUIRED "${fixture_name}"
      )
      tinkerhp_label_test("${name}"
        scope:integration claim:functional claim:numerical-regression
        method:differential)
    endfunction()
    function(tinkerhp_add_parallel_scientific_test name ranks precision reference key decomposition)
      set(_tinkerhp_parallel_scientific_features)
      if(precision STREQUAL "MIXED" AND TINKERHP_ENABLE_FPA)
        list(APPEND _tinkerhp_parallel_scientific_features --feature FPA)
      endif()
      add_test(
        NAME "${name}"
        COMMAND "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/tests/runtime_refactor_oracle.py" run
          --program dynamic
          --executable "$<TARGET_FILE:dynamic>"
          --build-report "${CMAKE_BINARY_DIR}/tinkerhp-build-report.txt"
          --precision "${precision}"
          ${_tinkerhp_parallel_scientific_features}
          --backend GPU
          --mpi-launcher mpirun
          --mpi-ranks "${ranks}"
          --mpi-layout "${decomposition}"
          --fixture-dir "${CMAKE_SOURCE_DIR}/examples"
          --key "${key}"
          --reference "${CMAKE_SOURCE_DIR}/examples/${reference}"
          --expected-decomposition "${decomposition}"
          --output-root "${_tinkerhp_parallel_scientific_output_root}"
      )
      set_tests_properties("${name}" PROPERTIES RUN_SERIAL TRUE TIMEOUT 180)
      tinkerhp_label_test("${name}"
        scope:integration claim:functional claim:numerical-regression
        method:conformance)
    endfunction()
    if(TINKERHP_PRECISION STREQUAL "DOUBLE")
      tinkerhp_add_parallel_dhfr_analyze_test(
        double DOUBLE dhfr_analyze_double_sequential)
      tinkerhp_add_parallel_minimize_watersmall_test(
        tinkerhp-parallel-scientific-minimize-watersmall-double DOUBLE
        "${CMAKE_BINARY_DIR}/scientific-test-contract/watersmall-minimize-double-sequential"
        watersmall_minimize_double_sequential)
      tinkerhp_add_pimd_layout_test(
        tinkerhp-parallel-scientific-pimd-double-layout-energy DOUBLE
        watersmall.double.2b.pimd.ref "${_tinkerhp_parallel_scientific_output_root}")
    endif()
    if(TINKERHP_PRECISION STREQUAL "MIXED")
      if(TINKERHP_ENABLE_FPA)
        tinkerhp_add_parallel_dhfr_analyze_test(
          fpa MIXED dhfr_analyze_fpa_sequential)
        tinkerhp_add_parallel_minimize_watersmall_test(
          tinkerhp-parallel-scientific-minimize-watersmall-fpa MIXED
          "${CMAKE_BINARY_DIR}/scientific-test-contract/watersmall-minimize-fpa-sequential"
          watersmall_minimize_fpa_sequential)
      endif()
      if(TINKERHP_ENABLE_FPA AND TINKERHP_ENABLE_COLVARS)
        tinkerhp_add_pimd_layout_test(
          tinkerhp-parallel-scientific-pimd-fpa-colvars-layout-energy MIXED
          watersmall.fpa.2b.pimd.ref "${_tinkerhp_parallel_scientific_output_root}"
          --feature FPA --feature COLVARS)
      endif()
      tinkerhp_add_pimd_restart_test(
        tinkerhp-parallel-functional-pimd-restart-complete-4 4 MIXED complete
        "${_tinkerhp_parallel_scientific_output_root}")
      tinkerhp_add_pimd_restart_test(
        tinkerhp-parallel-functional-pimd-restart-resumed-4 4 MIXED resumed
        "${_tinkerhp_parallel_scientific_output_root}")
      tinkerhp_add_pimd_restart_test(
        tinkerhp-parallel-functional-pimd-restart-missing-4 4 MIXED missing
        "${_tinkerhp_parallel_scientific_output_root}")
      tinkerhp_add_pimd_restart_test(
        tinkerhp-parallel-functional-pimd-restart-malformed-4 4 MIXED malformed
        "${_tinkerhp_parallel_scientific_output_root}")
      tinkerhp_add_pimd_watersmall_scientific_test(
        tinkerhp-parallel-scientific-pimd-watersmall-4 4 MIXED
        "${_tinkerhp_parallel_scientific_output_root}" 300)
      tinkerhp_add_pimd_piqtb_functional_test(
        tinkerhp-parallel-functional-pimd-piqtb-watersmall-4
        "${_tinkerhp_parallel_scientific_output_root}")
      if(TINKERHP_ENABLE_FPA)
        add_test(
          NAME tinkerhp-parallel-scientific-dynamic-bench6-fpa-2
          COMMAND "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/tests/runtime_refactor_oracle.py" run
            --profile dynamic
            --executable "$<TARGET_FILE:dynamic>"
            --build-report "${CMAKE_BINARY_DIR}/tinkerhp-build-report.txt"
            --precision MIXED
            --backend GPU
            --mpi-launcher mpirun
            --mpi-ranks 2
            --mpi-layout auto
            --feature FPA
            --dynamic-energy-threshold 1e-4
            --dynamic-parallel-initial-threshold 8.1e-7
            --dynamic-parallel-tail-threshold 4.3e-6
            --dynamic-energy-acceptance-ratio 0.43
            --dynamic-parallel-start-threshold 5e-8
            --fixture-dir "${CMAKE_SOURCE_DIR}/examples"
            --key bench6_r.key
            --primary-xyz bench6.xyz
            --restart bench6.dyn
            --restart-reference "${CMAKE_SOURCE_DIR}/examples/bench6.dyn"
            --format-reference "${CMAKE_SOURCE_DIR}/examples/bench6.dyn"
            --reference "${CMAKE_SOURCE_DIR}/examples/bench6_r.fpa.dyn.ref"
            --expected-decomposition auto
            --dynamic-temperature-target 300
            --output-root "${_tinkerhp_parallel_scientific_output_root}"
            -- bench6 50 2 20 2 300 -k bench6_r
        )
        tinkerhp_label_test(tinkerhp-parallel-scientific-dynamic-bench6-fpa-2 scope:integration claim:functional claim:numerical-regression method:conformance)
        set_tests_properties(tinkerhp-parallel-scientific-dynamic-bench6-fpa-2
          PROPERTIES RUN_SERIAL TRUE TIMEOUT 180)
      endif()
      tinkerhp_add_parallel_scientific_test(
        tinkerhp-parallel-scientific-dynamic-mixed-auto-2 2 MIXED
        dhfr2dyn_r.dyn.ref dhfr2dyn_r.key auto)
      tinkerhp_add_parallel_scientific_test(
        tinkerhp-parallel-scientific-dynamic-mixed-1d-4 4 MIXED
        dhfr2dyn_r.dyn.ref dhfr2dyn_r_1d.key explicit-1d)
    elseif(TINKERHP_PRECISION STREQUAL "DOUBLE")
      tinkerhp_add_parallel_scientific_test(
        tinkerhp-parallel-scientific-dynamic-double-auto-4 4 DOUBLE
        dhfr2dyn_r.double.dyn.ref dhfr2dyn_r.key auto)
    endif()
  endif()
  add_test(
    NAME tinkerhp-acceptance
    COMMAND "${CMAKE_COMMAND}"
      "-Ddynamic=$<TARGET_FILE:dynamic>"
      "-Dsource_dir=${CMAKE_SOURCE_DIR}"
      "-Droot_dir=${CMAKE_BINARY_DIR}/acceptance-test"
      -P "${CMAKE_SOURCE_DIR}/tests/run_acceptance.cmake"
  )
  tinkerhp_label_test(tinkerhp-acceptance scope:integration claim:functional method:conformance)
  set_tests_properties(tinkerhp-acceptance PROPERTIES
    ENVIRONMENT "OMPI_MCA_oob_tcp_if_include=lo;OMPI_MCA_btl=self,vader"
    RUN_SERIAL TRUE
    TIMEOUT 300
  )
endif()
