add_library(tinkerhp_2decomp_fft STATIC
  src/decomp_2d.f90
  src/io.f90
  src/glassman.f90
  src/fft_generic.f90
  src/fft_cufft.f90
)

set_source_files_properties(
  src/decomp_2d.f90
  src/io.f90
  src/glassman.f90
  src/fft_generic.f90
  src/fft_cufft.f90
  PROPERTIES Fortran_PREPROCESS ON
)

set_target_properties(tinkerhp_2decomp_fft PROPERTIES
  ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
  Fortran_MODULE_DIRECTORY "${TINKERHP_2DECOMP_MODULE_DIR}"
)

target_compile_definitions(tinkerhp_2decomp_fft PRIVATE
  ${TINKERHP_2DECOMP_PREPROCESS_DEFINITIONS}
)

target_compile_options(tinkerhp_2decomp_fft PRIVATE
  $<$<COMPILE_LANGUAGE:Fortran>:${TINKERHP_COMMON_FORTRAN_FLAGS}>
)

target_include_directories(tinkerhp_2decomp_fft PUBLIC
  "${TINKERHP_2DECOMP_MODULE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/include"
  "${CMAKE_CURRENT_SOURCE_DIR}/src"
)
