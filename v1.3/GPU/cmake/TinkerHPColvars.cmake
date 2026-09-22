set(COLVARS_HOME "${CMAKE_SOURCE_DIR}/colvars" CACHE PATH
  "Colvars source/build directory used when TINKERHP_ENABLE_COLVARS is enabled."
)
set(TCL_HOME "" CACHE PATH
  "Optional TCL installation prefix used when TINKERHP_ENABLE_COLVARS is enabled."
)
set(COLVARS_MAKE "make" CACHE STRING
  "Make command used to build the external Colvars library."
)
set(TINKERHP_MPIRUN "" CACHE FILEPATH
  "Optional MPI launcher used to detect the OpenMPI C++ compatibility library."
)

function(tinkerhp_detect_colvars_mpi_cxx_library)
  set(_tinkerhp_mpirun "${TINKERHP_MPIRUN}")
  if(_tinkerhp_mpirun STREQUAL "")
    get_filename_component(_tinkerhp_mpi_wrapper_dir "${CMAKE_Fortran_COMPILER}" DIRECTORY)
    unset(_tinkerhp_mpirun)
    unset(_tinkerhp_mpirun CACHE)
    find_program(_tinkerhp_mpirun
      NAMES mpirun
      HINTS "${_tinkerhp_mpi_wrapper_dir}"
      NO_DEFAULT_PATH
    )
  endif()

  set(_tinkerhp_mpi_version "")
  set(_tinkerhp_mpi_cxx_library "")
  if(_tinkerhp_mpirun)
    execute_process(
      COMMAND "${_tinkerhp_mpirun}" --version
      RESULT_VARIABLE _tinkerhp_mpirun_result
      OUTPUT_VARIABLE _tinkerhp_mpirun_stdout
      ERROR_VARIABLE _tinkerhp_mpirun_stderr
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_STRIP_TRAILING_WHITESPACE
    )
    if(_tinkerhp_mpirun_result EQUAL 0)
      set(_tinkerhp_mpirun_output "${_tinkerhp_mpirun_stdout}\n${_tinkerhp_mpirun_stderr}")
      string(REGEX MATCH "mpirun \\(Open MPI\\) ([0-9]+\\.[0-9]+\\.[0-9]+)" _tinkerhp_openmpi_match "${_tinkerhp_mpirun_output}")
      if(_tinkerhp_openmpi_match)
        set(_tinkerhp_mpi_version "${CMAKE_MATCH_1}")
        if(_tinkerhp_mpi_version VERSION_GREATER_EQUAL "3.0" AND _tinkerhp_mpi_version VERSION_LESS "4.0")
          get_filename_component(_tinkerhp_mpirun_realpath "${_tinkerhp_mpirun}" REALPATH)
          get_filename_component(_tinkerhp_mpi_bin_dir "${_tinkerhp_mpirun_realpath}" DIRECTORY)
          get_filename_component(_tinkerhp_mpi_prefix "${_tinkerhp_mpi_bin_dir}" DIRECTORY)
          unset(_tinkerhp_mpi_cxx_library)
          unset(_tinkerhp_mpi_cxx_library CACHE)
          find_library(_tinkerhp_mpi_cxx_library
            NAMES mpi_cxx
            PATHS
              "${_tinkerhp_mpi_prefix}/lib"
              "${_tinkerhp_mpi_prefix}/lib64"
            NO_DEFAULT_PATH
          )
        endif()
      endif()
    endif()
  endif()

  set(TINKERHP_COLVARS_MPI_VERSION "${_tinkerhp_mpi_version}" CACHE STRING
    "OpenMPI version detected for the Colvars link." FORCE
  )
  set(TINKERHP_COLVARS_MPI_CXX_LIBRARY "${_tinkerhp_mpi_cxx_library}" CACHE FILEPATH
    "OpenMPI C++ compatibility library used by the Colvars-enabled link." FORCE
  )
endfunction()

function(tinkerhp_configure_colvars)
  if(NOT TINKERHP_ENABLE_COLVARS)
    set(TINKERHP_COLVARS_TCL_INCLUDE_DIR "" CACHE PATH
      "TCL headers used by the Colvars-enabled Tinker-HP build." FORCE
    )
    set(TINKERHP_COLVARS_TCL_LIBRARY "" CACHE FILEPATH
      "TCL library used by the Colvars-enabled Tinker-HP link." FORCE
    )
    set(TINKERHP_COLVARS_MPI_VERSION "" CACHE STRING
      "OpenMPI version detected for the Colvars link." FORCE
    )
    set(TINKERHP_COLVARS_MPI_CXX_LIBRARY "" CACHE FILEPATH
      "OpenMPI C++ compatibility library used by the Colvars-enabled link." FORCE
    )
    return()
  endif()

  if(NOT EXISTS "${COLVARS_HOME}/Makefile")
    message(FATAL_ERROR
      "COLVARS_HOME=${COLVARS_HOME} does not contain the required Colvars Makefile."
    )
  endif()

  # Re-evaluate the library when a user switches between bundled and external Tcl.
  unset(_tinkerhp_tcl_library CACHE)
  if(TCL_HOME STREQUAL "")
    set(_tinkerhp_tcl_include_dir "${CMAKE_SOURCE_DIR}/tcl8.6.13/generic")
    find_library(_tinkerhp_tcl_library NAMES tcl8.6 tcl)
  else()
    set(_tinkerhp_tcl_include_dir "${TCL_HOME}/include")
    find_library(_tinkerhp_tcl_library
      NAMES tcl8.6 tcl
      PATHS
        "${TCL_HOME}/lib"
        "${TCL_HOME}/lib64"
      NO_DEFAULT_PATH
    )
  endif()

  if(NOT EXISTS "${_tinkerhp_tcl_include_dir}/tcl.h")
    message(FATAL_ERROR
      "TCL header not found: configure with -DTCL_HOME=/path/to/tcl or provide ${CMAKE_SOURCE_DIR}/tcl8.6.13/generic/tcl.h"
    )
  endif()

  if(NOT _tinkerhp_tcl_library)
    message(FATAL_ERROR
      "TCL NOT FIND: install libtcl8.6 or configure with -DTCL_HOME=/path/to/tcl"
    )
  endif()

  set(TINKERHP_COLVARS_TCL_LIBRARY "${_tinkerhp_tcl_library}" CACHE FILEPATH
    "TCL library used by the Colvars-enabled Tinker-HP link." FORCE
  )
  set(TINKERHP_COLVARS_TCL_INCLUDE_DIR "${_tinkerhp_tcl_include_dir}" CACHE PATH
    "TCL headers used by the Colvars-enabled Tinker-HP build." FORCE
  )
  tinkerhp_detect_colvars_mpi_cxx_library()

  set(_colvars_cxxflags ${ARGV})
  if(TINKERHP_LINK_MODE STREQUAL "SHARED")
    list(APPEND _colvars_cxxflags -fPIC)
  endif()
  string(JOIN " " _colvars_cxxflags ${_colvars_cxxflags})
  string(STRIP "${_colvars_cxxflags}" _colvars_cxxflags)
  set(_tinkerhp_colvars_make_arguments
    "TCL_SUPPORT=1"
    "CXX=${CMAKE_CXX_COMPILER}"
    "CXXFLAGS=${_colvars_cxxflags}"
  )
  if(NOT TCL_HOME STREQUAL "")
    list(APPEND _tinkerhp_colvars_make_arguments "TCL_HOME=${TCL_HOME}")
  endif()

  # Colvars is built in its source directory.  A static archive produced by a
  # prior non-PIC configuration is not suitable for a shared Tinker-HP core,
  # and make does not consider CXXFLAGS when deciding whether to rebuild it.
  # Reset those artifacts once for each shared CMake build tree before making
  # the archive with -fPIC.
  if(TINKERHP_LINK_MODE STREQUAL "SHARED")
    set(_tinkerhp_colvars_pic_reset_stamp
      "${CMAKE_BINARY_DIR}/CMakeFiles/tinkerhp-colvars-pic-reset.stamp"
    )
    add_custom_command(
      OUTPUT "${_tinkerhp_colvars_pic_reset_stamp}"
      COMMAND "${COLVARS_MAKE}" clean
      COMMAND "${CMAKE_COMMAND}" -E touch "${_tinkerhp_colvars_pic_reset_stamp}"
      WORKING_DIRECTORY "${COLVARS_HOME}"
      VERBATIM
      USES_TERMINAL
    )
    add_custom_target(tinkerhp_colvars_pic_reset
      DEPENDS "${_tinkerhp_colvars_pic_reset_stamp}"
    )
  endif()

  add_custom_target(tinkerhp_colvars
    BYPRODUCTS "${COLVARS_HOME}/libcolvars.a"
    COMMAND "${CMAKE_COMMAND}" -E make_directory "${COLVARS_HOME}/.deps"
    COMMAND "${COLVARS_MAKE}"
            ${_tinkerhp_colvars_make_arguments}
            libcolvars.a
    WORKING_DIRECTORY "${COLVARS_HOME}"
    VERBATIM
    USES_TERMINAL
  )
  if(TARGET tinkerhp_colvars_pic_reset)
    add_dependencies(tinkerhp_colvars tinkerhp_colvars_pic_reset)
  endif()

  add_custom_target(tinkerhp_colvars_clean
    COMMAND "${COLVARS_MAKE}" clean
    WORKING_DIRECTORY "${COLVARS_HOME}"
    VERBATIM
    USES_TERMINAL
  )

  add_library(tinkerhp_colvars_link INTERFACE)
  add_dependencies(tinkerhp_colvars_link tinkerhp_colvars)
  set(_tinkerhp_colvars_link_libraries
    "${COLVARS_HOME}/libcolvars.a"
    "${TINKERHP_COLVARS_TCL_LIBRARY}"
  )
  if(TINKERHP_COLVARS_MPI_CXX_LIBRARY)
    list(APPEND _tinkerhp_colvars_link_libraries "${TINKERHP_COLVARS_MPI_CXX_LIBRARY}")
  endif()
  target_link_libraries(tinkerhp_colvars_link INTERFACE
    ${_tinkerhp_colvars_link_libraries}
  )
endfunction()
