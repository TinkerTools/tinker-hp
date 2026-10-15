set(PLUMED_SOURCE_HOME "${CMAKE_SOURCE_DIR}/plumed2" CACHE PATH
  "PLUMED source/build directory containing configure and its Makefile."
)
set(PLUMED_HOME "${PLUMED_SOURCE_HOME}/plumed/Portland" CACHE PATH
  "PLUMED installation prefix used when TINKERHP_ENABLE_PLUMED is enabled."
)
set(PLUMED_MAKE "make" CACHE STRING
  "Make command used to build the external PLUMED installation."
)

function(tinkerhp_configure_plumed)
  if(NOT TINKERHP_ENABLE_PLUMED)
    set(TINKERHP_PLUMED_LIBRARY "" CACHE FILEPATH
      "Installed PLUMED library used by the Tinker-HP link." FORCE
    )
    set(TINKERHP_PLUMED_KERNEL_LIBRARY "" CACHE FILEPATH
      "Installed PLUMED kernel library used by the Tinker-HP link." FORCE
    )
    return()
  endif()

  if(NOT EXISTS "${PLUMED_SOURCE_HOME}/configure")
    message(FATAL_ERROR
      "PLUMED_SOURCE_HOME=${PLUMED_SOURCE_HOME} does not contain the required configure script."
    )
  endif()

  # PLUMED's configure writes these files into the source tree. A validate-only
  # configure restores them afterwards, so the validator leaves no trace there.
  # A real build keeps them, because the PLUMED build step needs config.status.
  set(_tinkerhp_plumed_written_files
    Makefile.conf sourceme.sh config.status config.log stamp-h)
  set(_tinkerhp_plumed_snapshot "${CMAKE_BINARY_DIR}/plumed-configure-snapshot")
  if(TINKERHP_INTERNAL_VALIDATE_ONLY)
    file(REMOVE_RECURSE "${_tinkerhp_plumed_snapshot}")
    file(MAKE_DIRECTORY "${_tinkerhp_plumed_snapshot}")
    foreach(_tinkerhp_plumed_file IN LISTS _tinkerhp_plumed_written_files)
      if(EXISTS "${PLUMED_SOURCE_HOME}/${_tinkerhp_plumed_file}")
        file(COPY "${PLUMED_SOURCE_HOME}/${_tinkerhp_plumed_file}"
          DESTINATION "${_tinkerhp_plumed_snapshot}")
      endif()
    endforeach()
  endif()

  message(STATUS "Configuring external PLUMED in ${PLUMED_SOURCE_HOME}")
  execute_process(
    COMMAND "${PLUMED_SOURCE_HOME}/configure"
            "--prefix=${PLUMED_HOME}"
            "FC=${CMAKE_Fortran_COMPILER}"
            "CXX=${CMAKE_CXX_COMPILER}"
            "CC=${CMAKE_C_COMPILER}"
    WORKING_DIRECTORY "${PLUMED_SOURCE_HOME}"
    RESULT_VARIABLE _tinkerhp_plumed_configure_result
    OUTPUT_VARIABLE _tinkerhp_plumed_configure_stdout
    ERROR_VARIABLE _tinkerhp_plumed_configure_stderr
  )
  if(TINKERHP_INTERNAL_VALIDATE_ONLY)
    foreach(_tinkerhp_plumed_file IN LISTS _tinkerhp_plumed_written_files)
      if(EXISTS "${_tinkerhp_plumed_snapshot}/${_tinkerhp_plumed_file}")
        # configure_file copies by content; file(COPY) can skip a file that is
        # newer at the destination, which the configure has just rewritten.
        configure_file("${_tinkerhp_plumed_snapshot}/${_tinkerhp_plumed_file}"
          "${PLUMED_SOURCE_HOME}/${_tinkerhp_plumed_file}" COPYONLY)
      elseif(EXISTS "${PLUMED_SOURCE_HOME}/${_tinkerhp_plumed_file}")
        file(REMOVE "${PLUMED_SOURCE_HOME}/${_tinkerhp_plumed_file}")
      endif()
    endforeach()
  endif()
  if(NOT _tinkerhp_plumed_configure_result EQUAL 0)
    message(FATAL_ERROR
      "PLUMED configure failed in ${PLUMED_SOURCE_HOME} (exit ${_tinkerhp_plumed_configure_result}).\n"
      "stdout:\n${_tinkerhp_plumed_configure_stdout}\n"
      "stderr:\n${_tinkerhp_plumed_configure_stderr}"
    )
  endif()

  set(_tinkerhp_plumed_library "${PLUMED_HOME}/lib/libplumed.so")
  set(_tinkerhp_plumed_kernel_library "${PLUMED_HOME}/lib/libplumedKernel.so")
  set(TINKERHP_PLUMED_LIBRARY "${_tinkerhp_plumed_library}" CACHE FILEPATH
    "Installed PLUMED library used by the Tinker-HP link." FORCE
  )
  set(TINKERHP_PLUMED_KERNEL_LIBRARY "${_tinkerhp_plumed_kernel_library}" CACHE FILEPATH
    "Installed PLUMED kernel library used by the Tinker-HP link." FORCE
  )

  add_custom_target(tinkerhp_plumed
    BYPRODUCTS
      "${_tinkerhp_plumed_library}"
      "${_tinkerhp_plumed_kernel_library}"
    COMMAND "${PLUMED_MAKE}" install VERBOSE=1
    WORKING_DIRECTORY "${PLUMED_SOURCE_HOME}"
    VERBATIM
    USES_TERMINAL
  )

  add_custom_target(tinkerhp_plumed_clean
    COMMAND "${PLUMED_MAKE}" clean
    WORKING_DIRECTORY "${PLUMED_SOURCE_HOME}"
    VERBATIM
    USES_TERMINAL
  )

  add_library(tinkerhp_plumed_link INTERFACE)
  add_dependencies(tinkerhp_plumed_link tinkerhp_plumed)
  target_link_libraries(tinkerhp_plumed_link INTERFACE
    "${_tinkerhp_plumed_library}"
    "${_tinkerhp_plumed_kernel_library}"
  )
endfunction()
