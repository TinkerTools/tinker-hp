include_guard(GLOBAL)

set(DEEPHP_PYTHON_INTERPRETER "" CACHE FILEPATH
  "Python interpreter for a Deep-HP development environment.")
set(DEEPHP_UV_EXECUTABLE "" CACHE FILEPATH
  "uv executable used to create the managed Deep-HP environment.")
set(DEEPHP_PROJECT_DIR "${CMAKE_SOURCE_DIR}/deep-hp" CACHE PATH
  "Deep-HP managed-environment project directory.")
mark_as_advanced(DEEPHP_UV_EXECUTABLE DEEPHP_PROJECT_DIR)

function(tinkerhp_configure_deep_hp)
  set(TINKERHP_DEEPHP_ENVIRONMENT_ORIGIN "disabled" PARENT_SCOPE)
  set(TINKERHP_DEEPHP_PYTHON "" PARENT_SCOPE)
  set(TINKERHP_DEEPHP_PYTHON_VERSION "" PARENT_SCOPE)
  set(TINKERHP_DEEPHP_PYTHON_LIBRARY "" PARENT_SCOPE)
  set(TINKERHP_DEEPHP_PYTHON_LIBRARY_DIR "" PARENT_SCOPE)
  set(TINKERHP_DEEPHP_CUDA_PROFILE "" PARENT_SCOPE)

  if(NOT TINKERHP_ENABLE_NN)
    return()
  endif()

  if(NOT EXISTS "${DEEPHP_PROJECT_DIR}/pyproject.toml")
    message(FATAL_ERROR "Deep-HP project directory does not contain pyproject.toml: ${DEEPHP_PROJECT_DIR}")
  endif()

  set(_managed_python "${DEEPHP_PROJECT_DIR}/.venv/bin/python")
  if(EXISTS "${_managed_python}")
    set(_python "${_managed_python}")
    set(_origin "managed-local")
    set(_minimum_version "3.11")
    set(_exact_managed_version TRUE)
  elseif(DEEPHP_PYTHON_INTERPRETER)
    if(NOT EXISTS "${DEEPHP_PYTHON_INTERPRETER}")
      message(FATAL_ERROR "DEEPHP_PYTHON_INTERPRETER does not exist: ${DEEPHP_PYTHON_INTERPRETER}")
    endif()
    set(_python "${DEEPHP_PYTHON_INTERPRETER}")
    set(_origin "explicit")
    set(_minimum_version "3.10")
    set(_exact_managed_version FALSE)
  else()
    if(NOT TINKERHP_DETECTED_CUDA_VERSION MATCHES "^([0-9]+)\\.")
      message(FATAL_ERROR "Cannot select a Deep-HP CUDA profile from '${TINKERHP_DETECTED_CUDA_VERSION}'.")
    endif()
    set(_cuda_major "${CMAKE_MATCH_1}")
    if(NOT _cuda_major STREQUAL "12" AND NOT _cuda_major STREQUAL "13")
      message(FATAL_ERROR
        "Deep-HP managed environment supports NVHPC CUDA major versions 12 and 13; detected ${TINKERHP_DETECTED_CUDA_VERSION}."
      )
    endif()
    if(DEEPHP_UV_EXECUTABLE)
      set(_uv "${DEEPHP_UV_EXECUTABLE}")
    else()
      find_program(_uv NAMES uv REQUIRED)
    endif()
    execute_process(
      COMMAND "${_uv}" sync --locked --extra "cuda${_cuda_major}"
      WORKING_DIRECTORY "${DEEPHP_PROJECT_DIR}"
      RESULT_VARIABLE _uv_result
    )
    if(NOT _uv_result EQUAL 0)
      message(FATAL_ERROR "Failed to create the managed Deep-HP environment with uv sync --locked --extra cuda${_cuda_major}.")
    endif()
    if(NOT EXISTS "${_managed_python}")
      message(FATAL_ERROR "uv completed without creating ${_managed_python}.")
    endif()
    set(_python "${_managed_python}")
    set(_origin "managed-created")
    set(_minimum_version "3.11")
    set(_exact_managed_version TRUE)
    set(_profile "cuda${_cuda_major}")
  endif()

  execute_process(
    COMMAND "${_python}" -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}')"
    RESULT_VARIABLE _version_result
    OUTPUT_VARIABLE _version
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  if(NOT _version_result EQUAL 0 OR NOT _version MATCHES "^[0-9]+\\.[0-9]+\\.[0-9]+$")
    message(FATAL_ERROR "Unable to determine the Python version of selected Deep-HP interpreter: ${_python}")
  endif()
  if(_exact_managed_version AND NOT _version MATCHES "^3\\.11\\.")
    message(FATAL_ERROR "Selected managed Deep-HP interpreter must use Python 3.11; ${_python} reports ${_version}.")
  endif()
  if(NOT _exact_managed_version AND _version VERSION_LESS _minimum_version)
    message(FATAL_ERROR "DEEPHP_PYTHON_INTERPRETER must use Python 3.10 or newer; ${_python} reports ${_version}.")
  endif()

  execute_process(
    COMMAND "${_python}" -c "import cffi, numpy, yaml, cupy, jax, fennol"
    RESULT_VARIABLE _modules_result
    ERROR_VARIABLE _modules_error
  )
  if(NOT _modules_result EQUAL 0)
    message(FATAL_ERROR "Selected Deep-HP interpreter is missing a required module (cffi, numpy, yaml, cupy, jax, fennol): ${_python}\n${_modules_error}")
  endif()

  # The generated CFFI plugin has a DT_NEEDED entry for this interpreter's
  # libpython.  Expose that dependency to CMake so the final Fortran link can
  # resolve the plugin's Python symbols instead of relying on an unfindable
  # indirect dependency.
  execute_process(
    COMMAND "${_python}" -c "import sysconfig; print(sysconfig.get_config_var('LIBDIR') or ''); print(sysconfig.get_config_var('LDLIBRARY') or sysconfig.get_config_var('LIBRARY') or '')"
    RESULT_VARIABLE _python_library_result
    OUTPUT_VARIABLE _python_library_details
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  if(_python_library_result EQUAL 0 AND NOT _python_library_details STREQUAL "")
    string(REPLACE "\n" ";" _python_library_details "${_python_library_details}")
    list(LENGTH _python_library_details _python_library_detail_count)
    if(_python_library_detail_count GREATER 1)
      list(GET _python_library_details 0 _python_library_dir)
      list(GET _python_library_details 1 _python_library_name)
      if(NOT _python_library_dir STREQUAL "" AND NOT _python_library_name STREQUAL "")
        set(_python_library "${_python_library_dir}/${_python_library_name}")
      endif()
    endif()
  endif()
  if(NOT _python_library OR NOT EXISTS "${_python_library}" OR _python_library MATCHES "\\.(a|lib)$")
    message(FATAL_ERROR
      "Selected Deep-HP interpreter does not provide a usable shared libpython: ${_python}\n"
      "sysconfig reported library '${_python_library}'. Set DEEPHP_PYTHON_INTERPRETER to a Python installation with a shared native library."
    )
  endif()

  if(NOT DEFINED _profile)
    set(_profile "external")
  endif()
  set(TINKERHP_DEEPHP_ENVIRONMENT_ORIGIN "${_origin}" PARENT_SCOPE)
  set(TINKERHP_DEEPHP_PYTHON "${_python}" PARENT_SCOPE)
  set(TINKERHP_DEEPHP_PYTHON_VERSION "${_version}" PARENT_SCOPE)
  set(TINKERHP_DEEPHP_PYTHON_LIBRARY "${_python_library}" PARENT_SCOPE)
  set(TINKERHP_DEEPHP_PYTHON_LIBRARY_DIR "${_python_library_dir}" PARENT_SCOPE)
  set(TINKERHP_DEEPHP_CUDA_PROFILE "${_profile}" PARENT_SCOPE)
endfunction()

function(tinkerhp_add_deep_hp_plugin)
  if(NOT TINKERHP_ENABLE_NN)
    return()
  endif()

  set(_plugin_dir "${CMAKE_BINARY_DIR}/deephp")
  set(_plugin "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/${CMAKE_SHARED_LIBRARY_PREFIX}mlinterface${CMAKE_SHARED_LIBRARY_SUFFIX}")
  set(_module "${_plugin_dir}/mlinterface.py")
  set(_header "${_plugin_dir}/mlinterface.h")
  set(_builder "${_plugin_dir}/mlbuilder.py")
  add_custom_command(
    OUTPUT "${_plugin}"
    COMMAND "${CMAKE_COMMAND}" -E make_directory "${_plugin_dir}" "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}"
    COMMAND "${CMAKE_COMMAND}" -E copy_if_different "${CMAKE_SOURCE_DIR}/source/mlinterface.py" "${_module}"
    COMMAND "${CMAKE_COMMAND}" -E copy_if_different "${CMAKE_SOURCE_DIR}/source/mlinterface.h" "${_header}"
    COMMAND "${CMAKE_COMMAND}" -E copy_if_different "${CMAKE_SOURCE_DIR}/source/mlbuilder.py" "${_builder}"
    COMMAND "${TINKERHP_DEEPHP_PYTHON}" "${_builder}" --module "${_module}" --header "${_header}" --output "${_plugin}" --python-library-dir "${TINKERHP_DEEPHP_PYTHON_LIBRARY_DIR}"
    DEPENDS "${CMAKE_SOURCE_DIR}/source/mlbuilder.py" "${CMAKE_SOURCE_DIR}/source/mlinterface.py" "${CMAKE_SOURCE_DIR}/source/mlinterface.h"
    VERBATIM
    COMMENT "Generating Deep-HP CFFI mlinterface plugin"
  )
  add_custom_target(tinkerhp_mlinterface_build DEPENDS "${_plugin}")
  add_library(tinkerhp_mlinterface SHARED IMPORTED GLOBAL)
  set_target_properties(tinkerhp_mlinterface PROPERTIES IMPORTED_LOCATION "${_plugin}")
  add_library(tinkerhp_deephp_python SHARED IMPORTED GLOBAL)
  set_target_properties(tinkerhp_deephp_python PROPERTIES
    IMPORTED_LOCATION "${TINKERHP_DEEPHP_PYTHON_LIBRARY}"
  )
  target_link_libraries(tinkerhp_mlinterface INTERFACE tinkerhp_deephp_python)
  add_dependencies(tinkerhp_mlinterface tinkerhp_mlinterface_build)
  set(TINKERHP_DEEPHP_PLUGIN "${_plugin}" PARENT_SCOPE)
endfunction()
