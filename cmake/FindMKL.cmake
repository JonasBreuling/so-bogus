#[=======================================================================[.rst:
FindMKL
-------

Finds the core Intel MKL (Math Kernel Library), used for bogus's optional
Block-Sparse-Row multiply bindings.

Respects the ``MKL_ROOT`` environment variable for non-standard installs.

Imported Targets
^^^^^^^^^^^^^^^^^

``MKL::MKL``
  The MKL core library, if found.

Result Variables
^^^^^^^^^^^^^^^^^

``MKL_FOUND``
  True if MKL was found.
``MKL_INCLUDE_DIRS``
  Include directory needed to use MKL.
``MKL_LIBRARIES``
  Library to link against for the MKL core.
#]=======================================================================]

find_library(MKL_LIBRARY
  NAMES mkl_core
  PATHS
    "${CMAKE_INSTALL_PREFIX}/lib/intel/em64t"
    "$ENV{MKL_ROOT}/lib/em64t"
)
find_path(MKL_INCLUDE_DIR
  NAMES mkl.h
  PATHS
    "${CMAKE_INSTALL_PREFIX}/include"
    "$ENV{MKL_ROOT}/include"
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MKL
  REQUIRED_VARS          MKL_LIBRARY MKL_INCLUDE_DIR
  REASON_FAILURE_MESSAGE "MKL library not found - set the MKL_ROOT environment variable"
)

if(MKL_FOUND AND NOT TARGET MKL::MKL)
  add_library(MKL::MKL UNKNOWN IMPORTED)
  set_target_properties(MKL::MKL PROPERTIES
    IMPORTED_LOCATION             "${MKL_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${MKL_INCLUDE_DIR}"
  )
endif()

# Legacy variable names, kept for callers that pre-date the MKL::MKL target.
set(MKL_LIBRARIES    "${MKL_LIBRARY}")
set(MKL_INCLUDE_DIRS "${MKL_INCLUDE_DIR}")

mark_as_advanced(MKL_LIBRARY MKL_INCLUDE_DIR)
