#[=======================================================================[.rst:
FindFCLib
---------

Finds the FCLib library (HDF5-based friction-contact problem format), used
by the optional ``FCLibLoader`` app.

Respects the ``FCLIB_ROOT`` environment variable for non-standard installs.

Imported Targets
^^^^^^^^^^^^^^^^^

``FCLib::FCLib``
  The FCLib library, if found.

Result Variables
^^^^^^^^^^^^^^^^^

``FCLib_FOUND``
  True if FCLib was found.
``FCLib_INCLUDE_DIRS``
  Include directory needed to use FCLib.
``FCLib_LIBRARIES``
  Library to link against for FCLib.
#]=======================================================================]

find_library(FCLib_LIBRARY
  NAMES fclib
  PATHS
    "${CMAKE_INSTALL_PREFIX}/lib"
    "$ENV{FCLIB_ROOT}/lib"
)
find_path(FCLib_INCLUDE_DIR
  NAMES fclib.h
  PATHS
    "${CMAKE_INSTALL_PREFIX}/include"
    "$ENV{FCLIB_ROOT}/include"
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FCLib
  REQUIRED_VARS          FCLib_LIBRARY FCLib_INCLUDE_DIR
  REASON_FAILURE_MESSAGE "FCLib library not found - set the FCLIB_ROOT environment variable"
)

if(FCLib_FOUND AND NOT TARGET FCLib::FCLib)
  add_library(FCLib::FCLib UNKNOWN IMPORTED)
  set_target_properties(FCLib::FCLib PROPERTIES
    IMPORTED_LOCATION             "${FCLib_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${FCLib_INCLUDE_DIR}"
  )
endif()

# Legacy variable names, kept for callers that pre-date the FCLib::FCLib target.
set(FCLib_LIBRARIES    "${FCLib_LIBRARY}")
set(FCLib_INCLUDE_DIRS "${FCLib_INCLUDE_DIR}")

mark_as_advanced(FCLib_LIBRARY FCLib_INCLUDE_DIR)
