cmake_minimum_required(VERSION 3.14)
 
project(bogus
  VERSION       1.4.2
  DESCRIPTION   "Sparse block-matrix library for Second Order Cone problems"
  HOMEPAGE_URL  "https://github.com/JonasBreuling/so-bogus"
  LANGUAGES     CXX
)
 
# ── Guard against in-source builds ───────────────────────────────────────────
if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(FATAL_ERROR
    "In-source builds are not supported.\n"
    "Please create a separate build directory and run cmake from there.")
endif()
 
# ── Default build type ───────────────────────────────────────────────────────
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    Debug Release RelWithDebInfo MinSizeRel)
  message(STATUS "Build type not set – defaulting to Release")
endif()
 
# ── Detect whether we are the top-level project ──────────────────────────────
# When consumed via FetchContent the options default to OFF so we don't
# build tests/apps/sobogus inside the consumer's project.
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  set(_BOGUS_IS_TOPLEVEL ON)
else()
  set(_BOGUS_IS_TOPLEVEL OFF)
endif()
 
# ── Options ──────────────────────────────────────────────────────────────────
option(BOGUS_BUILD_TESTS
  "Build the bogus unit-test suite"             ${_BOGUS_IS_TOPLEVEL})
option(BOGUS_BUILD_APPS
  "Build the example applications"              ${_BOGUS_IS_TOPLEVEL})
option(BOGUS_BUILD_SOBOGUS
  "Build the compiled So-bogus interface layer" ${_BOGUS_IS_TOPLEVEL})
option(BOGUS_BUILD_EXAMPLES
  "Build more bogus examples"             	  	${_BOGUS_IS_TOPLEVEL})
 
option(BOGUS_WITH_OPENMP
  "Parallelise via OpenMP"                      ON)
option(BOGUS_WITH_BOOST
  "Enable boost::serialization I/O support"     ON)
option(BOGUS_WITH_MKL
  "Enable Intel MKL Block-Sparse-Row multiply"  OFF)
 
# Dimension flags only affect the compiled So-bogus library
option(BOGUS_WITH_2D
  "Instantiate 2-D Second Order Cone solvers"   ON)
option(BOGUS_WITH_3D
  "Instantiate 3-D Second Order Cone solvers"   ON)
option(BOGUS_WITH_DYNAMIC
  "Instantiate dynamic-dimension SOC solvers"   OFF)
 
# ── Custom finders (FindMKL, FindFCLib) ──────────────────────────
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
# Keep BOGUS_ROOT so the legacy finders can resolve relative search paths
set(BOGUS_ROOT "${CMAKE_CURRENT_SOURCE_DIR}")
 
# ══════════════════════════════════════════════════════════════════════════════
#  Header-only INTERFACE target  bogus::bogus
# ══════════════════════════════════════════════════════════════════════════════
add_library(bogus INTERFACE)
add_library(bogus::bogus ALIAS bogus)
 
# Require C++17 from consumers
target_compile_features(bogus INTERFACE cxx_std_17)
 
# Build-tree includes point at src/; installed tree includes point at
# <prefix>/include/bogus/ so users always write  #include <Core/...>
target_include_directories(bogus INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:include/bogus>
)
 
# ── Eigen3 (required) ────────────────────────────────────────────────────────
find_package(Eigen3 3.4...5.0.0 QUIET CONFIG)
if(NOT Eigen3_FOUND)
    include(FetchContent)
    FetchContent_Declare(
        Eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG 5.0.0
        GIT_SHALLOW TRUE
    )
    FetchContent_MakeAvailable(Eigen)
else()
    message(STATUS "Found system Eigen3: ${EIGEN3_VERSION}")
endif()
# target_link_libraries(bogus INTERFACE Eigen3::Eigen)
# message(STATUS "bogus: Eigen3 ${EIGEN3_VERSION} – OK")

# IMPORTED targets (Eigen3::Eigen, OpenMP::OpenMP_CXX, Boost::serialization …)
# cannot appear in an install() export set.  Wrap them in BUILD_INTERFACE so
# they are used during the build but are NOT written into bogus-targets.cmake.
# bogus-config.cmake.in calls find_dependency() to recreate them on the
# consumer side and then re-attaches them via set_property().
target_link_libraries(bogus INTERFACE
  $<BUILD_INTERFACE:Eigen3::Eigen>)
message(STATUS "bogus: Eigen3 ${EIGEN3_VERSION} – OK")
# Record for substitution into bogus-config.cmake.in
set(BOGUS_EIGEN3_VERSION "${EIGEN3_VERSION}")
 
# ── OpenMP (optional) ─────────────────────────────────────────────────────────
if(BOGUS_WITH_OPENMP)
  find_package(OpenMP QUIET)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(bogus INTERFACE OpenMP::OpenMP_CXX)
    # Eigen's own parallelisation clashes with bogus's; disable it
    target_compile_definitions(bogus INTERFACE EIGEN_DONT_PARALLELIZE)
    message(STATUS "bogus: OpenMP ${OpenMP_CXX_VERSION} – enabled")
  else()
    message(STATUS "bogus: OpenMP not found – single-threaded build")
    # Suppress spurious "unknown pragma omp" warnings
    target_compile_options(bogus INTERFACE
      $<$<OR:
          $<CXX_COMPILER_ID:GNU>,
          $<CXX_COMPILER_ID:Clang>,
          $<CXX_COMPILER_ID:AppleClang>>:
        -Wno-unknown-pragmas>)
  endif()
endif()
 
# ── Boost serialization (optional) ───────────────────────────────────────────
if(BOGUS_WITH_BOOST)
  find_package(Boost 1.36 QUIET COMPONENTS serialization)
  if(Boost_FOUND)
    target_compile_definitions(bogus INTERFACE BOGUS_WITH_BOOST_SERIALIZATION)
    target_link_libraries(bogus INTERFACE Boost::serialization)
    message(STATUS "bogus: Boost::serialization ${Boost_VERSION} – enabled")
  else()
    message(STATUS "bogus: Boost not found – serialisation disabled")
  endif()
endif()
 
# ── Intel MKL (optional) ──────────────────────────────────────────────────────
if(BOGUS_WITH_MKL)
  find_package(MKL QUIET)     # uses cmake/FindMKL.cmake
  if(MKL_FOUND)
    target_compile_definitions(bogus INTERFACE BOGUS_WITH_MKL)
    target_include_directories(bogus INTERFACE "${MKL_INCLUDE_DIRS}")
    target_link_libraries(bogus INTERFACE
      "${MKL_LIBRARIES}"
      mkl_intel_lp64 mkl_intel_thread mkl_core iomp5)
    message(STATUS "bogus: Intel MKL BSR bindings – enabled")
  else()
    message(STATUS "bogus: MKL not found – BSR bindings disabled")
  endif()
endif()
 
# ══════════════════════════════════════════════════════════════════════════════
#  Sub-directories
# ══════════════════════════════════════════════════════════════════════════════
if(BOGUS_BUILD_SOBOGUS)
  add_subdirectory(src)
endif()
 
if(BOGUS_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
 
if(BOGUS_BUILD_APPS)
  add_subdirectory(apps)
endif()
 
if(BOGUS_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()
 
# ══════════════════════════════════════════════════════════════════════════════
#  Installation
# ══════════════════════════════════════════════════════════════════════════════
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
 
# Install the bogus::bogus INTERFACE target
install(TARGETS bogus
  EXPORT bogus-targets)
 
# Headers: Core and Extra land under  <prefix>/include/bogus/{Core,Extra}/
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src/Core"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/bogus"
  FILES_MATCHING PATTERN "*.hpp" PATTERN "*.h")
 
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src/Extra"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/bogus"
  FILES_MATCHING PATTERN "*.hpp")
 
# Version file
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/bogus-config-version.cmake"
  VERSION       ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion)
 
# Config file (fills in @PACKAGE_INIT@ and paths)
configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/bogus-config.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/bogus-config.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/bogus")
 
# Export the targets file
install(EXPORT bogus-targets
  FILE        bogus-targets.cmake
  NAMESPACE   bogus::
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/bogus")
 
# Install config + version files
install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/bogus-config.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/bogus-config-version.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/bogus")