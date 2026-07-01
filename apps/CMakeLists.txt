# ── Example applications ──────────────────────────────────────────────────────
# Both apps require the compiled So-bogus interface layer (bogus::sobogus).
# They are only built when BOGUS_BUILD_SOBOGUS is ON (the default when bogus
# is the top-level project).
# ─────────────────────────────────────────────────────────────────────────────

if(NOT TARGET bogus::sobogus)
  message(STATUS
    "bogus: sobogus not built – skipping apps. "
    "Enable BOGUS_BUILD_SOBOGUS to build them.")
  return()
endif()

# ── LoadFrictionProblem ───────────────────────────────────────────────────────
# Reads and solves serialised bogus friction-problem files.
# Requires Boost::serialization (enabled by BOGUS_WITH_BOOST).
if(TARGET Boost::serialization)
  add_executable(LoadFrictionProblem
    "${CMAKE_CURRENT_SOURCE_DIR}/LoadFrictionProblem.cpp"
  )
  target_link_libraries(LoadFrictionProblem PRIVATE bogus::sobogus)
  target_compile_options(LoadFrictionProblem PRIVATE
    $<$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>>:-Wall -Wextra -Wno-ignored-attributes>
    $<$<CXX_COMPILER_ID:MSVC>:/W3>
  )
  message(STATUS "bogus: app LoadFrictionProblem – enabled")
else()
  message(STATUS
    "bogus: app LoadFrictionProblem – skipped (needs Boost::serialization)")
endif()

# ── FCLibLoader ───────────────────────────────────────────────────────────────
# Reads FCLib HDF5 friction-contact problems and solves them.
# Requires the FCLib library (set FCLIB_ROOT env-var for non-standard installs).
find_package(FCLib QUIET)   # uses cmake/FindFCLib.cmake
if(FCLib_FOUND)
  add_executable(FCLibLoader
    "${CMAKE_CURRENT_SOURCE_DIR}/FCLibLoader.cpp"
  )
  target_link_libraries(FCLibLoader PRIVATE
    bogus::sobogus
    FCLib::FCLib
  )
  target_compile_options(FCLibLoader PRIVATE
    $<$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>>:-Wall -Wextra -Wno-ignored-attributes>
    $<$<CXX_COMPILER_ID:MSVC>:/W3>
  )
  message(STATUS "bogus: app FCLibLoader – enabled")
else()
  message(STATUS
    "bogus: app FCLibLoader – skipped (FCLib not found; set FCLIB_ROOT)")
endif()