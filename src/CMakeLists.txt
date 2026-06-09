# ── So-bogus compiled library ─────────────────────────────────────────────────
#
# Compiles only the three Interfaces .cpp files.  Everything else in src/ is
# header-only and is consumed through bogus::bogus.
#
# License: GPL-2.0-or-later
# ─────────────────────────────────────────────────────────────────────────────

set(SOBOGUS_SOURCES
  "${CMAKE_CURRENT_SOURCE_DIR}/Interfaces/FrictionProblem.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/Interfaces/MecheInterface.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/Interfaces/PrimalFrictionProblem.cpp"
)

add_library(sobogus SHARED ${SOBOGUS_SOURCES})
add_library(bogus::sobogus ALIAS sobogus)

# Inherit Eigen, OpenMP, Boost, MKL transitively from the INTERFACE target
target_link_libraries(sobogus PUBLIC bogus::bogus)

# SOC dimension instantiation (only affects this compiled unit)
target_compile_definitions(sobogus PRIVATE
  $<$<BOOL:${BOGUS_WITH_2D}>:BOGUS_INSTANTIATE_2D_SOC>
  $<$<BOOL:${BOGUS_WITH_3D}>:BOGUS_INSTANTIATE_3D_SOC>
  $<$<BOOL:${BOGUS_WITH_DYNAMIC}>:BOGUS_INSTANTIATE_DYNAMIC_SOC>
)

# Compiler warnings – applied only to this target, not propagated
# TODO: Fix this
# target_compile_options(sobogus PRIVATE
#   $<$<OR:
#       $<CXX_COMPILER_ID:GNU>,
#       $<CXX_COMPILER_ID:Clang>,
#       $<CXX_COMPILER_ID:AppleClang>>:
#     -Wall -Wextra -Wno-ignored-attributes>
#   $<$<CXX_COMPILER_ID:MSVC>:
#     /W3 /wd4305>
# )

set_target_properties(sobogus PROPERTIES
  VERSION                 ${PROJECT_VERSION}
  SOVERSION               ${PROJECT_VERSION_MAJOR}
  EXPORT_NAME             sobogus
  POSITION_INDEPENDENT_CODE ON
)

# ── Install ───────────────────────────────────────────────────────────────────
include(GNUInstallDirs)

install(TARGETS sobogus
  EXPORT   bogus-targets
  LIBRARY  DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  ARCHIVE  DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  RUNTIME  DESTINATION "${CMAKE_INSTALL_BINDIR}"
)

# Interfaces headers alongside Core/Extra under <prefix>/include/bogus/
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/Interfaces"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/bogus"
  FILES_MATCHING PATTERN "*.hpp"
)