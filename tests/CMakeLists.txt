# ── Tests ─────────────────────────────────────────────────────────────────────
# Pulls Google Test at configure time via FetchContent (no internet dependency
# at build time after the first configure).
# ─────────────────────────────────────────────────────────────────────────────

include(FetchContent)

# ── Google Test v1.14.0 ───────────────────────────────────────────────────────
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        v1.14.0
  GIT_SHALLOW    TRUE          # fetch only the tag, not full history
)

# Prevent gtest from hijacking our install() targets
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
# On Windows, link the CRT the same way the rest of the project does
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(googletest)

# ── Test sources ──────────────────────────────────────────────────────────────
# All tests exercise the header-only bogus::bogus target directly;
# they do not need the compiled sobogus library.
set(BOGUS_TEST_SOURCES
  ADMM.cpp
  BlockMV.cpp
  Flat.cpp
  GaussSeidel.cpp
  Krylov.cpp
  Map.cpp
  Mkl.cpp
  Polynomial.cpp
  ProductGaussSeidel.cpp
  ProjectedGradient.cpp
  Serialization.cpp
  SparseBlock.cpp
)

list(TRANSFORM BOGUS_TEST_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/")

# ── Test executable ───────────────────────────────────────────────────────────
add_executable(testbogus ${BOGUS_TEST_SOURCES})

target_link_libraries(testbogus PRIVATE
  bogus::bogus      # header-only Core + Extra, pulls Eigen / OpenMP / Boost
  GTest::gtest_main # runner + main()
)

# Compiler warnings for the test binary (not propagated)
target_compile_options(testbogus PRIVATE
  $<$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>>:-Wall -Wextra -Wno-ignored-attributes>
  $<$<CXX_COMPILER_ID:MSVC>:/W3 /wd4305>
)

# ── Register with CTest ───────────────────────────────────────────────────────
include(GoogleTest)
gtest_discover_tests(testbogus
  PROPERTIES TIMEOUT 120     # seconds per test case
)