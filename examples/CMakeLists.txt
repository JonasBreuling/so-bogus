add_executable(hello hello.cpp)
target_link_libraries(hello PRIVATE
  bogus::bogus      # header-only Core + Extra, pulls Eigen / OpenMP / Boost
)

add_executable(solver solver.cpp)
target_link_libraries(solver PRIVATE
  bogus::bogus      # header-only Core + Extra, pulls Eigen / OpenMP / Boost
)

add_executable(block_direct_solver block_direct_solver.cpp)
target_link_libraries(block_direct_solver PRIVATE
  bogus::bogus      # header-only Core + Extra, pulls Eigen / OpenMP / Boost
)
