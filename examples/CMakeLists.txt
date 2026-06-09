cmake_minimum_required(VERSION 3.16)

add_executable(hello hello.cpp)
# target_link_libraries(hello PRIVATE bogus)
# target_include_directories(hello PRIVATE
#     ../src/
# )
target_link_libraries(hello PRIVATE
  bogus::bogus      # header-only Core + Extra, pulls Eigen / OpenMP / Boost
)

add_executable(solver solver.cpp)
target_link_libraries(solver PRIVATE bogus)
# target_include_directories(solver PRIVATE
#     ../src/
# )
target_link_libraries(solver PRIVATE
  bogus::bogus      # header-only Core + Extra, pulls Eigen / OpenMP / Boost
)
