cmake_minimum_required(VERSION 3.16)

add_executable(hello hello.cpp)
target_link_libraries(hello PRIVATE bogus)
target_include_directories(hello PRIVATE
    ../src/
)
