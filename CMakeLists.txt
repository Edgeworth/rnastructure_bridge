cmake_minimum_required(VERSION 3.2)

project(miles_rnastructure)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w -std=c++14 -march=native -O3 -Wno-write-strings -fsched-spec-load -D NDEBUG")
file(GLOB SOURCE "src/*.cpp" "include/*.h")
include_directories(include)
add_library(miles_rnastructure ${SOURCE})
