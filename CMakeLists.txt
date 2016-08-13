cmake_minimum_required(VERSION 3.2)

project(rnastructure)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -O3 -Wno-write-strings -fsched-spec-load -D NDEBUG")
file(GLOB SOURCE "src/*.cpp" "src/*.h")
include_directories(src)
add_library(librnastructure ${SOURCE})
