cmake_minimum_required(VERSION 3.2)

project(miles_rnastructure)

include(CheckCCompilerFlag)
check_c_compiler_flag("-fsched-spec-load" HAVE_RNASTRUCTURE_FLAG)
if(HAVE_RNASTRUCTURE_FLAG)
  set(RNASTRUCTURE_FLAGS "-fsched-spec-load")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${RNASTRUCTURE_FLAGS} -w -std=c++14 -march=native \
  -O3 -Wno-write-strings -D NDEBUG")
file(GLOB_RECURSE SOURCE src "*.cpp")
file(GLOB_RECURSE HEADERS src "*.h")
include_directories(include)
add_library(miles_rnastructure ${SOURCE} ${HEADERS})
