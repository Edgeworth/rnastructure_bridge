cmake_minimum_required(VERSION 3.12)

project(rnastructure_bridge VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-fsched-spec-load" HAVE_RNASTRUCTURE_FLAG)
if(HAVE_RNASTRUCTURE_FLAG)
  set(RNASTRUCTURE_FLAGS "-fsched-spec-load")
endif()

file(GLOB_RECURSE SOURCE src "*.cpp")
file(GLOB_RECURSE HEADERS src "*.h")

add_library(rnastructure_bridge ${SOURCE} ${HEADERS})
target_include_directories(rnastructure_bridge PUBLIC include)
target_compile_options(rnastructure_bridge PRIVATE
  ${RNASTRUCTURE_FLAGS} -w -march=native -O3 -Wno-write-strings -DNDEBUG)
