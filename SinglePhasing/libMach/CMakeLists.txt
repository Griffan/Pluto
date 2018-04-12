cmake_minimum_required(VERSION 2.8.4)
project(libMach)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3")

file(GLOB SOURCE_FILES
    "*.h"
    "*.cpp"
)
add_library(Mach ${SOURCE_FILES})
target_include_directories(Mach PRIVATE ${LIBSTATGEN_INCLUDE_DIR})
add_dependencies(Mach LIBSTATGEN)
