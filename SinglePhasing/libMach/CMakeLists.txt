cmake_minimum_required(VERSION 2.8)
project(libMach)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3")

file(GLOB SOURCE_FILES
    "*.h"
    "*.cpp"
)
add_library(Mach ${SOURCE_FILES})
#set_target_properties(Mach PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
target_include_directories(Mach PRIVATE "../libStatGen/include")
