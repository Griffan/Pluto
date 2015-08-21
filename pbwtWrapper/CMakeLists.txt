cmake_minimum_required(VERSION 2.8)
project(PBWTWrapper)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
add_subdirectory(pbwt)
set(SOURCE_FILES PBWTWrapper.cpp)
add_library(PBWTWrapper ${SOURCE_FILES})
target_link_libraries(PBWTWrapper pbwt)