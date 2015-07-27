cmake_minimum_required(VERSION 3.2)
project(Pluto)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
add_subdirectory(pbwt)
add_subdirectory(SinglePhasing)
add_subdirectory(pbwtWrapper)
set(SOURCE_FILES main.cpp)
add_executable(Pluto ${SOURCE_FILES})