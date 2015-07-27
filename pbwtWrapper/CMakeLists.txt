cmake_minimum_required(VERSION 3.2)
project(PBWTWrapper)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES PBWTWrapper.cpp)
add_library(PBWTWrapper ${SOURCE_FILES})