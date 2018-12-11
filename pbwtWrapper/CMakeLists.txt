cmake_minimum_required(VERSION 2.8.4)
project(PBWTWrapper)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3  -std=c++11")
set(SOURCE_FILES PBWTWrapper.cpp KSTest.cpp DAG.cpp)

add_library(PBWTWrapper ${SOURCE_FILES})
add_dependencies(PBWTWrapper LIBSTATGEN)
target_include_directories(PBWTWrapper PRIVATE ${LIBSTATGEN_INCLUDE_DIR})



