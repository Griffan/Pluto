cmake_minimum_required(VERSION 2.8)
project(PBWTWrapper)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3  -std=c++11")
add_subdirectory(pbwt)
set(SOURCE_FILES PBWTWrapper.cpp ks.cpp)
add_library(PBWTWrapper ${SOURCE_FILES})
target_link_libraries(PBWTWrapper ${CMAKE_CURRENT_SOURCE_DIR}/../Debug/libpbwt.a ${CMAKE_CURRENT_SOURCE_DIR}/../Debug/libRmath.a)

add_library(PBWTViewer PBWTViewer.cpp PBWTViewer.h ks.cpp )
target_link_libraries(PBWTViewer PBWTWrapper ${CMAKE_CURRENT_SOURCE_DIR}/../Debug/libpbwt.a ${CMAKE_CURRENT_SOURCE_DIR}/../Debug/libRmath.a)


