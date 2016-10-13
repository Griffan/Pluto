cmake_minimum_required(VERSION 2.8)
project(PBWTWrapper)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3  -std=c++11")


add_subdirectory(pbwt)
set(SOURCE_FILES PBWTWrapper.cpp ks.cpp)

add_library(PBWTWrapper ${SOURCE_FILES})
target_link_libraries(PBWTWrapper pbwt)
target_include_directories(PBWTWrapper PRIVATE ${LIBRMATH_INCLUDE_DIR} ${LIBSTATGEN_INCLUDE_DIR})
add_dependencies(PBWTWrapper LIBRMATH)
#add_library(PBWTViewer PBWTViewer.cpp PBWTViewer.h ks.cpp )
#target_link_libraries(PBWTViewer PBWTWrapper pbwt ${CMAKE_CURRENT_SOURCE_DIR}/../Debug/libRmath.a)


