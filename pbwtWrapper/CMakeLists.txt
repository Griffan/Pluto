cmake_minimum_required(VERSION 2.8.4)
project(PBWTWrapper)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3  -std=c++11")


#add_subdirectory(pbwt)

set(SOURCE_FILES PBWTWrapper.cpp ks.cpp)

add_library(PBWTWrapper ${SOURCE_FILES})
add_dependencies(PBWTWrapper LIBRMATH)
#target_link_libraries(PBWTWrapper pbwt)
target_include_directories(PBWTWrapper PRIVATE ${LIBRMATH_INCLUDE_DIR})



