cmake_minimum_required(VERSION 2.8)
project(Pluto)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -Wall -O3 -std=c++11 -Wno-char-subscripts -Wno-unused-function -Wno-unused-variable -Wno-tautological-compare -Wno-logical-op-parentheses -Wno-overloaded-virtual")
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
#SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
#SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")


add_subdirectory(pbwtWrapper)
add_subdirectory(SinglePhasing)
#add_subdirectory(TestUnit)

set(SOURCE_FILES main.cpp)
add_executable(Pluto ${SOURCE_FILES})
target_link_libraries(Pluto  PBWTWrapper thunderVCF)
