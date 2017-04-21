cmake_minimum_required(VERSION 2.8)
project(Pluto)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -Wall -g -O3 -std=c++11 -Wextra")
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
#SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
#SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")

add_definitions(-D__ZLIB_AVAILABLE__)

#find_library(LIBRMATH NAMES libRmath Rmath)
#find_path(LIBRMATH_INCLUDE NAMES Rmath.h)
if(NOT LIBRMATH OR NOT LIBRMATH_INCLUDE)
    include(ExternalProject)
    ExternalProject_Add(LIBRMATH
            PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/build
            STAMP_DIR ${CMAKE_CURRENT_SOURCE_DIR}/build/stamp

            GIT_REPOSITORY https://github.com/Griffan/Rmath.git
            GIT_TAG "master"

            UPDATE_COMMAND ""

            SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/pbwtWrapper/libRmath
            CONFIGURE_COMMAND ""
            BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/pbwtWrapper/libRmath
            BUILD_COMMAND make
            INSTALL_COMMAND ""
            )

    set(LIBRMATH_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/pbwtWrapper/libRmath/Rmath)
    set(LIBRMATH_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/pbwtWrapper/libRmath)
else()
    message(STATUS "libRmath found: " ${LIBRMATH})
    message(STATUS "libRmath_include found: " ${LIBRMATH_INCLUDE})
    set(LIBRMATH_INCLUDE_DIR ${LIBRMATH_INCLUDE})
    get_filename_component(LIBRMATH_PARENT_DIR ${LIBRMATH} DIRECTORY)
    set(LIBRMATH_LIB_DIR ${LIBRMATH_PARENT_DIR})
endif()





find_library(LIBSTATGEN NAMES libStatGen StatGen)
if(NOT LIBSTATGEN)
    include(ExternalProject)
    ExternalProject_Add(LIBSTATGEN
            PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/build
            STAMP_DIR ${CMAKE_CURRENT_SOURCE_DIR}/build/stamp
            GIT_REPOSITORY https://github.com/statgen/libStatGen.git
            GIT_TAG "master"

            UPDATE_COMMAND ""

            SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/SinglePhasing/libStatGen
            CONFIGURE_COMMAND ""
            BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/SinglePhasing/libStatGen
            BUILD_COMMAND make opt
            INSTALL_COMMAND ""
            )
    set(LIBSTATGEN_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/SinglePhasing/libStatGen/include)
    set(LIBSTATGEN_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/SinglePhasing/libStatGen)
else()
    message(STATUS "libStatGen found: " ${LIBSTATGEN})

    get_filename_component(LIBSTATGEN_INCLUDE_PARENT_DIR ${LIBSTATGEN} DIRECTORY)
    set(LIBRMATH_INCLUDE_DIR ${LIBRMATH_INCLUDE_PARENT_DIR}/include)
    get_filename_component(LIBSTATGEN_PARENT_DIR ${LIBSTATGEN} DIRECTORY)
    set(LIBRMATH_LIB_DIR ${LIBSTATGEN_PARENT_DIR})
endif()


message(STATUS "LIBRMATH_INCLUDE_DIR: " ${LIBRMATH_INCLUDE_DIR})
message(STATUS "LIBSTATGEN_INCLUDE_DIR: " ${LIBSTATGEN_INCLUDE_DIR})

add_subdirectory(pbwtWrapper)
add_subdirectory(SinglePhasing)
add_subdirectory(TestUnit)

set(SOURCE_FILES main.cpp)
add_executable(Pluto ${SOURCE_FILES})
target_link_libraries(Pluto  PBWTWrapper thunderVCF)
