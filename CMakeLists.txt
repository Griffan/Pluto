cmake_minimum_required(VERSION 2.8.4)
project(Pluto)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -Wall -g -O3 -Wextra -Wno-char-subscripts -Wno-unused-parameter -Wno-unused-value -fno-omit-frame-pointer -mno-omit-leaf-frame-pointer ")
set (CMAKE_CXX_STANDARD 11)
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
#SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
#SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")

add_definitions(-D__ZLIB_AVAILABLE__)
add_definitions(-DBUILD_TESTING=OFF)

find_library(LIBSTATGEN NAMES libStatGen StatGen)
if(NOT LIBSTATGEN)
    message(STATUS "libStatGen not found")
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
            BUILD_COMMAND make opt CFLAGS="-fpermissive"# CXX="/usr/local/bin/g++-7" CC="/usr/local/bin/g++-7"
            INSTALL_COMMAND ""
            )
    set(LIBSTATGEN_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/SinglePhasing/libStatGen/include/)
    set(LIBSTATGEN_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/SinglePhasing/libStatGen/)
    set(LIBSTATGEN {LIBSTATGEN_LIB_DIR}/libStatGen.a)
else()
    message(STATUS "libStatGen found: " ${LIBSTATGEN})
    get_filename_component(LIBSTATGEN_INCLUDE_PARENT_DIR ${LIBSTATGEN} DIRECTORY)
    set(LIBSTATGEN_INCLUDE_DIR ${LIBSTATGEN_INCLUDE_PARENT_DIR}/include/)
    get_filename_component(LIBSTATGEN_PARENT_DIR ${LIBSTATGEN} DIRECTORY)
    set(LIBSTATGEN_LIB_DIR ${LIBSTATGEN_PARENT_DIR})
endif()
message(STATUS "LIBSTATGEN_INCLUDE_DIR: " ${LIBSTATGEN_INCLUDE_DIR})

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/abseil-cpp/)

add_subdirectory(abseil-cpp)
add_subdirectory(SinglePhasing)
add_subdirectory(pbwtWrapper)

#add_subdirectory(TestUnit)


set(SOURCE_FILES main.cpp)
add_executable(Pluto ${SOURCE_FILES})
target_link_libraries(Pluto  PBWTWrapper thunderVCF absl::base absl::flat_hash_map absl::flat_hash_set ${ZLIB_LIBRARIES})