cmake_minimum_required(VERSION 2.8)
project(SinglePhasing)

#set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
#set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib")
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/build")
#set(EXECUTABLE_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/build/")
#set(LIBRARY_OUTPUT_PATH "${PROJECT_SOURCE_DIR}")

find_library(LIBSTATGEN libStatGen)
if(NOT LIBSTATGEN)
    include(ExternalProject)
    ExternalProject_Add(LIBSTATGEN
        PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/build
        STAMP_DIR ${CMAKE_CURRENT_SOURCE_DIR}/build/stamp
        GIT_REPOSITORY https://github.com/statgen/libStatGen.git
        GIT_TAG "master"

        UPDATE_COMMAND ""

        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/libStatGen
        CONFIGURE_COMMAND ""
        BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/libStatGen
        BUILD_COMMAND make opt
        INSTALL_COMMAND ""
    )
    set(LibStatGen_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/libStatGen/include)
    set(LibStatGen_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/libStatGen/)
endif()

add_subdirectory(libMach)
add_subdirectory(libVcf)
add_subdirectory(thunderVCF)


