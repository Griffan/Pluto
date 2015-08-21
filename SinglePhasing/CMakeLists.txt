cmake_minimum_required(VERSION 2.8)
project(SinglePhasing)

#set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
#set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib")
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/build")
#set(EXECUTABLE_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/build/")
#set(LIBRARY_OUTPUT_PATH "${PROJECT_SOURCE_DIR}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")


message(STATUS "dir is " ${PROJECT_SOURCE_DIR})
#set(SOURCE_FILES src/FASTQuick.cpp src/BwtIndexer.cpp src/BwtMapper.cpp src/PopulationIdentifier.cpp src/RefBuilder.cpp src/StatCollector.cpp)

#add_subdirectory(libStatGen)
add_subdirectory(libMach)
LINK_DIRECTORIES(libStatGen)
add_subdirectory(libVcf)
add_subdirectory(thunderVCF)

#include_directories()

#add_executable(thunderVCF ${SOURCE_FILES})
#target_link_libraries(thunderVCF bwa mpu libStatGen.a)
#target_include_directories(thunderVCF PRIVATE "/Users/fanzhang/BioLib/libStatGen/include")
