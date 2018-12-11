cmake_minimum_required(VERSION 2.8.4)
project(SinglePhasing)

#set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
#set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib")
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/build")
#set(EXECUTABLE_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/build/")
#set(LIBRARY_OUTPUT_PATH "${PROJECT_SOURCE_DIR}")

add_subdirectory(libMach)
add_subdirectory(libVcf)
add_subdirectory(thunderVCF)


