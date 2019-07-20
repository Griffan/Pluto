cmake_minimum_required(VERSION 2.8)
project(thunderVCF)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -std=c++11 -O3")

file(GLOB SOURCE_FILES
    "*.h"
    "*.cpp"
)
find_path(ZLIB_INCLUDE_DIRS zlib.h HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)
if(NOT ZLIB_INCLUDE_DIRS )
    message(FATAL_ERROR "libz ZLIB_INCLUDE_DIRS not found")
endif()

message(STATUS ${HTS_INCLUDE_DIRS})

find_library(ZLIB_LIBRARIES z HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)
if(NOT ZLIB_LIBRARIES)
    message(FATAL_ERROR "libz ZLIB_LIBRARIES not found")
endif()

add_library(thunderVCF ${SOURCE_FILES} FwdBwdLocalParameter.cpp FwdBwdLocalParameter.h)
add_dependencies(thunderVCF LIBSTATGEN)
target_link_libraries(thunderVCF Mach Vcf  PBWTWrapper  ${ZLIB_LIBRARIES} ${LIBSTATGEN_LIB_DIR}/libStatGen.a)
target_include_directories(thunderVCF PRIVATE  ${LIBSTATGEN_INCLUDE_DIR})
target_include_directories(thunderVCF PRIVATE "../libMach/")
target_include_directories(thunderVCF PRIVATE "../libVcf/")
target_include_directories(thunderVCF PRIVATE "../../pbwtWrapper/")
target_include_directories(thunderVCF PRIVATE ${LIBRMATH_INCLUDE_DIR} ${ZLIB_INCLUDE_DIRS})
