cmake_minimum_required(VERSION 2.8)
project(thunderVCF)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg -g -Wall -std=c++11 -O3")

file(GLOB SOURCE_FILES
    "*.h"
    "*.cpp"
)

add_library(thunderVCF ${SOURCE_FILES})
add_dependencies(thunderVCF LIBSTATGEN)
target_link_libraries(thunderVCF Mach Vcf  PBWTWrapper  z ${LIBSTATGEN_LIB_DIR}/libStatGen.a)
target_include_directories(thunderVCF PRIVATE  ${LIBSTATGEN_INCLUDE_DIR})
target_include_directories(thunderVCF PRIVATE "../libMach/")
target_include_directories(thunderVCF PRIVATE "../libVcf/")
target_include_directories(thunderVCF PRIVATE "../../pbwtWrapper/")
target_include_directories(thunderVCF PRIVATE ${LIBRMATH_INCLUDE_DIR})
