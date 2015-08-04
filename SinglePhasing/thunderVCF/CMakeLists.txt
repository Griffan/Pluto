cmake_minimum_required(VERSION 3.2)
project(thunderVCF)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -lz -O3")

file(GLOB SOURCE_FILES
    "*.h"
    "*.cpp"
)
LINK_DIRECTORIES("../libStatGen")

#add_executable(thunderVCF ${SOURCE_FILES})
add_library(thunderVCF ${SOURCE_FILES})
target_link_libraries(thunderVCF Mach Vcf StatGen)
target_include_directories(thunderVCF PRIVATE "/Users/fanzhang/BioLib/libStatGen/include")
target_include_directories(thunderVCF PRIVATE "../libMach/")
target_include_directories(thunderVCF PRIVATE "../libVcf/")
