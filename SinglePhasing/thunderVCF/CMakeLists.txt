cmake_minimum_required(VERSION 2.8)
project(thunderVCF)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -std=c++11 -O3")

file(GLOB SOURCE_FILES
    "*.h"
    "*.cpp"
)
LINK_DIRECTORIES("../libStatGen")

#add_executable(thunderVCF ${SOURCE_FILES})
add_library(thunderVCF ${SOURCE_FILES})
target_link_libraries(thunderVCF Mach Vcf StatGen PBWTWrapper z)
target_include_directories(thunderVCF PRIVATE "../libStatGen/include")
target_include_directories(thunderVCF PRIVATE "../libMach/")
target_include_directories(thunderVCF PRIVATE "../libVcf/")
target_include_directories(thunderVCF PRIVATE "../../pbwtWrapper/")
