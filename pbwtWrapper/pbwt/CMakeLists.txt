cmake_minimum_required(VERSION 2.8)
project(pbwt)
#set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "/Users/fanzhang/BioLib")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}  -O3")
include_directories("/net/wonderland/home/fanzhang/WorkingSpace/install-packages/htslib")
file(GLOB SOURCE_FILES
    "*.h"
    "*.c"
)
LINK_DIRECTORIES(/net/wonderland/home/fanzhang/WorkingSpace/lib/)
#FIND_LIBRARY(PBWT_LIBRARY NAMES "hts" PATHS /Users/fanzhang/BioLib/lib)
#add_executable(pbwt ${SOURCE_FILES})
add_library(pbwt ${SOURCE_FILES})
#set_target_properties(Vcf PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/libVcf")

#FIND_LIBRARY(APRUTIL NAMES "aprutil-1" PATHS ${PROJECT_SOURCE_DIR}/tools/apr/libs NO_DEFAULT_PATH)
TARGET_LINK_LIBRARIES(pbwt hts "/usr/lib/x86_64-linux-gnu/libz.a")
