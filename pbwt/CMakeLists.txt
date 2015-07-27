cmake_minimum_required(VERSION 3.2)
project(pbwt)
#set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "/Users/fanzhang/BioLib")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}  -O3")
include_directories("/Users/fanzhang/BioLib/include/")
file(GLOB SOURCE_FILES
    "*.h"
    "*.c"
)
LINK_DIRECTORIES(/Users/fanzhang/BioLib/lib)
#FIND_LIBRARY(PBWT_LIBRARY NAMES "hts" PATHS /Users/fanzhang/BioLib/lib)
add_executable(pbwt ${SOURCE_FILES})
#set_target_properties(Vcf PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/libVcf")

#FIND_LIBRARY(APRUTIL NAMES "aprutil-1" PATHS ${PROJECT_SOURCE_DIR}/tools/apr/libs NO_DEFAULT_PATH)
TARGET_LINK_LIBRARIES(pbwt z hts)