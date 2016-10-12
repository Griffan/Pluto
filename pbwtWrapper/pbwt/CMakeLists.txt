cmake_minimum_required(VERSION 2.8)
project(pbwt)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}  -O3")

INCLUDE_DIRECTORIES(/Users/fanzhang/BioLib/include)
file(GLOB SOURCE_FILES
    "*.h"
    "*.c"
)

FIND_LIBRARY(LIBHTS hts HINTS /Users/fanzhang/BioLib/lib)

add_library(pbwt ${SOURCE_FILES})

TARGET_LINK_LIBRARIES(pbwt ${LIBHTS} z)
