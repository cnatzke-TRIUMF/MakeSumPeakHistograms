#------------------------------------------------------------------
#Setup the project
#------------------------------------------------------------------
cmake_minimum_required(VERSION 3.2)

project(SumPeakHistograms VERSION 0.1)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

# configure header file
configure_file(SumPeakHistograms.h.in SumPeakHistograms.h)

set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER "g++")

set(CMAKE_CXX_FLAGS "-Wall -O3 ${CMAKE_CXX_FLAGS} ${GRSI_CONFIG}")

# Connect ROOT to project
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

#link_directories($ENV{GRSISYS}/lib $ENV{GRSISYS}/GRSIData/lib)
#Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT CONFIG REQUIRED)
include(${ROOT_USE_FILE})

# Connect GRSISORT headers
set(GRSI_INCLUDE_DIRS $ENV{GRSISYS}/include $ENV{GRSISYS}/GRSIData/include)

# Grsisort libraries for linking
execute_process(COMMAND grsi-config --all-libs --GRSIData-libs --root
   OUTPUT_VARIABLE GRSI_CONFIG
   OUTPUT_STRIP_TRAILING_WHITESPACE
)
# Message for debugging
#message(STATUS "Found Grsisort libraries: ${GRSI_CONFIG}")


# Adding src files
file(GLOB SOURCES ${PROJECT_SOURCE_DIR}/src/*.cpp)

# Naming main executable
add_executable(SumPeakHistograms ${SOURCES})

add_custom_command(TARGET SumPeakHistograms PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${PROJECT_SOURCE_DIR}/misc/ $<TARGET_FILE_DIR:SumPeakHistograms>)

# linking libraries
target_link_libraries(SumPeakHistograms PUBLIC
   ${GRSI_CONFIG}
)

# add the binary tree to the search path for include files so that we will find header files
target_include_directories(SumPeakHistograms PUBLIC
    "${PROJECT_BINARY_DIR}"
    "${GRSI_INCLUDE_DIRS}"
    include
    )

# add install targets
install(TARGETS SumPeakHistograms DESTINATION "${PROJECT_BINARY_DIR}/bin")
install(FILES "${PROJECT_BINARY_DIR}/SumPeakHistograms.h"
   DESTINATION "${PROJECT_BINARY_DIR}/include"
)
