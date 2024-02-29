# Ian Wark July 6 2017

cmake_minimum_required(VERSION 3.1)
project(HFold)

# enable c++11
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} \ -std=c++17 \ -O3")
# only continue if we know we already have boost
# create simfold library
add_subdirectory(simfold)

set(SOURCE
    W_final.cpp 
    pseudo_loop.cpp 
    V_final.cpp 
    VM_final.cpp
    h_common.cpp 
    hfold_validation.cpp 
    HFold.cpp
    cmdline.cc
    Result.cc
)

# pass SIMFOLD_HOME to binary
set(SIMFOLD_HOME ${CMAKE_SOURCE_DIR}/simfold)
add_definitions("-DSIMFOLD_HOME=\"${SIMFOLD_HOME}\"")

# set path for executable to be created
# create executable
add_executable(HFold ${SOURCE})

# enable C++11
target_compile_features(HFold PRIVATE cxx_nullptr)

# link to simfold
target_link_libraries(HFold LINK_PUBLIC simfold)
