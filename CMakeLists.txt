# Ian Wark July 6 2017

cmake_minimum_required(VERSION 2.8)
project(HFold)

# only continue if we know we already have boost
# create simfold library
add_subdirectory(simfold)

set(SOURCE
    W_final.h W_final.cpp 
    pseudo_loop.h pseudo_loop.cpp 
    V_final.h V_final.cpp 
    VM_final.h VM_final.cpp
    h_common.h h_common.cpp 
    hfold_validation.h hfold_validation.cpp 
    hfold.h HFold.cpp
)

# enable c++11
add_compile_options(-std=c++11)

# set path for executable to be created
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR})
# create executable
add_executable(HFold ${SOURCE})

# link to simfold
target_link_libraries(HFold LINK_PUBLIC simfold)
