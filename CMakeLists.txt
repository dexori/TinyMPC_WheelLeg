cmake_minimum_required(VERSION 3.16)
project(TinyMPCWheelLeg VERSION 0.1.0 LANGUAGES C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

option(BUILD_MPC "Build the extracted H7 C++ MPC task (requires Eigen)" ON)
option(BUILD_COLLISION "Build the independent collision observer" ON)
include(CTest)

if(BUILD_MPC)
    add_subdirectory(mpc)
endif()

if(BUILD_COLLISION)
    add_subdirectory(collision_detection)
endif()

