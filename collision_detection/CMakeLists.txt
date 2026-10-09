cmake_minimum_required(VERSION 3.16)
project(wheelleg_collision_observer LANGUAGES C)

add_library(wheelleg_collision STATIC
    mcu/collision_detect.c mcu/wheelleg_mcg.c)
target_include_directories(wheelleg_collision PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/mcu)
target_compile_features(wheelleg_collision PUBLIC c_std_99)
if(NOT MSVC)
    target_compile_options(wheelleg_collision PRIVATE -Wall -Wextra -Wpedantic)
    target_link_libraries(wheelleg_collision PUBLIC m)
endif()

include(CTest)
if(BUILD_TESTING)
    add_executable(test_collision tests/test_collision.c)
    target_link_libraries(test_collision PRIVATE wheelleg_collision)
    add_test(NAME collision_observer COMMAND test_collision)
endif()
