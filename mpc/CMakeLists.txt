cmake_minimum_required(VERSION 3.16)
project(wheelleg_mpc_module LANGUAGES CXX)
include(CTest)

# 只保留原 H7 的 C++ 实现；测试由 mcu 子目录注册。
add_subdirectory(mcu)
