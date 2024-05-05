cmake_minimum_required(VERSION 3.27)
project(1_solve_ss)

set(CMAKE_CXX_STANDARD 11)

file(GLOB _SRCS "*.cpp")
message(STATUS "srcs: ${_SRCS}")

add_executable(${PROJECT_NAME} ${_SRCS})
