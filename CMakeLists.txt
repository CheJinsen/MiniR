cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

project(main)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)

add_compile_options(-Wall -g -O2)

add_subdirectory(src)