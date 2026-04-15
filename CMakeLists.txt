cmake_minimum_required(VERSION 3.20)
project(EF LANGUAGES CXX)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


add_executable(ef_test test.cpp)

target_include_directories(ef_test PRIVATE ${CMAKE_SOURCE_DIR}/include)

target_compile_options(ef_test PRIVATE
    -O3
    -Wall
    -Wextra
    -march=native
)