cmake_minimum_required(VERSION 3.14)
project(NumberAnalysis C)

set(CMAKE_C_STANDARD 99)

add_executable(NumberAnalysis main.c)
add_executable(final final/main.c final/csr_sparse.h)