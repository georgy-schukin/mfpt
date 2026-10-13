cmake_minimum_required(VERSION 3.10)

project(MFPT LANGUAGES CXX Fortran)

add_subdirectory(fortran)
add_subdirectory(cpp)
add_subdirectory(mpi)