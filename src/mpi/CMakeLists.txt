cmake_minimum_required(VERSION 3.0...3.5)

project(MFPT_MPI CXX)

find_package(MPI REQUIRED COMPONENTS CXX)

add_executable(mfpt_mpi
    block_decomp.cpp
    comm_util.cpp
    mfpt_mpi.cpp
    output.cpp
)

target_link_libraries(mfpt_mpi PRIVATE MPI::MPI_CXX)

install(TARGETS mfpt_mpi DESTINATION mpi)