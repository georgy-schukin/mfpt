cmake_minimum_required(VERSION 3.0...3.5)

project(MFPT_MPI CXX)

find_package(MPI REQUIRED COMPONENTS CXX)

add_executable(mfpt_mpi
    block_decomp.cpp
    distributed_array2d.cpp
    comm_util.cpp
    mfpt_mpi.cpp
)

add_executable(mfpt2_mpi
    block_decomp.cpp
    distributed_array2d.cpp
    comm_util.cpp
    mfpt2_mpi.cpp
)

target_link_libraries(mfpt_mpi PRIVATE MPI::MPI_CXX)
target_link_libraries(mfpt2_mpi PRIVATE MPI::MPI_CXX)

install(TARGETS mfpt_mpi DESTINATION mpi)
install(TARGETS mfpt2_mpi DESTINATION mpi2)