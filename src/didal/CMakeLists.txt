cmake_minimum_required(VERSION 3.0...3.5)

project(MFPT_DDL CXX)

find_package(MPI REQUIRED COMPONENTS CXX)

add_executable(mfpt_ddl
    mfpt_ddl.cpp
    output.cpp
)

target_link_libraries(mfpt_ddl PRIVATE MPI::MPI_CXX)
target_link_libraries(mfpt_ddl didal)

install(TARGETS mfpt_ddl DESTINATION didal)