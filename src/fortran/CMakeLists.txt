cmake_minimum_required(VERSION 3.0...3.5)

project(MFPT_F Fortran)

add_executable(mfpt_f
    brbz003c.f
)

install(TARGETS mfpt_f DESTINATION fortran)
