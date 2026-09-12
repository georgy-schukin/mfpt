cmake_minimum_required(VERSION 3.0...3.5)

project(MFPT2_F Fortran)

add_executable(mfpt2_f
    brbz003f.f
)

install(TARGETS mfpt2_f DESTINATION fortran2)
