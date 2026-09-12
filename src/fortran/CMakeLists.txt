cmake_minimum_required(VERSION 3.0...3.5)

project(MFPT_F Fortran)

add_executable(mfpt_f
    brbz003c.f
)

add_executable(mfpt2_f
    brbz003f.f
)

install(TARGETS mfpt_f DESTINATION fortran)
install(TARGETS mfpt2_f DESTINATION fortran2)
