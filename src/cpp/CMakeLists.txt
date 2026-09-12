cmake_minimum_required(VERSION 3.0...3.5)

project(MFPT_CPP CXX)

add_executable(mfpt
    mfpt.cpp
)

add_executable(mfpt2
    mfpt2.cpp
)

install(TARGETS mfpt DESTINATION cpp)
install(TARGETS mfpt2 DESTINATION cpp2)
