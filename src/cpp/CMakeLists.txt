cmake_minimum_required(VERSION 3.0...3.5)

project(MFPT_CPP CXX)

add_executable(mfpt
    mfpt.cpp
)

install(TARGETS mfpt DESTINATION cpp)
