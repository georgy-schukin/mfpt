TEMPLATE = app
CONFIG += console c++14
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXX = mpicxx
QMAKE_CC = mpicc
QMAKE_LINK = mpicxx

TARGET = mfpt_mpi

SOURCES +=   \
    block_decomp.cpp \
    comm_util.cpp \
    mfpt_mpi.cpp \
    output.cpp

HEADERS +=   \
    block_decomp.h \
    comm_util.h \
    common.h \
    output.h \
    shadowed_array2d.h
