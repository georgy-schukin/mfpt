TEMPLATE = app
CONFIG += console c++14
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXX = mpicxx
QMAKE_CC = mpicc
QMAKE_LINK = mpicxx

TARGET = mfpt2_mpi

exists(local.pri) {
    include(local.pri)
}

SOURCES +=   \
    block_decomp.cpp \
    comm_util.cpp \
    distributed_array2d.cpp \
    mfpt2_mpi.cpp

HEADERS +=   \
    ../common/output.h \
    ../common/timer.h \
    block_decomp.h \
    comm_util.h \
    defs.h \
    distributed_array2d.h \
    shadowed_array2d.h
