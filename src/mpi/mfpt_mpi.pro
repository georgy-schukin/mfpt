TEMPLATE = app
CONFIG += console c++14
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXX = mpicxx
QMAKE_CC = mpicc
QMAKE_LINK = mpicxx

TARGET = mfpt

SOURCES +=   \
    block_decomp.cpp \
    mfpt_mpi.cpp

HEADERS +=   \
    array2d.h \
    block_decomp.h
