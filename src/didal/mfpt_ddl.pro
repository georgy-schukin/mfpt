TEMPLATE = app
CONFIG += console c++14
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXX = mpicxx
QMAKE_CC = mpicc
QMAKE_LINK = mpicxx

TARGET = mfpt_ddl

SOURCES +=   \
    mfpt_ddl.cpp \
    output.cpp

HEADERS +=   \
    common.h \
    output.h \
    shadowed_array2d.h

exists(local.pri) {
    include(local.pri)
}

