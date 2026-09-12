TEMPLATE = app
CONFIG += console c++14
CONFIG -= app_bundle
CONFIG -= qt

TARGET = mfpt2

SOURCES +=   \
    mfpt2.cpp

HEADERS +=   \
    ../common/output.h \
    ../common/timer.h \
    array2d.h \
    defs.h
