TEMPLATE = app
CONFIG += console c++14
CONFIG -= app_bundle
CONFIG -= qt

TARGET = mfpt

SOURCES +=   \
    mfpt.cpp

HEADERS +=   \
    ../common/output.h \
    ../common/timer.h \
    array2d.h \
    defs.h
