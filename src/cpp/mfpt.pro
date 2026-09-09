TEMPLATE = app
CONFIG += console c++14
CONFIG -= app_bundle
CONFIG -= qt

TARGET = mfpt

SOURCES +=   \
    mfpt.cpp \
    output.cpp

HEADERS +=   \
    array2d.h \
    defs.h \
    output.h \
    timer.h
