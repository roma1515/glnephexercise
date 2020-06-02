TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        src/diffs.cpp \
        src/glnsvpos.cpp \
        src/main.cpp \
        src/rungekutta.cpp

HEADERS += \
    include/libglnsvpos/diffs.h \
    include/libglnsvpos/glnsvpos.h \
    include/libglnsvpos/rungekutta.h
