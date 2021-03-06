TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += c++11
QMAKE_CXXFLAGS += -g -std=c++11 -O3 -ipo

SOURCES += main.cpp \
    mcintegrator.cpp \
    vec3.cpp \
    random.cpp \
    mcresult.cpp \
    wavefunction.cpp

include(deployment.pri)
qtcAddDeployment()

HEADERS += \
    mcintegrator.h \
    vec3.h \
    random.h \
    mcresult.h \
    wavefunction.h \
    mcintegratorparameters.h

