CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
FILE(TO_CMAKE_PATH "$ENV{IMAGINEPP_ROOT}" d)
IF(NOT EXISTS "${d}")
    MESSAGE(FATAL_ERROR "Error: environment variable IMAGINEPP_ROOT=" "${d}")
ENDIF(NOT EXISTS "${d}")
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${d}/CMake")
FIND_PACKAGE(Imagine)

PROJECT(Panorama)
add_executable(Fundamental
        Fundamental.cpp
        Imagine/SIFT_VL.cpp
        Imagine/vl/generic.c Imagine/vl/host.c Imagine/vl/imop.c Imagine/vl/sift.c)
ImagineUseModules(Fundamental LinAlg Graphics Images)
