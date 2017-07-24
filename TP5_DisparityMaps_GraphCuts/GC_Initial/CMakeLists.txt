PROJECT(Exo4)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} $ENV{IMAGINEPP_ROOT}/CMake)
FIND_PACKAGE(Imagine)

ImagineAddExecutable(exampleGC exampleGC.cpp maxflow/graph.cpp)
ImagineUseModules(exampleGC Images LinAlg Common Graphics)

ImagineAddExecutable(stereoGCtemplate stereoGCtemplate.cpp maxflow/graph.cpp)
ImagineUseModules(stereoGCtemplate Images)
