cmake_minimum_required (VERSION 2.8)
project (Karplus)

include_directories(${SC_PATH}/include/plugin_interface)
include_directories(${SC_PATH}/include/common)
include_directories(${SC_PATH}/external_libraries/libsndfile/)


set(CMAKE_SHARED_MODULE_PREFIX "")
if(APPLE OR WIN32)
set(CMAKE_SHARED_MODULE_SUFFIX ".scx")
endif()

add_library(Karplus MODULE Karplus.cpp)

add_library(Karplus_supernova MODULE Karplus.cpp)
set_property(TARGET Karplus_supernova APPEND PROPERTY COMPILE_DEFINITIONS  "SUPERNOVA")
include_directories(${SC_PATH}/external_libraries
					${SC_PATH}/external_libraries/boost-lockfree
					${SC_PATH}/external_libraries/nova-simd
					${SC_PATH}/external_libraries/nova-tt
					)
include_directories(${SC_PATH}/external_libraries/boost)
include_directories(${SC_PATH}/external_libraries/nova-simd)
add_definitions("-DNOVA_SIMD")

INSTALL(TARGETS Karplus Karplus_supernova LIBRARY DESTINATION "C:/SupercolliderRepos/plugins/mingw48/source/")