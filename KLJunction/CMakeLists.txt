cmake_minimum_required (VERSION 2.8)
set(PROYECTO "KLJunction")
project (${PROYECTO})

add_library(${PROYECTO} MODULE ${PROYECTO}.cpp)

if (SUPERNOVA)
    add_library(${PROYECTO}_supernova MODULE ${PROYECTO}.cpp)
    set_property(TARGET ${PROYECTO}_supernova APPEND PROPERTY COMPILE_DEFINITIONS  "SUPERNOVA")
endif()