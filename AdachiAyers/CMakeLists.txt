cmake_minimum_required (VERSION 3.13)
set(PROYECTO "AdachiAyers")
project (${PROYECTO})

add_library(${PROYECTO} MODULE ${PROYECTO}.cpp)

if (SUPERNOVA)
    add_library(${PROYECTO}_supernova MODULE ${PROYECTO}.cpp)
    set_property(TARGET ${PROYECTO}_supernova APPEND PROPERTY COMPILE_DEFINITIONS  "SUPERNOVA")
endif()

install(DIRECTORY ${CMAKE_SOURCE_DIR}/${PROYECTO}/impulses/ DESTINATION ${INSTALL_DESTINATION_DISTRO}/${PROYECTO}/impulses/ PATTERN "*")