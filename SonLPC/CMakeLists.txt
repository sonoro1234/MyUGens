cmake_minimum_required (VERSION 2.8)
set(PROYECTO "SonLPC")
project (${PROYECTO})

add_library(${PROYECTO} MODULE ${PROYECTO}.cpp)

if(APPLE)
    target_link_libraries(${PROYECTO} "-framework Accelerate")
endif()

if (SUPERNOVA)
    add_library(${PROYECTO}_supernova MODULE ${PROYECTO}.cpp)
    set_property(TARGET ${PROYECTO}_supernova APPEND PROPERTY COMPILE_DEFINITIONS  "SUPERNOVA")
    if(APPLE)
        target_link_libraries(${PROYECTO}_supernova "-framework Accelerate")
    endif()
endif()