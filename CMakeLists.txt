cmake_minimum_required (VERSION 3.13) #for subdirectory targets
project (MyUGens)

include_directories(${SC_PATH}/include/plugin_interface)
include_directories(${SC_PATH}/include/common)
include_directories(${SC_PATH}/common)
include_directories(${SC_PATH}/external_libraries/libsndfile/)

include_directories(${SC_PATH}/external_libraries/nova-simd)
add_definitions(-DNOVA_SIMD)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_SHARED_MODULE_PREFIX "")
if(APPLE OR WIN32)
	set(CMAKE_SHARED_MODULE_SUFFIX ".scx")
endif()
if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
	set(CMAKE_COMPILER_IS_CLANG 1)
endif()
if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_CLANG)
    #add_definitions(-fvisibility=hidden)

    include (CheckCCompilerFlag)
    include (CheckCXXCompilerFlag)

    CHECK_C_COMPILER_FLAG(-msse HAS_SSE)
    CHECK_CXX_COMPILER_FLAG(-msse HAS_CXX_SSE)

    if (HAS_SSE)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse")
    endif()
    if (HAS_CXX_SSE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse")
    endif()

    CHECK_C_COMPILER_FLAG(-msse2 HAS_SSE2)
    CHECK_CXX_COMPILER_FLAG(-msse2 HAS_CXX_SSE2)

    if (HAS_SSE2)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse2")
    endif()
    if (HAS_CXX_SSE2)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse2")
    endif()

    CHECK_C_COMPILER_FLAG(-mfpmath=sse HAS_FPMATH_SSE)
    CHECK_CXX_COMPILER_FLAG(-mfpmath=sse HAS_CXX_FPMATH_SSE)

    if (HAS_FPMATH_SSE)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mfpmath=sse")
    endif()
    if (HAS_CXX_FPMATH_SSE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mfpmath=sse")
    endif()


	#if(CPP11)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
		if(CMAKE_COMPILER_IS_CLANG)
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
		endif()
	#endif()
endif()
if(MINGW)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mstackrealign")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mstackrealign")
endif()


if (SUPERNOVA)
    include_directories(${SC_PATH}/external_libraries/nova-tt)
    # actually just boost.atomic
    include_directories(${SC_PATH}/external_libraries/boost)
    include_directories(${SC_PATH}/external_libraries/boost_lockfree)
    include_directories(${SC_PATH}/external_libraries/boost-lockfree)
endif()


#install paths
if (CMAKE_SYSTEM_NAME MATCHES "Linux")
    set(INSTALL_DESTINATION "lib/SuperCollider/plugins")
    set(INSTALL_DESTINATION_DISTRO "share/SuperCollider/Extensions/Myplugins")
else()
    set(INSTALL_DESTINATION "Myplugins")
    set(INSTALL_DESTINATION_DISTRO "Myplugins")
endif()

include_directories(DWGReverb/dwglib/)

include(cmake/functions.cmake)

add_subdirectory(DWGClarinet)
InstallMyUgens(DWGClarinet)

add_subdirectory(DWGFlute)
InstallMyUgens(DWGFlute)

add_subdirectory(KLJunction)
InstallMyUgens(KLJunction)

add_subdirectory(SonLPC)
InstallMyUgens(SonLPC)

add_subdirectory(Karplus)
InstallMyUgens(Karplus)

add_subdirectory(IIRf)
InstallMyUgens(IIRf)

add_subdirectory(AdachiAyers)
InstallMyUgens(AdachiAyers)

add_subdirectory(PitchTracker)
InstallMyUgens(PitchTracker)

add_subdirectory(PluckSynth)
InstallMyUgens(PluckSynth)

add_subdirectory(MembraneV)
InstallMyUgens(MembraneV)

add_subdirectory(MyPlucked)
InstallMyUgens(MyPlucked)

##add_subdirectory(MyBowed)
##InstallMyUgens(MyBowed)

add_subdirectory(DWGReverb)