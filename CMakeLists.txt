################################################################################
# General Informations
################################################################################

cmake_minimum_required(VERSION 3.1)
project(SaddlePointsElasticKnots)

# CMP0063: Honor visibility properties for all target types.
if(POLICY CMP0063)
    cmake_policy(SET CMP0063 NEW)
endif()

################################################################################

set(THIRD_PARTY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/ElasticKnots/3rdparty)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# set(ALLOW_DUPLICATE_CUSTOM_TARGETS TRUE)  # eliminates 'add_library cannot create target' error, but gives The install of the tbb target requires changing an RPATH from the build tree, but this is not supported with the Ninja generator unless on an ELF-based or XCOFF-based platform

set(CMAKE_CXX_FLAGS_RELWITHASSERT        "-O3" CACHE STRING "Flags used during RelWithAssertions builds" FORCE)
set(CMAKE_C_FLAGS_RELWITHASSERT          "-O3" CACHE STRING "Flags used during RelWithAssertions builds" FORCE)
set(CMAKE_EXE_LINKER_FLAGS_RELWITHASSERT    "" CACHE STRING "Flags used during RelWithAssertions builds" FORCE)
set(CMAKE_SHARED_LINKER_FLAGS_RELWITHASSERT "" CACHE STRING "Flags used during RelWithAssertions builds" FORCE)
set(CMAKE_MODULE_LINKER_FLAGS_RELWITHASSERT "" CACHE STRING "Flags used during RelWithAssertions builds" FORCE)
set(CMAKE_STATIC_LINKER_FLAGS_RELWITHASSERT "" CACHE STRING "Flags used during RelWithAssertions builds" FORCE)
mark_as_advanced(
    CMAKE_CXX_FLAGS_RELWITHASSERT
    CMAKE_C_FLAGS_RELWITHASSERT
    CMAKE_EXE_LINKER_FLAGS_RELWITHASSERT
    CMAKE_SHARED_LINKER_FLAGS_RELWITHASSERT
    CMAKE_MODULE_LINKER_FLAGS_RELWITHASSERT
    CMAKE_STATIC_LINKER_FLAGS_RELWITHASSERT
    )

set(CMAKE_CXX_FLAGS_RELEASENATIVE        "-O3 -march=native -DNDEBUG" CACHE STRING "Flags used during native release builds" FORCE)
set(CMAKE_C_FLAGS_RELEASENATIVE          "-O3 -march=native -DNDEBUG" CACHE STRING "Flags used during native release builds" FORCE)
set(CMAKE_EXE_LINKER_FLAGS_RELEASENATIVE    "" CACHE STRING "Flags used during native release builds" FORCE)
set(CMAKE_SHARED_LINKER_FLAGS_RELEASENATIVE "" CACHE STRING "Flags used during native release builds" FORCE)
set(CMAKE_MODULE_LINKER_FLAGS_RELEASENATIVE "" CACHE STRING "Flags used during native release builds" FORCE)
set(CMAKE_STATIC_LINKER_FLAGS_RELEASENATIVE "" CACHE STRING "Flags used during native release builds" FORCE)
mark_as_advanced(
    CMAKE_CXX_FLAGS_RELEASENATIVE
    CMAKE_C_FLAGS_RELEASENATIVE
    CMAKE_EXE_LINKER_FLAGS_RELEASENATIVE
    CMAKE_SHARED_LINKER_FLAGS_RELEASENATIVE
    CMAKE_MODULE_LINKER_FLAGS_RELEASENATIVE
    CMAKE_STATIC_LINKER_FLAGS_RELEASENATIVE
    )

set(CMAKE_CXX_FLAGS_DEBUGNATIVE        "-g -march=native" CACHE STRING "Flags used during native debug builds" FORCE)
set(CMAKE_C_FLAGS_DEBUGNATIVE          "-g -march=native" CACHE STRING "Flags used during native debug builds" FORCE)
set(CMAKE_EXE_LINKER_FLAGS_DEBUGNATIVE    "" CACHE STRING "Flags used during native debug builds" FORCE)
set(CMAKE_SHARED_LINKER_FLAGS_DEBUGNATIVE "" CACHE STRING "Flags used during native debug builds" FORCE)
set(CMAKE_MODULE_LINKER_FLAGS_DEBUGNATIVE "" CACHE STRING "Flags used during native debug builds" FORCE)
set(CMAKE_STATIC_LINKER_FLAGS_DEBUGNATIVE "" CACHE STRING "Flags used during native debug builds" FORCE)
mark_as_advanced(
    CMAKE_CXX_FLAGS_DEBUGNATIVE
    CMAKE_C_FLAGS_DEBUGNATIVE
    CMAKE_EXE_LINKER_FLAGS_DEBUGNATIVE
    CMAKE_SHARED_LINKER_FLAGS_DEBUGNATIVE
    CMAKE_MODULE_LINKER_FLAGS_DEBUGNATIVE
    CMAKE_STATIC_LINKER_FLAGS_DEBUGNATIVE
    )


if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to RelWithAssert")
    set(CMAKE_BUILD_TYPE "RelWithAssert")
endif()

set(CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}" CACHE STRING "Choose the type of build. Options are: None Debug Release RelWithDebInfo MinSizeRel RelWithAssert" FORCE)

# Make sure warnings/errors are still colorized when using Ninja for building.
if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    add_definitions(-fdiagnostics-color=always)
endif()

# Export compile flags(used for autocompletion of the C++ code)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

# CMake plugin for vscode
include(CMakeToolsHelpers OPTIONAL)

# Enable more warnings
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -Wno-comment")

# We need -fPIC when compiling our libraries and our dependencies for
# the python bindings to link.
set(CMAKE_POSITION_INDEPENDENT_CODE ON)



################################################################################
# Python Bindings
################################################################################
#add_subdirectory(python_bindings)

################################################################################
# More binaries
################################################################################
# add_subdirectory(tests)  # TODO: write unit tests
#[[
#Disable all python bindings

set(MESHFEM_WITH_PYTHON OFF CACHE BOOL "" FORCE)
set(MESHFEM_WITH_IGL OFF CACHE BOOL "" FORCE)       # optional, but avoids extra issues
set(ELASTICROD_WITH_PYTHON OFF CACHE BOOL "" FORCE)
set(MESHFEM_WITH_PYTHON OFF CACHE BOOL "" FORCE)
set(MESHFEM_BUILD_PYTHON_BINDINGS OFF CACHE BOOL "" FORCE)
set(PYBIND11_FINDPYTHON OFF CACHE BOOL "" FORCE)
set(PYBIND11_INSTALL OFF CACHE BOOL "" FORCE)
]]


set(LIBIGL_WITH_OPENGL_GLFW OFF CACHE BOOL "" FORCE)
set(LIBIGL_WITH_OPENGL_GLFW_IMGUI OFF CACHE BOOL "" FORCE)
#set(POLYSCOPE_USE_SYSTEM_GLAD OFF CACHE BOOL "Use system glad")
add_subdirectory(3rdparty/polyscope) #Polyscope

add_subdirectory(3rdparty/ElasticKnots)



if(NOT MSVC)
    add_compile_options(-Wno-error=changes-meaning)
endif()


# MeshFEM target
target_include_directories(ElasticKnots
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/ElasticKnots/3rdparty/ElasticRods/3rdparty/MeshFEM/3rdparty/imgui
    ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/ElasticKnots/3rdparty/ElasticRods/3rdparty/MeshFEM/3rdparty
    ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/ElasticKnots/3rdparty/ElasticRods/3rdparty/libigl/include/igl/opengl
    ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/ElasticKnots/3rdparty/ElasticRods/3rdparty/libigl/include/igl
)

# Polyscope target
target_include_directories(polyscope
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/polyscope/deps/imgui
    ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/polyscope/deps/imgui/imgui
    ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/polyscope/deps/imgui/implot
    ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/polyscope/src/render/opengl
    ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/polyscope/src/render
    ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/polyscope/include/polyscope/render/opengl
    ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/polyscope/include/polyscope/render
)

add_executable(test1
    test1.cpp
    helpers.cpp
    KnotVisualizer.cpp
)

add_executable(FindSaddle
    FindSaddlePoint.cpp
    helpers.cpp
    KnotVisualizer.cpp
)


target_link_libraries(test1 PRIVATE polyscope ElasticKnots)
target_link_libraries(FindSaddle PRIVATE polyscope ElasticKnots)

