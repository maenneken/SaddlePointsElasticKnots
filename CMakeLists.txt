cmake_minimum_required(VERSION 3.20)
project(KnotSimulation LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Optional: default build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# -------------------------------
# 1. Add 3rdparty libraries
# -------------------------------
#TODO coment out Python Bindings in CMakefile of ElasticRods
add_subdirectory(3rdparty/polyscope)

if(NOT MSVC)
    add_compile_options(-Wno-error=changes-meaning)
endif()
set(LIBIGL_WITH_OPENGL_GLFW OFF CACHE BOOL "" FORCE)
set(LIBIGL_WITH_OPENGL_GLFW_IMGUI OFF CACHE BOOL "" FORCE)

set(MESHFEM_WITH_PYTHON OFF CACHE BOOL "" FORCE)
set(MESHFEM_WITH_IGL OFF CACHE BOOL "" FORCE)       # optional, but avoids extra issues
set(ELASTICROD_WITH_PYTHON OFF CACHE BOOL "" FORCE)
set(MESHFEM_WITH_PYTHON OFF CACHE BOOL "" FORCE)
set(MESHFEM_BUILD_PYTHON_BINDINGS OFF CACHE BOOL "" FORCE)
set(PYBIND11_FINDPYTHON OFF CACHE BOOL "" FORCE)
set(PYBIND11_INSTALL OFF CACHE BOOL "" FORCE)



add_subdirectory(3rdparty/ElasticKnots/3rdparty/ElasticRods)

target_include_directories(ElasticRods PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/ElasticKnots/3rdparty/ElasticRods/3rdparty/libigl/external
    ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/ElasticKnots/3rdparty/ElasticRods/3rdparty/libigl/external/imgui
)


#add_subdirectory(3rdparty/ElasticKnots/3rdparty/ElasticRods/3rdparty/MeshFEM)

#target_include_directories(ElasticRods PRIVATE 3rdparty/ElasticKnots/3rdparty/ElasticRods/3rdparty/libigl/external)




# -------------------------------
# 3. Create main executable
# -------------------------------
add_executable(test1
    test1.cpp
    helpers.cpp
    KnotVisualizer.cpp
)

# Suppress warnings
if(MSVC)
    target_compile_options(test1 PRIVATE /W0)
else()
    target_compile_options(test1 PRIVATE -w)
endif()

# -------------------------------
# 4. Link libraries
# -------------------------------

target_link_libraries(test1 PRIVATE polyscope ElasticRods MeshFEM)

