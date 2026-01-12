cmake_minimum_required(VERSION 3.16)
project(dsd_parallel LANGUAGES CXX)

# ========================
# C++ standard
# ========================
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ========================
# Default build type
# ========================
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()

# ========================
# OpenMP support
# ========================
# On macOS with Apple Clang, we need to manually set OpenMP paths
if(APPLE)
    set(OPENMP_INCLUDE_DIR "/opt/homebrew/opt/libomp/include")
    set(OPENMP_LIB_DIR "/opt/homebrew/opt/libomp/lib")
    # Apple Clang requires -Xpreprocessor -fopenmp
    set(OpenMP_CXX_FLAGS "-Xpreprocessor;-fopenmp")
    set(OpenMP_CXX_LIB_NAMES "omp")
    set(OpenMP_omp_LIBRARY "${OPENMP_LIB_DIR}/libomp.dylib")
    include_directories(${OPENMP_INCLUDE_DIR})
    link_directories(${OPENMP_LIB_DIR})
    message(STATUS "OpenMP configured for macOS: ${OPENMP_INCLUDE_DIR}")
else()
    find_package(OpenMP REQUIRED)
    if(OpenMP_CXX_FOUND)
        message(STATUS "OpenMP found")
        message(STATUS "OpenMP C++ flags: ${OpenMP_CXX_FLAGS}")
    else()
        message(FATAL_ERROR "OpenMP not found")
    endif()
endif()

# ========================
# Include headers
# ========================
include_directories(inc)

# ========================
# Executable
# ========================
add_executable(dsd_parallel
    main.cpp
    src/graph.cpp
    src/common.cpp
    src/helpers.cpp
)

# ========================
# Link OpenMP
# ========================
if(APPLE)
    target_link_libraries(dsd_parallel PUBLIC omp)
    target_compile_options(dsd_parallel PRIVATE ${OpenMP_CXX_FLAGS})
else()
    target_link_libraries(dsd_parallel PUBLIC OpenMP::OpenMP_CXX)
endif()

# ========================
# Compile flags
# ========================
target_compile_options(dsd_parallel PRIVATE
    -Wall
    -g
    -O0
)
