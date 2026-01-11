cmake_minimum_required(VERSION 3.16)
project(dsd LANGUAGES CXX)

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
# Include headers
# ========================
include_directories(inc)

# ========================
# Executable
# ========================
add_executable(dsd
    main.cpp
    src/graph.cpp
    src/common.cpp
    src/helpers.cpp
)

# ========================
# Compile flags
# ========================
target_compile_options(dsd PRIVATE
    -Wall
    -g
    -O0
)
