cmake_minimum_required(VERSION 3.16)
project(dsd LANGUAGES CXX)

# ========================
# C++ Standard
# ========================
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ========================
# Default to Release
# ========================
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# ========================
# Source files
# ========================
add_executable(dsd
    main.cpp
    src/graph.cpp
    src/common.cpp
    src/helpers.cpp
)

# ========================
# Include directories
# ========================
target_include_directories(dsd PRIVATE inc)

# ========================
# Compile options (PRODUCTION)
# ========================
target_compile_options(dsd PRIVATE
    -O3
    -march=native
    -ffast-math
    -funroll-loops
    -DNDEBUG
    -Wall
    -Wextra
)

# ========================
# Link options
# ========================
target_link_options(dsd PRIVATE
    -O3
    -march=native
    -flto
)

# ========================
# Interprocedural optimization (LTO)
# ========================
include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_supported OUTPUT error)
if(ipo_supported)
    set_property(TARGET dsd PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()
