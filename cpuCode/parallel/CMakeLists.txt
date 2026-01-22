cmake_minimum_required(VERSION 3.16)
project(dsd_parallel LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set default build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Find OpenMP
find_package(OpenMP REQUIRED)

# Executable
add_executable(dsd_parallel
    main.cpp
    src/graph.cpp
    src/common.cpp
    src/helpers.cpp
)

# Link OpenMP
target_link_libraries(dsd_parallel PRIVATE OpenMP::OpenMP_CXX)

# Optimizations
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    target_compile_options(dsd_parallel PRIVATE
        -O3
        -march=native
        -DNDEBUG
    )
    target_link_options(dsd_parallel PRIVATE
        -flto
    )
else()  # Debug
    target_compile_options(dsd_parallel PRIVATE
        -O0
        -g
        -Wall
        -Wextra
        -Wshadow
        -Wconversion
        -fsanitize=address
        -fsanitize=undefined
    )
endif()