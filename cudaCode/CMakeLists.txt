cmake_minimum_required(VERSION 3.10)

# Set the project name
project(CDS LANGUAGES CXX CUDA)

# Set the CUDA architecture
set(CMAKE_CUDA_ARCHITECTURES 75 80 86)

# Add executable
add_executable(main 
    src/graph.cpp 
    src/common.cpp 
    src/gpuMemoryAllocation.cu 
    src/helpers.cu 
    main.cu
)

# Enable extended lambda (if needed)
target_compile_options(main PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:--extended-lambda>)

# Additional include directories (if needed)
target_include_directories(main PRIVATE ${CMAKE_SOURCE_DIR}/src)

# Link libraries (if needed)
# target_link_libraries(main PRIVATE some_library)