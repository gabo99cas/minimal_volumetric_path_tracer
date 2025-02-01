cmake_minimum_required(VERSION 3.10)
project(RTProject LANGUAGES CXX)

# Set the C++ compiler to clang++
set(CMAKE_CXX_COMPILER "clang++")

# Set C++ standard and other compile flags
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "-Xclang -fopenmp -mcpu=apple-m1 -march=native -O3")
set(CMAKE_EXE_LINKER_FLAGS "-L/opt/homebrew/opt/libomp/lib -lomp")

# Include directories for OpenMP
include_directories(/opt/homebrew/opt/libomp/include)

# Add the includes directory to the search path
include_directories(${CMAKE_SOURCE_DIR}/include)

# Add executable
add_executable(rt src/rt.cpp
        include/Vector.cpp
        include/Vector.h
        include/Ray.cpp
        include/Ray.h
        include/Sphere.cpp
        include/Sphere.h
        include/mathUtilities.h
        include/pathTracingUtilities.h
        include/samplingFunctions.h
        include/microFacetUtilities.h
        include/misSamplingFunctions.h
        include/shadeMethods.h
        include/rayMarchingUtilities.h
        include/rayMarchingMethods.h
        include/volumetricBasicFunctions.h
        include/vptShadeMethods.h
        include/vptSamplingFunctions.h
)

# Clean command (if using custom scripts for cleaning)
set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "rt")
