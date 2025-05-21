cmake_minimum_required(VERSION 3.10)
project(RTProject LANGUAGES CXX)

# Usar clang++ de Homebrew (LLVM)
set(CMAKE_CXX_COMPILER "/opt/homebrew/opt/llvm/bin/clang++")

# Detectar y usar el SDK de macOS actual para evitar errores con 'cmath'
execute_process(
        COMMAND xcrun --sdk macosx --show-sdk-path
        OUTPUT_VARIABLE MACOSX_SDK
        OUTPUT_STRIP_TRAILING_WHITESPACE
)
set(CMAKE_OSX_SYSROOT ${MACOSX_SDK})

# Estándar de C++ y optimizaciones
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "-Xclang -fopenmp -mcpu=apple-m1 -march=native -O3")

# Linkeo con OpenMP de Homebrew
set(CMAKE_EXE_LINKER_FLAGS "-L/opt/homebrew/opt/libomp/lib -lomp")

# Rutas de include y lib para LLVM y libomp (Homebrew)
include_directories(/opt/homebrew/opt/libomp/include)
include_directories(/opt/homebrew/opt/llvm/include)
link_directories(/opt/homebrew/opt/llvm/lib)

# Incluir tu código
include_directories(${CMAKE_SOURCE_DIR}/include)

# Ejecutable
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

# Limpieza
set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "rt")
