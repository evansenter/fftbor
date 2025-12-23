cmake_minimum_required(VERSION 3.16)

project(FFTBOR1D
    VERSION 1.0.0
    LANGUAGES CXX
    DESCRIPTION "FFT-based RNA partition function calculator"
)

# Options
option(FFTBOR_BUILD_TESTS "Build the test suite" ON)
option(FFTBOR_USE_OPENMP "Use OpenMP for parallelization" ON)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)

# Add custom cmake modules path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake_modules")

# Compiler warnings
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall -Wextra -Wno-write-strings)
    if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
        add_compile_options(-O3 -funroll-loops -ffast-math)
    endif()
endif()

# Find dependencies
find_package(FFTW REQUIRED)

if(FFTBOR_USE_OPENMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        message(STATUS "OpenMP found, enabling parallel execution")
    else()
        message(WARNING "OpenMP not found, parallel execution disabled")
    endif()
endif()

# Define the library with core functionality (for testing)
add_library(fftbor_core STATIC
    delta.cpp
    misc.cpp
    parameter_parser.cpp
    globals.cpp
)

target_include_directories(fftbor_core PUBLIC
    ${CMAKE_SOURCE_DIR}
    ${FFTW_INCLUDES}
)

target_link_libraries(fftbor_core PUBLIC
    ${FFTW_LIBRARIES}
    m
)

if(OpenMP_CXX_FOUND)
    target_link_libraries(fftbor_core PUBLIC OpenMP::OpenMP_CXX)
endif()

# Main executable
add_executable(FFTbor main.cpp)
target_link_libraries(FFTbor PRIVATE fftbor_core)

# Installation
install(TARGETS FFTbor DESTINATION bin)
install(FILES rna_turner1999.par rna_turner2004.par DESTINATION bin)

# Testing
if(FFTBOR_BUILD_TESTS)
    enable_testing()

    # Fetch GoogleTest
    include(FetchContent)
    FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG v1.14.0
    )
    # For Windows: Prevent overriding the parent project's compiler/linker settings
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)

    # Add test executable
    add_executable(fftbor_tests
        tests/test_main.cpp
        tests/test_parameter_parser.cpp
        tests/test_misc.cpp
        tests/test_energy.cpp
        tests/test_integration.cpp
    )

    target_link_libraries(fftbor_tests PRIVATE
        fftbor_core
        GTest::gtest
        GTest::gtest_main
    )

    target_include_directories(fftbor_tests PRIVATE
        ${CMAKE_SOURCE_DIR}
    )

    # Copy parameter files to bin directory for tests
    file(COPY ${CMAKE_SOURCE_DIR}/rna_turner2004.par DESTINATION ${CMAKE_SOURCE_DIR}/bin)
    file(COPY ${CMAKE_SOURCE_DIR}/rna_turner1999.par DESTINATION ${CMAKE_SOURCE_DIR}/bin)

    include(GoogleTest)
    gtest_discover_tests(fftbor_tests WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
endif()
