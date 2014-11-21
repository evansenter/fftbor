include(CheckCXXCompilerFlag)

cmake_minimum_required(VERSION 2.6)

if(POLICY CMP0042)
  cmake_policy(SET CMP0042 OLD)
endif()

project(FFTBOR1D)

if (CMAKE_BUILD_TYPE)
  string(TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE_UPPER)
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake_modules")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${FFTBOR1D_SOURCE_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${FFTBOR1D_SOURCE_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${FFTBOR1D_SOURCE_DIR}/lib)

macro(cxx_compiler_flag_test flag_name flag fallback_flag required)
  check_cxx_compiler_flag(${flag} ${flag_name})
  if("${flag_name}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flag}")
  elseif(NOT "${fallback_flag}" STREQUAL "")
    cxx_compiler_flag_test(${flag_name}_FALLBACK ${fallback_flag} "" ${required})
  elseif("${required}")
    message(FATAL_ERROR "${flag} support is required by ${CMAKE_CXX_COMPILER}")
  endif()
endmacro(cxx_compiler_flag_test)

cxx_compiler_flag_test(CXX_STD_CXX98 "-std=c++98" "" 1)
cxx_compiler_flag_test(CXX_PEDANTIC "-pedantic" "" 0)
cxx_compiler_flag_test(CXX_WALL "-Wall" "" 0)
cxx_compiler_flag_test(CXX_WEXTRA "-Wextra" "" 0)
cxx_compiler_flag_test(CXX_WNO_WRITE_STR "-Wno-write-strings" "" 0)
cxx_compiler_flag_test(CXX_WVLA "-Wno-vla" "" 0)
cxx_compiler_flag_test(CXX_MARCH "-march=native" "" 0)
cxx_compiler_flag_test(CXX_CLANG_ASSEMBLER "-Wa,-q" "" 0)

if (CMAKE_BUILD_TYPE_UPPER STREQUAL "DEBUG")
  cxx_compiler_flag_test(CXX_DEBUG "-g" "" 1)
elseif(CMAKE_BUILD_TYPE_UPPER STREQUAL "SAFE")
  cxx_compiler_flag_test(CXX_O2 "-O2" "" 0)
else()
  cxx_compiler_flag_test(CXX_FUNROLL_LOOPS "-funroll-loops" "" 0)
  cxx_compiler_flag_test(CXX_OFAST "-Ofast" "-O3" 1)
  cxx_compiler_flag_test(C_FAST_MATH "-ffast-math" "" 0)
endif()

INSTALL(FILES rna_turner1999.par rna_turner2004.par DESTINATION bin)

find_package(OpenMP REQUIRED)
find_package(FFTW REQUIRED)
find_library(RNA libRNA.a PATHS)

include_directories(${FFTW_INCLUDES})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

add_executable(FFTbor main.cpp delta.cpp misc.cpp)

target_link_libraries(FFTbor ${FFTW_LIBRARIES})
target_link_libraries(FFTbor ${RNA})
target_link_libraries(FFTbor m)

install(TARGETS FFTbor DESTINATION bin)
