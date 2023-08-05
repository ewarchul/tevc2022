cmake_minimum_required(VERSION 3.20)

project(
  ppmf
  VERSION 0.1
  LANGUAGES CXX
  DESCRIPTION "Source code for PPMF experiments")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Aliases
set(ROOT_DIR ${CMAKE_SOURCE_DIR})
set(SOURCE_DIR ${ROOT_DIR}/src)
set(CMAKE_MODULE_PATH ${ROOT_DIR}/cmake ${ROOT_DIR}/cmake/packages)

# Artifatcs directories
set(ARTIFACT_DIRECTORY ${CMAKE_SOURCE_DIR}/artifacts CACHE STRING "artifacts")
set(CMAKE_BINARY_DIR ${ARTIFACT_DIRECTORY}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/exe)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/archive)

# Dependencies
find_package(spdlog REQUIRED)
find_package(fmt 10.0.0 REQUIRED)
find_package(range-v3 REQUIRED)
find_package(Eigen3 3.0.0 REQUIRED)
find_package(blaze REQUIRED)
find_package(libcmaes REQUIRED)
find_package(LAPACK REQUIRED)

# Directory with experiments source code
add_subdirectory(src)
