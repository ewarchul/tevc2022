cmake_minimum_required(VERSION 3.20)

file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS ${SOURCE_DIR}/*.cpp
     ${SOURCE_DIR}/**/*.cpp)
file(GLOB_RECURSE HEADERS CONFIGURE_DEPENDS ${SOURCE_DIR}/*.h
     ${SOURCE_DIR}/**/*.h)

set(EXECUTABLE "ppmf")
set(DEPENDENCIES spdlog::spdlog fmt::fmt LAPACK::LAPACK blaze::blaze
                 libcmaes::cmaes)

add_executable(${EXECUTABLE} ${SOURCES} ${HEADERS})
add_dependencies(${EXECUTABLE} ${DEPENDENCIES})
target_link_libraries(${EXECUTABLE} ${DEPENDENCIES})
target_compile_options(${EXECUTABLE} PRIVATE -O0)
set_property(TARGET ${EXECUTABLE} PROPERTY CXX_STANDARD 20)
