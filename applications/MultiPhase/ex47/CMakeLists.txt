cmake_minimum_required(VERSION 3.15)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  project(MultiPhase_ex47 LANGUAGES CXX)
endif()

get_filename_component(EX47_REPOSITORY_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../../.." ABSOLUTE)
add_executable(MultiPhase_ex47 ex47.cpp "${EX47_REPOSITORY_ROOT}/src/fe/Line.cpp" "${EX47_REPOSITORY_ROOT}/src/fe/Quadrilateral.cpp" "${EX47_REPOSITORY_ROOT}/src/fe/Triangle.cpp" "${EX47_REPOSITORY_ROOT}/src/fe/Hexaedron.cpp" "${EX47_REPOSITORY_ROOT}/src/fe/Tetrahedral.cpp" "${EX47_REPOSITORY_ROOT}/src/fe/Wedge.cpp")
target_include_directories(MultiPhase_ex47 PRIVATE "${EX47_REPOSITORY_ROOT}/src/fe" "${EX47_REPOSITORY_ROOT}/src/enums" "${EX47_REPOSITORY_ROOT}/src/quadrature/cutfem_ls")

target_compile_features(MultiPhase_ex47 PRIVATE cxx_std_14)
find_path(EX47_BOOST_INCLUDE_DIR boost/multiprecision/cpp_bin_float.hpp HINTS "${BOOST_INCLUDEDIR}" "${BOOST_ROOT}" ${Boost_INCLUDE_DIRS} ENV BOOST_INCLUDEDIR ENV BOOST_ROOT)
if(NOT EX47_BOOST_INCLUDE_DIR)
  message(FATAL_ERROR "ex47 reference checks need Boost headers; set EX47_BOOST_INCLUDE_DIR.")
endif()
target_include_directories(MultiPhase_ex47 PRIVATE "${EX47_BOOST_INCLUDE_DIR}")

option(EX47_TRACK_ALLOCATIONS "Count heap allocations in the ex47 quadrature benchmark" OFF)
if(EX47_TRACK_ALLOCATIONS)
  target_compile_definitions(MultiPhase_ex47 PRIVATE EX47_TRACK_ALLOCATIONS)
endif()

option(EX47_SANITIZE "Build the ex47 main with address and undefined-behavior sanitizers" OFF)
if(EX47_SANITIZE)
  target_compile_options(MultiPhase_ex47 PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer)
  target_link_libraries(MultiPhase_ex47 -fsanitize=address,undefined)
endif()
