cmake_minimum_required(VERSION 3.12)

include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)

project(vgrna-project-scripts)

find_package(OpenMP REQUIRED)

ExternalProject_Add(libvgio-proj
  GIT_REPOSITORY "https://github.com/vgteam/libvgio.git"
  GIT_TAG "5f9e7cbd873ad11c82dcfddb0488e1b769d6dd9b"
  SOURCE_DIR "${PROJECT_SOURCE_DIR}/external/libvgio"
  BUILD_IN_SOURCE TRUE
  UPDATE_COMMAND ""
  CONFIGURE_COMMAND cmake -DCMAKE_INSTALL_PREFIX=${PROJECT_SOURCE_DIR}/external/libvgio .
  INSTALL_COMMAND ""
)

ExternalProject_Add(seqlib-proj
  GIT_REPOSITORY "https://github.com/walaj/SeqLib.git"
  GIT_TAG "08771285f60f3434c796cd4126ea02c577632ddb"
  SOURCE_DIR "${PROJECT_SOURCE_DIR}/external/seqlib"
  BUILD_IN_SOURCE TRUE
  UPDATE_COMMAND ""
  CONFIGURE_COMMAND ./configure
  INSTALL_COMMAND ""
)

include_directories(
  ${PROJECT_SOURCE_DIR}/external/libvgio
  ${PROJECT_SOURCE_DIR}/external/libvgio/include
  ${PROJECT_SOURCE_DIR}/external/seqlib
  ${PROJECT_SOURCE_DIR}/external/seqlib/htslib
)

link_directories(
  ${PROJECT_SOURCE_DIR}/external/libvgio
  ${PROJECT_SOURCE_DIR}/external/seqlib/src
)

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

  set(CMAKE_CXX_FLAGS "--std=c++14 -Xpreprocessor -fopenmp -g -O3 -Wall")

elseif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")

  set(CMAKE_CXX_FLAGS "--std=c++14 -fopenmp -g -O3")

endif()

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

add_executable(calc_allele_read_coverage src/calc_allele_read_coverage.cpp)
add_dependencies(calc_allele_read_coverage seqlib-proj)
target_link_libraries(calc_allele_read_coverage seqlib hts)

