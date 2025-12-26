cmake_minimum_required(VERSION 3.16)
project(STIR2PETSIRD VERSION 0.3.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(GNUInstallDirs)

find_package(STIR REQUIRED)
find_package(PETSIRD CONFIG REQUIRED)

enable_testing()


add_library(stir2petsird
  STIR_PETSIRD_convertor.cpp
)

add_library(STIR2PETSIRD::stir2petsird ALIAS stir2petsird)

target_include_directories(stir2petsird
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(stir2petsird
  PUBLIC
    PETSIRD::petsird
    ${STIR_LIBRARIES}
)

target_compile_features(stir2petsird PUBLIC cxx_std_17)

add_executable(STIR_PETSIRD_convertor
  STIR_PETSIRD_convertor
)

target_link_libraries(STIR_PETSIRD_convertor
  PRIVATE
    STIR2PETSIRD::stir2petsird
)

# -----------------------------
# Installation
# -----------------------------
install(
  TARGETS stir2petsird STIR_PETSIRD_convertor
  EXPORT STIR2PETSIRDTargets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(
  DIRECTORY 
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# -----------------------------
# CMake package export
# -----------------------------
install(
  EXPORT STIR2PETSIRDTargets
  FILE STIR2PETSIRDTargets-noconfig.cmake
  NAMESPACE STIR2PETSIRD::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/STIR2PETSIRD
)

include(CMakePackageConfigHelpers)

configure_package_config_file(
  cmake/STIR2PETSIRDConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/STIR2PETSIRDConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/STIR2PETSIRD
)

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/STIR2PETSIRDConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

install(
  FILES
    ${CMAKE_CURRENT_BINARY_DIR}/STIR2PETSIRDConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/STIR2PETSIRDConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/STIR2PETSIRD
)

add_executable(test_stir2petsird tests/test_stir2petsird.cxx)
target_include_directories(test_stir2petsird
  PRIVATE
    ${STIR_INCLUDE_DIRS}
)
target_link_libraries(test_stir2petsird PRIVATE ${STIR_LIBRARIES} stir2petsird)
add_test(NAME test_stir2petsird COMMAND test_stir2petsird)