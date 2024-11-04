cmake_minimum_required(VERSION 3.12.0)
project(STIR2PETSIRD VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 17)

if(WIN32)
  add_compile_options(/W3 /WX)
else()
  add_compile_options(-Wall -Wextra -pedantic)
endif()

find_package(STIR REQUIRED)
if (STIR_VERSION VERSION_LESS 5.2.0) # probably will need 6.0 really
  message(FATAL_ERROR "Your STIR version is very old (${STIR_VERSION}). Use 6.0 or later")
endif()

add_subdirectory(../PETSIRD/cpp/generated PETSIRD_generated)
include_directories(../PETSIRD/cpp/generated)

add_executable(STIR_PETSIRD_convertor STIR_PETSIRD_convertor.cpp ${STIR_REGISTRIES})
target_link_libraries(STIR_PETSIRD_convertor PUBLIC petsird_generated ${STIR_LIBRARIES})
install(TARGETS STIR_PETSIRD_convertor DESTINATION bin)

