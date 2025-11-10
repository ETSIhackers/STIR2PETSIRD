cmake_minimum_required(VERSION 3.12.0)
project(STIR2PETSIRD VERSION 0.3.0)

set(CMAKE_CXX_STANDARD 17)

if(WIN32)
  add_compile_options(/W3 /WX)
else()
  add_compile_options(-Wall) # -Wextra -pedantic)
endif()

find_package(STIR REQUIRED)
if (STIR_VERSION VERSION_LESS 6.0.0) # probably will need 6.0 really
  message(FATAL_ERROR "Your STIR version is very old (${STIR_VERSION}). Use 6.0 or later")
endif()

set(PETSIRD_dir ../PETSIRD/cpp/)
add_subdirectory(${PETSIRD_dir} petsird)

add_executable(STIR_PETSIRD_convertor STIR_PETSIRD_convertor.cpp ${STIR_REGISTRIES})
target_link_libraries(STIR_PETSIRD_convertor PUBLIC petsird_helpers ${STIR_LIBRARIES})
install(TARGETS STIR_PETSIRD_convertor DESTINATION bin)

