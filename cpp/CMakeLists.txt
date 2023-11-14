cmake_minimum_required(VERSION 3.12.0) # older would work, but could give warnings on policy CMP0074
project(STIR2PETSIRD VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 17)

if(WIN32)
  add_compile_options(/W3 /WX)
else()
  add_compile_options(-Wall -Wextra -pedantic)
endif()

find_package(STIR)

# Example lines for a new executable
add_executable(STIR_PETSIRD_convertor STIR_PETSIRD_convertor.cpp)

target_link_libraries(STIR_PETSIRD_convertor PUBLIC prd_generated ${STIR_LIBRARIES})
install(TARGETS STIR_PETSIRD_convertor DESTINATION bin)

add_subdirectory(../PETSIRD/cpp/generated PETSIRD_generated)
