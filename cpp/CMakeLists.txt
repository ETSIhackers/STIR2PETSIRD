cmake_minimum_required(VERSION 3.12.0) # older would work, but could give warnings on policy CMP0074
project(PETSIRDUseCaseTemplate VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 17)

if(WIN32)
  add_compile_options(/W3 /WX)
else()
  add_compile_options(-Wall -Wextra -pedantic)
endif()

# Example lines for a new executable
# add_executable(my_prog my_prog.cpp)
# target_link_libraries(my_prog PUBLIC prd_generated)
# install(TARGETS my_prog DESTINATION bin)

add_subdirectory(../PETSIRD/cpp/generated PETSIRD_generated)
