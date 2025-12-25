cmake_minimum_required(VERSION 3.12.0)
project(STIR2PETSIRD VERSION 0.4.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(WIN32)
  add_compile_options(/W3 /WX)
else()
  add_compile_options(-Wall) # -Wextra -pedantic)
endif()

find_package(STIR CONFIG REQUIRED)
if (STIR_VERSION VERSION_LESS 6.0.0) # probably will need 6.0 really
  message(FATAL_ERROR "Your STIR version is very old (${STIR_VERSION}). Use 6.0 or later")
endif()

find_package(PETSIRD CONFIG REQUIRED)

add_executable(STIR_PETSIRD_convertor STIR_PETSIRD_convertor.cpp)

target_link_libraries(STIR_PETSIRD_convertor PUBLIC ${STIR_LIBRARIES} PETSIRD::petsird)
target_include_directories(STIR_PETSIRD_convertor
  PUBLIC
    ${STIR_INCLUDE_DIRS}
)

install(TARGETS STIR_PETSIRD_convertor RUNTIME DESTINATION bin)

