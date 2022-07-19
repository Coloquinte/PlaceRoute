
cmake_minimum_required(VERSION 3.10)
PROJECT(COLOQUINTE)

add_subdirectory(thirdparty/xtl)
add_subdirectory(thirdparty/xtensor)
add_subdirectory(thirdparty/kahypar)

FIND_PACKAGE(Boost REQUIRED COMPONENTS system filesystem iostreams program_options unit_test_framework)
FIND_PACKAGE(Eigen3 3.3 REQUIRED NO_MODULE)

INCLUDE_DIRECTORIES(
    ${COLOQUINTE_SOURCE_DIR}/src
)

SET(SOURCES
  src/coloquinte.cpp
  src/common.cpp
  src/place_global/wirelength_model.cpp
  src/place_global/gradient_descent.cpp
  src/place_global/density_legalizer.cpp
  src/place_global/density_grid.cpp
  src/place_global/place_global.cpp
  src/place_detailed/legalizer.cpp
  src/place_detailed/place_detailed.cpp
)

add_library(coloquinte SHARED ${SOURCES})
target_compile_features(coloquinte PRIVATE cxx_std_17)
target_include_directories(coloquinte PRIVATE ${xtensor_INCLUDE_DIRS})
target_link_libraries(coloquinte
  ${Boost_LIBRARIES}
  xtensor)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    enable_testing()
    add_subdirectory(test)
endif()
