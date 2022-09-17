
cmake_minimum_required(VERSION 3.10)
PROJECT(COLOQUINTE)

add_subdirectory(thirdparty/lemon)
add_subdirectory(thirdparty/pybind11)

SET(Boost_USE_STATIC_LIBS ON)
FIND_PACKAGE(Boost REQUIRED COMPONENTS system filesystem iostreams program_options unit_test_framework)
FIND_PACKAGE(Eigen3 3.3 REQUIRED NO_MODULE)

INCLUDE_DIRECTORIES(
    ${COLOQUINTE_SOURCE_DIR}/src
    ${COLOQUINTE_SOURCE_DIR}/thirdparty/lemon/lemon
    ${COLOQUINTE_SOURCE_DIR}/thirdparty/pybind11
)

SET(SOURCES
  src/coloquinte.cpp
  src/common.cpp
  src/place_global/net_model.cpp
  src/place_global/density_legalizer.cpp
  src/place_global/density_grid.cpp
  src/place_global/place_global.cpp
  src/place_detailed/legalizer.cpp
  src/place_detailed/row_legalizer.cpp
  src/place_detailed/place_detailed.cpp
  src/place_detailed/detailed_placement.cpp
  src/place_detailed/incr_net_model.cpp
  src/place_detailed/row_neighbourhood.cpp
)

add_library(coloquinte STATIC ${SOURCES})
set_target_properties(coloquinte PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
target_compile_features(coloquinte PRIVATE cxx_std_17)
target_link_libraries(coloquinte
  ${Boost_LIBRARIES}
  lemon
)
add_subdirectory(pycoloquinte)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    enable_testing()
    add_subdirectory(test)
endif()
