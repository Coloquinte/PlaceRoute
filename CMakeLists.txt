
cmake_minimum_required(VERSION 3.10)
PROJECT(COLOQUINTE)
SET(CMAKE_CXX_STANDARD 17)

OPTION(BUILD_PYCOLOQUINTE "Build Coloquinte Python package")

IF(BUILD_PYCOLOQUINTE)
    SET(Boost_USE_STATIC_LIBS ON)
ENDIF()

FIND_PACKAGE(Boost REQUIRED COMPONENTS system filesystem iostreams program_options unit_test_framework)
FIND_PACKAGE(LEMON NAMES LEMON lemon REQUIRED)
FIND_PACKAGE(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(Threads REQUIRED)

INCLUDE_DIRECTORIES(
    ${COLOQUINTE_SOURCE_DIR}/src
    ${LEMON_INCLUDE_DIRS}
    ${COLOQUINTE_SOURCE_DIR}/thirdparty/pybind11
)

SET(SOURCES
  src/coloquinte.cpp
  src/parameters.cpp
  src/place_global/net_model.cpp
  src/place_global/density_legalizer.cpp
  src/place_global/density_grid.cpp
  src/place_global/place_global.cpp
  src/place_detailed/legalizer.cpp
  src/place_detailed/abacus_legalizer.cpp
  src/place_detailed/tetris_legalizer.cpp
  src/place_detailed/row_legalizer.cpp
  src/place_detailed/place_detailed.cpp
  src/place_detailed/detailed_placement.cpp
  src/place_detailed/incr_net_model.cpp
  src/place_detailed/row_neighbourhood.cpp
  src/place_global/transportation.cpp
  src/place_global/transportation_1d.cpp
)

IF(BUILD_PYCOLOQUINTE)
    add_library(coloquinte STATIC ${SOURCES})
ELSE()
    add_library(coloquinte SHARED ${SOURCES})
ENDIF()

set_target_properties(coloquinte PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
target_compile_features(coloquinte PRIVATE cxx_std_17)
target_link_libraries(coloquinte
  PRIVATE
  ${Boost_LIBRARIES}
  Threads::Threads
)

IF(BUILD_PYCOLOQUINTE)
    add_subdirectory(thirdparty/pybind11)
    add_subdirectory(pycoloquinte)
ENDIF(BUILD_PYCOLOQUINTE)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    enable_testing()
    add_subdirectory(test)
endif()

INSTALL(
  TARGETS coloquinte
  LIBRARY DESTINATION lib
  COMPONENT library
)
INSTALL(
  FILES ${COLOQUINTE_SOURCE_DIR}/src/coloquinte.hpp
  DESTINATION include/coloquinte
  COMPONENT headers
)
