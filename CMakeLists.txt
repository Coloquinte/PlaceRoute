
cmake_minimum_required(VERSION 3.10)
PROJECT(COLOQUINTE)

#add_subdirectory(thirdparty/kahypar)

FIND_PACKAGE(Boost REQUIRED COMPONENTS system filesystem iostreams program_options unit_test_framework)
FIND_PACKAGE(Eigen3 3.3 REQUIRED NO_MODULE)

INCLUDE_DIRECTORIES(
    ${COLOQUINTE_SOURCE_DIR}/src
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
)

add_library(coloquinte SHARED ${SOURCES})
target_compile_features(coloquinte PRIVATE cxx_std_17)
target_link_libraries(coloquinte
  ${Boost_LIBRARIES}
)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    enable_testing()
    add_subdirectory(test)
endif()
