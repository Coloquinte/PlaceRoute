
cmake_minimum_required(VERSION 2.8)
PROJECT(COLOQUINTE)

FIND_PACKAGE(Boost REQUIRED COMPONENTS system filesystem iostreams program_options unit_test_framework)
FIND_PACKAGE(XTensor REQUIRED)

INCLUDE_DIRECTORIES(
    ${COLOQUINTE_SOURCE_DIR}/src
)

SET(SOURCES
  src/coloquinte.cpp
  src/place_global/topology.cpp
  src/place_global/descent.cpp
)

add_library(coloquinte SHARED ${SOURCES})
target_include_directories(coloquinte PUBLIC ${xtensor_INCLUDE_DIRS})
target_link_libraries(coloquinte
  ${Boost_LIBRARIES}
  xtensor)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    enable_testing()
    add_subdirectory(test)
endif()
