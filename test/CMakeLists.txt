
FIND_PACKAGE(Boost REQUIRED COMPONENTS unit_test_framework)

SET(TESTS
    test_net_model
    test_density_legalizer
    test_density_grid
    test_row_legalizer
    test_legalizer
    test_incr_net_model
    test_detailed_placement
    test_row_neighbourhood
    test_transportation
    test_interior_point_transportation
)

FOREACH(TEST IN LISTS TESTS)
    add_executable(${TEST} ${TEST}.cpp)
    target_include_directories(${TEST} PRIVATE ${Boost_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR}/include/)
    target_link_libraries(${TEST} ${Boost_LIBRARIES} coloquinte)
    add_test(${TEST} ${TEST})
ENDFOREACH(TEST)


