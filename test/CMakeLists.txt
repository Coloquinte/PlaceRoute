
FIND_PACKAGE(Boost REQUIRED COMPONENTS unit_test_framework)

SET(TESTS
    test_wirelength_model
    test_density_legalizer
)

FOREACH(TEST IN LISTS TESTS)
    add_executable(${TEST} ${TEST}.cpp)
    target_include_directories(${TEST} PRIVATE ${Boost_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR}/include/)
    target_link_libraries(${TEST} ${Boost_LIBRARIES} coloquinte)
    add_test(${TEST} ${TEST})
ENDFOREACH(TEST)


