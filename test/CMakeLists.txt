
FIND_PACKAGE(Boost REQUIRED COMPONENTS unit_test_framework)

SET(TESTS
    net_topology
)

FOREACH(TEST IN LISTS TESTS)
    add_executable(${TEST}_test ${TEST}.cpp)
    target_include_directories(${TEST}_test PRIVATE ${Boost_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR}/include/)
    target_link_libraries(${TEST}_test ${Boost_LIBRARIES} coloquinte)
    add_test(${TEST} ${TEST}_test)
ENDFOREACH(TEST)


