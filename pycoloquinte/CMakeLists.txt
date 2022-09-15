
pybind11_add_module(pycoloquinte module.cpp)
target_link_libraries(pycoloquinte PRIVATE coloquinte ${Boost_LIBRARIES})
