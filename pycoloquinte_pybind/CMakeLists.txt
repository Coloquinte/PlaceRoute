
pybind11_add_module(pycoloquinte_pybind module.cpp)
target_link_libraries(pycoloquinte_pybind PRIVATE coloquinte ${Boost_LIBRARIES})
