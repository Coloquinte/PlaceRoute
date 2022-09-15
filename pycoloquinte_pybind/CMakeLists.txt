
pybind11_add_module(coloquinte_pybind module.cpp)
target_link_libraries(coloquinte_pybind PRIVATE coloquinte ${Boost_LIBRARIES})
