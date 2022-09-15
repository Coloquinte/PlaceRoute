
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "coloquinte.hpp"
#include "place_detailed/place_detailed.hpp"
#include "place_global/place_global.hpp"

using namespace coloquinte;

namespace py = pybind11;

PYBIND11_MODULE(pycoloquinte, m) {
  m.doc() = R"pbdoc(
        Coloquinte VLSI placer
        -----------------------
        .. currentmodule:: pycoloquinte
        .. autosummary::
           :toctree: _generate
           Circuit
    )pbdoc";

  py::enum_<LegalizationModel>(m, "LegalizationModel")
      .value("L1", LegalizationModel::L1)
      .value("L2", LegalizationModel::L2)
      .value("LInf", LegalizationModel::LInf)
      .value("L2Squared", LegalizationModel::L2Squared)
      .export_values();

  py::enum_<CellOrientation>(m, "CellOrientation")
      .value("N", CellOrientation::N, "North")
      .value("S", CellOrientation::S, "South")
      .value("W", CellOrientation::W, "West")
      .value("E", CellOrientation::E, "East")
      .value("FN", CellOrientation::FN, "Flipped + North")
      .value("FS", CellOrientation::FS, "Flipped + South")
      .value("FW", CellOrientation::FW, "Flipped + West")
      .value("FE", CellOrientation::FE, "Flipped + East")
      .export_values();

  py::class_<Circuit>(m, "Circuit")
      .def(py::init<int>(), R"pbdoc(
Construct a circuit.

:param int nb_cells: Number of cells
)pbdoc",
           py::arg("nb_cells"))
      .def_property_readonly("nb_cells", &Circuit::nbCells, "Number of cells")
      .def_property_readonly("nb_nets", &Circuit::nbNets, "Number of nets")
      .def_property_readonly("nb_rows", &Circuit::nbRows,
                             "Number of standard cell rows")
      .def_property_readonly("nb_pins", &Circuit::nbPins,
                             "Total number of pins")
      .def_property("cell_x", &Circuit::setCellX, &Circuit::setCellX,
                    "X position of the cells")
      .def_property("cell_y", &Circuit::setCellY, &Circuit::setCellY,
                    "Y position of the cells")
      .def_property("cell_width", &Circuit::cellWidth, &Circuit::setCellWidth,
                    "Width of the cells")
      .def_property("cell_height", &Circuit::cellHeight,
                    &Circuit::setCellHeight, "Height of the cells")
      .def_property("cell_fixed", &Circuit::cellFixed, &Circuit::setCellFixed,
                    "Fixed status flag of the cells")
      .def_property("cell_orientation", &Circuit::cellOrientation,
                    &Circuit::setCellOrientation, "Orientation of the cells")
      .def_property("rows", &Circuit::rows,
                    &Circuit::setRows, "Standard cell rows")
      .def("add_net", &Circuit::addNet, "Add a net to the circuit")
      .def("hpwl", &Circuit::hpwl, "Compute the half-perimeter wirelength")
      .def("check", &Circuit::check, "Check the datastructure");
}
