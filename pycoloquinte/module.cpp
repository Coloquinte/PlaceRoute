
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "coloquinte.hpp"

using namespace coloquinte;

namespace py = pybind11;

PYBIND11_MODULE(coloquinte_pybind, m) {
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

  py::class_<Rectangle>(m, "Rectangle")
      .def(py::init<int, int, int, int>(), R"pbdoc(
        Construct a rectangle
)pbdoc",
           py::arg("min_x"), py::arg("max_x"), py::arg("min_y"),
           py::arg("max_y"))
      .def_readwrite("min_x", &Rectangle::minX)
      .def_readwrite("max_x", &Rectangle::maxX)
      .def_readwrite("min_y", &Rectangle::minY)
      .def_readwrite("max_y", &Rectangle::maxY)
      .def("__str__", &Rectangle::toString)
      .def("__repr__", &Rectangle::toString);

  py::class_<GlobalPlacerParameters>(m, "GlobalPlacerParameters")
      .def(py::init<int>(), R"pbdoc(
Construct the parameters

:param int effort: Effort level
)pbdoc",
           py::arg("effort"))
      .def_readwrite("max_nb_steps", &GlobalPlacerParameters::maxNbSteps)
      .def_readwrite("gap_tolerance", &GlobalPlacerParameters::gapTolerance)
      .def_readwrite("penalty_cutoff_distance",
                     &GlobalPlacerParameters::penaltyCutoffDistance)
      .def_readwrite("initial_penalty", &GlobalPlacerParameters::initialPenalty)
      .def_readwrite("penalty_update_factor",
                     &GlobalPlacerParameters::penaltyUpdateFactor)
      .def_readwrite("approximation_distance",
                     &GlobalPlacerParameters::approximationDistance)
      .def_readwrite("max_nb_conjugate_gradient_steps",
                     &GlobalPlacerParameters::maxNbConjugateGradientSteps)
      .def_readwrite("conjugate_gradient_error_tolerance",
                     &GlobalPlacerParameters::conjugateGradientErrorTolerance);

  py::class_<DetailedPlacerParameters>(m, "DetailedPlacerParameters")
      .def(py::init<int>(), R"pbdoc(
Construct the parameters

:param int effort: Effort level
)pbdoc",
           py::arg("effort"))
      .def_readwrite("nb_passes", &DetailedPlacerParameters::nbPasses)
      .def_readwrite("local_search_nb_neighbours",
                     &DetailedPlacerParameters::localSearchNbNeighbours)
      .def_readwrite("local_search_nb_rows",
                     &DetailedPlacerParameters::localSearchNbRows)
      .def_readwrite("shift_nb_rows", &DetailedPlacerParameters::shiftNbRows);

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
      .def_property("cell_x", &Circuit::cellX, &Circuit::setCellX,
                    "X position of the cells")
      .def_property("cell_y", &Circuit::cellY, &Circuit::setCellY,
                    "Y position of the cells")
      .def_property("cell_width", &Circuit::cellWidth, &Circuit::setCellWidth,
                    "Width of the cells")
      .def_property("cell_height", &Circuit::cellHeight,
                    &Circuit::setCellHeight, "Height of the cells")
      .def_property("cell_fixed", &Circuit::cellFixed, &Circuit::setCellFixed,
                    "Fixed status flag of the cells")
      .def_property("cell_orientation", &Circuit::cellOrientation,
                    &Circuit::setCellOrientation, "Orientation of the cells")
      .def_property("rows", &Circuit::rows, &Circuit::setRows,
                    "Standard cell rows")
      .def("add_net", &Circuit::addNet, "Add a net to the circuit")
      .def("hpwl", &Circuit::hpwl, "Compute the half-perimeter wirelength")
      .def("place", &Circuit::place,
           "Run the whole placement algorithm (global and detailed)")
      .def(
          "place_global",
          [](Circuit &circuit, int effort) {
            py::gil_scoped_release release;
            circuit.placeGlobal(effort);
          },
          "Run the global placement algorithm")
      .def(
          "place_global",
          [](Circuit &circuit, const GlobalPlacerParameters &params) {
            py::gil_scoped_release release;
            circuit.placeGlobal(params);
          },
          "Run the global placement algorithm")
      .def(
          "place_detailed",
          [](Circuit &circuit, int effort) {
            py::gil_scoped_release release;
            circuit.placeDetailed(effort);
          },
          "Run the detailed placement algorithm")
      .def(
          "place_detailed",
          [](Circuit &circuit, const DetailedPlacerParameters &params) {
            py::gil_scoped_release release;
            circuit.placeDetailed(params);
          },
          "Run the detailed placement algorithm")
      .def("check", &Circuit::check, "Check the datastructure")
      .def("__str__", &Circuit::toString)
      .def("__repr__", &Circuit::toString);
}
