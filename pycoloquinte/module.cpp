
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

  py::enum_<LegalizationModel>(m, "LegalizationModel")
      .value("L1", LegalizationModel::L1)
      .value("L2", LegalizationModel::L2)
      .value("LInf", LegalizationModel::LInf)
      .value("L1Squared", LegalizationModel::L1Squared)
      .value("L2Squared", LegalizationModel::L2Squared)
      .value("LInfSquared", LegalizationModel::LInfSquared)
      .export_values();

  py::enum_<NetModelOption>(m, "NetModel")
      .value("BoundToBound", NetModelOption::BoundToBound)
      .value("Star", NetModelOption::Star)
      .value("Clique", NetModelOption::Clique)
      .value("LightStar", NetModelOption::LightStar)
      .export_values();

  py::class_<GlobalPlacerParameters>(m, "GlobalPlacerParameters")
      .def(py::init<int, int>(), R"pbdoc(
Construct the parameters

:param int effort: Effort level
:param int seed: Random seed
)pbdoc",
           py::arg("effort") = 3, py::arg("seed") = -1)
      .def_readwrite("max_nb_steps", &GlobalPlacerParameters::maxNbSteps)
      .def_readwrite("nb_initial_steps",
                     &GlobalPlacerParameters::nbInitialSteps)
      .def_readwrite("gap_tolerance", &GlobalPlacerParameters::gapTolerance)
      .def_readwrite("penalty_cutoff_distance",
                     &GlobalPlacerParameters::penaltyCutoffDistance)
      .def_readwrite("initial_penalty", &GlobalPlacerParameters::initialPenalty)
      .def_readwrite("penalty_update_factor",
                     &GlobalPlacerParameters::penaltyUpdateFactor)
      .def_readwrite("net_model", &GlobalPlacerParameters::netModel)
      .def_readwrite("approximation_distance",
                     &GlobalPlacerParameters::approximationDistance)
      .def_readwrite("max_nb_conjugate_gradient_steps",
                     &GlobalPlacerParameters::maxNbConjugateGradientSteps)
      .def_readwrite("conjugate_gradient_error_tolerance",
                     &GlobalPlacerParameters::conjugateGradientErrorTolerance)
      .def_readwrite("rough_legalization_cost_model",
                     &GlobalPlacerParameters::roughLegalizationCostModel)
      .def_readwrite("rough_legalization_nb_steps",
                     &GlobalPlacerParameters::roughLegalizationNbSteps)
      .def_readwrite("rough_legalization_bin_size",
                     &GlobalPlacerParameters::roughLegalizationBinSize)
      .def_readwrite("export_weighting",
                     &GlobalPlacerParameters::exportWeighting)
      .def_readwrite("seed", &GlobalPlacerParameters::seed)
      .def("check", &GlobalPlacerParameters::check)
      .def("__str__", &GlobalPlacerParameters::toString)
      .def("__repr__", &GlobalPlacerParameters::toString);

  py::class_<DetailedPlacerParameters>(m, "DetailedPlacerParameters")
      .def(py::init<int, int>(), R"pbdoc(
Construct the parameters

:param int effort: Effort level
:param int seed: Random seed
)pbdoc",
           py::arg("effort") = 3, py::arg("seed") = -1)
      .def_readwrite("nb_passes", &DetailedPlacerParameters::nbPasses)
      .def_readwrite("local_search_nb_neighbours",
                     &DetailedPlacerParameters::localSearchNbNeighbours)
      .def_readwrite("local_search_nb_rows",
                     &DetailedPlacerParameters::localSearchNbRows)
      .def_readwrite("shift_nb_rows", &DetailedPlacerParameters::shiftNbRows)
      .def_readwrite("shift_max_nb_cells",
                     &DetailedPlacerParameters::shiftMaxNbCells)
      .def_readwrite("legalization_cost_model",
                     &DetailedPlacerParameters::legalizationCostModel)
      .def_readwrite("seed", &DetailedPlacerParameters::seed)
      .def("check", &DetailedPlacerParameters::check)
      .def("__str__", &DetailedPlacerParameters::toString)
      .def("__repr__", &DetailedPlacerParameters::toString);

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
      .def_property("cell_is_fixed", &Circuit::cellIsFixed,
                    &Circuit::setCellIsFixed, "Fixed status flag of the cells")
      .def_property("cell_is_obstruction", &Circuit::cellIsObstruction,
                    &Circuit::setCellIsObstruction,
                    "Obstruction status flag of the cells")
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
          "legalize",
          [](Circuit &circuit, int effort) {
            py::gil_scoped_release release;
            circuit.legalize(effort);
          },
          "Run the detailed placement algorithm")
      .def(
          "legalize",
          [](Circuit &circuit, const DetailedPlacerParameters &params) {
            py::gil_scoped_release release;
            circuit.legalize(params);
          },
          "Run the detailed placement algorithm")
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

  py::class_<GlobalRoutingPin>(m, "GlobalRoutingPin")
      .def(py::init<int, int, int>(), R"pbdoc(
        Construct a routing pin
)pbdoc",
           py::arg("x"), py::arg("y"), py::arg("z"))
      .def_readwrite("x", &GlobalRoutingPin::x)
      .def_readwrite("y", &GlobalRoutingPin::y)
      .def_readwrite("z", &GlobalRoutingPin::z);

  py::class_<GlobalRoutingSegment>(m, "GlobalRoutingSegment")
      .def(py::init<GlobalRoutingPin, GlobalRoutingPin>(), R"pbdoc(
        Construct a routing segment
)pbdoc",
           py::arg("a"), py::arg("b"))
      .def_property_readonly("is_horizontal",
                             &GlobalRoutingSegment::isHorizontal)
      .def_property_readonly("is_vertical", &GlobalRoutingSegment::isVertical)
      .def_property_readonly("is_via", &GlobalRoutingSegment::isVia)
      .def_property_readonly("length", &GlobalRoutingSegment::length)
      .def_readwrite("a", &GlobalRoutingSegment::a)
      .def_readwrite("b", &GlobalRoutingSegment::b);

  py::class_<GlobalRoutingProblem>(m, "GlobalRoutingProblem")
      .def(py::init<int, int, int>(), R"pbdoc(
        Initialize a routing problem
)pbdoc",
           py::arg("width"), py::arg("height"), py::arg("nb_layers"))
      .def_property_readonly("width", &GlobalRoutingProblem::width,
                             "Width of the grid")
      .def_property_readonly("height", &GlobalRoutingProblem::height,
                             "Height of the grid")
      .def_property_readonly("nb_layers", &GlobalRoutingProblem::nbLayers,
                             "Number of layers of the grid")
      .def("set_horizontal_capacity",
           &GlobalRoutingProblem::setHorizontalCapacity,
           "Set the horizontal capacity for a whole layer")
      .def("set_vertical_capacity", &GlobalRoutingProblem::setVerticalCapacity,
           "Set the vertical capacity for a whole layer")
      .def("set_via_capacity", &GlobalRoutingProblem::setViaCapacity,
           "Set the via capacity globally")
      .def("get_capacity", &GlobalRoutingProblem::capacity,
           "Get the capacity between neighbours")
      .def("set_capacity", &GlobalRoutingProblem::setCapacity,
           "Set the capacity between neighbours")
      .def("add_net", &GlobalRoutingProblem::addNet, "Add a net to the problem")
      .def("get_pins", &GlobalRoutingProblem::pins, "Get the pins of a net")
      .def("get_routing", &GlobalRoutingProblem::routing,
           "Get the routing of a net")
      .def("check", &GlobalRoutingProblem::check, "Check the datastructure")
      .def("__str__", &GlobalRoutingProblem::toString)
      .def("__repr__", &GlobalRoutingProblem::toString);
}
