
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <utility>

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

  py::enum_<CellRowPolarity>(m, "CellRowPolarity")
      .value("SAME", CellRowPolarity::SAME, "Same")
      .value("OPPOSITE", CellRowPolarity::OPPOSITE, "Opposite")
      .value("ANY", CellRowPolarity::ANY, "Any")
      .value("NW", CellRowPolarity::ANY, "North or West")
      .value("SE", CellRowPolarity::ANY, "South or East")
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
      .def_property_readonly("height", &Rectangle::height, "Height")
      .def_property_readonly("width", &Rectangle::width, "Width")
      .def("__str__", &Rectangle::toString)
      .def("__repr__", &Rectangle::toString);

  py::class_<Row, Rectangle>(m, "Row")
      .def(py::init<Rectangle, CellOrientation>(), R"pbdoc(
        Construct a row
)pbdoc",
           py::arg("area"), py::arg("orientation"))
      .def_readwrite("orientation", &Row::orientation);

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

  py::enum_<PlacementStep>(m, "PlacementStep")
      .value("LowerBound", PlacementStep::LowerBound)
      .value("UpperBound", PlacementStep::UpperBound)
      .value("Detailed", PlacementStep::Detailed)
      .export_values();

  py::class_<ColoquinteParameters>(m, "ColoquinteParameters")
      .def(py::init<int, int>(), R"pbdoc(
Construct the parameters

:param int effort: Effort level
)pbdoc",
           py::arg("effort") = 3, py::arg("seed") = -1)

      .def_readwrite("global", &ColoquinteParameters::global)
      .def_readwrite("legalization", &ColoquinteParameters::legalization)
      .def_readwrite("detailed", &ColoquinteParameters::detailed)
      .def("check", &ColoquinteParameters::check)
      .def("__str__", &ColoquinteParameters::toString)
      .def("__repr__", &ColoquinteParameters::toString);

  py::class_<ContinuousModelParameters>(m, "ContinuousModelParameters")
      .def(py::init<int>(), R"pbdoc(
Construct the parameters

:param int effort: Effort level
)pbdoc",
           py::arg("effort") = 3)
      .def_readwrite("net_model", &ContinuousModelParameters::netModel)
      .def_readwrite("approximation_distance",
                     &ContinuousModelParameters::approximationDistance)
      .def_readwrite(
          "approximation_distance_update_factor",
          &ContinuousModelParameters::approximationDistanceUpdateFactor)
      .def_readwrite("max_nb_conjugate_gradient_steps",
                     &ContinuousModelParameters::maxNbConjugateGradientSteps)
      .def_readwrite(
          "conjugate_gradient_error_tolerance",
          &ContinuousModelParameters::conjugateGradientErrorTolerance)
      .def("check", &ContinuousModelParameters::check)
      .def("__str__", &ContinuousModelParameters::toString)
      .def("__repr__", &ContinuousModelParameters::toString);

  py::class_<RoughLegalizationParameters>(m, "RoughLegalizationParameters")
      .def(py::init<int>(), R"pbdoc(
Construct the parameters

:param int effort: Effort level
)pbdoc",
           py::arg("effort") = 3)
      .def_readwrite("cost_model", &RoughLegalizationParameters::costModel)
      .def_readwrite("nb_steps", &RoughLegalizationParameters::nbSteps)
      .def_readwrite("bin_size", &RoughLegalizationParameters::binSize)
      .def_readwrite("line_reopt_size",
                     &RoughLegalizationParameters::lineReoptSize)
      .def_readwrite("line_reopt_overlap",
                     &RoughLegalizationParameters::lineReoptOverlap)
      .def_readwrite("diag_reopt_size",
                     &RoughLegalizationParameters::diagReoptSize)
      .def_readwrite("diag_reopt_overlap",
                     &RoughLegalizationParameters::diagReoptOverlap)
      .def_readwrite("square_reopt_size",
                     &RoughLegalizationParameters::squareReoptSize)
      .def_readwrite("square_reopt_overlap",
                     &RoughLegalizationParameters::squareReoptOverlap)
      .def_readwrite("target_blending",
                     &RoughLegalizationParameters::targetBlending)
      .def_readwrite("unidimensional_transport",
                     &RoughLegalizationParameters::unidimensionalTransport)
      .def_readwrite("quadratic_penalty",
                     &RoughLegalizationParameters::quadraticPenalty)
      .def_readwrite("side_margin", &RoughLegalizationParameters::sideMargin)
      .def_readwrite("coarsening_limit",
                     &RoughLegalizationParameters::coarseningLimit)
      .def("check", &RoughLegalizationParameters::check)
      .def("__str__", &RoughLegalizationParameters::toString)
      .def("__repr__", &RoughLegalizationParameters::toString);

  py::class_<PenaltyParameters>(m, "PenaltyParameters")
      .def(py::init<int>(), R"pbdoc(
Construct the parameters

:param int effort: Effort level
)pbdoc",
           py::arg("effort") = 3)
      .def_readwrite("cutoff_distance", &PenaltyParameters::cutoffDistance)
      .def_readwrite("cutoff_distance_update_factor",
                     &PenaltyParameters::cutoffDistanceUpdateFactor)
      .def_readwrite("area_exponent", &PenaltyParameters::areaExponent)
      .def_readwrite("initial_value", &PenaltyParameters::initialValue)
      .def_readwrite("update_factor", &PenaltyParameters::updateFactor)
      .def_readwrite("target_blending", &PenaltyParameters::targetBlending)
      .def("check", &PenaltyParameters::check)
      .def("__str__", &PenaltyParameters::toString)
      .def("__repr__", &PenaltyParameters::toString);

  py::class_<GlobalPlacerParameters>(m, "GlobalPlacerParameters")
      .def(py::init<int>(), R"pbdoc(
Construct the parameters

:param int effort: Effort level
)pbdoc",
           py::arg("effort") = 3)
      .def_readwrite("rough_legalization",
                     &GlobalPlacerParameters::roughLegalization)
      .def_readwrite("continuous_model",
                     &GlobalPlacerParameters::continuousModel)
      .def_readwrite("penalty", &GlobalPlacerParameters::penalty)
      .def_readwrite("max_nb_steps", &GlobalPlacerParameters::maxNbSteps)
      .def_readwrite("nb_initial_steps",
                     &GlobalPlacerParameters::nbInitialSteps)
      .def_readwrite("nb_steps_before_rough_legalization",
                     &GlobalPlacerParameters::nbStepsBeforeRoughLegalization)
      .def_readwrite("gap_tolerance", &GlobalPlacerParameters::gapTolerance)
      .def_readwrite("distance_tolerance",
                     &GlobalPlacerParameters::distanceTolerance)
      .def_readwrite("export_blending", &GlobalPlacerParameters::exportBlending)
      .def_readwrite("noise", &GlobalPlacerParameters::noise)
      .def("check", &GlobalPlacerParameters::check)
      .def("__str__", &GlobalPlacerParameters::toString)
      .def("__repr__", &GlobalPlacerParameters::toString);

  py::class_<LegalizationParameters>(m, "LegalizationParameters")
      .def(py::init<int>(), R"pbdoc(
Construct the parameters

:param int effort: Effort level
)pbdoc",
           py::arg("effort") = 3)
      .def_readwrite("cost_model", &LegalizationParameters::costModel)
      .def_readwrite("ordering_height", &LegalizationParameters::orderingHeight)
      .def_readwrite("ordering_width", &LegalizationParameters::orderingWidth)
      .def_readwrite("ordering_y", &LegalizationParameters::orderingY)
      .def("check", &LegalizationParameters::check)
      .def("__str__", &LegalizationParameters::toString)
      .def("__repr__", &LegalizationParameters::toString);

  py::class_<DetailedPlacerParameters>(m, "DetailedPlacerParameters")
      .def(py::init<int>(), R"pbdoc(
Construct the parameters

:param int effort: Effort level
)pbdoc",
           py::arg("effort") = 3)
      .def_readwrite("nb_passes", &DetailedPlacerParameters::nbPasses)
      .def_readwrite("local_search_nb_neighbours",
                     &DetailedPlacerParameters::localSearchNbNeighbours)
      .def_readwrite("local_search_nb_rows",
                     &DetailedPlacerParameters::localSearchNbRows)
      .def_readwrite("reordering_nb_rows",
                     &DetailedPlacerParameters::reorderingNbRows)
      .def_readwrite("reordering_max_nb_cells",
                     &DetailedPlacerParameters::reorderingMaxNbCells)
      .def_readwrite("shift_nb_rows", &DetailedPlacerParameters::shiftNbRows)
      .def_readwrite("shift_max_nb_cells",
                     &DetailedPlacerParameters::shiftMaxNbCells)
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

      .def_property("cell_width", &Circuit::cellWidth, &Circuit::setCellWidth,
                    "Width of the cells")
      .def_property("cell_height", &Circuit::cellHeight,
                    &Circuit::setCellHeight, "Height of the cells")
      .def_property("cell_is_fixed", &Circuit::cellIsFixed,
                    &Circuit::setCellIsFixed, "Fixed status flag of the cells")
      .def_property("cell_is_obstruction", &Circuit::cellIsObstruction,
                    &Circuit::setCellIsObstruction,
                    "Obstruction status flag of the cells")
      .def_property(
          "cell_row_polarity", &Circuit::cellRowPolarity,
          &Circuit::setCellRowPolarity,
          "Polarity of the cell with respect to the standard cell rows")
      .def_property("cell_x", &Circuit::cellX, &Circuit::setCellX,
                    "X position of the cells")
      .def_property("cell_y", &Circuit::cellY, &Circuit::setCellY,
                    "Y position of the cells")
      .def_property("cell_orientation", &Circuit::cellOrientation,
                    &Circuit::setCellOrientation, "Orientation of the cells")
      .def_property_readonly("cell_placement", &Circuit::cellPlacement,
                             "Place occupied by the cells")
      .def_property_readonly("placement_area", &Circuit::computePlacementArea,
                             "Bounding box of the placement area")
      .def_property_readonly("row_height", &Circuit::rowHeight,
                             "Standard-cell row height")
      .def_property("rows", &Circuit::rows, &Circuit::setRows,
                    "Standard cell rows")
      .def("add_net", &Circuit::addNet, "Add a net to the circuit")
      .def("hpwl", &Circuit::hpwl, "Compute the half-perimeter wirelength")
      .def("setup_rows", &Circuit::setupRows,
           "Setup the rows from a placement area")
      .def("place", &Circuit::place,
           "Run the whole placement algorithm (global and detailed)")
      .def(
          "place_global",
          [](Circuit &circuit, const ColoquinteParameters &params,
             std::optional<PlacementCallback> callback) {
            py::gil_scoped_release release;
            circuit.placeGlobal(params, std::move(callback));
          },
          "Run the global placement algorithm")
      .def(
          "legalize",
          [](Circuit &circuit, const ColoquinteParameters &params,
             std::optional<PlacementCallback> callback) {
            py::gil_scoped_release release;
            circuit.legalize(params, std::move(callback));
          },
          "Run the detailed placement algorithm")
      .def(
          "place_detailed",
          [](Circuit &circuit, const ColoquinteParameters &params,
             std::optional<PlacementCallback> callback) {
            py::gil_scoped_release release;
            circuit.placeDetailed(params, std::move(callback));
          },
          "Run the detailed placement algorithm")
      .def("expand_cells_to_density", &Circuit::expandCellsToDensity,
           py::arg("target_density"), py::arg("row_side_margin") = 0.0,
           "Expand the standard cells to reach the target density")
      .def("expand_cells_by_factor", &Circuit::expandCellsByFactor,
           py::arg("expansion_factor"), py::arg("max_density") = 1.0,
           py::arg("row_side_margin") = 0.0,
           "Expand the standard cells by an individual factor")
      .def("check", &Circuit::check, "Check the datastructure")
      .def("report", &Circuit::report)
      .def("__str__", &Circuit::toString)
      .def("__repr__", &Circuit::toString);
}
