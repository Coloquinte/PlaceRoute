
#include "coloquinte.hpp"

#include <boost/polygon/polygon.hpp>
#include <iostream>
#include <sstream>

#include "coloquinte.hpp"
#include "place_detailed/place_detailed.hpp"
#include "place_global/density_legalizer.hpp"
#include "place_global/net_model.hpp"
#include "place_global/place_global.hpp"

namespace bpl = boost::polygon;
namespace coloquinte {

std::string toString(LegalizationModel model) {
  switch (model) {
    case LegalizationModel::L1:
      return "L1";
    case LegalizationModel::L2:
      return "L2";
    case LegalizationModel::LInf:
      return "LInf";
    case LegalizationModel::L1Squared:
      return "L1Squared";
    case LegalizationModel::L2Squared:
      return "L2Squared";
    case LegalizationModel::LInfSquared:
      return "LInfSquared";
    default:
      return "UnknownLegalizationModel";
  }
}

std::string toString(NetModelOption model) {
  switch (model) {
    case NetModelOption::BoundToBound:
      return "BoundToBound";
    case NetModelOption::Star:
      return "Star";
    default:
      return "UnknownNetModel";
  }
}

std::string Rectangle::toString() const {
  std::stringstream ss;
  ss << "Rectangle " << minX << ".." << maxX << " x " << minY << ".." << maxY;
  return ss.str();
}

std::string GlobalPlacerParameters::toString() const {
  std::stringstream ss;
  ss << "Global placer params:"
     << "\n\tMax nb steps: " << maxNbSteps
     << "\n\tGap tolerance: " << gapTolerance
     << "\n\tPenalty cutoff distance: " << penaltyCutoffDistance
     << "\n\tInitial penalty: " << initialPenalty
     << "\n\tPenalty update factor: " << penaltyUpdateFactor
     << "\n\tNet model: " << coloquinte::toString(netModel)
     << "\n\tApproximation distance: " << approximationDistance
     << "\n\tMax nb CG steps: " << maxNbConjugateGradientSteps
     << "\n\tCG error tolerance: " << conjugateGradientErrorTolerance
     << "\n\tRough legalization cost model: "
     << coloquinte::toString(roughLegalizationCostModel)
     << "\n\tNb rough legalization steps: " << nbRoughLegalizationSteps
     << "\n\tRough legalization bin size: " << roughLegalizationBinSize
     << std::endl;
  return ss.str();
}

std::string DetailedPlacerParameters::toString() const {
  std::stringstream ss;
  ss << "Detailed placer params:"
     << "\n\tNb passes: " << nbPasses
     << "\n\tLocal search nb neighbours: " << localSearchNbNeighbours
     << "\n\tLocal search nb rows: " << localSearchNbRows
     << "\n\tShift nb rows: " << shiftNbRows
     << "\n\tShift max nb cells: " << shiftMaxNbCells
     << "\n\tLegalization cost model: "
     << coloquinte::toString(legalizationCostModel) << std::endl;
  return ss.str();
}

Circuit::Circuit(int nbCells) {
  cellWidth_.resize(nbCells);
  cellHeight_.resize(nbCells);
  cellIsFixed_.resize(nbCells, false);
  cellIsObstruction_.resize(nbCells, true);
  cellX_.resize(nbCells);
  cellY_.resize(nbCells);
  cellOrientation_.resize(nbCells);
  netLimits_.push_back(0);
  check();
}

void Circuit::addNet(const std::vector<int> &cells,
                     const std::vector<int> &xOffsets,
                     const std::vector<int> &yOffsets) {
  if (cells.size() != xOffsets.size() || cells.size() != yOffsets.size()) {
    throw std::runtime_error("Inconsistent number of pins for the net");
  }
  if (cells.empty()) return;
  netLimits_.push_back(netLimits_.back() + cells.size());
  pinCells_.insert(pinCells_.end(), cells.begin(), cells.end());
  pinXOffsets_.insert(pinXOffsets_.end(), xOffsets.begin(), xOffsets.end());
  pinYOffsets_.insert(pinYOffsets_.end(), yOffsets.begin(), yOffsets.end());
}

void Circuit::setNets(const std::vector<int> &limits,
                      const std::vector<int> &cells,
                      const std::vector<int> &xOffsets,
                      const std::vector<int> &yOffsets) {
  assert(!limits.empty());
  assert(limits.front() == 0);
  assert(limits.back() == cells.size());
  assert(limits.back() == xOffsets.size());
  assert(limits.back() == yOffsets.size());
  netLimits_ = limits;
  pinCells_ = cells;
  pinXOffsets_ = xOffsets;
  pinYOffsets_ = yOffsets;
}

void Circuit::setCellX(const std::vector<int> &x) {
  if (x.size() != nbCells()) {
    throw std::runtime_error(
        "Number of elements is not the same as the number of cells of the "
        "circuit");
  }
  cellX_ = x;
}

void Circuit::setCellY(const std::vector<int> &y) {
  if (y.size() != nbCells()) {
    throw std::runtime_error(
        "Number of elements is not the same as the number of cells of the "
        "circuit");
  }
  cellY_ = y;
}

void Circuit::setCellIsFixed(const std::vector<bool> &f) {
  if (f.size() != nbCells()) {
    throw std::runtime_error(
        "Number of elements is not the same as the number of cells of the "
        "circuit");
  }
  cellIsFixed_ = f;
}

void Circuit::setCellIsObstruction(const std::vector<bool> &f) {
  if (f.size() != nbCells()) {
    throw std::runtime_error(
        "Number of elements is not the same as the number of cells of the "
        "circuit");
  }
  cellIsObstruction_ = f;
}

void Circuit::setCellOrientation(const std::vector<CellOrientation> &orient) {
  if (orient.size() != nbCells()) {
    throw std::runtime_error(
        "Number of elements is not the same as the number of cells of the "
        "circuit");
  }
  cellOrientation_ = orient;
}

void Circuit::setCellWidth(const std::vector<int> &widths) {
  if (widths.size() != nbCells()) {
    throw std::runtime_error(
        "Number of elements is not the same as the number of cells of the "
        "circuit");
  }
  cellWidth_ = widths;
}

void Circuit::setCellHeight(const std::vector<int> &heights) {
  if (heights.size() != nbCells()) {
    throw std::runtime_error(
        "Number of elements is not the same as the number of cells of the "
        "circuit");
  }
  cellHeight_ = heights;
}

int Circuit::placedWidth(int cell) const {
  CellOrientation orient = orientation(cell);
  bool turned = orient == CellOrientation::E || orient == CellOrientation::W ||
                orient == CellOrientation::FW || orient == CellOrientation::FE;
  return turned ? cellHeight_[cell] : cellWidth_[cell];
}

int Circuit::placedHeight(int cell) const {
  CellOrientation orient = orientation(cell);
  bool turned = orient == CellOrientation::E || orient == CellOrientation::W ||
                orient == CellOrientation::FW || orient == CellOrientation::FE;
  return turned ? cellWidth_[cell] : cellHeight_[cell];
}

int Circuit::pinXOffset(int net, int i) const {
  int cell = pinCell(net, i);
  CellOrientation orient = orientation(cell);
  bool turned = orient == CellOrientation::E || orient == CellOrientation::W ||
                orient == CellOrientation::FW || orient == CellOrientation::FE;
  int offs = turned ? pinYOffsets_[netLimits_[net] + i]
                    : pinXOffsets_[netLimits_[net] + i];
  bool flipped = orient == CellOrientation::S || orient == CellOrientation::W ||
                 orient == CellOrientation::FN || orient == CellOrientation::FE;
  return flipped ? placedWidth(cell) - offs : offs;
}

int Circuit::pinYOffset(int net, int i) const {
  int cell = pinCell(net, i);
  CellOrientation orient = orientation(cell);
  bool turned = orient == CellOrientation::E || orient == CellOrientation::W ||
                orient == CellOrientation::FW || orient == CellOrientation::FE;
  int offs = turned ? pinXOffsets_[netLimits_[net] + i]
                    : pinYOffsets_[netLimits_[net] + i];
  bool flipped = orient == CellOrientation::S || orient == CellOrientation::E ||
                 orient == CellOrientation::FS || orient == CellOrientation::FE;
  return flipped ? placedHeight(cell) - offs : offs;
}

long long Circuit::hpwl() const {
  long long ret = 0;
  for (int net = 0; net < nbNets(); ++net) {
    if (nbPinsNet(net) == 0) {
      continue;
    }
    int minX = std::numeric_limits<int>::max();
    int maxX = std::numeric_limits<int>::min();
    int minY = std::numeric_limits<int>::max();
    int maxY = std::numeric_limits<int>::min();
    for (int pin = 0; pin < nbPinsNet(net); ++pin) {
      int cell = pinCell(net, pin);
      int px = x(cell) + pinXOffset(net, pin);
      int py = y(cell) + pinYOffset(net, pin);
      minX = std::min(px, minX);
      maxX = std::max(px, maxX);
      minY = std::min(py, minY);
      maxY = std::max(py, maxY);
    }
    ret += (maxX - minX);
    ret += (maxY - minY);
  }
  return ret;
}

Rectangle Circuit::computePlacementArea() const {
  int minX = std::numeric_limits<int>::max();
  int maxX = std::numeric_limits<int>::min();
  int minY = std::numeric_limits<int>::max();
  int maxY = std::numeric_limits<int>::min();
  if (rows_.empty()) {
    return Rectangle(0, 0, 0, 0);
  }
  for (Rectangle row : rows_) {
    minX = std::min(row.minX, minX);
    maxX = std::max(row.maxX, maxX);
    minY = std::min(row.minY, minY);
    maxY = std::max(row.maxY, maxY);
  }
  return Rectangle(minX, maxX, minY, maxY);
}

std::vector<Rectangle> Circuit::computeRows() const {
  std::vector<Rectangle> obstacles;
  for (int i = 0; i < nbCells(); ++i) {
    if (!isFixed(i)) continue;
    if (!isObstruction(i)) continue;
    obstacles.emplace_back(x(i), x(i) + placedWidth(i), y(i),
                           y(i) + placedHeight(i));
  }
  std::vector<Rectangle> ret;
  for (Rectangle row : rows_) {
    bpl::polygon_90_set_data<int> row_set;
    row_set.insert(
        bpl::rectangle_data<int>(row.minX, row.minY, row.maxX, row.maxY));
    for (Rectangle r : obstacles) {
      row_set.insert(bpl::rectangle_data<int>(r.minX, r.minY, r.maxX, r.maxY),
                     true);
    }

    std::vector<bpl::rectangle_data<int> > diff;
    bpl::get_rectangles(diff, row_set);
    for (auto r : diff) {
      Rectangle newRow(bpl::xl(r), bpl::xh(r), bpl::yl(r), bpl::yh(r));
      // Filter out partially covered rows
      if (newRow.height() == row.height()) {
        ret.push_back(newRow);
      }
    }
  }

  return ret;
}

std::string Circuit::toString() const {
  std::stringstream ss;
  ss << "Circuit with " << nbCells() << " cells, " << nbNets() << " nets and "
     << nbPins() << " pins";
  return ss.str();
}

void Circuit::check() const {
  assert(cellWidth_.size() == nbCells());
  assert(cellHeight_.size() == nbCells());
  assert(cellIsFixed_.size() == nbCells());
  assert(cellIsObstruction_.size() == nbCells());
  assert(cellX_.size() == nbCells());
  assert(cellY_.size() == nbCells());
  assert(cellOrientation_.size() == nbCells());
  assert(!netLimits_.empty());
  assert(netLimits_.front() == 0);
  assert(pinCells_.size() == nbPins());
  assert(pinXOffsets_.size() == nbPins());
  assert(pinYOffsets_.size() == nbPins());
}

void Circuit::placeGlobal(const GlobalPlacerParameters &params) {
  GlobalPlacer::place(*this, params);
}

void Circuit::placeDetailed(const DetailedPlacerParameters &params) {
  DetailedPlacer::place(*this, params);
}

extern "C" {
int place_ispd(int nb_cells, int nb_nets, int *cell_widths, int *cell_heights,
               char *cell_fixed, int *net_limits, int *pin_cells,
               int *pin_x_offsets, int *pin_y_offsets, int *cell_x, int *cell_y,
               int *cell_orientation, int nb_rows, int *row_min_x,
               int *row_max_x, int *row_min_y, int *row_max_y, int effort) {
  Circuit circuit(nb_cells);

  try {
    circuit.setCellWidth(std::vector<int>(cell_widths, cell_widths + nb_cells));
    circuit.setCellHeight(
        std::vector<int>(cell_heights, cell_heights + nb_cells));
    std::vector<bool> cell_fixed_vec;
    for (char *f = cell_fixed; f != cell_fixed + nb_cells; ++f) {
      cell_fixed_vec.push_back(*f);
    }
    circuit.setCellIsFixed(cell_fixed_vec);
    int nb_pins = net_limits[nb_nets];
    circuit.setNets(std::vector<int>(net_limits, net_limits + nb_nets + 1),
                    std::vector<int>(pin_cells, pin_cells + nb_pins),
                    std::vector<int>(pin_x_offsets, pin_x_offsets + nb_pins),
                    std::vector<int>(pin_y_offsets, pin_y_offsets + nb_pins));
    circuit.setCellX(std::vector<int>(cell_x, cell_x + nb_cells));
    circuit.setCellY(std::vector<int>(cell_y, cell_y + nb_cells));
    std::vector<CellOrientation> orient;
    for (int i = 0; i < nb_cells; ++i) {
      orient.push_back(static_cast<CellOrientation>(cell_orientation[i]));
    }
    circuit.setCellOrientation(orient);
    std::vector<Rectangle> rows;
    for (int i = 0; i < nb_rows; ++i) {
      rows.emplace_back(row_min_x[i], row_max_x[i], row_min_y[i], row_max_y[i]);
    }
    circuit.setRows(rows);
  } catch (const std::exception &e) {
    std::cout << "initialization terminated with exception: " << e.what()
              << std::endl;
    return -1;
  }

  try {
    std::cout << "Placing circuit with " << circuit.nbCells() << " cells, "
              << circuit.nbNets() << " nets and " << circuit.nbPins()
              << " pins." << std::endl;
    circuit.place(effort);
  } catch (const std::exception &e) {
    std::cout << "placement terminated with exception: " << e.what()
              << std::endl;
    return -1;
  }

  for (int i = 0; i < nb_cells; ++i) {
    cell_x[i] = circuit.cellX_[i];
    cell_y[i] = circuit.cellY_[i];
  }

  // Return with success.
  return 0;
}
}
}  // namespace coloquinte
