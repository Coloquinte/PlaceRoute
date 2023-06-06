
#include "coloquinte.hpp"

#include <boost/polygon/polygon.hpp>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <unordered_set>
#include <utility>

#include "place_detailed/place_detailed.hpp"
#include "place_global/density_legalizer.hpp"
#include "place_global/net_model.hpp"
#include "place_global/place_global.hpp"
#include "utils/norm.hpp"

namespace bpl = boost::polygon;
namespace coloquinte {

Circuit::Circuit(int nbCells) {
  cellWidth_.resize(nbCells);
  cellHeight_.resize(nbCells);
  cellIsFixed_.resize(nbCells, false);
  cellIsObstruction_.resize(nbCells, true);
  cellRowPolarity_.resize(nbCells, CellRowPolarity::ANY);
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
  if (cells.empty()) {
    return;
  }
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

void Circuit::setCellRowPolarity(const std::vector<CellRowPolarity> &orient) {
  if (orient.size() != nbCells()) {
    throw std::runtime_error(
        "Number of elements is not the same as the number of cells of the "
        "circuit");
  }
  cellRowPolarity_ = orient;
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

void Circuit::setSolution(const PlacementSolution &sol) {
  if (sol.size() != nbCells()) {
    throw std::runtime_error(
        "Number of elements is not the same as the number of cells of the "
        "circuit");
  }
  cellX_.clear();
  cellY_.clear();
  cellOrientation_.clear();
  for (auto [p, o] : sol) {
    cellX_.push_back(p.x);
    cellY_.push_back(p.y);
    cellOrientation_.push_back(o);
  }
}

int Circuit::placedWidth(int cell) const {
  CellOrientation orient = orientation(cell);
  return isTurn(orient) ? cellHeight_[cell] : cellWidth_[cell];
}

int Circuit::placedHeight(int cell) const {
  CellOrientation orient = orientation(cell);
  return isTurn(orient) ? cellWidth_[cell] : cellHeight_[cell];
}

int Circuit::pinXOffset(int net, int i) const {
  int cell = pinCell(net, i);
  CellOrientation orient = orientation(cell);
  int offs = isTurn(orient) ? pinYOffsets_[netLimits_[net] + i]
                            : pinXOffsets_[netLimits_[net] + i];
  bool flipped = orient == CellOrientation::S || orient == CellOrientation::W ||
                 orient == CellOrientation::FN || orient == CellOrientation::FE;
  return flipped ? placedWidth(cell) - offs : offs;
}

int Circuit::pinYOffset(int net, int i) const {
  int cell = pinCell(net, i);
  CellOrientation orient = orientation(cell);
  int offs = isTurn(orient) ? pinXOffsets_[netLimits_[net] + i]
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

std::vector<Row> Circuit::computeRows(
    const std::vector<Rectangle> &additionalObstacles) const {
  std::vector<Rectangle> obstacles = additionalObstacles;
  for (int i = 0; i < nbCells(); ++i) {
    if (!isFixed(i)) {
      continue;
    }
    if (!isObstruction(i)) {
      continue;
    }
    obstacles.emplace_back(placement(i));
  }
  // Use boost::polygon ro compute the difference of each row to every other
  std::vector<Row> ret;
  for (Row row : rows_) {
    auto o = row.freespace(obstacles);
    ret.insert(ret.end(), o.begin(), o.end());
  }
  return ret;
}

std::vector<Row> Row::freespace(const std::vector<Rectangle> &obstacles) const {
  std::vector<Row> ret;
  bpl::polygon_90_set_data<int> row_set;
  row_set.insert(bpl::rectangle_data<int>(minX, minY, maxX, maxY));
  for (Rectangle r : obstacles) {
    row_set.insert(bpl::rectangle_data<int>(r.minX, r.minY, r.maxX, r.maxY),
                   true);
  }

  std::vector<bpl::rectangle_data<int> > diff;
  bpl::get_rectangles(diff, row_set);
  for (const auto &r : diff) {
    Rectangle newRow(bpl::xl(r), bpl::xh(r), bpl::yl(r), bpl::yh(r));
    // Filter out partially covered rows
    if (newRow.height() == height()) {
      ret.emplace_back(newRow, orientation);
    }
  }
  return ret;
}

int Circuit::rowHeight() const {
  if (nbRows() == 0) {
    throw std::runtime_error(
        "Cannot compute row height as no row has been defined");
  }
  int ret = rows_[0].height();
  for (Rectangle row : rows_) {
    if (row.height() != ret) {
      throw std::runtime_error(
          "The circuit contains rows of different heights");
    }
  }
  return ret;
}

void Circuit::setupRows(Rectangle placementArea, int rowHeight,
                        bool alternatingOrientation) {
  if (rowHeight <= 0) {
    throw std::runtime_error("Row height for row creation must be positive");
  }
  rows_.clear();
  bool orient = true;
  for (int y = placementArea.minY; y + rowHeight < placementArea.maxY;
       y += rowHeight) {
    CellOrientation dir = orient ? CellOrientation::N : CellOrientation::FS;
    if (alternatingOrientation) {
      orient = !orient;
    }
    rows_.emplace_back(placementArea.minX, placementArea.maxX, y, y + rowHeight,
                       dir);
  }
}

std::string Circuit::toString() const {
  std::stringstream ss;
  ss << "Circuit with " << nbCells() << " cells, " << nbNets() << " nets and "
     << nbPins() << " pins";
  return ss.str();
}

void Circuit::check() const {
  if (cellWidth_.size() != nbCells()) {
    throw std::runtime_error("Size mismatch");
  }
  if (cellHeight_.size() != nbCells()) {
    throw std::runtime_error("Size mismatch");
  }
  if (cellIsFixed_.size() != nbCells()) {
    throw std::runtime_error("Size mismatch");
  }
  if (cellIsObstruction_.size() != nbCells()) {
    throw std::runtime_error("Size mismatch");
  }
  if (cellX_.size() != nbCells()) {
    throw std::runtime_error("Size mismatch");
  }
  if (cellY_.size() != nbCells()) {
    throw std::runtime_error("Size mismatch");
  }
  if (cellOrientation_.size() != nbCells()) {
    throw std::runtime_error("Size mismatch");
  }
  if (netLimits_.empty()) {
    throw std::runtime_error("Size mismatch");
  }
  if (netLimits_.front() != 0) {
    throw std::runtime_error("Size mismatch");
  }
  if (pinCells_.size() != nbPins()) {
    throw std::runtime_error("Size mismatch");
  }
  if (pinXOffsets_.size() != nbPins()) {
    throw std::runtime_error("Size mismatch");
  }
  if (pinYOffsets_.size() != nbPins()) {
    throw std::runtime_error("Size mismatch");
  }
}

std::string Circuit::report() const {
  int stdCellHeight = rowHeight();
  int nbMacros = 0;
  int nbSingleRowCells = 0;
  int nbMultiRowCells = 0;
  int nbUnalignedCells = 0;
  int nbPlaceableMacros = 0;
  long long singleRowCellArea = 0;
  long long multiRowCellArea = 0;
  long long unalignedCellArea = 0;
  long long placeableMacroArea = 0;
  for (int i = 0; i < nbCells(); ++i) {
    if (isFixed(i)) {
      if (isObstruction(i)) {
        ++nbMacros;
      }
    } else if (cellHeight_[i] > 4 * stdCellHeight) {
      ++nbPlaceableMacros;
      placeableMacroArea += area(i);
    } else if (cellHeight_[i] % stdCellHeight != 0) {
      ++nbUnalignedCells;
      unalignedCellArea += area(i);
    } else if (cellHeight_[i] == stdCellHeight) {
      ++nbSingleRowCells;
      singleRowCellArea += area(i);
    } else {
      ++nbMultiRowCells;
      multiRowCellArea += area(i);
    }
  }
  int nbPlaceableCells =
      nbSingleRowCells + nbMultiRowCells + nbUnalignedCells + nbPlaceableMacros;
  long long placeableCellArea = singleRowCellArea + multiRowCellArea +
                                unalignedCellArea + placeableMacroArea;
  int nbFixedPins = 0;
  for (int i = 0; i < nbNets(); ++i) {
    for (int j = 0; j < nbPinsNet(i); ++j) {
      int c = pinCell(i, j);
      if (isFixed(c)) {
        ++nbFixedPins;
      }
    }
  }
  int nbMoveablePins = nbPins() - nbFixedPins;
  long long availableArea = 0;
  for (Rectangle row : computeRows()) {
    availableArea += row.area();
  }
  long long totalArea = 0;
  for (Rectangle row : rows_) {
    totalArea += row.area();
  }
  std::stringstream ss;
  ss << std::fixed << std::setprecision(1);
  ss << "Circuit report:\n";

  /**
   * Report cells
   */
  ss << "\t" << nbPlaceableCells << " cells (" << nbSingleRowCells
     << " single-row cells";
  if (nbMultiRowCells > 0) {
    ss << ", " << nbMultiRowCells << " multi-row cells";
  }
  if (nbUnalignedCells > 0) {
    ss << ", " << nbUnalignedCells << " unaligned cells";
  }
  if (nbPlaceableMacros > 0) {
    ss << ", " << nbPlaceableMacros << " placeable macro blocks";
  } else {
    ss << ", no placeable macro block";
  }
  ss << ")\n";

  /**
   * Report nets
   */
  ss << "\t" << nbNets() << " nets, " << nbMoveablePins << " pins + "
     << nbFixedPins << " fixed"
     << ", " << (float)nbPins() / nbNets() << " pins/net\n";

  /**
   * Report rows
   */
  std::unordered_set<CellOrientation> orients;
  for (Row r : rows_) {
    orients.insert(r.orientation);
  }
  ss << "\t" << nbRows() << " rows ";
  if (orients.size() <= 1) {
    ss << "(single orientation)\n";
  } else {
    ss << "(" << orients.size() << " orientations)\n";
  }

  /**
   * Report area breakdown
   */
  ss << "\t";
  if (totalArea != availableArea) {
    ss << 100.0 * (totalArea - availableArea) / totalArea
       << "% fixed macro blocks, ";
  }
  ss << 100.0 * placeableCellArea / availableArea << "% density";
  if (multiRowCellArea > 0 || placeableMacroArea > 0) {
    ss << " (" << 100.0 * singleRowCellArea / availableArea
       << "% single-row cells";
    if (multiRowCellArea > 0) {
      ss << ", " << 100.0 * multiRowCellArea / availableArea
         << "% multi-row cells";
    }
    if (unalignedCellArea > 0) {
      ss << ", " << 100.0 * unalignedCellArea / availableArea
         << "% unaligned cells";
    }
    if (placeableMacroArea > 0) {
      ss << ", " << 100.0 * placeableMacroArea / availableArea
         << "% placeable macro blocks";
    }
    ss << ")";
  }
  ss << "\n";

  // Analyze and report possible issues
  if (nbUnalignedCells > 0) {
    ss << "WARNING: some cells are not aligned on a standard cell heights, "
       << "which is not supported.\n";
  }
  if (nbPlaceableMacros > 0) {
    ss << std::endl;
    ss << "WARNING: movable macros are present; this is supported but would "
          "benefit from quality improvements.\n";
  }

  bool badRows = false;
  for (CellOrientation o : {CellOrientation::E, CellOrientation::W,
                            CellOrientation::FE, CellOrientation::FW}) {
    if (orients.count(o)) {
      badRows = true;
    }
  }
  if (badRows) {
    ss << "WARNING: some rows have orientations other than N/S. This is not "
          "supported.\n";
  }

  bool badRowHeight = false;
  bool badRowOrient = false;
  for (int i = 0; i + 1 < rows_.size(); ++i) {
    if (rows_[i].height() != rows_[i + 1].height()) {
      badRowHeight = true;
    }
    if (rows_[i].minY == rows_[i + 1].minY &&
        rows_[i].orientation != rows_[i + 1].orientation) {
      badRowOrient = true;
    }
  }
  if (badRowHeight) {
    ss << "WARNING: some rows have different heights. This is not supported.\n";
  }
  if (badRowOrient) {
    ss << "WARNING: some rows at the same y have different orientations. This "
          "is not supported.\n";
  }

  int oddRowWithBadPolarity = 0;
  int evenRowWithBadPolarity = 0;
  int standardCellWithBadPolarity = 0;
  for (int i = 0; i < nbCells(); ++i) {
    if (isFixed(i)) {
      continue;
    }
    if (cellHeight_[i] > 4 * stdCellHeight) {
      continue;
    }
    CellRowPolarity pol = cellRowPolarity_[i];
    if (pol == CellRowPolarity::ANY) {
      ++standardCellWithBadPolarity;
      continue;
    }
    int nbRows = cellHeight_[i] / stdCellHeight;
    if (nbRows % 2 != 0) {
      if (pol != CellRowPolarity::OPPOSITE && pol != CellRowPolarity::SAME) {
        ++oddRowWithBadPolarity;
      }
    } else {
      if (pol != CellRowPolarity::NW && pol != CellRowPolarity::SE) {
        ++evenRowWithBadPolarity;
      }
    }
  }

  if (standardCellWithBadPolarity) {
    ss << "WARNING: " << standardCellWithBadPolarity
       << " small cells have no polarity specified with respect to rows.\n";
  }
  if (oddRowWithBadPolarity) {
    ss << "WARNING: " << oddRowWithBadPolarity
       << " cells occupy an odd number of rows but are only allowed some rows "
          "(NW or SE polarity).\n";
  }
  if (evenRowWithBadPolarity) {
    ss << "WARNING: " << evenRowWithBadPolarity
       << " cells occupy an eveb number of rows but are allowed to go to any "
          "row (SAME or OPPOSITE polarity).\n";
  }

  return ss.str();
}

void Circuit::placeGlobal(const ColoquinteParameters &params,
                          const std::optional<PlacementCallback> &callback) {
  GlobalPlacer::place(*this, params, callback);
}

void Circuit::legalize(const ColoquinteParameters &params,
                       const std::optional<PlacementCallback> &callback) {
  DetailedPlacer::legalize(*this, params, callback);
}

void Circuit::placeDetailed(const ColoquinteParameters &params,
                            const std::optional<PlacementCallback> &callback) {
  DetailedPlacer::place(*this, params, callback);
}

float Circuit::meanDisruption(const PlacementSolution &a,
                              const PlacementSolution &b,
                              LegalizationModel costModel) {
  auto dist = allDistances(a, b, costModel);
  float disp = 0.0f;
  float totArea = 0.0f;
  for (int i = 0; i < nbCells(); ++i) {
    disp += area(i) * dist[i];
    totArea += area(i);
  }

  disp /= totArea;
  if (costModel == LegalizationModel::L1Squared ||
      costModel == LegalizationModel::L2Squared ||
      costModel == LegalizationModel::LInfSquared) {
    disp = std::sqrt(disp);
  }
  return disp;
}

float Circuit::rmsDisruption(const PlacementSolution &a,
                             const PlacementSolution &b,
                             LegalizationModel costModel) {
  auto dist = allDistances(a, b, costModel);
  float disp = 0.0f;
  float totArea = 0.0f;
  for (int i = 0; i < nbCells(); ++i) {
    disp += area(i) * dist[i] * dist[i];
    totArea += area(i);
  }
  disp = std::sqrt(disp / totArea);
  if (costModel == LegalizationModel::L1Squared ||
      costModel == LegalizationModel::L2Squared ||
      costModel == LegalizationModel::LInfSquared) {
    disp = std::sqrt(disp);
  }
  return disp;
}
float Circuit::maxDisruption(const PlacementSolution &a,
                             const PlacementSolution &b,
                             LegalizationModel costModel) {
  auto dist = allDistances(a, b, costModel);
  return *std::max_element(dist.begin(), dist.end());
}

std::vector<float> Circuit::allDistances(const PlacementSolution &a,
                                const PlacementSolution &b,
                                LegalizationModel costModel) {
  if (a.size() != nbCells() || b.size() != nbCells()) {
    throw std::runtime_error("Solution size doesn't match number of cells");
  }
  std::vector<float> ret;
  for (int i = 0; i < nbCells(); ++i) {
    Point pa = a[i].position;
    Point pb = b[i].position;
    ret.push_back(norm((float) (pa.x - pb.x), (float) (pa.y - pb.y), costModel));
  }
  return ret;
}
}  // namespace coloquinte
