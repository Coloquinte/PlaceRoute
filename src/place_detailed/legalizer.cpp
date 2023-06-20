
#include "place_detailed/legalizer.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>

#include "place_detailed/abacus_legalizer.hpp"
#include "place_detailed/tetris_legalizer.hpp"
#include "utils/norm.hpp"

namespace coloquinte {
Legalizer Legalizer::fromIspdCircuit(const Circuit &circuit) {
  // Represent fixed cells with -1 width so they are not considered
  int rowHeight = circuit.rowHeight();
  std::vector<int> widths;
  std::vector<int> heights;
  std::vector<int> x;
  std::vector<int> y;
  std::vector<CellRowPolarity> polarities;
  std::vector<CellOrientation> orient;
  for (int i = 0; i < circuit.nbCells(); ++i) {
    if (circuit.cellIsFixed_[i]) {
      continue;
    }
    widths.push_back(circuit.placedWidth(i));
    heights.push_back(circuit.placedHeight(i));
    polarities.push_back(circuit.cellRowPolarity_[i]);
    x.push_back(circuit.cellX_[i]);
    y.push_back(circuit.cellY_[i]);
    orient.push_back(circuit.cellOrientation_[i]);
  }
  return Legalizer(circuit.computeRows(), widths, heights, polarities, x, y,
                   orient);
}

LegalizerBase::LegalizerBase(
    const std::vector<Row> &rows, const std::vector<int> &width,
    const std::vector<int> &height,
    const std::vector<CellRowPolarity> &polarities,
    const std::vector<int> &targetX, const std::vector<int> &targetY,
    const std::vector<CellOrientation> &targetOrientation)
    : cellWidth_(width),
      cellHeight_(height),
      cellRowPolarity_(polarities),
      cellTargetX_(targetX),
      cellTargetY_(targetY),
      cellTargetOrientation_(targetOrientation) {
  assert(width.size() == height.size());
  assert(width.size() == polarities.size());
  assert(width.size() == targetX.size());
  assert(width.size() == targetY.size());
  assert(width.size() == targetOrientation.size());
  // Sort the rows
  rows_ = rows;
  std::stable_sort(rows_.begin(), rows_.end(), [](Row a, Row b) -> bool {
    return a.minY < b.minY || (a.minY == b.minY && a.minX < b.minX);
  });
  cellToX_ = cellTargetX_;
  cellToY_ = cellTargetY_;
  cellToOrientation_ = cellTargetOrientation_;
  cellIsPlaced_.assign(width.size(), false);
}

void LegalizerBase::check() const {
  if (cellWidth_.size() != nbCells()) {
    throw std::runtime_error("Number of cell widths does not match");
  }
  if (cellHeight_.size() != nbCells()) {
    throw std::runtime_error("Number of cell heights does not match");
  }
  if (cellTargetX_.size() != nbCells()) {
    throw std::runtime_error("Number of cell x targets does not match");
  }
  if (cellTargetY_.size() != nbCells()) {
    throw std::runtime_error("Number of cell y targets does not match");
  }
  if (cellTargetOrientation_.size() != nbCells()) {
    throw std::runtime_error(
        "Number of cell orientation targets does not match");
  }
  if (cellToX_.size() != nbCells()) {
    throw std::runtime_error("Number of cell x positions does not match");
  }
  if (cellToY_.size() != nbCells()) {
    throw std::runtime_error("Number of cell y positions does not match");
  }
  if (cellToOrientation_.size() != nbCells()) {
    throw std::runtime_error("Number of cell orientations does not match");
  }
  for (Row r : rows_) {
    if (r.height() != rowHeight()) {
      throw std::runtime_error("Rows have different heights");
    }
  }
}

void LegalizerBase::checkAllPlaced() const {
  for (int i = 0; i < nbCells(); ++i) {
    if (!isPlaced(i)) {
      throw std::runtime_error("Not all cells have been placed");
    }
  }
}

void LegalizerBase::importLegalization(const LegalizerBase &leg,
                                       const std::vector<int> &cells) {
  std::vector<int> x = leg.cellLegalX();
  std::vector<int> y = leg.cellLegalY();
  std::vector<CellOrientation> o = leg.cellLegalOrientation();
  assert(cells.size() == leg.nbCells());
  for (int i = 0; i < cells.size(); ++i) {
    int c = cells[i];
    if (leg.cellIsPlaced_[i]) {
      cellToX_[c] = x[i];
      cellToY_[c] = y[i];
      cellToOrientation_[c] = o[i];
      cellIsPlaced_[c] = true;
    }
  }
}

std::vector<int> LegalizerBase::computeCellOrder(float weightX,
                                                 float weightWidth,
                                                 float weightY,
                                                 float weightHeight) const {
  // Sort the cells by target X coordinate
  std::vector<std::pair<float, int> > sortedCells;
  for (int i = 0; i < nbCells(); ++i) {
    float val = weightX * cellTargetX_[i] + weightWidth * cellWidth_[i] +
                weightY * cellTargetY_[i] + weightHeight * cellHeight_[i];
    sortedCells.emplace_back(val, i);
  }
  std::stable_sort(sortedCells.begin(), sortedCells.end());
  std::vector<int> cells;
  cells.reserve(sortedCells.size());

  for (auto p : sortedCells) {
    cells.push_back(p.second);
  }
  return cells;
}

int LegalizerBase::closestRow(int y) const {
  auto it = std::lower_bound(rows_.begin(), rows_.end(), y,
                             [](Rectangle r, int v) { return r.minY < v; });
  if (it == rows_.end()) {
    return nbRows() - 1;
  }
  if (it == rows_.begin()) {
    return 0;
  }
  int row = it - rows_.begin();
  if (row > 0 && rows_[row].minY - y > y - rows_[row - 1].minY) {
    return row - 1;
  }
  return row;
}

std::vector<Row> LegalizerBase::remainingRows() const {
  std::vector<Rectangle> obstacles;
  for (int i = 0; i < nbCells(); ++i) {
    if (!isPlaced(i)) {
      continue;
    }
    obstacles.emplace_back(cellToX_[i], cellToX_[i] + cellWidth_[i],
                           cellToY_[i], cellToY_[i] + cellHeight_[i]);
  }
  std::vector<Row> ret;
  for (Row r : rows_) {
    auto o = r.freespace(obstacles);
    ret.insert(ret.end(), o.begin(), o.end());
  }
  return ret;
}

const std::vector<int> &LegalizerBase::cellLegalX() const { return cellToX_; }

const std::vector<int> &LegalizerBase::cellLegalY() const { return cellToY_; }

const std::vector<CellOrientation> &LegalizerBase::cellLegalOrientation()
    const {
  return cellToOrientation_;
}

void Legalizer::exportPlacement(Circuit &circuit) {
  std::vector<int> cellX = cellLegalX();
  std::vector<int> cellY = cellLegalY();
  std::vector<CellOrientation> cellOrient = cellLegalOrientation();
  int j = 0;
  for (int i = 0; i < circuit.nbCells(); ++i) {
    if (circuit.cellIsFixed_[i]) {
      continue;
    }
    if (j >= nbCells()) {
      throw std::runtime_error("Circuit does not match legalizer for export");
    }
    if (isPlaced(j)) {
      circuit.cellX_[i] = cellX[j];
      circuit.cellY_[i] = cellY[j];
      circuit.cellOrientation_[i] = cellOrient[j];
    }
    ++j;
  }
}

std::vector<float> LegalizerBase::allDistances(LegalizationModel model) const {
  std::vector<int> cellX = cellLegalX();
  std::vector<int> cellY = cellLegalY();
  std::vector<int> targetX = cellTargetX_;
  std::vector<int> targetY = cellTargetY_;
  std::vector<float> distances;
  distances.reserve(nbCells());
  for (int i = 0; i < nbCells(); ++i) {
    float dx = targetX[i] - cellX[i];
    float dy = targetY[i] - cellY[i];
    distances.push_back(norm(dx, dy, model));
  }
  return distances;
}

float LegalizerBase::meanDistance(LegalizationModel model) const {
  std::vector<float> dist = allDistances(model);
  float disp = 0.0f;
  for (int i = 0; i < nbCells(); ++i) {
    disp += dist[i] * cellWidth_[i] * cellHeight_[i];
  }
  disp /= totalCellArea();
  if (model == LegalizationModel::L1Squared ||
      model == LegalizationModel::L2Squared ||
      model == LegalizationModel::LInfSquared) {
    disp = std::sqrt(disp);
  }
  return disp;
}

float LegalizerBase::rmsDistance(LegalizationModel model) const {
  std::vector<float> dist = allDistances(model);
  float disp = 0.0f;
  for (int i = 0; i < nbCells(); ++i) {
    disp += dist[i] * dist[i] * cellWidth_[i] * cellHeight_[i];
  }
  disp = std::sqrt(disp / totalCellArea());
  if (model == LegalizationModel::L1Squared ||
      model == LegalizationModel::L2Squared ||
      model == LegalizationModel::LInfSquared) {
    disp = std::sqrt(disp);
  }
  return disp;
}

float LegalizerBase::maxDistance(LegalizationModel model) const {
  std::vector<float> dist = allDistances(model);
  return *std::max_element(dist.begin(), dist.end());
}

long long LegalizerBase::totalCellArea() const {
  long long ret = 0;
  for (int c = 0; c < nbCells(); ++c) {
    ret += static_cast<long long>(cellWidth_[c]) *
           static_cast<long long>(cellHeight_[c]);
  }
  return ret;
}

CellOrientation LegalizerBase::getOrientation(int cell, int row) const {
  CellRowPolarity pol = cellRowPolarity_[cell];
  CellOrientation rowOrientation = rows_[row].orientation;
  CellOrientation orient = cellOrientationInRow(pol, rowOrientation);
  if (orient == CellOrientation::UNKNOWN) {
    // Keep the same orientation
    return cellTargetOrientation_[cell];
  }
  return orient;
}

int LegalizerBase::rowHeight() const {
  if (rows_.empty()) {
    throw std::runtime_error("No row present");
  }
  return rows_.front().height();
}

Legalizer::Legalizer(const std::vector<Row> &rows,
                     const std::vector<int> &width,
                     const std::vector<int> &height,
                     const std::vector<CellRowPolarity> &polarities,
                     const std::vector<int> &targetX,
                     const std::vector<int> &targetY,
                     const std::vector<CellOrientation> &targetOrientation)
    : LegalizerBase(rows, width, height, polarities, targetX, targetY,
                    targetOrientation) {}

void Legalizer::run(const ColoquinteParameters &params) {
  std::vector<int> cellOrder = computeCellOrder(
      1.0, params.legalization.orderingWidth, params.legalization.orderingY,
      params.legalization.orderingHeight);
  // Run the Tetris legalizer on the macros and large cells
  runTetris(cellOrder);
  // Run the Abacus legalizer on the remaining cells
  runAbacus(cellOrder);
  // Check that everything is legalized
  checkAllPlaced();
}

void Legalizer::runTetris(const std::vector<int> &cells) {
  std::vector<Row> r = remainingRows();
  std::vector<int> w, h, x, y, remainingCells;
  std::vector<CellRowPolarity> p;
  std::vector<CellOrientation> o;
  for (int c : cells) {
    if (isPlaced(c)) {
      continue;
    }
    if (cellHeight_[c] <= rowHeight()) {
      continue;
    }
    remainingCells.push_back(c);
    w.push_back(cellWidth_[c]);
    h.push_back(cellHeight_[c]);
    p.push_back(cellRowPolarity_[c]);
    x.push_back(cellTargetX_[c]);
    y.push_back(cellTargetY_[c]);
    o.push_back(cellTargetOrientation_[c]);
  }

  TetrisLegalizer leg(r, w, h, p, x, y, o);
  leg.run();
  importLegalization(leg, remainingCells);
}

void Legalizer::runAbacus(const std::vector<int> &cells) {
  std::vector<Row> r = remainingRows();
  std::vector<int> w, h, x, y, remainingCells;
  std::vector<CellRowPolarity> p;
  std::vector<CellOrientation> o;
  for (int c : cells) {
    if (isPlaced(c)) {
      continue;
    }
    if (cellHeight_[c] != rowHeight()) {
      continue;
    }
    remainingCells.push_back(c);
    w.push_back(cellWidth_[c]);
    h.push_back(cellHeight_[c]);
    p.push_back(cellRowPolarity_[c]);
    x.push_back(cellTargetX_[c]);
    y.push_back(cellTargetY_[c]);
    o.push_back(cellTargetOrientation_[c]);
  }

  AbacusLegalizer leg(r, w, h, p, x, y, o);
  leg.run();
  importLegalization(leg, remainingCells);
}
}  // namespace coloquinte