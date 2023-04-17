
#include "place_detailed/legalizer.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>

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
  for (int i = 0; i < circuit.nbCells(); ++i) {
    if (circuit.cellIsFixed_[i]) {
      continue;
    }
    if (circuit.cellHeight_[i] != rowHeight) {
      throw std::runtime_error(
          "Some placeable cells have a height that is different from the row "
          "height");
    }
    widths.push_back(circuit.cellWidth_[i]);
    heights.push_back(circuit.cellHeight_[i]);
    polarities.push_back(circuit.cellRowPolarity_[i]);
    x.push_back(circuit.cellX_[i]);
    y.push_back(circuit.cellY_[i]);
  }
  return Legalizer(circuit.computeRows(), widths, heights, polarities, x, y);
}

Legalizer::Legalizer(const std::vector<Row> &rows,
                     const std::vector<int> &width,
                     const std::vector<int> &height,
                     const std::vector<CellRowPolarity> &polarities,
                     const std::vector<int> &targetX,
                     const std::vector<int> &targetY)
    : cellWidth_(width),
      cellHeight_(height),
      cellRowPolarity_(polarities),
      cellTargetX_(targetX),
      cellTargetY_(targetY) {
  assert(width.size() == height.size());
  assert(width.size() == polarities.size());
  assert(width.size() == targetX.size());
  assert(width.size() == targetY.size());
  // Sort the rows
  rows_ = rows;
  std::stable_sort(rows_.begin(), rows_.end(), [](Row a, Row b) -> bool {
    return a.minY < b.minY || (a.minY == b.minY && a.minX < b.minX);
  });
  for (const Row &row : rows_) {
    rowLegalizers_.emplace_back(row.minX, row.maxX);
  }
  rowToCells_.resize(rows_.size());
  cellToX_ = cellTargetX_;
  cellToY_ = cellTargetY_;
  cellToRow_.assign(width.size(), -1);
}

void Legalizer::run(const ColoquinteParameters &params) {
  std::vector<int> cellOrder = computeCellOrder(
      1.0, params.legalization.orderingWidth, params.legalization.orderingY,
      params.legalization.orderingHeight);

  for (int c : cellOrder) {
    placeCellOptimally(c, params.legalization.costModel);
  }
  for (int i = 0; i < nbRows(); ++i) {
    std::vector<int> pl = rowLegalizers_[i].getPlacement();
    assert(pl.size() == rowToCells_[i].size());
    for (int j = 0; j < pl.size(); ++j) {
      int cell = rowToCells_[i][j];
      cellToX_[cell] = pl[j];
      cellToY_[cell] = rows_[i].minY;
      cellToRow_[cell] = i;
    }
  }
  for (RowLegalizer &leg : rowLegalizers_) {
    leg.clear();
  }
  check();
}

void Legalizer::check() const {
  if (cellWidth_.size() != nbCells()) {
    throw std::runtime_error("Number of cell widths does not match");
  }
  if (cellTargetX_.size() != nbCells()) {
    throw std::runtime_error("Number of cell x targets does not match");
  }
  if (cellTargetY_.size() != nbCells()) {
    throw std::runtime_error("Number of cell y targets does not match");
  }
  if (cellToX_.size() != nbCells()) {
    throw std::runtime_error("Number of cell x positions does not match");
  }
  if (cellToY_.size() != nbCells()) {
    throw std::runtime_error("Number of cell y positions does not match");
  }
  if (cellToRow_.size() != nbCells()) {
    throw std::runtime_error("Number of cell row positions does not match");
  }
  if (rowToCells_.size() != nbRows()) {
    throw std::runtime_error("Number of row cells does not match");
  }
  if (rowLegalizers_.size() != nbRows()) {
    throw std::runtime_error("Number of row legalizers does not match");
  }
  for (int i = 0; i < nbRows(); ++i) {
    for (int c : rowToCells_[i]) {
      if (cellToRow_[c] != i) {
        throw std::runtime_error(
            "Cell allocation does not match row allocation");
      }
      if (cellToX_[c] < rows_[i].minX) {
        throw std::runtime_error("Cell placed before the row");
      }
      if (cellToX_[c] + cellWidth_[c] > rows_[i].maxX) {
        throw std::runtime_error("Cell placed after the row");
      }
    }
  }
  for (int i = 0; i < nbRows(); ++i) {
    for (int j = 0; j + 1 < rowToCells_[i].size(); ++j) {
      int c1 = rowToCells_[i][j];
      int c2 = rowToCells_[i][j + 1];
      if (cellToX_[c1] + cellWidth_[c1] > cellToX_[c2]) {
        throw std::runtime_error("Cell overlap detected");
      }
    }
  }
}

bool Legalizer::placeCellOptimally(int cell, LegalizationModel costModel) {
  /**
   * Simple algorithm that tries close row first and stops early if no
   * improvement can be found
   */
  int targetX = cellTargetX_[cell];
  int targetY = cellTargetY_[cell];
  int bestX = 0;
  int bestRow = -1;
  long long bestDist = std::numeric_limits<long long>::max();

  auto tryPlace = [&](int row) {
    long long yDist =
        cellWidth_[cell] * norm(0, rows_[row].minY - targetY, costModel);
    if (bestRow != -1 && yDist > bestDist) {
      // Not possible to do better since the rows are sorted
      return true;
    }
    // Find the best position for the cell
    auto [ok, xDist] = placeCellOptimally(cell, row);
    // TODO: extend this to non-L1 cases
    long long dist = xDist + yDist;
    if (!ok) {
      // Not possible to place in this row, but cannot stop yet
      return false;
    }
    if (bestRow == -1 || dist < bestDist) {
      bestRow = row;
      bestDist = dist;
    }
    // Cannot stop yet
    return false;
  };

  // Try promising candidates first
  int initialRow = closestRow(targetY);
  for (int row = initialRow; row < nbRows(); ++row) {
    bool canStop = tryPlace(row);
    if (canStop) {
      break;
    }
  }
  for (int row = initialRow - 1; row >= 0; --row) {
    bool canStop = tryPlace(row);
    if (canStop) {
      break;
    }
  }

  if (bestRow == -1) {
    throw std::runtime_error(
        "Unable to place a cell with the greedy legalization algorithm");
  }
  rowLegalizers_[bestRow].push(cellWidth_[cell], targetX);
  rowToCells_[bestRow].push_back(cell);
  cellToRow_[cell] = bestRow;
  return true;
}

std::pair<bool, long long> Legalizer::placeCellOptimally(int cell, int row) {
  if (rowLegalizers_[row].remainingSpace() < cellWidth_[cell]) {
    return std::make_pair(false, 0);
  }
  int dist = rowLegalizers_[row].getCost(cellWidth_[cell], cellTargetX_[cell]);
  return std::make_pair(true, dist);
}

std::vector<int> Legalizer::computeCellOrder(float weightX, float weightWidth,
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

int Legalizer::closestRow(int y) const {
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

std::vector<int> Legalizer::cellLegalX() const {
  std::vector<int> ret(nbCells());
  for (int r = 0; r < nbRows(); ++r) {
    for (int c : rowToCells_[r]) {
      ret[c] = cellToX_[c];
    }
  }
  return ret;
}

std::vector<int> Legalizer::cellLegalY() const {
  std::vector<int> ret(nbCells());
  for (int r = 0; r < nbRows(); ++r) {
    for (int i = 0; i < rowToCells_[r].size(); ++i) {
      int c = rowToCells_[r][i];
      ret[c] = rows_[r].minY;
    }
  }
  return ret;
}

std::vector<CellOrientation> Legalizer::cellLegalOrientation() const {
  std::vector<CellOrientation> ret(nbCells());
  for (int r = 0; r < nbRows(); ++r) {
    for (int c : rowToCells_[r]) {
      ret[c] = rows_[r].orientation;
    }
  }
  return ret;
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
    circuit.cellX_[i] = cellX[j];
    circuit.cellY_[i] = cellY[j];
    circuit.cellOrientation_[i] = cellOrient[j];
    ++j;
  }
}

std::vector<float> Legalizer::allDistances(LegalizationModel model) const {
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

float Legalizer::meanDistance(LegalizationModel model) const {
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

float Legalizer::rmsDistance(LegalizationModel model) const {
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

float Legalizer::maxDistance(LegalizationModel model) const {
  std::vector<float> dist = allDistances(model);
  return *std::max_element(dist.begin(), dist.end());
}

long long Legalizer::totalCellArea() const {
  long long ret = 0;
  for (int c = 0; c < nbCells(); ++c) {
    ret += static_cast<long long>(cellWidth_[c]) *
           static_cast<long long>(cellHeight_[c]);
  }
  return ret;
}

void Legalizer::runTetris(const std::vector<int> &cells) {}
}  // namespace coloquinte