
#include "place_detailed/legalizer.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>

namespace coloquinte {
Legalizer Legalizer::fromIspdCircuit(const Circuit &circuit) {
  // Represent fixed cells with -1 width so they are not considered
  std::vector<int> widths = circuit.cellWidth_;
  for (int i = 0; i < circuit.nbCells(); ++i) {
    if (circuit.cellFixed_[i]) {
      widths[i] = -1;
    }
  }
  return Legalizer(circuit.computeRows(), widths, circuit.cellX_, circuit.cellY_);
}

Legalizer::Legalizer(const std::vector<Rectangle> &rows,
                     const std::vector<int> &width,
                     const std::vector<int> &targetX,
                     const std::vector<int> &targetY)
    : costModel_(LegalizationModel::L1),
      cellWidth_(width),
      cellTargetX_(targetX),
      cellTargetY_(targetY) {
  assert(width.size() == targetX.size());
  assert(width.size() == targetY.size());
  // Sort the rows
  rows_ = rows;
  std::stable_sort(
      rows_.begin(), rows_.end(), [](Rectangle a, Rectangle b) -> bool {
        return a.minY < b.minY || (a.minY == b.minY && a.minX < b.minX);
      });
  for (Rectangle row : rows) {
    rowLegalizers_.emplace_back(row.minX, row.maxX);
  }
  rowToCells_.resize(rows.size());
  cellToX_ = cellTargetX_;
  cellToY_ = cellTargetY_;
  cellToRow_.assign(width.size(), -1);
}

void Legalizer::run() {
  std::vector<int> cellOrder = computeCellOrder(1.0, 0.5, 0.01);

  for (int c : cellOrder) {
    placeCellOptimally(c);
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

bool Legalizer::placeCellOptimally(int cell) {
  /**
   * Simple algorithm that tries close row first and stops early if no
   * improvement can be found
   */
  if (isIgnored(cell)) {
    return true;
  }
  int targetX = cellTargetX_[cell];
  int targetY = cellTargetY_[cell];
  int bestX = 0;
  int bestRow = -1;
  long long bestDist = std::numeric_limits<long long>::max();

  auto tryPlace = [&](int row) {
    int yDist = norm(0, rows_[row].minY - targetY, costModel_);
    if (bestRow != -1 && yDist > bestDist) {
      // Not possible to do better since the rows are sorted
      return true;
    }
    // Find the best position for the cell
    auto [ok, dist] = placeCellOptimally(cell, row);
    dist += yDist;
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
        "Unable to place a cell with the naive legalization algorithm");
  }
  rowLegalizers_[bestRow].push(cellWidth_[cell], targetX);
  rowToCells_[bestRow].push_back(cell);
  cellToRow_[cell] = bestRow;
  return true;
}

std::pair<bool, int> Legalizer::placeCellOptimally(int cell, int row) {
  if (rowLegalizers_[row].remainingSpace() < cellWidth_[cell]) {
    return std::make_pair(false, 0);
  }
  int dist = rowLegalizers_[row].getCost(cellWidth_[cell], cellTargetX_[cell]);
  return std::make_pair(true, dist);
}

std::vector<int> Legalizer::computeCellOrder(float weightX, float weightWidth,
                                             float weightY) const {
  // Sort the cells by target X coordinate
  std::vector<std::pair<float, int> > sortedCells;
  for (int i = 0; i < nbCells(); ++i) {
    float val = weightX * cellTargetX_[i] + weightWidth * cellWidth_[i] +
                weightY * cellTargetY_[i];
    sortedCells.emplace_back(val, i);
  }
  std::stable_sort(sortedCells.begin(), sortedCells.end());
  std::vector<int> cells;
  for (auto p : sortedCells) {
    cells.push_back(p.second);
  }
  return cells;
}

int Legalizer::closestRow(int y) const {
  auto it = std::lower_bound(rows_.begin(), rows_.end(), y,
                             [](Rectangle r, int v) { return r.minY < v; });
  if (it == rows_.end()) return nbRows() - 1;
  if (it == rows_.begin()) return 0;
  int row = it - rows_.begin();
  assert(row >= 1);
  if (rows_[row].minX - y > y - rows_[row - 1].minX) {
    return row - 1;
  } else {
    return row;
  }
}

std::vector<int> Legalizer::cellLegalX() const {
  std::vector<int> ret(nbCells());
  for (int r = 0; r < nbRows(); ++r) {
    for (int i = 0; i < rowToCells_[r].size(); ++i) {
      int c = rowToCells_[r][i];
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

void Legalizer::exportPlacement(Circuit &circuit) {
  std::vector<int> cellX = cellLegalX();
  std::vector<int> cellY = cellLegalY();
  for (int i = 0; i < circuit.nbCells(); ++i) {
    if (isIgnored(i)) continue;
    circuit.cellX_[i] = cellX[i];
    circuit.cellY_[i] = cellY[i];
  }
}

void Legalizer::report(bool verbose) const {
  std::cout << "Legalizer with " << nbCells() << " cells on " << nbRows()
            << " rows" << std::endl;
  std::cout << "Mean dist: " << meanDistance(costModel_) << std::endl;
  std::cout << "RMS dist " << rmsDistance(costModel_) << std::endl;
  std::cout << "Max dist " << maxDistance(costModel_) << std::endl;
  if (!verbose) return;
  for (int row = 0; row < nbRows(); ++row) {
    std::cout << "Row " << rows_[row].minY << ", " << rows_[row].minX << " to "
              << rows_[row].maxX << ": ";
    for (int c : rowToCells_[row]) {
      int x = cellToX_[c];
      int w = cellWidth_[c];
      std::cout << c << " (" << x << " - " << x + w << ", target "
                << cellTargetX_[c] << "," << cellTargetY_[c] << ") ";
    }
    std::cout << std::endl;
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
    if (isIgnored(i)) distances.push_back(0.0f);
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
    disp += cellWidth_[i] * dist[i];
  }
  return disp / totalCellWidth();
}

float Legalizer::rmsDistance(LegalizationModel model) const {
  std::vector<float> dist = allDistances(model);
  float disp = 0.0f;
  for (int i = 0; i < nbCells(); ++i) {
    disp += cellWidth_[i] * dist[i] * dist[i];
  }
  return std::sqrt(disp / totalCellWidth());
}

float Legalizer::maxDistance(LegalizationModel model) const {
  std::vector<float> dist = allDistances(model);
  return *std::max_element(dist.begin(), dist.end());
}

int Legalizer::totalCellWidth() const {
  int ret = 0;
  for (int c = 0; c < nbCells(); ++c) {
    if (isIgnored(c)) continue;
    ret += cellWidth_[c];
  }
  return ret;
}
}  // namespace coloquinte