
#include "place_global/density_grid.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

DensityGrid::DensityGrid(Rectangle area)
    : DensityGrid(std::vector<Rectangle>({area})) {}

DensityGrid::DensityGrid(std::vector<Rectangle> regions,
                         std::vector<Rectangle> obstacles) {
  regions_ = regions;
  obstacles_ = obstacles;
  placementArea_ = computePlacementArea();
  updateBinsToNumber(1, 1);
  check();
}

DensityGrid DensityGrid::fromIspdCircuit(const Circuit &circuit,
                                         float sizeFactor) {
  DensityGrid ret(circuit.rows);
  int minCellHeight = std::numeric_limits<int>::max();
  for (int i = 0; i < circuit.nbCells(); ++i) {
    int height = circuit.cellHeights[i];
    if (height > 0) {
      minCellHeight = std::min(height, minCellHeight);
    }
  }
  ret.updateBinsToSize(sizeFactor * minCellHeight);
  return ret;
}

void DensityGrid::updateBinsToNumber(int binsX, int binsY) {
  binLimitX_.clear();
  binLimitY_.clear();
  binX_.clear();
  binY_.clear();
  for (int i = 0; i < binsX + 1; ++i) {
    binLimitX_.push_back(
        placementArea_.minX +
        (i * (placementArea_.maxX - placementArea_.minX) / binsX));
  }
  assert(binLimitX_.size() == binsX + 1);
  for (int i = 0; i < binsY + 1; ++i) {
    binLimitY_.push_back(
        placementArea_.minY +
        (i * (placementArea_.maxY - placementArea_.minY) / binsY));
  }
  assert(binLimitY_.size() == binsY + 1);
  for (int i = 0; i < binsX; ++i) {
    binX_.push_back(0.5f * (binLimitX_[i] + binLimitX_[i + 1]));
  }
  for (int i = 0; i < binsY; ++i) {
    binY_.push_back(0.5f * (binLimitY_[i] + binLimitY_[i + 1]));
  }
  // TODO: actually compute bin capacity from rows
  binCapacity_.assign(binsX, std::vector<long long>(binsY, 0));
  for (int i = 0; i < binsX; ++i) {
    for (int j = 0; j < binsY; ++j) {
      long long w = binLimitX_[i + 1] - binLimitX_[i];
      long long h = binLimitY_[j + 1] - binLimitY_[j];
      assert(w >= 0);
      assert(h >= 0);
      binCapacity_[i][j] = w * h;
    }
  }
}

void DensityGrid::updateBinsToSize(int maxXSize, int maxYSize) {
  int binsX = placementArea_.width() / maxXSize;
  int binsY = placementArea_.height() / maxYSize;
  updateBinsToNumber(binsX, binsY);
}

void DensityGrid::check() const {
  assert(binCapacity_.size() == nbBinsX());
  for (const auto &bc : binCapacity_) {
    assert(bc.size() == nbBinsY());
  }
  assert(binX_.size() == nbBinsX());
  assert(binY_.size() == nbBinsY());
  assert(binLimitX_.size() == nbBinsX() + 1);
  assert(binLimitY_.size() == nbBinsY() + 1);
}

long long DensityGrid::totalCapacity() const {
  long long ret = 0;
  for (int i = 0; i < nbBinsX(); ++i) {
    for (int j = 0; j < nbBinsY(); ++j) {
      ret += binCapacity_[i][j];
    }
  }
  return ret;
}

float DensityGrid::groupCenterX(DensityGrid::BinGroup g) const {
  float capa = 0.0;
  float coord = 0.0;
  for (int i = g.minXCoord; i < g.maxXCoord; ++i) {
    for (int j = g.minYCoord; j < g.maxYCoord; ++j) {
      capa += binCapacity(i, j);
      coord += binX(i) * binCapacity(i, j);
    }
  }
  if (capa > 0) {
    return coord / capa;
  } else {
    return binX((g.minXCoord + g.maxXCoord) / 2);
  }
}

float DensityGrid::groupCenterY(DensityGrid::BinGroup g) const {
  float capa = 0.0;
  float coord = 0.0;
  for (int i = g.minXCoord; i < g.maxXCoord; ++i) {
    for (int j = g.minYCoord; j < g.maxYCoord; ++j) {
      capa += binCapacity(i, j);
      coord += binY(i) * binCapacity(i, j);
    }
  }
  if (capa > 0) {
    return coord / capa;
  } else {
    return binY((g.minYCoord + g.maxYCoord) / 2);
  }
}

Rectangle DensityGrid::computePlacementArea() const {
  int minX = std::numeric_limits<int>::max();
  int maxX = std::numeric_limits<int>::min();
  int minY = std::numeric_limits<int>::max();
  int maxY = std::numeric_limits<int>::min();
  if (regions_.empty()) {
    return Rectangle(0, 0, 0, 0);
  }
  for (Rectangle row : regions_) {
    minX = std::min(row.minX, minX);
    maxX = std::max(row.maxX, maxX);
    minY = std::min(row.minY, minY);
    maxY = std::max(row.maxY, maxY);
  }
  return Rectangle(minX, maxX, minY, maxY);
}
