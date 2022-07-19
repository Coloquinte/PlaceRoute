
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
  for (int i = 0; i < nbBinsX(); ++i) {
    assert(binLimitX_[i] <= binLimitX_[i + 1]);
  }
  for (int i = 0; i < nbBinsY(); ++i) {
    assert(binLimitY_[i] <= binLimitY_[i + 1]);
  }
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

long long DensityGrid::binCapacity(BinGroup g) const {
  long long ret = 0;
  for (int i = g.minXCoord; i < g.maxXCoord; ++i) {
    for (int j = g.minYCoord; j < g.maxYCoord; ++j) {
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

long long DensityPlacement::totalDemand() const {
  long long ret = 0;
  for (int demand : cellDemand_) {
    ret += demand;
  }
  return ret;
}

long long DensityPlacement::totalOverflow() const {
  long long ret = 0;
  for (int i = 0; i < nbBinsX(); ++i) {
    for (int j = 0; j < nbBinsY(); ++j) {
      ret += std::max(binUsage(i, j) - binCapacity(i, j), 0LL);
    }
  }
  return ret;
}

float DensityPlacement::overflowRatio() const {
  float demand = totalDemand();
  if (demand > 0) {
    return totalOverflow() / demand;
  } else {
    return 0.0;
  }
}

long long DensityPlacement::binUsage(int x, int y) const {
  long long usage = 0;
  for (int c : binCells(x, y)) {
    usage += cellDemand(c);
  }
  return usage;
}

std::vector<float> DensityPlacement::simpleCoordX() const {
  std::vector<float> ret(nbCells(), 0.0f);
  for (int i = 0; i < nbBinsX(); ++i) {
    for (int j = 0; j < nbBinsY(); ++j) {
      for (int c : binCells_[i][j]) {
        ret[c] = binX(i);
      }
    }
  }
  return ret;
}

std::vector<float> DensityPlacement::simpleCoordY() const {
  std::vector<float> ret(nbCells(), 0.0f);
  for (int i = 0; i < nbBinsY(); ++i) {
    for (int j = 0; j < nbBinsX(); ++j) {
      for (int c : binCells_[i][j]) {
        ret[c] = binY(i);
      }
    }
  }
  return ret;
}

void DensityPlacement::check() const { DensityGrid::check(); }

HierarchicalDensityState::HierarchicalDensityState(DensityGrid grid,
                                                   std::vector<int> cellDemand)
    : DensityGrid(grid), cellDemand_(cellDemand) {
  xLimits_.push_back(0);
  xLimits_.push_back(nbBinsX());
  yLimits_.push_back(0);
  yLimits_.push_back(nbBinsY());
  std::vector<int> allCells;
  for (int c = 0; c < nbCells(); ++c) {
    allCells.push_back(c);
  }
  binCells_.emplace_back();
  binCells_.back().push_back(allCells);
  check();
}

HierarchicalDensityState::HierarchicalDensityState(DensityPlacement placement)
    : DensityGrid(placement) {
  for (int i = 0; i <= nbBinsX(); ++i) {
    xLimits_.push_back(i);
  }
  for (int i = 0; i <= nbBinsY(); ++i) {
    yLimits_.push_back(i);
  }
  binCells_ = placement.binCells_;
}

long long HierarchicalDensityState::hBinUsage(int x, int y) const {
  long long usage = 0;
  for (int c : binCells(x, y)) {
    usage += cellDemand(c);
  }
  return usage;
}

namespace {
/**
 * @brief Update the limits between bins after a split is done, and return the
 * number of bins split
 */
std::vector<int> hierarchicalSplitHelper(std::vector<int> &limits) {
  int maxSz = 1;
  for (int i = 0; i + 1 < limits.size(); ++i) {
    int sz = limits[i + 1] - limits[i];
    maxSz = std::max(maxSz, sz);
  }
  // Minimum number of bins to split: only split the larger ones
  int minSplitSize = std::max(2, maxSz / 2 + 1);
  std::vector<int> ret;
  std::vector<int> newLimits;
  newLimits.push_back(limits.front());
  for (int i = 0; i + 1 < limits.size(); ++i) {
    int sz = limits[i + 1] - limits[i];
    if (sz < minSplitSize) {
      ret.push_back(1);
      newLimits.push_back(limits[i+1]);
    }
    else {
      ret.push_back(2);
      newLimits.push_back((limits[i] + limits[i+1]) / 2);
      newLimits.push_back(limits[i+1]);
    }
  }
  limits = newLimits;
  return ret;
}

std::vector<int> numberOfSplitToAssociation(const std::vector<int> &nb) {
  std::vector<int> ret;
  ret.push_back(0);
  for (int n : nb) {
    ret.push_back(ret.back() + n);
  }
  return ret;
}
}  // namespace

std::vector<int> HierarchicalDensityState::splitX() {
  std::vector<int> nbSplit = hierarchicalSplitHelper(xLimits_);
  std::vector<std::vector<std::vector<int> > > newBinCells;
  for (int i = 0; i < binCells_.size(); ++i) {
    newBinCells.push_back(binCells_[i]);
    if (nbSplit[i] > 1) {
      newBinCells.emplace_back();
    }
  }
  binCells_ = newBinCells;
  check();
  return numberOfSplitToAssociation(nbSplit);
}

std::vector<int> HierarchicalDensityState::splitY() {
  std::vector<int> nbSplit = hierarchicalSplitHelper(yLimits_);
  std::vector<std::vector<std::vector<int> > > newBinCells(binCells_.size());
  for (int i = 0; i < binCells_.size(); ++i) {
    for (int j = 0; j < binCells_[i].size(); ++j) {
      newBinCells[i].push_back(binCells_[i][j]);
      if (nbSplit[j] > 1) {
        newBinCells[i].emplace_back();
      }
    }
  }
  binCells_ = newBinCells;
  check();
  return numberOfSplitToAssociation(nbSplit);
}

void HierarchicalDensityState::check() const {
  for (int i = 0; i < hNbBinsX(); ++i) {
    assert(xLimits_[i] < xLimits_[i + 1]);
  }
  for (int i = 0; i < hNbBinsY(); ++i) {
    assert(yLimits_[i] < yLimits_[i + 1]);
  }
  assert(binCells_.size() == nbBinxX());
  for (auto &bc : binCells_) {
    assert(bc.size() == nbBinxY());
  }
  std::vector<int> nbPlaced(nbCells(), 0);
  for (auto &bc : binCells_) {
    for (auto &bcc : bc) {
      for (int c : bcc) {
        assert(c >= 0);
        assert(c < nbCells());
        nbPlaced[c]++;
      }
    }
  }
  for (int nb : nbPlaced) {
    assert(nb == 1);
  }
}