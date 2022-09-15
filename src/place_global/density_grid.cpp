
#include "place_global/density_grid.hpp"

#include <algorithm>
#include <boost/polygon/polygon.hpp>
#include <cmath>
#include <numeric>

namespace bpl = boost::polygon;

namespace coloquinte {
DensityGrid::DensityGrid(int binSize, Rectangle area)
    : DensityGrid(binSize, std::vector<Rectangle>({area})) {}

DensityGrid::DensityGrid(int binSize, std::vector<Rectangle> regions,
                         std::vector<Rectangle> obstacles) {
  placementArea_ = computePlacementArea(regions);
  updateBinsToSize(binSize);
  std::vector<Rectangle> actualRegions =
      computeActualRegions(regions, obstacles);
  updateBinCapacity(actualRegions);
  check();
}

DensityGrid DensityGrid::fromIspdCircuit(const Circuit &circuit,
                                         float sizeFactor) {
  int minCellHeight = std::numeric_limits<int>::max();
  for (int i = 0; i < circuit.nbCells(); ++i) {
    int height = circuit.cellHeight_[i];
    if (height > 0) {
      minCellHeight = std::min(height, minCellHeight);
    }
  }
  std::vector<Rectangle> obstacles;
  for (int i = 0; i < circuit.nbCells(); ++i) {
    if (!circuit.fixed(i)) continue;
    int x = circuit.cellX_[i];
    int y = circuit.cellY_[i];
    obstacles.emplace_back(x, x + circuit.cellWidth_[i], y,
                           y + circuit.cellHeight_[i]);
  }
  return DensityGrid(sizeFactor * minCellHeight, circuit.rows, obstacles);
}

DensityGrid::DensityGrid(std::vector<int> xLimits, std::vector<int> yLimits,
                         std::vector<std::vector<long long> > binCapacity) {
  binLimitX_ = xLimits;
  binLimitY_ = yLimits;
  binCapacity_ = binCapacity;
  updateBinCenters();
}

std::vector<Rectangle> DensityGrid::computeActualRegions(
    const std::vector<Rectangle> &regions,
    const std::vector<Rectangle> &obstacles) {
  bpl::polygon_90_set_data<int> region_set;
  for (Rectangle r : regions) {
    region_set.insert(bpl::rectangle_data<int>(r.minX, r.minY, r.maxX, r.maxY));
  }
  for (Rectangle r : obstacles) {
    region_set.insert(bpl::rectangle_data<int>(r.minX, r.minY, r.maxX, r.maxY),
                      true);
  }
  std::vector<bpl::rectangle_data<int> > diff;
  bpl::get_rectangles(diff, region_set);
  std::vector<Rectangle> ret;
  for (auto r : diff) {
    ret.emplace_back(bpl::xl(r), bpl::xh(r), bpl::yl(r), bpl::yh(r));
  }
  return ret;
}

void DensityGrid::updateBinCenters() {
  int binsX = binLimitX_.size() - 1;
  for (int i = 0; i < binsX; ++i) {
    binX_.push_back(0.5f * (binLimitX_[i] + binLimitX_[i + 1]));
  }
  int binsY = binLimitY_.size() - 1;
  for (int i = 0; i < binsY; ++i) {
    binY_.push_back(0.5f * (binLimitY_[i] + binLimitY_[i + 1]));
  }
}

void DensityGrid::updateBinCapacity() {
  int binsX = binLimitX_.size() - 1;
  int binsY = binLimitY_.size() - 1;
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

void DensityGrid::updateBinCapacity(const std::vector<Rectangle> &regions) {
  for (int i = 0; i < nbBinsX(); ++i) {
    for (int j = 0; j < nbBinsY(); ++j) {
      binCapacity_[i][j] = 0;
    }
  }
  for (Rectangle reg : regions) {
    // TODO: this should do a binary search first to limit the number of bins
    // that are accessed
    for (int i = 0; i < nbBinsX(); ++i) {
      for (int j = 0; j < nbBinsY(); ++j) {
        Rectangle binReg = region(i, j);
        if (reg.intersects(binReg)) {
          binCapacity_[i][j] += Rectangle::intersection(reg, binReg).area();
        }
      }
    }
  }
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
  updateBinCenters();
  updateBinCapacity();
}

void DensityGrid::updateBinsToSize(int maxSize) {
  int binsX = std::max(1, placementArea_.width() / maxSize);
  int binsY = std::max(1, placementArea_.height() / maxSize);
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
      coord += binY(j) * binCapacity(i, j);
    }
  }
  if (capa > 0) {
    return coord / capa;
  } else {
    return binY((g.minYCoord + g.maxYCoord) / 2);
  }
}

Rectangle DensityGrid::computePlacementArea(
    const std::vector<Rectangle> &regions) {
  int minX = std::numeric_limits<int>::max();
  int maxX = std::numeric_limits<int>::min();
  int minY = std::numeric_limits<int>::max();
  int maxY = std::numeric_limits<int>::min();
  if (regions.empty()) {
    return Rectangle(0, 0, 0, 0);
  }
  for (Rectangle row : regions) {
    minX = std::min(row.minX, minX);
    maxX = std::max(row.maxX, maxX);
    minY = std::min(row.minY, minY);
    maxY = std::max(row.maxY, maxY);
  }
  return Rectangle(minX, maxX, minY, maxY);
}

HierarchicalDensityPlacement HierarchicalDensityPlacement::fromIspdCircuit(
    const Circuit &circuit, float sizeFactor) {
  DensityGrid grid = DensityGrid::fromIspdCircuit(circuit, sizeFactor);
  std::vector<int> demands;
  for (int i = 0; i < circuit.nbCells(); ++i) {
    if (circuit.fixed(i)) {
      demands.push_back(0LL);
    } else {
      demands.push_back(circuit.area(i));
    }
  }
  return HierarchicalDensityPlacement(grid, demands);
}

int HierarchicalDensityPlacement::nbNonEmptyCells() const {
  int ret = 0;
  for (int i = 0; i < nbCells(); ++i) {
    if (cellDemand(i) > 0) ret++;
  }
  return ret;
}

long long HierarchicalDensityPlacement::totalDemand() const {
  long long ret = 0;
  for (int demand : cellDemand_) {
    ret += demand;
  }
  return ret;
}

long long HierarchicalDensityPlacement::totalOverflow() const {
  long long ret = 0;
  for (int i = 0; i < nbBinsX(); ++i) {
    for (int j = 0; j < nbBinsY(); ++j) {
      ret += std::max(binUsage(i, j) - binCapacity(i, j), 0LL);
    }
  }
  return ret;
}

float HierarchicalDensityPlacement::overflowRatio() const {
  float demand = totalDemand();
  if (demand > 0) {
    return totalOverflow() / demand;
  } else {
    return 0.0;
  }
}

std::vector<float> HierarchicalDensityPlacement::simpleCoordX() const {
  std::vector<float> ret(nbCells(), 0.0f);
  for (int i = 0; i < nbBinsX(); ++i) {
    for (int j = 0; j < nbBinsY(); ++j) {
      float coord = binX(i, j);
      for (int c : binCells(i, j)) {
        ret[c] = coord;
      }
    }
  }
  return ret;
}

std::vector<float> HierarchicalDensityPlacement::simpleCoordY() const {
  std::vector<float> ret(nbCells(), 0.0f);
  for (int i = 0; i < nbBinsX(); ++i) {
    for (int j = 0; j < nbBinsY(); ++j) {
      float coord = binY(i, j);
      for (int c : binCells(i, j)) {
        ret[c] = coord;
      }
    }
  }
  return ret;
}

namespace {
std::vector<float> spreadCells(const std::vector<float> &targets,
                               const std::vector<float> &demands,
                               float minCoord, float maxCoord) {
  assert(targets.size() == demands.size());
  std::vector<std::pair<float, int> > order;
  for (int i = 0; i < targets.size(); ++i) {
    order.emplace_back(targets[i], i);
  }
  std::sort(order.begin(), order.end());
  float invTotalDemand =
      1.0f / std::accumulate(demands.begin(), demands.end(), 0.0f);
  float dem = 0.0f;
  std::vector<float> coords(order.size(), 0.0f);
  for (int i = 0; i < order.size(); ++i) {
    int c = order[i].second;
    dem += 0.5f * demands[c] * invTotalDemand;
    coords[c] = dem * maxCoord + (1.0f - dem) * minCoord;
    dem += 0.5f * demands[c] * invTotalDemand;
  }
  return coords;
}
}  // namespace

std::vector<float> HierarchicalDensityPlacement::spreadCoordX(
    const std::vector<float> &target) const {
  std::vector<float> ret(nbCells(), 0.0f);
  for (int i = 0; i < nbBinsX(); ++i) {
    for (int j = 0; j < nbBinsY(); ++j) {
      std::vector<float> binTargets;
      std::vector<float> binDemands;
      for (int c : binCells(i, j)) {
        binTargets.push_back(target[c]);
        binDemands.push_back(cellDemand(c));
      }
      std::vector<float> coords =
          spreadCells(binTargets, binDemands, binLimitX(i), binLimitX(i + 1));
      for (int k = 0; k < coords.size(); ++k) {
        ret[binCells(i, j)[k]] = coords[k];
      }
    }
  }
  return ret;
}

std::vector<float> HierarchicalDensityPlacement::spreadCoordY(
    const std::vector<float> &target) const {
  std::vector<float> ret(nbCells(), 0.0f);
  for (int i = 0; i < nbBinsX(); ++i) {
    for (int j = 0; j < nbBinsY(); ++j) {
      std::vector<float> binTargets;
      std::vector<float> binDemands;
      for (int c : binCells(i, j)) {
        binTargets.push_back(target[c]);
        binDemands.push_back(cellDemand(c));
      }
      std::vector<float> coords =
          spreadCells(binTargets, binDemands, binLimitY(j), binLimitY(j + 1));
      for (int k = 0; k < coords.size(); ++k) {
        ret[binCells(i, j)[k]] = coords[k];
      }
    }
  }
  return ret;
}

HierarchicalDensityPlacement::HierarchicalDensityPlacement(DensityGrid grid,
                                                           int nbCells)
    : HierarchicalDensityPlacement(grid, std::vector<int>(nbCells, 0)) {}

HierarchicalDensityPlacement::HierarchicalDensityPlacement(
    DensityGrid grid, std::vector<int> cellDemand)
    : grid_(grid), cellDemand_(cellDemand) {
  setupHierarchy();
  levelX_ = nbLevelX() - 1;
  levelY_ = nbLevelY() - 1;
  std::vector<int> allCells;
  for (int c = 0; c < nbCells(); ++c) {
    allCells.push_back(c);
  }
  binCells_.emplace_back();
  binCells_.back().push_back(allCells);
  updateCellToBin();
  check();
}

int HierarchicalDensityPlacement::findBinByX(int coord) const {
  int mn = 0;
  int mx = nbBinsX();
  while (mx > mn + 1) {
    int mid = (mx + mn) / 2;
    if (binLimitX(mid) > coord) {
      mx = mid;
    } else {
      mn = mid;
    }
  }
  assert(binLimitX(mn) <= coord || mn == 0);
  assert(binLimitX(mn + 1) > coord || mn == nbBinsX() - 1);
  return mn;
}

int HierarchicalDensityPlacement::findBinByY(int coord) const {
  int mn = 0;
  int mx = nbBinsY();
  while (mx > mn + 1) {
    int mid = (mx + mn) / 2;
    if (binLimitY(mid) > coord) {
      mx = mid;
    } else {
      mn = mid;
    }
  }
  assert(binLimitY(mn) <= coord || mn == 0);
  assert(binLimitY(mn + 1) > coord || mn == nbBinsY() - 1);
  return mn;
}

long long HierarchicalDensityPlacement::binUsage(int x, int y) const {
  long long usage = 0;
  for (int c : binCells(x, y)) {
    usage += cellDemand(c);
  }
  return usage;
}

void HierarchicalDensityPlacement::setBinCells(int x, int y,
                                               std::vector<int> cells) {
  assert(x < nbBinsX());
  assert(y < nbBinsY());
  for (int c : cells) {
    cellBinX_[c] = x;
    cellBinY_[c] = y;
  }
  binCells_[x][y] = cells;
}

namespace {
bool canRefine(const std::vector<int> &limits) {
  for (int i = 0; i + 1 < limits.size(); ++i) {
    if (limits[i + 1] - limits[i] > 1) return true;
  }
  return false;
}
void refine(const std::vector<int> &oldLimits, std::vector<int> &limits,
            std::vector<int> &parents) {
  limits.push_back(0);
  int minSplitSize = 2;
  for (int i = 0; i + 1 < oldLimits.size(); ++i) {
    int e = oldLimits[i + 1];
    int b = oldLimits[i];
    if (e - b >= minSplitSize) {
      limits.push_back((e + b) / 2);
      parents.push_back(i);
    }
    limits.push_back(e);
    parents.push_back(i);
  }
}
void setupHierarchyHelper(int nbBins, std::vector<std::vector<int> > &limits,
                          std::vector<std::vector<int> > &parents) {
  limits.clear();
  parents.clear();
  limits.push_back(std::vector<int>({0, nbBins}));
  parents.push_back(std::vector<int>({0}));
  while (canRefine(limits.back())) {
    std::vector<int> nextLimits, nextParents;
    refine(limits.back(), nextLimits, nextParents);
    limits.push_back(nextLimits);
    parents.push_back(nextParents);
  }
  std::reverse(limits.begin(), limits.end());
  std::reverse(parents.begin(), parents.end());
}
}  // namespace

void HierarchicalDensityPlacement::setupHierarchy() {
  setupHierarchyHelper(grid_.nbBinsX(), xLimits_, parentX_);
  setupHierarchyHelper(grid_.nbBinsY(), yLimits_, parentY_);
}

void HierarchicalDensityPlacement::updateCellToBin() {
  cellBinX_.assign(nbCells(), -1);
  cellBinY_.assign(nbCells(), -1);
  for (int i = 0; i < nbBinsX(); ++i) {
    for (int j = 0; j < nbBinsY(); ++j) {
      for (int c : binCells_[i][j]) {
        cellBinX_[c] = i;
        cellBinY_[c] = j;
      }
    }
  }
}

void HierarchicalDensityPlacement::coarsenX() {
  assert(levelX_ + 1 < nbLevelX());
  std::vector<std::vector<std::vector<int> > > newCells;
  newCells.assign(nbBinsX(levelX_ + 1),
                  std::vector<std::vector<int> >(nbBinsY()));
  for (int i = 0; i < nbBinsX(); ++i) {
    for (int j = 0; j < nbBinsY(); ++j) {
      int p = parentX(i);
      for (int c : binCells(i, j)) {
        newCells[p][j].push_back(c);
      }
    }
  }
  binCells_ = newCells;
  levelX_++;
  updateCellToBin();
  check();
}

void HierarchicalDensityPlacement::coarsenY() {
  assert(levelY_ + 1 < nbLevelY());
  std::vector<std::vector<std::vector<int> > > newCells;
  newCells.assign(nbBinsX(),
                  std::vector<std::vector<int> >(nbBinsY(levelY_ + 1)));
  for (int i = 0; i < nbBinsX(); ++i) {
    for (int j = 0; j < nbBinsY(); ++j) {
      int p = parentY(j);
      for (int c : binCells(i, j)) {
        newCells[i][p].push_back(c);
      }
    }
  }
  binCells_ = newCells;
  levelY_++;
  updateCellToBin();
  check();
}

void HierarchicalDensityPlacement::coarsenFully() {
  while (levelX() + 1 < nbLevelX()) {
    coarsenX();
  }
  while (levelY() + 1 < nbLevelY()) {
    coarsenY();
  }
}

void HierarchicalDensityPlacement::refineFully() {
  while (levelX() > 0) {
    refineX();
  }
  while (levelY() > 0) {
    refineY();
  }
}

void HierarchicalDensityPlacement::refineX() {
  assert(levelX_ >= 1);
  levelX_--;
  std::vector<std::vector<std::vector<int> > > newCells;
  newCells.assign(nbBinsX(), std::vector<std::vector<int> >(nbBinsY()));
  for (int i = 0; i < nbBinsX(); ++i) {
    for (int j = 0; j < nbBinsY(); ++j) {
      // Only assign the cells to the first child bin
      if (i != 0 && parentX(i) == parentX(i - 1)) {
        continue;
      }
      int p = parentX(i);
      for (int c : binCells_[p][j]) {
        newCells[i][j].push_back(c);
      }
    }
  }
  binCells_ = newCells;
  updateCellToBin();
  check();
}

void HierarchicalDensityPlacement::refineY() {
  assert(levelY_ >= 1);
  levelY_--;
  std::vector<std::vector<std::vector<int> > > newCells;
  newCells.assign(nbBinsX(), std::vector<std::vector<int> >(nbBinsY()));
  for (int i = 0; i < nbBinsX(); ++i) {
    for (int j = 0; j < nbBinsY(); ++j) {
      // Only assign the cells to the first child bin
      if (j != 0 && parentY(j) == parentY(j - 1)) {
        continue;
      }
      int p = parentY(j);
      for (int c : binCells_[i][p]) {
        newCells[i][j].push_back(c);
      }
    }
  }
  binCells_ = newCells;
  updateCellToBin();
  check();
}

void HierarchicalDensityPlacement::check() const {
  // Limits of the bins consistent
  for (auto &l : xLimits_) {
    assert(l.size() >= 1);
    assert(l.front() == 0);
    assert(l.back() == grid_.nbBinsX());
    for (int i = 0; i + 1 < l.size(); ++i) {
      assert(l[i] < l[i + 1]);
    }
  }
  for (auto &l : yLimits_) {
    assert(l.size() >= 1);
    assert(l.front() == 0);
    assert(l.back() == grid_.nbBinsY());
    for (int i = 0; i + 1 < l.size(); ++i) {
      assert(l[i] < l[i + 1]);
    }
  }
  // Number of levels consistent
  assert(xLimits_.size() == parentX_.size());
  for (int i = 0; i < nbLevelX(); ++i) {
    assert(xLimits_[i].size() == parentX_[i].size() + 1);
  }
  assert(yLimits_.size() == parentY_.size());
  for (int i = 0; i < nbLevelY(); ++i) {
    assert(yLimits_[i].size() == parentY_[i].size() + 1);
  }
  // Size of the allocation consistent
  assert(binCells_.size() == nbBinsX());
  for (auto &bc : binCells_) {
    assert(bc.size() == nbBinsY());
  }
  // All cells placed once and consistent
  std::vector<int> placeX(nbCells(), -1);
  std::vector<int> placeY(nbCells(), -1);
  for (int i = 0; i < nbBinsX(); ++i) {
    for (int j = 0; j < nbBinsY(); ++j) {
      for (int c : binCells_[i][j]) {
        assert(c >= 0);
        assert(c < nbCells());
        assert(placeX[c] == -1);
        assert(placeY[c] == -1);
        placeX[c] = i;
        placeY[c] = j;
      }
    }
  }
  for (int c = 0; c < nbCells(); ++c) {
    assert(placeX[c] != -1);
    assert(placeY[c] != -1);
    assert(cellBinX(c) == placeX[c]);
    assert(cellBinY(c) == placeY[c]);
  }
  // Two ways to compute the demand consistent
  long long usage = 0;
  for (int i = 0; i < nbBinsX(); ++i) {
    for (int j = 0; j < nbBinsY(); ++j) {
      usage += binUsage(i, j);
    }
  }
  assert(usage == totalDemand());
  // Two ways to compute the capacity consistent
  long long capacity = 0;
  for (int i = 0; i < nbBinsX(); ++i) {
    for (int j = 0; j < nbBinsY(); ++j) {
      capacity += binCapacity(i, j);
    }
  }
  assert(capacity == totalCapacity());
}
}