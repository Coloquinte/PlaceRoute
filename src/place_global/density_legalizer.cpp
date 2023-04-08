
#include "place_global/density_legalizer.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <numeric>
#include <utility>

#include "place_global/transportation.hpp"
#include "place_global/transportation_1d.hpp"
#include "utils/norm.hpp"

namespace coloquinte {
DensityLegalizer::Parameters::Parameters() {
  costModel = LegalizationModel::L1;
  nbSteps = 1;
  reoptimizationLength = 2;
  reoptimizationSquareSize = 1;
  quadraticPenaltyFactor = 0.0;
  coarseningLimit = 1.0;
  unidimensionalTransport = false;
}

DensityLegalizer::DensityLegalizer(DensityGrid grid,
                                   std::vector<int> cellDemand,
                                   const Parameters &params)
    : HierarchicalDensityPlacement(std::move(grid), std::move(cellDemand)) {
  cellTargetX_.assign(nbCells(), 0.0f);
  cellTargetY_.assign(nbCells(), 0.0f);
}

DensityLegalizer::DensityLegalizer(HierarchicalDensityPlacement pl,
                                   const Parameters &params)
    : HierarchicalDensityPlacement(std::move(pl)) {
  cellTargetX_.assign(nbCells(), 0.0f);
  cellTargetY_.assign(nbCells(), 0.0f);
}

DensityLegalizer DensityLegalizer::fromIspdCircuit(const Circuit &circuit,
                                                   float sizeFactor,
                                                   float sideMargin) {
  return DensityLegalizer(HierarchicalDensityPlacement::fromIspdCircuit(
      circuit, sizeFactor, sideMargin));
}

std::vector<float> DensityLegalizer::allDistances() const {
  std::vector<float> cellX = simpleCoordX();
  std::vector<float> cellY = simpleCoordY();
  std::vector<float> distances;
  distances.reserve(nbCells());
  for (int i = 0; i < nbCells(); ++i) {
    float dx = cellTargetX(i) - cellX[i];
    float dy = cellTargetY(i) - cellY[i];
    distances.push_back(distance(dx, dy));
  }
  return distances;
}

float DensityLegalizer::meanDistance() const {
  std::vector<float> dist = allDistances();
  float disp = 0.0f;
  for (int i = 0; i < nbCells(); ++i) {
    disp += cellDemand(i) * dist[i];
  }
  disp /= totalDemand();
  if (params_.costModel == LegalizationModel::L1Squared ||
      params_.costModel == LegalizationModel::L2Squared ||
      params_.costModel == LegalizationModel::LInfSquared) {
    disp = std::sqrt(disp);
  }
  return disp;
}

float DensityLegalizer::rmsDistance() const {
  std::vector<float> dist = allDistances();
  float disp = 0.0f;
  for (int i = 0; i < nbCells(); ++i) {
    disp += cellDemand(i) * dist[i] * dist[i];
  }
  return std::sqrt(disp / totalDemand());
}

float DensityLegalizer::maxDistance() const {
  std::vector<float> dist = allDistances();
  return *std::max_element(dist.begin(), dist.end());
}

inline float DensityLegalizer::distance(float x, float y) const {
  float d = norm(x, y, params_.costModel);
  float val = d * (1.0f + (float)params_.quadraticPenaltyFactor * d);
  return val;
}

void DensityLegalizer::check() const { HierarchicalDensityPlacement::check(); }

void DensityLegalizer::report(bool verbose) const {
  std::cout << "Total demand " << totalDemand() << std::endl;
  std::cout << "Total capacity " << totalCapacity() << std::endl;
  Rectangle a = grid_.placementArea();
  std::cout << "Area (" << a.minX << ", " << a.maxX << ") x (" << a.minY << ", "
            << a.maxY << ")" << std::endl;
  std::cout << "Bins " << nbBinsX() << " x " << nbBinsY() << std::endl;
  std::cout << "Overflow: " << overflowRatio() << std::endl;
  std::cout << "Mean dist: " << meanDistance() << std::endl;
  std::cout << "RMS dist " << rmsDistance() << std::endl;
  std::cout << "Max dist " << maxDistance() << std::endl;
  if (!verbose) {
    return;
  }
  std::cout << std::endl;
  for (int j = 0; j < nbBinsY(); ++j) {
    for (int i = 0; i < nbBinsX(); ++i) {
      std::cout << "\t" << binUsage(i, j) << "/" << binCapacity(i, j);
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

std::vector<std::pair<float, int> > DensityLegalizer::computeCellCosts(
    float cx1, float cy1, float cx2, float cy2,
    const std::vector<int> &cells) const {
  // TODO: always add a secondary objective using squared distance to improve
  // stability
  std::vector<std::pair<float, int> > cellCosts;
  for (int c : cells) {
    float x = cellTargetX_[c];
    float y = cellTargetY_[c];
    float cost = distance(x - cx1, y - cy1) - distance(x - cx2, y - cy2);
    cellCosts.emplace_back(cost, c);
  }
  std::sort(cellCosts.begin(), cellCosts.end(),
            [](std::pair<float, int> a, std::pair<float, int> b) {
              return a.first < b.first;
            });
  return cellCosts;
}

std::pair<std::vector<int>, std::vector<int> > DensityLegalizer::doSplit(
    const std::vector<std::pair<float, int> > &cellCosts, int ind) const {
  std::pair<std::vector<int>, std::vector<int> > ret;
  for (int i = 0; i < ind; ++i) {
    ret.first.push_back(cellCosts[i].second);
  }
  for (int i = ind; i < cellCosts.size(); ++i) {
    ret.second.push_back(cellCosts[i].second);
  }
  return ret;
}

int DensityLegalizer::findIdealSplitPos(
    const std::vector<std::pair<float, int> > &cellCosts) const {
  // Find the ideal split position
  int splitInd = 0;
  for (; splitInd < cellCosts.size(); ++splitInd) {
    if (cellCosts[splitInd].first > 0.0) {
      break;
    }
  }
  return splitInd;
}

int DensityLegalizer::findConstrainedSplitPos(
    const std::vector<std::pair<float, int> > &cellCosts, int targetPos,
    long long capa1, long long capa2) const {
  long long demand1 = 0;
  long long demand2 = 0;
  for (int i = 0; i < targetPos; ++i) {
    demand1 += cellDemand_[cellCosts[i].second];
  }
  for (int i = targetPos; i < cellCosts.size(); ++i) {
    demand2 += cellDemand_[cellCosts[i].second];
  }
  int splitPos = targetPos;
  // Remove from the left if overflowed
  while (splitPos > 0 && demand1 - capa1 > 0 && capa2 > 0) {
    int dem = cellDemand_[cellCosts[splitPos - 1].second];
    if (capa1 > 0 && demand1 - capa1 < demand2 - capa2 + dem) {
      // Stop if overflow would be bigger on the other side
      break;
    }
    demand1 -= dem;
    demand2 += dem;
    --splitPos;
  }
  // Remove from the right if overflowed
  while (splitPos < cellCosts.size() && demand2 - capa2 > 0 && capa1 > 0) {
    int dem = cellDemand_[cellCosts[splitPos].second];
    if (capa2 > 0 && demand2 - capa2 < demand1 - capa1 + dem) {
      // Stop if overflow would be bigger on the other side
      break;
    }
    demand2 -= dem;
    demand1 += dem;
    ++splitPos;
  }
  return splitPos;
}

void DensityLegalizer::rebisect(int x1, int y1, int x2, int y2) {
  if (x1 == x2 && y1 == y2) {
    return;
  }

  // Get all cells
  std::vector<int> cells;
  cells.insert(cells.end(), binCells_[x1][y1].begin(), binCells_[x1][y1].end());
  cells.insert(cells.end(), binCells_[x2][y2].begin(), binCells_[x2][y2].end());

  std::vector<std::pair<float, int> > cellCosts = computeCellCosts(
      binX(x1, y1), binY(x1, y1), binX(x2, y2), binY(x2, y2), cells);
  int idealSplitPos = findIdealSplitPos(cellCosts);
  int splitPos = findConstrainedSplitPos(
      cellCosts, idealSplitPos, binCapacity(x1, y1), binCapacity(x2, y2));
  auto b = doSplit(cellCosts, splitPos);
  setBinCells(x1, y1, b.first);
  setBinCells(x2, y2, b.second);
}

void DensityLegalizer::reoptimize(
    const std::vector<std::pair<int, int> > &binCandidates) {
  if (binCandidates.size() == 2) {
    auto [x1, y1] = binCandidates[0];
    auto [x2, y2] = binCandidates[1];
    rebisect(x1, y1, x2, y2);
    return;
  }

  std::vector<int> cells;
  std::vector<int> assignment;
  std::vector<std::pair<int, int> > bins;
  for (auto [x, y] : binCandidates) {
    if (binCapacity(x, y) > 0) {
      bins.emplace_back(x, y);
    }
    // Allocate to the last non-empty bin
    int binCnt = bins.empty() ? 0 : bins.size() - 1;
    for (int c : binCells_[x][y]) {
      cells.push_back(c);
      assignment.push_back(binCnt);
    }
  }
  if (bins.empty()) {
    return;
  }

  for (auto [x, y] : binCandidates) {
    binCells_[x][y].clear();
  }

  if (bins.size() == 1) {
    auto [x, y] = bins[0];
    setBinCells(x, y, cells);
    return;
  }

  // Build a transportation problem
  std::vector<long long> demands;
  for (int c : cells) {
    demands.push_back(cellDemand(c));
  }

  std::vector<long long> capacities;
  for (auto [x, y] : bins) {
    capacities.push_back(binCapacity(x, y));
  }

  std::vector<std::vector<float> > costs;
  for (auto [x, y] : bins) {
    float bx = binX(x, y);
    float by = binY(x, y);
    std::vector<float> binCosts;
    for (int c : cells) {
      float cx = cellTargetX(c);
      float cy = cellTargetY(c);
      float cost = distance(bx - cx, by - cy);
      binCosts.push_back(cost);
    }
    costs.push_back(binCosts);
  }

  TransportationProblem solver(capacities, demands, costs);
  solver.increaseCapacity();
  solver.solve();
  assignment = solver.toAssignment();

  // Reallocate the cells
  std::vector<std::vector<int> > binCells(bins.size());
  for (int i = 0; i < cells.size(); ++i) {
    binCells[assignment[i]].push_back(cells[i]);
  }
  for (int b = 0; b < bins.size(); ++b) {
    auto [x, y] = bins[b];
    setBinCells(x, y, binCells[b]);
  }
}

void DensityLegalizer::refine() {
  // Refine both if they are at the same level, otherwise refine only the
  // coarsest level
  bool doX = levelX() >= levelY();
  bool doY = levelY() >= levelX();

  if (doX && doY && params_.reoptimizationSquareSize >= 2) {
    refineX();
    refineY();
    improveSquareNeighbours();
    improveSquareNeighbours(false, false);
  } else {
    if (doX) {
      refineX();
      improveXNeighbours();
      improveXNeighbours(false);
    }
    if (doY) {
      refineY();
      improveYNeighbours();
      improveYNeighbours(false);
    }
  }
}

void DensityLegalizer::improve() {
  for (int i = 0; i < params_.nbSteps; ++i) {
    improveSquare();
    improveXY();
    improveDiagonals();
  }
}

void DensityLegalizer::improveXNeighbours(bool sameParent) {
  for (int i = 0; i + 1 < nbBinsX(); ++i) {
    if ((parentX(i) == parentX(i + 1)) != sameParent) {
      continue;
    }
    for (int j = 0; j < nbBinsY(); ++j) {
      rebisect(i, j, i + 1, j);
    }
  }
}

void DensityLegalizer::improveYNeighbours(bool sameParent) {
  for (int j = 0; j + 1 < nbBinsY(); ++j) {
    if ((parentY(j) == parentY(j + 1)) != sameParent) {
      continue;
    }
    for (int i = 0; i < nbBinsX(); ++i) {
      rebisect(i, j, i, j + 1);
    }
  }
}

void DensityLegalizer::improveSquareNeighbours(bool sameParentX,
                                               bool sameParentY) {
  for (int i = 0; i < nbBinsX(); ++i) {
    if (i + 1 < nbBinsX() && (parentX(i) == parentX(i + 1)) != sameParentX) {
      continue;
    }
    for (int j = 0; j < nbBinsY(); ++j) {
      if (j + 1 < nbBinsY() && (parentY(j) == parentY(j + 1)) != sameParentY) {
        continue;
      }
      improveRectangle(i, j, 2, 2);
    }
  }
}

void DensityLegalizer::improveXY() {
  if (params_.unidimensionalTransport) {
    improveX();
    improveY();
  } else {
    int nb = params_.reoptimizationLength;
    if (nb == 1) return;

    int mid = (nb + 1) / 2;
    int stride = 2 * mid;

    // X (first)
    improveRectangles(nb, 1, stride, 1, 0, 0);
    // Y (first)
    improveRectangles(1, nb, 1, stride, 0, 0);
    // X (second, overlapping)
    improveRectangles(nb, 1, stride, 1, mid, 0);
    // Y (second, overlapping)
    improveRectangles(1, nb, 1, stride, 0, mid);
  }
}

void DensityLegalizer::improveX() {
  for (int j = 0; j < nbBinsY(); ++j) {
    std::vector<int> cells;
    std::vector<long long> u;
    std::vector<long long> v;
    std::vector<long long> s;
    std::vector<long long> d;
    for (int i = 0; i < nbBinsX(); ++i) {
      v.push_back(binX(i, j));
      d.push_back(binCapacity(i, j));
      for (int c : binCells(i, j)) {
        cells.push_back(c);
        u.push_back(cellTargetX(c));
        s.push_back(cellDemand(c));
      }
    }
    Transportation1d pb(u, v, s, d);
    pb.balanceDemand();
    std::vector<int> assignment = pb.assign();
    std::vector<std::vector<int> > binCells(nbBinsX());
    for (int i = 0; i < cells.size(); ++i) {
      binCells[assignment[i]].push_back(cells[i]);
    }
    for (int i = 0; i < nbBinsX(); ++i) {
      setBinCells(i, j, binCells[i]);
    }
  }
  check();
}

void DensityLegalizer::improveY() {
  for (int i = 0; i < nbBinsX(); ++i) {
    std::vector<int> cells;
    std::vector<long long> u;
    std::vector<long long> v;
    std::vector<long long> s;
    std::vector<long long> d;
    for (int j = 0; j < nbBinsY(); ++j) {
      v.push_back(binY(i, j));
      d.push_back(binCapacity(i, j));
      for (int c : binCells(i, j)) {
        cells.push_back(c);
        u.push_back(cellTargetY(c));
        s.push_back(cellDemand(c));
      }
    }
    Transportation1d pb(u, v, s, d);
    pb.balanceDemand();
    std::vector<int> assignment = pb.assign();
    std::vector<std::vector<int> > binCells(nbBinsY());
    for (int i = 0; i < cells.size(); ++i) {
      binCells[assignment[i]].push_back(cells[i]);
    }
    for (int j = 0; j < nbBinsY(); ++j) {
      setBinCells(i, j, binCells[j]);
    }
  }
  check();
}

void DensityLegalizer::improveDiagonals() {
  int nb = params_.reoptimizationLength;
  if (nb == 1) return;

  int mid = (nb + 1) / 2;
  int stride = 2 * mid;

  // X + Y direction (first)
  improveDiagonalRectangles(1, nb, 1, stride, 0, 0);
  // X - Y direction (first)
  improveDiagonalRectangles(nb, 1, stride, 1, 0, 0);
  // X + Y direction (second, overlapping)
  improveDiagonalRectangles(1, nb, 1, stride, 0, mid);
  // X - Y direction (second, overlapping)
  improveDiagonalRectangles(nb, 1, stride, 1, mid, 0);
}

void DensityLegalizer::improveSquare() {
  int nb = params_.reoptimizationSquareSize;
  if (nb == 1) return;

  int mid = (nb + 1) / 2;
  int stride = 2 * mid;

  improveRectangles(nb, nb, stride, stride, 0, 0);
  improveRectangles(nb, nb, stride, stride, mid, mid);
  improveRectangles(nb, nb, stride, stride, 0, mid);
  improveRectangles(nb, nb, stride, stride, mid, 0);
}

void DensityLegalizer::improveRectangles(int width, int height, int strideX,
                                         int strideY, int startX, int startY) {
  assert(width >= 1 && height >= 1 && strideX >= 1 && strideY >= 1);
  if (width * height == 1) return;
  for (int i = startX; i < nbBinsX(); i += strideX) {
    for (int j = startY; j < nbBinsY(); j += strideY) {
      improveRectangle(i, j, width, height);
    }
  }
}

void DensityLegalizer::improveRectangle(int i, int j, int width, int height) {
  assert(width >= 1 && height >= 1);
  std::vector<std::pair<int, int> > bins;
  for (int k = i; k < nbBinsX() && k < i + width; ++k) {
    for (int l = j; l < nbBinsY() && l < j + height; ++l) {
      bins.emplace_back(k, l);
    }
  }
  reoptimize(bins);
}

void DensityLegalizer::improveDiagonalRectangles(int xmySize, int xpySize,
                                                 int strideX, int strideY,
                                                 int startX, int startY) {
  assert(xmySize >= 1 && xpySize >= 1 && strideX >= 1 && strideY >= 1);
  if (xmySize * xpySize == 1) return;
  for (int i = startX; i < nbBinsX(); i += strideX) {
    for (int j = startY; j < nbBinsY(); j += strideY) {
      std::vector<std::pair<int, int> > bins;
      for (int k = 0; k < xmySize; ++k) {
        for (int l = 0; l < xpySize; ++l) {
          int x = i + k + l;
          int y = j - k + l;
          if (x < 0 || x >= nbBinsX()) {
            continue;
          }
          if (y < 0 || y >= nbBinsY()) {
            continue;
          }
          bins.emplace_back(x, y);
        }
      }
      reoptimize(bins);
    }
  }
}

void DensityLegalizer::run() {
  runCoarsening();
  runRefinement();
}

void DensityLegalizer::runCoarsening() {
  double dist = params_.coarseningLimit * meanDistance();
  while (true) {
    bool doX = levelX() + 1 < nbLevelX();
    bool doY = levelY() + 1 < nbLevelY();
    double distX = placementArea().width() / (double)nbBinsX();
    double distY = placementArea().height() / (double)nbBinsY();
    doX &= distX <= dist;
    doY &= distY <= dist;
    if (!doX && !doY) {
      // Not possible to coarsen
      break;
    }
    if (doX) {
      coarsenX();
    }
    if (doY) {
      coarsenY();
    }
  }
  while (levelX() + 1 < nbLevelX()) {
    coarsenX();
  }
  while (levelY() + 1 < nbLevelY()) {
    coarsenY();
  }
}

void DensityLegalizer::runRefinement() {
  while (levelX() > 0 || levelY() > 0) {
    refine();
    improve();
  }
}
}  // namespace coloquinte
