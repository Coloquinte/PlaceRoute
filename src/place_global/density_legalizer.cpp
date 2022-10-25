
#include "place_global/density_legalizer.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <numeric>
#include <utility>

#include "place_global/transportation.hpp"
#include "utils/norm.hpp"

namespace coloquinte {
DensityLegalizer::Parameters::Parameters() {
  costModel = LegalizationModel::L1;
  nbSteps = 1;
  reoptimizationLength = 2;
  reoptimizationSquareSize = 1;
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

std::vector<float> DensityLegalizer::allDistances(
    LegalizationModel model) const {
  std::vector<float> cellX = simpleCoordX();
  std::vector<float> cellY = simpleCoordY();
  std::vector<float> distances;
  distances.reserve(nbCells());
  for (int i = 0; i < nbCells(); ++i) {
    float dx = cellTargetX(i) - cellX[i];
    float dy = cellTargetY(i) - cellY[i];
    distances.push_back(norm(dx, dy, model));
  }
  return distances;
}

float DensityLegalizer::meanDistance(LegalizationModel model) const {
  std::vector<float> dist = allDistances(model);
  float disp = 0.0f;
  for (int i = 0; i < nbCells(); ++i) {
    disp += cellDemand(i) * dist[i];
  }
  return disp / totalDemand();
}

float DensityLegalizer::rmsDistance(LegalizationModel model) const {
  std::vector<float> dist = allDistances(model);
  float disp = 0.0f;
  for (int i = 0; i < nbCells(); ++i) {
    disp += cellDemand(i) * dist[i] * dist[i];
  }
  return std::sqrt(disp / totalDemand());
}

float DensityLegalizer::maxDistance(LegalizationModel model) const {
  std::vector<float> dist = allDistances(model);
  return *std::max_element(dist.begin(), dist.end());
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
  std::cout << "Mean dist: " << meanDistance(params_.costModel) << std::endl;
  std::cout << "RMS dist " << rmsDistance(params_.costModel) << std::endl;
  std::cout << "Max dist " << maxDistance(params_.costModel) << std::endl;
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
    float cost = norm(x - cx1, y - cy1, params_.costModel) -
                 norm(x - cx2, y - cy2, params_.costModel);
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
  std::vector<std::pair<int, int> > bins;
  for (auto [x, y] : binCandidates) {
    if (binCapacity(x, y) > 0) {
      bins.emplace_back(x, y);
    }
  }
  if (bins.size() <= 1) {
    return;
  }

  if (bins.size() == 2) {
    auto [x1, y1] = bins[0];
    auto [x2, y2] = bins[1];
    rebisect(x1, y1, x2, y2);
    return;
  }

  std::vector<int> cells;
  for (auto [x, y] : bins) {
    cells.insert(cells.end(), binCells_[x][y].begin(), binCells_[x][y].end());
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
  long long totalDemand = std::accumulate(demands.begin(), demands.end(), 0LL);
  long long totalCapa =
      std::accumulate(capacities.begin(), capacities.end(), 0LL);
  long long missingCapa = std::max(totalDemand - totalCapa, 0LL);
  for (int i = 0; i < bins.size(); ++i) {
    capacities[i] += (missingCapa + bins.size() - 1) / bins.size();
  }

  std::vector<std::vector<float> > costs;
  for (auto [x, y] : bins) {
    float bx = binX(x, y);
    float by = binY(x, y);
    std::vector<float> binCosts;
    for (int c : cells) {
      float cx = cellTargetX(c);
      float cy = cellTargetY(c);
      float cost = norm(bx - cx, by - cy, params_.costModel);
      binCosts.push_back(cost);
    }
    costs.push_back(binCosts);
  }

  // Solve
  std::vector<std::vector<long long> > allocation =
      solveTransportation(capacities, demands, costs);

  // Reallocate the cells
  std::vector<std::vector<int> > binCells(bins.size());
  for (int i = 0; i < cells.size(); ++i) {
    int cell = cells[i];
    int bestBin = 0;
    int highestAlloc = 0;
    for (int b = 0; b < bins.size(); ++b) {
      if (allocation[b][i] > highestAlloc) {
        highestAlloc = allocation[b][i];
        bestBin = b;
      }
    }
    binCells[bestBin].push_back(cell);
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
  if (doX) {
    refineX();
    improveXNeighbours();
  }
  if (doY) {
    refineY();
    improveYNeighbours();
  }
}

void DensityLegalizer::improve() {
  for (int i = 0; i < params_.nbSteps; ++i) {
    improveSquare();
    improveX();
    improveY();
    improveDiag();
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

void DensityLegalizer::improveX() {
  int nb = params_.reoptimizationLength;
  assert(nb >= 2);
  for (int j = 0; j < nbBinsY(); ++j) {
    for (int i = 0; i < nbBinsX(); i += nb / 2) {
      std::vector<std::pair<int, int> > bins;
      for (int k = 0; i + k < nbBinsX() && k < nb; ++k) {
        bins.emplace_back(i + k, j);
      }
      reoptimize(bins);
    }
  }
}

void DensityLegalizer::improveY() {
  int nb = params_.reoptimizationLength;
  assert(nb >= 2);
  for (int i = 0; i < nbBinsX(); ++i) {
    for (int j = 0; j < nbBinsY(); j += nb / 2) {
      std::vector<std::pair<int, int> > bins;
      for (int k = 0; j + k < nbBinsY() && k < nb; ++k) {
        bins.emplace_back(i, j + k);
      }
      reoptimize(bins);
    }
  }
}

void DensityLegalizer::improveDiag() {
  int nb = params_.reoptimizationLength;
  assert(nb >= 2);
  for (int i = 0; i + 1 < nbBinsX(); i += nb / 2) {
    for (int j = 0; j + 1 < nbBinsY(); j += nb / 2) {
      std::vector<std::pair<int, int> > bins;
      for (int k = 0; i + k < nbBinsX() && j + k < nbBinsY() && k < nb; ++k) {
        bins.emplace_back(i + k, j + k);
      }
      reoptimize(bins);
    }
  }
}

void DensityLegalizer::improveSquare() {
  int nb = params_.reoptimizationSquareSize;
  if (nb == 1) return;

  for (int i = 0; i < nbBinsX(); i += nb / 2) {
    for (int j = 0; j < nbBinsY(); j += nb / 2) {
      std::vector<std::pair<int, int> > bins;
      for (int k = i; k < nbBinsX() && k < i + nb; ++k) {
        for (int l = j; l < nbBinsY() && l < j + nb; ++l) {
          bins.emplace_back(k, l);
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

void DensityLegalizer::runCoarsening() { coarsenFully(); }

void DensityLegalizer::runRefinement() {
  while (levelX() > 0 || levelY() > 0) {
    refine();
    improve();
  }
}
}  // namespace coloquinte