
#include "place_global/density_legalizer.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>

DensityLegalizer::DensityLegalizer(DensityGrid grid,
                                   std::vector<int> cellDemand)
    : HierarchicalDensityPlacement(grid, cellDemand) {
  costModel_ = LegalizationModel::L1;
  cellTargetX_.assign(nbCells(), 0.0f);
  cellTargetY_.assign(nbCells(), 0.0f);
}

DensityLegalizer::DensityLegalizer(HierarchicalDensityPlacement pl)
    : HierarchicalDensityPlacement(pl) {
  costModel_ = LegalizationModel::L1;
  cellTargetX_.assign(nbCells(), 0.0f);
  cellTargetY_.assign(nbCells(), 0.0f);
}

DensityLegalizer DensityLegalizer::fromIspdCircuit(const Circuit &circuit,
                                                   float sizeFactor) {
  return DensityLegalizer(
      HierarchicalDensityPlacement::fromIspdCircuit(circuit, sizeFactor));
}

void DensityLegalizer::exportPlacement(Circuit &circuit) {
  std::vector<float> cellX = simpleCoordX();
  std::vector<float> cellY = simpleCoordY();
  assert (nbCells() == circuit.nbCells());
  for (int i = 0; i < circuit.nbCells(); ++i) {
    if (circuit.fixed(i)) continue;
    circuit.cellX[i] = std::round(cellX[i]);
    circuit.cellY[i] = std::round(cellY[i]);
  }
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
  std::cout << "Mean dist: " << meanDistance(costModel_) << std::endl;
  std::cout << "RMS dist " << rmsDistance(costModel_) << std::endl;
  std::cout << "Max dist " << maxDistance(costModel_) << std::endl;
  if (!verbose) return;
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
    float cx1, float cy1, float cx2, float cy2, std::vector<int> cells) const {
  // TODO: always add a secondary objective using squared distance
  std::vector<std::pair<float, int> > cellCosts;
  for (int c : cells) {
    float x = cellTargetX_[c];
    float y = cellTargetY_[c];
    float cost =
        norm(x - cx1, y - cy1, costModel_) - norm(x - cx2, y - cy2, costModel_);
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
  while (splitPos > 0 && demand1 - capa1 > 0) {
    int dem = cellDemand_[cellCosts[splitPos - 1].second];
    if (demand1 - capa1 < demand2 - capa2 + dem) {
      // Stop if overflow would be bigger on the other side
      break;
    }
    demand1 -= dem;
    demand2 += dem;
    --splitPos;
  }
  // Remove from the right if overflowed
  while (splitPos < cellCosts.size() && demand2 - capa2 > 0) {
    int dem = cellDemand_[cellCosts[splitPos].second];
    if (demand2 - capa2 < demand1 - capa1 + dem) {
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
  if (x1 == x2 && y1 == y2) return;

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

void DensityLegalizer::refine() {
  if (levelX() >= levelY()) {
    refineX();
    improveX(true);
  } else {
    refineY();
    improveY(true);
  }
}

void DensityLegalizer::improve() {
  improveX(false);
  improveY(false);
  improveX(true);
  improveY(true);
}

void DensityLegalizer::improveX(bool sameParent) {
  for (int i = 0; i + 1 < nbBinsX(); ++i) {
    if ((parentX(i) == parentX(i + 1)) != sameParent) {
      continue;
    }
    for (int j = 0; j < nbBinsY(); ++j) {
      rebisect(i, j, i + 1, j);
    }
  }
}

void DensityLegalizer::improveY(bool sameParent) {
  for (int j = 0; j + 1 < nbBinsY(); ++j) {
    if ((parentY(j) == parentY(j + 1)) != sameParent) {
      continue;
    }
    for (int i = 0; i < nbBinsX(); ++i) {
      rebisect(i, j, i, j + 1);
    }
  }
}

void DensityLegalizer::run() {
  coarsenFully();
  while (levelX() > 0 || levelY() > 0) {
    refine();
    improve();
  }
}