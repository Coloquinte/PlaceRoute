#include "place_detailed.hpp"

#include <iostream>

#include "legalizer.hpp"

namespace coloquinte {
void DetailedPlacer::place(Circuit &circuit, int effort) {
  Legalizer leg = Legalizer::fromIspdCircuit(circuit);
  leg.run();
  leg.exportPlacement(circuit);
  DetailedPlacer pl(circuit);
  pl.runSwaps(3);
  pl.placement_.exportPlacement(circuit);
}

DetailedPlacer::DetailedPlacer(const Circuit &circuit)
    : placement_(DetailedPlacement::fromIspdCircuit(circuit)),
      xtopo_(IncrNetModel::xTopology(circuit)),
      ytopo_(IncrNetModel::yTopology(circuit)) {}

bool DetailedPlacer::trySwap(int c1, int c2) {
  if (!placement_.canSwap(c1, c2)) {
    return false;
  }
  long long oldValue = value();
  int x1 = placement_.cellX(c1);
  int x2 = placement_.cellX(c2);
  placement_.swap(c1, c2);
  updateCellPos(c1);
  updateCellPos(c2);
  if (value() < oldValue) {
    return true;
  }
  placement_.swapAt(c1, c2, x1, x2);
  updateCellPos(c1);
  updateCellPos(c2);
  return false;
}

bool DetailedPlacer::tryInsert(int c, int row, int pred) {
  if (!placement_.canInsert(c, row, pred)) {
    return false;
  }
  int oldRow = placement_.cellRow(c);
  int oldPred = placement_.cellPred(c);
  int oldX = placement_.cellX(c);
  long long oldValue = value();
  placement_.insert(c, row, pred);
  updateCellPos(c);
  if (value() < oldValue) {
    return true;
  }
  placement_.insertAt(c, oldRow, oldPred, oldX);
  updateCellPos(c);
  return false;
}

void DetailedPlacer::runSwaps(int nbNeighbours) {
  for (int i = 0; i < placement_.nbRows(); ++i) {
    runSwapsOneRow(i, nbNeighbours);
  }
  for (int i = 0; i + 1 < placement_.nbRows(); ++i) {
    runSwapsTwoRows(i, i + 1, nbNeighbours);
  }
  for (int i = 0; i + 1 < placement_.nbRows(); ++i) {
    runSwapsTwoRows(i + 1, i, nbNeighbours);
  }
}

void DetailedPlacer::runSwapsOneRow(int row, int nbNeighbours) {
  for (int c = placement_.rowFirstCell(row); c != -1;
       c = placement_.cellNext(c)) {
    for (int p = c, i = 0; (p != -1) && (i <= nbNeighbours);
         p = placement_.cellPred(p), ++i) {
      if (trySwap(p, c)) {
        // Continue before the cell
        p = c;
      }
    }
    for (int p = c, i = 0; (p != -1) && (i <= nbNeighbours);
         p = placement_.cellNext(p), ++i) {
      if (trySwap(p, c)) {
        // Continue after the cell
        p = c;
      }
    }
  }
}

void DetailedPlacer::runSwapsTwoRows(int r1, int r2, int nbNeighbours) {
  int c1 = placement_.rowFirstCell(r1);
  int c2 = placement_.rowFirstCell(r2);
  // TODO
}

void DetailedPlacer::updateCellPos(int c) {
  xtopo_.updateCellPos(c, placement_.cellX(c));
  ytopo_.updateCellPos(c, placement_.cellY(c));
}

void DetailedPlacer::check() const {
  placement_.check();
  xtopo_.check();
  ytopo_.check();
}
}  // namespace coloquinte
