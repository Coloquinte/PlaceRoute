#include "place_detailed.hpp"

#include <iostream>

#include "legalizer.hpp"

namespace coloquinte {
void DetailedPlacer::place(Circuit &circuit, int effort) {
  std::cout << "Wirelength before legalization: " << circuit.hpwl()
            << std::endl;
  Legalizer leg = Legalizer::fromIspdCircuit(circuit);
  leg.run();
  leg.exportPlacement(circuit);
  std::cout << "Wirelength after legalization: " << circuit.hpwl() << std::endl;
  DetailedPlacer pl(circuit);
  pl.check();
  for (int i = 0; i < effort / 3 + 1; ++i) {
    pl.runSwaps(effort / 2 + 1, effort / 2 + 1);
  }
  pl.check();
  pl.placement_.exportPlacement(circuit);
  std::cout << "Wirelength after detailed placement: " << circuit.hpwl()
            << std::endl;
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

void DetailedPlacer::runSwaps(int nbRows, int nbNeighbours) {
  // Optimize each row internally
  for (int i = 0; i < placement_.nbRows(); ++i) {
    runSwapsOneRow(i, nbNeighbours);
  }
  // Optimize each row with neighbours after it
  for (int d = 1; d <= nbRows; ++d) {
    for (int i = 0; i + d < placement_.nbRows(); ++i) {
      runSwapsTwoRows(i, i + d, nbNeighbours);
    }
  }
  // Optimize each row with neighbours before it;
  // We do both for symmetry and to allow large cell movements
  for (int d = 1; d <= nbRows; ++d) {
    for (int i = placement_.nbRows() - 1; i - d >= 0; --i) {
      runSwapsTwoRows(i, i - d, nbNeighbours);
    }
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
  int closest = placement_.rowFirstCell(r2);
  for (int c = placement_.rowFirstCell(r1); c != -1;
       c = placement_.cellNext(c)) {
    // Update the closest cell
    while (true) {
      int nextC = placement_.cellNext(closest);
      if (nextC == -1) break;
      if (placement_.cellX(nextC) > placement_.cellX(c)) break;
      closest = nextC;
    }
    for (int p = closest, i = 0; (p != -1) && (i <= nbNeighbours);
         p = placement_.cellPred(p), ++i) {
      if (trySwap(p, c)) {
        std::swap(c, p);
      }
    }
    for (int p = closest, i = 0; (p != -1) && (i <= nbNeighbours);
         p = placement_.cellNext(p), ++i) {
      if (trySwap(p, c)) {
        std::swap(c, p);
      }
    }
  }
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
