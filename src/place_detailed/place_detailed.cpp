#include "place_detailed.hpp"

#include <lemon/network_simplex.h>
#include <lemon/smart_graph.h>

#include <iostream>
#include <unordered_set>

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
  pl.runShifts(effort + 2);
  for (int i = 0; i < effort / 3 + 1; ++i) {
    pl.runSwaps(effort / 2 + 1, effort / 2 + 1);
  }
  pl.runShifts(effort + 2);
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
    runInsertsOneRow(i, nbNeighbours);
  }
  // Optimize each row with neighbours after it
  for (int d = 1; d <= nbRows; ++d) {
    for (int i = 0; i + d < placement_.nbRows(); ++i) {
      runSwapsTwoRows(i, i + d, nbNeighbours);
      runInsertsTwoRows(i, i + d, nbNeighbours);
    }
  }
  // Optimize each row with neighbours before it;
  // We do both for symmetry and to allow large cell movements
  for (int d = 1; d <= nbRows; ++d) {
    for (int i = placement_.nbRows() - 1; i - d >= 0; --i) {
      runSwapsTwoRows(i, i - d, nbNeighbours);
      runInsertsTwoRows(i, i - d, nbNeighbours);
    }
  }
}

void DetailedPlacer::runSwapsOneRow(int row, int nbNeighbours) {
  std::vector<int> cells = placement_.rowCells(row);
  for (int i = 0; i < cells.size(); ++i) {
    int b = std::max(0, i - nbNeighbours);
    int e = std::min((int)cells.size() - 1, i + nbNeighbours);
    for (int j = b; j <= e; ++j) {
      if (trySwap(cells[i], cells[j])) {
        std::swap(cells[i], cells[j]);
      }
    }
  }
}

void DetailedPlacer::runInsertsOneRow(int row, int nbNeighbours) {
  std::vector<int> cells = placement_.rowCells(row);
  // Consider insertion before the first cell
  cells.insert(cells.begin(), -1);
  for (int i = 1; i < cells.size(); ++i) {
    int b = std::max(0, i - nbNeighbours);
    int e = std::min((int)cells.size() - 1, i + nbNeighbours);
    for (int j = b; j <= e; ++j) {
      tryInsert(cells[i], row, cells[j]);
    }
  }
}

void DetailedPlacer::runSwapsTwoRows(int r1, int r2, int nbNeighbours) {
  std::vector<int> cells1 = placement_.rowCells(r1);
  std::vector<int> cells2 = placement_.rowCells(r2);
  std::vector<int> closestIndex = computeClosestIndexInRow(cells1, cells2);
  for (int i = 0; i < cells1.size(); ++i) {
    int closest = closestIndex[i];
    int b = std::max(0, closest - nbNeighbours);
    int e = std::min((int)cells2.size() - 1, closest + nbNeighbours);
    for (int j = b; j <= e; ++j) {
      trySwap(cells1[i], cells2[j]);
    }
  }
}

void DetailedPlacer::runInsertsTwoRows(int r1, int r2, int nbNeighbours) {
  std::vector<int> cells1 = placement_.rowCells(r1);
  std::vector<int> cells2 = placement_.rowCells(r2);
  cells2.insert(cells2.begin(), -1);
  std::vector<int> closestIndex = computeClosestIndexInRow(cells1, cells2);
  for (int i = 0; i < cells1.size(); ++i) {
    int closest = closestIndex[i];
    int b = std::max(0, closest - nbNeighbours);
    int e = std::min((int)cells2.size() - 1, closest + nbNeighbours);
    for (int j = b; j <= e; ++j) {
      tryInsert(cells1[i], r2, cells2[j]);
    }
  }
}

std::vector<int> DetailedPlacer::computeClosestIndexInRow(
    const std::vector<int> &row1Cells,
    const std::vector<int> &row2Cells) const {
  if (row2Cells.empty()) {
    return std::vector<int>(row1Cells.size(), 0);
  }
  std::vector<int> ret;
  int closest = 0;
  for (int i = 0; i < row1Cells.size(); ++i) {
    int x = placement_.cellX(row1Cells[i]);
    while (true) {
      int c = row2Cells[closest];
      if (closest == row2Cells.size() - 1) {
        break;
      }
      if (c != -1 && placement_.cellX(c) < x) {
        break;
      }
      ++closest;
    }
    ret.push_back(closest);
  }
  return ret;
}

void DetailedPlacer::runShifts(int nbRows) {
  if (nbRows < 2) return;
  for (int r = 0; r < placement_.nbRows(); r += nbRows / 2) {
    std::vector<int> cells;
    for (int i = r; i < std::min(r + nbRows, placement_.nbRows()); ++i) {
      for (int c : placement_.rowCells(i)) {
        cells.push_back(c);
      }
    }
    optimizeShift(cells);
  }
  placement_.check();
}

void DetailedPlacer::optimizeShift(const std::vector<int> &cells) {
  std::unordered_set<int> cell_set(cells.begin(), cells.end());
  if (cell_set.size() != cells.size()) {
    throw std::runtime_error("The given cells are not unique");
  }
  std::unordered_set<int> net_set;
  for (int c : cells) {
    for (int i = 0; i < xtopo_.nbCellPins(c); ++i) {
      int net = xtopo_.pinNet(c, i);
      net_set.insert(net);
    }
  }
  std::vector<int> nets(net_set.begin(), net_set.end());
  std::sort(nets.begin(), nets.end());
  // Create a graph with the cells as using namespace lemon;
  using namespace lemon;
  DIGRAPH_TYPEDEFS(SmartDigraph);

  // Create a graph with the cells and bounds of the nets as node
  SmartDigraph g;

  std::unordered_map<int, Node> cell_nodes(cells.size());
  for (int c : cells) {
    cell_nodes[c] = g.addNode();
  }
  std::unordered_map<int, Node> Lnet_nodes, Unet_nodes;
  for (int net : nets) {
    Lnet_nodes[net] = g.addNode();
    Unet_nodes[net] = g.addNode();
  }

  // Two nodes for position constraints
  Node fixed = g.addNode();
  typedef std::pair<SmartDigraph::Arc, int> arc_pair;
  typedef std::pair<SmartDigraph::Node, int> node_pair;

  // The arcs corresponding to constraints of the original problem
  std::vector<arc_pair> constraint_arcs;

  // Now we add every positional constraint, which becomes an arc in the
  // min-cost flow problem
  for (int c : cells) {
    int pred = placement_.cellPred(c);
    int next = placement_.cellNext(c);
    if (cell_set.count(next)) {
      // Two movable cells
      auto A = g.addArc(cell_nodes[next], cell_nodes[c]);
      constraint_arcs.push_back(arc_pair(A, -placement_.cellWidth(c)));
    }
    if (cell_set.count(pred) == 0) {
      // Predecessor fixed
      int boundary = placement_.boundaryBefore(c);
      auto A = g.addArc(cell_nodes[c], fixed);
      constraint_arcs.push_back(arc_pair(A, -boundary));
    }
    if (cell_set.count(next) == 0) {
      // Successor fixed
      int boundary = placement_.boundaryAfter(c);
      auto A = g.addArc(fixed, cell_nodes[c]);
      constraint_arcs.push_back(
          arc_pair(A, boundary - placement_.cellWidth(c)));
    }
  }

  // One arc for every pin of every net: arcs too
  for (int net : nets) {
    for (int i = 0; i < xtopo_.nbNetPins(net); ++i) {
      int c = xtopo_.pinCell(net, i);
      int pin_offs = xtopo_.netPinOffset(net, i);
      if (cell_set.count(c)) {
        Arc Al = g.addArc(cell_nodes[c], Lnet_nodes[net]);
        constraint_arcs.push_back(arc_pair(Al, pin_offs));
        Arc Ar = g.addArc(Unet_nodes[net], cell_nodes[c]);
        constraint_arcs.push_back(arc_pair(Ar, -pin_offs));
      } else {  // Fixed offset
        auto Al = g.addArc(fixed, Lnet_nodes[net]);
        constraint_arcs.push_back(arc_pair(Al, placement_.cellX(c) + pin_offs));
        auto Ar = g.addArc(Unet_nodes[net], fixed);
        constraint_arcs.push_back(
            arc_pair(Ar, -placement_.cellX(c) - pin_offs));
      }
    }
  }

  // Then the only capacitated arcs: the ones for the nets
  std::vector<node_pair> net_supplies;
  for (int net : nets) {
    net_supplies.push_back(node_pair(Unet_nodes[net], 1));
    net_supplies.push_back(node_pair(Lnet_nodes[net], -1));
  }

  // Create the maps to have cost and capacity for the arcs
  IntArcMap cost(g, 0);
  IntArcMap capacity(g, nets.size());
  IntNodeMap supply(g, 0);

  for (arc_pair A : constraint_arcs) {
    cost[A.first] = A.second;
  }

  for (node_pair N : net_supplies) {
    supply[N.first] = N.second;
  }

  // Then we (hope the solver can) solve it
  NetworkSimplex<SmartDigraph> ns(g);
  ns.supplyMap(supply).costMap(cost);
  auto res = ns.run();
  if (res != ns.OPTIMAL) {
    throw std::runtime_error("Could not solve the network flow optimallt");
  }

  // And we get the new positions as the dual values of the current solution
  // (compared to the fixed pin)
  for (int c : cells) {
    int pos = ns.potential(cell_nodes[c]) - ns.potential(fixed);
    placement_.cellX_[c] = pos;
    xtopo_.updateCellPos(c, pos);
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
  for (int c = 0; c < placement_.nbCells(); ++c) {
    if (placement_.isPlaced(c)) {
      if (xtopo_.cellPos(c) != placement_.cellX(c)) {
        throw std::runtime_error("X placement is not consistant");
      }
      if (ytopo_.cellPos(c) != placement_.cellY(c)) {
        throw std::runtime_error("Y placement is not consistant");
      }
    }
  }
}
}  // namespace coloquinte
