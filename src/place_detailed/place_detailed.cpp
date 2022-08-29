#include "place_detailed.hpp"

#include <lemon/network_simplex.h>
#include <lemon/smart_graph.h>

#include <iostream>
#include <unordered_set>

#include "legalizer.hpp"
#include "row_neighbourhood.hpp"

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

void DetailedPlacer::doSwap(int c1, int c2) {
  assert(placement_.canSwap(c1, c2));
  placement_.swap(c1, c2);
  updateCellPos(c1);
  updateCellPos(c2);
}

void DetailedPlacer::doInsert(int c, int row, int pred) {
  assert(placement_.canInsert(c, row, pred));
  placement_.insert(c, row, pred);
  updateCellPos(c);
}

void DetailedPlacer::runSwaps(int nbRows, int nbNeighbours) {
  // Optimize each row internally
  for (int i = 0; i < placement_.nbRows(); ++i) {
    runSwapsOneRow(i, nbNeighbours);
  }
  RowNeighbourhood rowsNeighbours(placement_.rows(), nbNeighbours);
  //  Optimize each row with neighbours after it
  for (int i = 0; i < placement_.nbRows(); ++i) {
    for (int j : rowsNeighbours.rowsAbove(i))
      runSwapsTwoRowsAmplify(i, j, nbNeighbours);
    for (int j : rowsNeighbours.rowsRight(i))
      runSwapsTwoRowsAmplify(i, j, nbNeighbours);
  }
  // Optimize each row with neighbours before it;
  // We do both for symmetry and to allow large cell movements
  for (int i = placement_.nbRows() - 1; i >= 1; --i) {
    for (int j : rowsNeighbours.rowsBelow(i))
      runSwapsTwoRowsAmplify(i, j, nbNeighbours);
    for (int j : rowsNeighbours.rowsLeft(i))
      runSwapsTwoRowsAmplify(i, j, nbNeighbours);
  }
}

void DetailedPlacer::runInserts(int nbRows, int nbNeighbours) {
  // Optimize each row internally
  for (int i = 0; i < placement_.nbRows(); ++i) {
    runInsertsOneRow(i, nbNeighbours);
  }
  // Optimize each row with neighbours after it
  for (int d = 1; d <= nbRows; ++d) {
    for (int i = 0; i + d < placement_.nbRows(); ++i) {
      runInsertsTwoRows(i, i + d, nbNeighbours);
    }
  }
  // Optimize each row with neighbours before it;
  // We do both for symmetry and to allow large cell movements
  for (int d = 1; d <= nbRows; ++d) {
    for (int i = placement_.nbRows() - 1; i - d >= 0; --i) {
      runInsertsTwoRows(i, i - d, nbNeighbours);
    }
  }
}

void DetailedPlacer::runSwapsOneRow(int row, int nbNeighbours) {
  std::vector<int> cells = placement_.rowCells(row);
  for (int i = 0; i < cells.size(); ++i) {
    int b = std::max(0, i - nbNeighbours);
    int e = std::min((int)cells.size(), i + nbNeighbours + 1);
    std::vector<int> candidates(cells.begin() + b, cells.begin() + e);
    bestSwap(cells[i], candidates);
  }
}

void DetailedPlacer::runInsertsOneRow(int row, int nbNeighbours) {
  std::vector<int> cells = placement_.rowCells(row);
  // Consider insertion before the first cell
  cells.insert(cells.begin(), -1);
  for (int i = 1; i < cells.size(); ++i) {
    int b = std::max(0, i - nbNeighbours);
    int e = std::min((int)cells.size(), i + nbNeighbours + 1);
    std::vector<int> candidates(cells.begin() + b, cells.begin() + e);
    bestInsert(cells[i], row, candidates);
  }
}

void DetailedPlacer::runSwapsTwoRows(int r1, int r2, int nbNeighbours) {
  std::vector<int> cells1 = placement_.rowCells(r1);
  std::vector<int> cells2 = placement_.rowCells(r2);
  std::vector<int> closestIndex = computeClosestIndexInRow(cells1, cells2);
  for (int i = 0; i < cells1.size(); ++i) {
    int closest = closestIndex[i];
    int b = std::max(0, closest - nbNeighbours);
    int e = std::min((int)cells2.size(), closest + nbNeighbours + 1);
    std::vector<int> candidates(cells2.begin() + b, cells2.begin() + e);
    bestSwap(cells1[i], candidates);
  }
}

void DetailedPlacer::runSwapsTwoRowsAmplify(int r1, int r2, int nbNeighbours) {
  int from = placement_.rowFirstCell(r2);
  for (int c = placement_.rowFirstCell(r1); c != -1;
       c = placement_.cellNext(c)) {
    while (bestSwapUpdate(c, from, nbNeighbours))
      ;
    from = findCellAfter(c, from);
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
    int e = std::min((int)cells2.size(), closest + nbNeighbours + 1);
    std::vector<int> candidates(cells2.begin() + b, cells2.begin() + e);
    bestInsert(cells1[i], r2, candidates);
  }
}

bool DetailedPlacer::bestSwap(int c, const std::vector<int> &candidates) {
  long long bestValue = value();
  bool found = false;
  int bestCandidate = -1;
  for (int candidate : candidates) {
    auto [feasible, val] = valueOnSwap(c, candidate);
    if (feasible && val < bestValue) {
      found = true;
      bestCandidate = candidate;
    }
  }
  if (found) {
    doSwap(c, bestCandidate);
  }
  return found;
}

bool DetailedPlacer::bestInsert(int c, int row,
                                const std::vector<int> &candidates) {
  long long bestValue = value();
  bool found = false;
  int bestCandidate = -1;
  for (int candidate : candidates) {
    auto [feasible, val] = valueOnInsert(c, row, candidate);
    if (feasible && val < bestValue) {
      found = true;
      bestCandidate = candidate;
    }
  }
  if (found) {
    doInsert(c, row, bestCandidate);
  }
  return found;
}

bool DetailedPlacer::bestSwapUpdate(int &c, int &from, int nbNeighbours) {
  long long bestValue = value();
  bool found = false;
  int bestCandidate = -1;
  for (int candidate = from; candidate != -1;
       candidate = placement_.cellNext(candidate)) {
    auto [feasible, val] = valueOnSwap(c, candidate);
    if (feasible && val < bestValue) {
      found = true;
      bestCandidate = candidate;
    }
  }
  for (int candidate = from; candidate != -1;
       candidate = placement_.cellPred(candidate)) {
    auto [feasible, val] = valueOnSwap(c, candidate);
    if (feasible && val < bestValue) {
      found = true;
      bestCandidate = candidate;
    }
  }
  if (found) {
    doSwap(c, bestCandidate);
    if (bestCandidate == from) from = c;
    c = bestCandidate;
  }
  return found;
}

std::pair<bool, long long> DetailedPlacer::valueOnSwap(int c1, int c2) {
  if (!placement_.canSwap(c1, c2)) {
    return std::make_pair(false, std::numeric_limits<long long>::max());
  }
  auto [newP1, newP2] = placement_.positionsOnSwap(c1, c2);
  auto oldP1 = placement_.cellPos(c1);
  auto oldP2 = placement_.cellPos(c2);
  updateCellPos(c1, newP1);
  updateCellPos(c2, newP2);
  long long newValue = value();
  updateCellPos(c1, oldP1);
  updateCellPos(c2, oldP2);
  return std::make_pair(true, newValue);
}

std::pair<bool, long long> DetailedPlacer::valueOnInsert(int c, int row,
                                                         int pred) {
  if (!placement_.canInsert(c, row, pred)) {
    return std::make_pair(false, std::numeric_limits<long long>::max());
  }
  auto newP = placement_.positionOnInsert(c, row, pred);
  auto oldP = placement_.cellPos(c);
  updateCellPos(c, newP);
  long long newValue = value();
  updateCellPos(c, oldP);
  return std::make_pair(true, newValue);
}

int DetailedPlacer::findCellAfter(int target, int fromCell) const {
  int c = fromCell;
  if (c == -1) return c;
  while (true) {
    int nextC = placement_.cellNext(c);
    if (nextC == -1) break;
    if (placement_.cellX(nextC) > placement_.cellX(target)) break;
    c = nextC;
  }
  return c;
}

int DetailedPlacer::findCellBefore(int target, int fromCell) const {
  int c = fromCell;
  if (c == -1) return c;
  while (true) {
    int nextC = placement_.cellPred(c);
    if (nextC == -1) break;
    if (placement_.cellX(nextC) + placement_.cellWidth(nextC) <
        placement_.cellX(target) + placement_.cellWidth(target))
      break;
    c = nextC;
  }
  return c;
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
  updateCellPos(c, placement_.cellPos(c));
}

void DetailedPlacer::updateCellPos(int c, Point pos) {
  xtopo_.updateCellPos(c, pos.x);
  ytopo_.updateCellPos(c, pos.y);
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
