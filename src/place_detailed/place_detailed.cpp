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
        // Continue iterating on the same row, from p
        if (closest == p) closest = c;
        std::swap(c, p);
      }
    }
    for (int p = closest, i = 0; (p != -1) && (i <= nbNeighbours);
         p = placement_.cellNext(p), ++i) {
      if (trySwap(p, c)) {
        // Continue iterating on the same row, from p
        if (closest == p) closest = c;
        std::swap(c, p);
      }
    }
  }
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
}
}  // namespace coloquinte
