#include "place_detailed.hpp"

#include <lemon/network_simplex.h>
#include <lemon/smart_graph.h>

#include <chrono>
#include <iomanip>
#include <iostream>
#include <unordered_set>

#include "legalizer.hpp"
#include "row_neighbourhood.hpp"

namespace coloquinte {

void DetailedPlacerParameters::check() const {
  if (nbPasses < 0) {
    throw std::runtime_error(
        "Number of detailed placement passes must be non-negative");
  }
  if (localSearchNbNeighbours < 0) {
    throw std::runtime_error(
        "Number of detailed placement neighbour cells must be non-negative");
  }
  if (localSearchNbRows < 0) {
    throw std::runtime_error(
        "Number of detailed placement rows must be non-negative");
  }
  if (shiftNbRows < 0) {
    throw std::runtime_error(
        "Number of detailed placement shift rows must be non-negative");
  }
  if (shiftMaxNbCells < 0) {
    throw std::runtime_error(
        "Number of detailed placement shift cells must be non-negative");
  }
  if (legalizationCostModel != LegalizationModel::L1) {
    throw std::runtime_error("Only L1 legalization model is supported");
  }
  if (legalizationOrderingWidth > 2.0 || legalizationOrderingWidth < -1.0) {
    throw std::runtime_error(
        "Legalization ordering width should be small (0 < ... < 1)");
  }
  if (legalizationOrderingY > 0.2 || legalizationOrderingY < -0.2) {
    throw std::runtime_error(
        "Legalization ordering y should be small (-0.1 < ... < 0.1)");
  }
}

void DetailedPlacer::legalize(Circuit &circuit,
                              const DetailedPlacerParameters &params,
                              std::optional<PlacementCallback> callback) {
  params.check();
  std::cout << "Legalization starting (WL " << circuit.hpwl() << ")"
            << std::endl;
  auto startTime = std::chrono::steady_clock::now();
  Legalizer leg = Legalizer::fromIspdCircuit(circuit);
  leg.run(params);
  auto endTime = std::chrono::steady_clock::now();
  std::chrono::duration<float> duration = endTime - startTime;
  float distance = leg.meanDistance(params.legalizationCostModel);
  leg.exportPlacement(circuit);
  std::cout << "Legalization done (WL " << circuit.hpwl();
  std::cout << std::fixed << std::setprecision(1) << ", dist " << distance;
  std::cout << std::fixed << std::setprecision(2) << ") in " << duration.count()
            << "s" << std::endl;
  if (callback.has_value()) {
    callback.value()(PlacementStep::Detailed);
  }
}

void DetailedPlacer::place(Circuit &circuit,
                           const DetailedPlacerParameters &params,
                           std::optional<PlacementCallback> callback) {
  legalize(circuit, params, callback);
  params.check();
  std::cout << "Detailed placement starting" << std::endl;
  auto startTime = std::chrono::steady_clock::now();
  DetailedPlacer pl(circuit, params);
  pl.callback_ = callback;
  pl.check();
  pl.run();
  pl.check();
  auto endTime = std::chrono::steady_clock::now();
  pl.exportPlacement(circuit);
  std::chrono::duration<float> duration = endTime - startTime;
  std::cout << "Detailed placement done (WL " << circuit.hpwl();
  std::cout << std::fixed << std::setprecision(2) << ") in " << duration.count()
            << "s" << std::endl;
}

DetailedPlacer::DetailedPlacer(Circuit &circuit,
                               const DetailedPlacerParameters &params)
    : placement_(DetailedPlacement::fromIspdCircuit(circuit)),
      xtopo_(IncrNetModel::xTopology(circuit)),
      ytopo_(IncrNetModel::yTopology(circuit)),
      params_(params),
      circuit_(circuit) {}

void DetailedPlacer::run() {
  for (int i = 1; i <= params_.nbPasses; ++i) {
    runSwaps(params_.localSearchNbNeighbours, params_.localSearchNbRows);
    auto swapValue = value();
    runShifts(params_.shiftNbRows, params_.shiftMaxNbCells);
    auto shiftValue = value();
    std::cout << "#" << i << ":\tSwaps " << swapValue << "\tShifts "
              << shiftValue << std::endl;
    callback();
  }
}

void DetailedPlacer::exportPlacement(Circuit &circuit) {
  placement_.exportPlacement(circuit);
}

void DetailedPlacer::callback() {
  if (!callback_.has_value()) return;
  exportPlacement(circuit_);
  callback_.value()(PlacementStep::Detailed);
}

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
  RowNeighbourhood rowsNeighbours(placement_.rows(), nbRows);
  //  Optimize each row with neighbours after it
  for (int i = 0; i < placement_.nbRows(); ++i) {
    for (int j : rowsNeighbours.rowsAbove(i)) {
      runSwapsTwoRowsAmplify(i, j, nbNeighbours);
    }
    for (int j : rowsNeighbours.rowsRight(i)) {
      runSwapsTwoRowsAmplify(i, j, nbNeighbours);
    }
  }
  // Optimize each row with neighbours before it;
  // We do both for symmetry and to allow large cell movements
  for (int i = placement_.nbRows() - 1; i >= 1; --i) {
    for (int j : rowsNeighbours.rowsBelow(i)) {
      runSwapsTwoRowsAmplify(i, j, nbNeighbours);
    }
    for (int j : rowsNeighbours.rowsLeft(i)) {
      runSwapsTwoRowsAmplify(i, j, nbNeighbours);
    }
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
    while (bestSwapUpdate(c, from, nbNeighbours)) {
      ;
    }
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
  for (int candidate = from, count = 0; candidate != -1 && count < nbNeighbours;
       candidate = placement_.cellNext(candidate), ++count) {
    auto [feasible, val] = valueOnSwap(c, candidate);
    if (feasible && val < bestValue) {
      found = true;
      bestCandidate = candidate;
    }
  }
  for (int candidate = from, count = 0; candidate != -1 && count < nbNeighbours;
       candidate = placement_.cellPred(candidate), ++count) {
    auto [feasible, val] = valueOnSwap(c, candidate);
    if (feasible && val < bestValue) {
      found = true;
      bestCandidate = candidate;
    }
  }
  if (found) {
    doSwap(c, bestCandidate);
    if (bestCandidate == from) {
      from = c;
    }
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
  if (c == -1) {
    return c;
  }
  while (true) {
    int nextC = placement_.cellNext(c);
    if (nextC == -1) {
      break;
    }
    if (placement_.cellX(nextC) > placement_.cellX(target)) {
      break;
    }
    c = nextC;
  }
  return c;
}

int DetailedPlacer::findCellBefore(int target, int fromCell) const {
  int c = fromCell;
  if (c == -1) {
    return c;
  }
  while (true) {
    int nextC = placement_.cellPred(c);
    if (nextC == -1) {
      break;
    }
    if (placement_.cellX(nextC) + placement_.cellWidth(nextC) <
        placement_.cellX(target) + placement_.cellWidth(target)) {
      break;
    }
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
  for (int row1Cell : row1Cells) {
    int x = placement_.cellX(row1Cell);
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

void DetailedPlacer::runShifts(int nbRows, int maxNbCells) {
  if (nbRows < 2) {
    return;
  }
  RowNeighbourhood rowsNeighbours(placement_.rows(), nbRows / 2);

  for (int r = 0; r < placement_.nbRows(); r += nbRows / 2) {
    std::vector<int> rows;
    rows.push_back(r);
    for (int i : rowsNeighbours.rowsBelow(r)) {
      rows.push_back(i);
    }
    for (int i : rowsNeighbours.rowsAbove(r)) {
      rows.push_back(i);
    }
    runShiftsOnRows(rows, maxNbCells);
  }
  placement_.check();
}

void DetailedPlacer::runShiftsOnRows(const std::vector<int> &rows,
                                     int maxNbCells) {
  // Gather all cells and sort them by x coordinate
  std::vector<std::pair<int, int> > sortedCells;
  for (int row : rows) {
    for (int c : placement_.rowCells(row)) {
      sortedCells.emplace_back(placement_.cellX(c), c);
    }
  }
  std::stable_sort(sortedCells.begin(), sortedCells.end());
  std::vector<int> cells;
  cells.reserve(sortedCells.size());

  for (auto p : sortedCells) {
    cells.push_back(p.second);
  }

  // Only select part of the cells so we never solve an optimization problem
  // bigger than maxNbCells
  int overlap = std::min(maxNbCells / 2, 10);
  for (int start = 0; start < cells.size(); start += maxNbCells - overlap) {
    int end = std::min(start + maxNbCells, (int)cells.size());
    std::vector<int> subproblem(cells.begin() + start, cells.begin() + end);
    runShiftsOnCells(subproblem);
  }
}

void DetailedPlacer::runShiftsOnCells(const std::vector<int> &cells) {
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
  using node_pair = std::pair<SmartDigraph::Node, int>;

  // The arcs corresponding to constraints of the original problem
  std::vector<arc_pair> constraint_arcs;

  // Now we add every positional constraint, which becomes an arc in the
  // min-cost flow problem
  for (int c : cells) {
    int pred = placement_.cellPred(c);
    int next = placement_.cellNext(c);
    if (cell_set.count(next) != 0u) {
      // Two movable cells
      auto A = g.addArc(cell_nodes[next], cell_nodes[c]);
      constraint_arcs.emplace_back(A, -placement_.cellWidth(c));
    }
    if (cell_set.count(pred) == 0) {
      // Predecessor fixed
      int boundary = placement_.boundaryBefore(c);
      auto A = g.addArc(cell_nodes[c], fixed);
      constraint_arcs.emplace_back(A, -boundary);
    }
    if (cell_set.count(next) == 0) {
      // Successor fixed
      int boundary = placement_.boundaryAfter(c);
      auto A = g.addArc(fixed, cell_nodes[c]);
      constraint_arcs.emplace_back(A, boundary - placement_.cellWidth(c));
    }
  }

  // One arc for every pin of every net: arcs too
  for (int net : nets) {
    for (int i = 0; i < xtopo_.nbNetPins(net); ++i) {
      int c = xtopo_.pinCell(net, i);
      int pin_offs = xtopo_.netPinOffset(net, i);
      if (cell_set.count(c) != 0u) {
        Arc Al = g.addArc(cell_nodes[c], Lnet_nodes[net]);
        constraint_arcs.emplace_back(Al, pin_offs);
        Arc Ar = g.addArc(Unet_nodes[net], cell_nodes[c]);
        constraint_arcs.emplace_back(Ar, -pin_offs);
      } else {  // Fixed offset
        auto Al = g.addArc(fixed, Lnet_nodes[net]);
        constraint_arcs.emplace_back(Al, xtopo_.cellPos(c) + pin_offs);
        auto Ar = g.addArc(Unet_nodes[net], fixed);
        constraint_arcs.emplace_back(Ar, -xtopo_.cellPos(c) - pin_offs);
      }
    }
  }

  // Then the only capacitated arcs: the ones for the nets
  std::vector<node_pair> net_supplies;
  for (int net : nets) {
    net_supplies.emplace_back(Unet_nodes[net], 1);
    net_supplies.emplace_back(Lnet_nodes[net], -1);
  }

  // Create the maps to have cost and capacity for the arcs
  IntArcMap cost(g, 0);
  IntNodeMap supply(g, 0);

  for (arc_pair A : constraint_arcs) {
    cost[A.first] = A.second;
  }

  for (node_pair N : net_supplies) {
    supply[N.first] = N.second;
  }

  // Then we (hope the solver can) solve it
  NetworkSimplex<SmartDigraph> ns(g);
  ns.supplyMap(supply);
  ns.costMap(cost);

  // May be necessary to specify a capacity for algorithms other than
  // NetworkSimplex IntArcMap capacity(g, nets.size()); ns.upperMap(capacity);

  auto res = ns.run();
  if (res != ns.OPTIMAL) {
    throw std::runtime_error("Could not solve the network flow optimally");
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
