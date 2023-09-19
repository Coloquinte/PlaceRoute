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

void DetailedPlacer::legalize(
    Circuit &circuit, const ColoquinteParameters &params,
    const std::optional<PlacementCallback> &callback) {
  params.check();
  std::cout << "Legalization starting (WL " << circuit.hpwl() << ")"
            << std::endl;
  auto startTime = std::chrono::steady_clock::now();
  Legalizer leg = Legalizer::fromIspdCircuit(circuit);
  leg.run(params);
  auto endTime = std::chrono::steady_clock::now();
  std::chrono::duration<float> duration = endTime - startTime;
  float distance = leg.meanDistance(params.legalization.costModel);
  leg.exportPlacement(circuit);
  std::cout << "Legalization done (WL " << circuit.hpwl();
  std::cout << std::fixed << std::setprecision(1) << ", dist " << distance;
  std::cout << std::fixed << std::setprecision(2) << ") in " << duration.count()
            << "s" << std::endl;
  if (callback.has_value()) {
    callback.value()(PlacementStep::Detailed);
  }
}

void DetailedPlacer::place(Circuit &circuit, const ColoquinteParameters &params,
                           const std::optional<PlacementCallback> &callback) {
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
                               const ColoquinteParameters &params)
    : placement_(DetailedPlacement::fromIspdCircuit(circuit)),
      xtopo_(IncrNetModel::xTopology(circuit)),
      ytopo_(IncrNetModel::yTopology(circuit)),
      params_(params),
      circuit_(circuit) {}

void DetailedPlacer::run() {
  for (int i = 1; i <= params_.detailed.nbPasses; ++i) {
    std::cout << "#" << i << ":" << std::flush;
    runSwaps(params_.detailed.localSearchNbNeighbours,
             params_.detailed.localSearchNbRows);
    auto swapValue = value();
    std::cout << "\tSwaps " << value() << std::flush;
    callback();
    if (params_.detailed.shiftMaxNbCells >= 2) {
      runShifts(params_.detailed.shiftNbRows, params_.detailed.shiftMaxNbCells);
      auto shiftValue = value();
      std::cout << "\tShifts " << shiftValue << std::flush;
      callback();
    }
    if (params_.detailed.reorderingMaxNbCells >= 2) {
      runReordering(params_.detailed.reorderingNbRows,
                    params_.detailed.reorderingMaxNbCells);
      auto reordValue = value();
      std::cout << "\tReordering " << reordValue << std::flush;
      callback();
    }
    std::cout << std::endl;
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
  std::vector<int> cells = placement_.rowCells(rows);

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

struct ReorderingRegion {
  int row;
  int minPos;
  int maxPos;
  int cellPred;
  int cellNext;

  int width() const { return maxPos - minPos; }
};

class RowReordering {
 public:
  /**
   * @brief Initialize the datastructure
   */
  RowReordering(DetailedPlacement &placement, IncrNetModel &xtopo,
                IncrNetModel &ytopo);

  /**
   * @brief Number of registered regions
   */
  size_t nbRegions() const { return regions_.size(); }

  /**
   * @brief Number of registered cells
   */
  size_t nbCells() const { return cells_.size(); }

  /**
   * @brief Add a region and its cells to be optimized
   */
  void addRow(int row, int cellPred, int cellNext);

  /**
   * @brief Add cells to be optimized at once
   */
  void addCells(const std::vector<int> &cells);

  /**
   * @brief Total cell width allocated to a given region
   */
  int allocatedWidth(int region) const;

  /**
   * @brief Run the full optimization
   */
  void run();
  void runRegionChoice(int cellInd);
  void runOrdering(int rowInd);
  long long value() const { return xtopo_.value() + ytopo_.value(); }

  void writeback();

  /**
   * @brief Check the consistency of the datastructure
   */
  void check() const;

 private:
  DetailedPlacement &placement_;
  IncrNetModel &xtopo_;
  IncrNetModel &ytopo_;
  std::vector<ReorderingRegion> regions_;
  std::vector<int> cells_;

  std::vector<std::vector<int> > order_;
  std::vector<std::vector<int> > positions_;

  long long bestVal_;
  std::vector<std::vector<int> > bestOrder_;
  std::vector<std::vector<int> > bestPositions_;
  bool improvement_;
};

RowReordering::RowReordering(DetailedPlacement &placement, IncrNetModel &xtopo,
                             IncrNetModel &ytopo)
    : placement_(placement), xtopo_(xtopo), ytopo_(ytopo) {
  bestVal_ = std::numeric_limits<long long>::max();
  improvement_ = false;
}

void RowReordering::addRow(int row, int cellPred, int cellNext) {
  ReorderingRegion newRow;
  newRow.row = row;
  newRow.cellPred = cellPred;
  newRow.cellNext = cellNext;
  assert(cellPred != cellNext);
  assert(cellPred == -1 || placement_.cellRow(cellPred) == row);
  assert(cellNext == -1 || placement_.cellRow(cellNext) == row);
  newRow.minPos = placement_.boundaryAfter(row, cellPred);
  newRow.maxPos = placement_.boundaryBefore(row, cellNext);
  regions_.push_back(newRow);
  std::vector<int> cells = placement_.cellsBetween(row, cellPred, cellNext);
  // Register cells
  for (int c : cells) {
    cells_.push_back(c);
  }
  // Add to best known order and position
  bestOrder_.push_back(cells);
  std::vector<int> positions;
  for (int c : cells) {
    positions.push_back(placement_.cellX(c));
  }
  bestPositions_.push_back(positions);
  order_.emplace_back();
  positions_.emplace_back();
}

void RowReordering::addCells(const std::vector<int> &cells) {
  std::unordered_set<int> cell_set(cells.begin(), cells.end());
  for (int c : cells) {
    if (cell_set.count(placement_.cellPred(c))) {
      // Not the start of a range of cells
      continue;
    }
    int cellPred = placement_.cellPred(c);
    int cellNext = c;
    while (cell_set.count(cellNext)) {
      cellNext = placement_.cellNext(cellNext);
    }
    addRow(placement_.cellRow(c), cellPred, cellNext);
  }
}

int RowReordering::allocatedWidth(int region) const {
  int ret = 0;
  for (int c : order_[region]) {
    ret += placement_.cellWidth(c);
  }
  return ret;
}

void RowReordering::run() {
  bestVal_ = xtopo_.value() + ytopo_.value();
  // Reverse-sorted to interact well with next_permutation
  std::sort(cells_.begin(), cells_.end(), std::greater<int>());
  runRegionChoice(cells_.size() - 1);
  writeback();
}

void RowReordering::runRegionChoice(int cellInd) {
  assert(cellInd < nbCells());
  if (cellInd < 0) {
    // Leaf case: run ordering optimization
    runOrdering(nbRegions() - 1);
  } else {
    for (int i = 0; i < nbRegions(); ++i) {
      order_[i].push_back(cells_[cellInd]);
      if (allocatedWidth(i) <= regions_[i].width()) {
        // Only if there is enough space left in the row
        ytopo_.updateCellPos(cells_[cellInd], placement_.rowY(regions_[i].row));
        runRegionChoice(cellInd - 1);
      }
      // Assumes cells are reverse-sorted: order_ is always sorted after
      // next_permutation calls
      assert(order_[i].back() == cells_[cellInd]);
      order_[i].pop_back();
    }
  }
}

void RowReordering::runOrdering(int rowInd) {
  if (rowInd < 0) {
    // Leaf case: evaluate
    long long value = xtopo_.value() + ytopo_.value();
    if (value < bestVal_) {
      bestVal_ = value;
      bestOrder_ = order_;
      bestPositions_ = positions_;
      improvement_ = true;
    }
  } else {
    // Iterate on all possible orderings
    while (
        std::next_permutation(order_[rowInd].begin(), order_[rowInd].end())) {
      // Setup the positions
      positions_[rowInd].clear();
      int predPos = regions_[rowInd].minPos;
      for (int c : order_[rowInd]) {
        positions_[rowInd].push_back(predPos);
        xtopo_.updateCellPos(c, predPos);
        predPos += placement_.cellWidth(c);
      }
      runOrdering(rowInd - 1);
    }
  }
}

void RowReordering::writeback() {
  if (improvement_) {
    // Save the new placement
    for (int c : cells_) {
      placement_.unplace(c);
    }
    for (int i = 0; i < nbRegions(); ++i) {
      int pred = regions_[i].cellPred;
      for (int j = 0; j < bestOrder_[i].size(); ++j) {
        int c = bestOrder_[i][j];
        placement_.place(c, regions_[i].row, pred, bestPositions_[i][j]);
        xtopo_.updateCellPos(c, placement_.cellX(c));
        ytopo_.updateCellPos(c, placement_.cellY(c));
        pred = c;
      }
    }
  } else {
    // Placement order is not modified, but we need to reset cell positions in
    // the net topologies
    for (int c : cells_) {
      xtopo_.updateCellPos(c, placement_.cellX(c));
      ytopo_.updateCellPos(c, placement_.cellY(c));
    }
  }
}

void RowReordering::check() const {
  if (nbRegions() != order_.size())
    throw std::runtime_error("Inconsistent number of regions");
  if (nbRegions() != positions_.size())
    throw std::runtime_error("Inconsistent number of regions");
  if (nbRegions() != bestOrder_.size())
    throw std::runtime_error("Inconsistent number of regions");
  if (nbRegions() != bestPositions_.size())
    throw std::runtime_error("Inconsistent number of regions");
  for (int i = 0; i < nbRegions(); ++i) {
    if (bestOrder_[i].size() != bestPositions_[i].size())
      throw std::runtime_error("Inconsistent number of cells in regions");
  }
  std::unordered_set<int> cell_set(cells_.begin(), cells_.end());
  if (cell_set.size() != cells_.size()) {
    throw std::runtime_error("The given cells are not unique");
  }
  for (const ReorderingRegion &region : regions_) {
    if (cell_set.count(region.cellPred)) {
      throw std::runtime_error("Some of the given regions are adjacent");
    }
    if (cell_set.count(region.cellNext)) {
      throw std::runtime_error("Some of the given regions are adjacent");
    }
  }
}

void DetailedPlacer::runReordering(int maxNbRows, int maxNbCells) {
  if (maxNbCells < 2) {
    return;
  }

  RowNeighbourhood rowsNeighbours(placement_.rows(), maxNbRows - 1);

  for (int row = 0; row < placement_.nbRows(); ++row) {
    std::vector<int> rows = {row};
    for (int r : rowsNeighbours.rowsAbove(row)) {
      rows.push_back(r);
    }
    runReorderingOnRows(rows, maxNbCells);
  }
  check();
}

void DetailedPlacer::runReorderingOnRows(const std::vector<int> &rows,
                                         int maxNbCells) {
  std::vector<int> cells = placement_.rowCells(rows);

  // Only select part of the cells so we never solve an optimization problem
  // bigger than maxNbCells
  int overlap = std::min(maxNbCells / 2, 10);
  for (int start = 0; start < cells.size(); start += maxNbCells - overlap) {
    int end = std::min(start + maxNbCells, (int)cells.size());
    std::vector<int> subproblem(cells.begin() + start, cells.begin() + end);
    runReorderingOnCells(subproblem);
  }
}

void DetailedPlacer::runReorderingOnCells(const std::vector<int> &cells) {
  RowReordering reord(placement_, xtopo_, ytopo_);
  reord.addCells(cells);
  reord.run();
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
