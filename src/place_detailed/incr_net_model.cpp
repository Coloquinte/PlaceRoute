#include "place_detailed/incr_net_model.hpp"

#include <limits>
#include <numeric>
#include <stdexcept>
#include <unordered_map>

namespace coloquinte {

IncrNetModelBuilder::IncrNetModelBuilder(int nbCells) : nbCells_(nbCells) {
  netLimits_.push_back(0);
}

void IncrNetModelBuilder::addNet(const std::vector<int> &cells,
                                 const std::vector<int> &pinOffsets) {
  assert(cells.size() == pinOffsets.size());
  if (cells.size() <= 1) {
    return;
  }
  netLimits_.push_back(netLimits_.back() + cells.size());
  netCells_.insert(netCells_.end(), cells.begin(), cells.end());
  netPinOffsets_.insert(netPinOffsets_.end(), pinOffsets.begin(),
                        pinOffsets.end());
}

IncrNetModel IncrNetModel::xTopology(const Circuit &circuit) {
  std::vector<int> cells;
  cells.reserve(circuit.nbCells());

  for (int c = 0; c < circuit.nbCells(); ++c) {
    cells.push_back(c);
  }
  return xTopology(circuit, cells);
}

IncrNetModel IncrNetModel::yTopology(const Circuit &circuit) {
  std::vector<int> cells;
  cells.reserve(circuit.nbCells());

  for (int c = 0; c < circuit.nbCells(); ++c) {
    cells.push_back(c);
  }
  return yTopology(circuit, cells);
}

namespace {
std::unordered_map<int, int> buildCellMapping(const Circuit &circuit,
                                              const std::vector<int> &cells) {
  std::unordered_map<int, int> cellMap;
  for (size_t i = 0; i < cells.size(); ++i) {
    int c = cells[i];
    assert(c >= 0 && c < circuit.nbCells());
    cellMap[c] = (int) i;
  }
  assert(cellMap.size() == cells.size());
  return cellMap;
}
}  // namespace

IncrNetModel IncrNetModel::xTopology(const Circuit &circuit,
                                     const std::vector<int> &cells) {
  std::unordered_map<int, int> cellMap = buildCellMapping(circuit, cells);

  int fixedCell = cells.size();
  std::vector<int> cellX;
  cellX.reserve(cells.size());

  for (int c : cells) {
    cellX.push_back(circuit.x(c));
  }
  cellX.push_back(0);

  IncrNetModelBuilder ret(cellX.size());
  for (int i = 0; i < circuit.nbNets(); ++i) {
    std::vector<int> cells;
    std::vector<int> offsets;
    bool hasFixed = false;
    int minFixed = std::numeric_limits<int>::max();
    int maxFixed = std::numeric_limits<int>::min();
    for (int j = 0; j < circuit.nbPinsNet(i); ++j) {
      int cell = circuit.pinCell(i, j);
      int offset = circuit.pinXOffset(i, j);
      if (cellMap.count(cell) != 0u) {
        cells.push_back(cellMap[cell]);
        offsets.push_back(offset);
      } else {
        int pos = circuit.x(cell) + offset;
        minFixed = std::min(pos, minFixed);
        maxFixed = std::max(pos, maxFixed);
        hasFixed = true;
      }
    }
    if (hasFixed) {
      cells.push_back(fixedCell);
      offsets.push_back(minFixed);
      if (minFixed != maxFixed) {
        cells.push_back(fixedCell);
        offsets.push_back(maxFixed);
      }
    }
    ret.addNet(cells, offsets);
  }
  return ret.build(cellX);
}

IncrNetModel IncrNetModel::yTopology(const Circuit &circuit,
                                     const std::vector<int> &cells) {
  std::unordered_map<int, int> cellMap = buildCellMapping(circuit, cells);

  int fixedCell = cells.size();
  std::vector<int> cellY;
  cellY.reserve(cells.size());

  for (int c : cells) {
    cellY.push_back(circuit.y(c));
  }
  cellY.push_back(0);

  IncrNetModelBuilder ret(cellY.size());
  for (int i = 0; i < circuit.nbNets(); ++i) {
    std::vector<int> cells;
    std::vector<int> offsets;
    bool hasFixed = false;
    int minFixed = std::numeric_limits<int>::max();
    int maxFixed = std::numeric_limits<int>::min();
    for (int j = 0; j < circuit.nbPinsNet(i); ++j) {
      int cell = circuit.pinCell(i, j);
      int offset = circuit.pinYOffset(i, j);
      if (cellMap.count(cell) != 0u) {
        cells.push_back(cellMap[cell]);
        offsets.push_back(offset);
      } else {
        int pos = circuit.y(cell) + offset;
        minFixed = std::min(pos, minFixed);
        maxFixed = std::max(pos, maxFixed);
        hasFixed = true;
      }
    }
    if (hasFixed) {
      cells.push_back(fixedCell);
      offsets.push_back(minFixed);
      if (minFixed != maxFixed) {
        cells.push_back(fixedCell);
        offsets.push_back(maxFixed);
      }
    }
    ret.addNet(cells, offsets);
  }
  return ret.build(cellY);
}

void IncrNetModel::exportPlacementX(Circuit &circuit) const {
  assert(circuit.nbCells() == nbCells());
  for (int i = 0; i < circuit.nbCells(); ++i) {
    if (!circuit.isFixed(i)) {
      circuit.cellX_[i] = cellPos_[i];
    }
  }
}

void IncrNetModel::exportPlacementY(Circuit &circuit) const {
  assert(circuit.nbCells() == nbCells());
  for (int i = 0; i < circuit.nbCells(); ++i) {
    if (!circuit.isFixed(i)) {
      circuit.cellY_[i] = cellPos_[i];
    }
  }
}

IncrNetModel IncrNetModelBuilder::build() const {
  std::vector<int> pos(nbCells(), 0);
  return build(pos);
}

IncrNetModel IncrNetModelBuilder::build(const std::vector<int> &pos) const {
  assert((int) pos.size() == nbCells());
  IncrNetModel ret;
  ret.netLimits_ = netLimits_;
  ret.netCells_ = netCells_;
  ret.netPinOffsets_ = netPinOffsets_;
  ret.cellPos_ = pos;
  ret.finalize();
  return ret;
}

void IncrNetModel::finalize() {
  cellLimits_.assign(nbCells() + 1, 0);
  cellNets_.assign(nbPins(), -1);
  cellPinOffsets_.assign(nbPins(), 0);
  // Compute the cell limits
  for (int i = 0; i < nbNets(); ++i) {
    for (int j = 0; j < nbNetPins(i); ++j) {
      ++cellLimits_[pinCell(i, j) + 1];
    }
  }
  std::partial_sum(cellLimits_.begin(), cellLimits_.end(), cellLimits_.begin());
  // Setup the nets and offsets
  std::vector<int> curIndex = cellLimits_;
  for (int i = 0; i < nbNets(); ++i) {
    for (int j = 0; j < nbNetPins(i); ++j) {
      int cell = pinCell(i, j);
      int ind = curIndex[cell]++;
      cellNets_[ind] = i;
      cellPinOffsets_[ind] = netPinOffset(i, j);
    }
  }

  // Setup the cost
  netMinMaxPos_ = computeNetMinMaxPos();
  value_ = computeValue();
}

std::pair<int, int> IncrNetModel::computeNetMinMaxPos(int net) const {
  int minPos = std::numeric_limits<int>::max();
  int maxPos = std::numeric_limits<int>::min();
  for (int j = 0; j < nbNetPins(net); ++j) {
    int c = pinCell(net, j);
    int pinPos = cellPos_[c] + netPinOffset(net, j);
    minPos = std::min(pinPos, minPos);
    maxPos = std::max(pinPos, maxPos);
  }
  return std::make_pair(minPos, maxPos);
}

std::vector<std::pair<int, int>> IncrNetModel::computeNetMinMaxPos() const {
  std::vector<std::pair<int, int>> ret(nbNets());
  for (int net = 0; net < nbNets(); ++net) {
    ret[net] = computeNetMinMaxPos(net);
  }
  return ret;
}

long long IncrNetModel::computeValue() const {
  long long ret = 0;
  for (int net = 0; net < nbNets(); ++net) {
    auto minMaxPos = netMinMaxPos_[net];
    ret += minMaxPos.second - minMaxPos.first;
  }
  return ret;
}

void IncrNetModel::updateCellPos(int cell, int pos) {
  // TODO: optimize performance: do not recompute every connected net
  cellPos_[cell] = pos;
  for (int i = 0; i < nbCellPins(cell); ++i) {
    int net = pinNet(cell, i);
    recomputeNet(net);
  }
}

void IncrNetModel::recomputeNet(int net) {
  auto newMinMaxPos = computeNetMinMaxPos(net);
  auto oldMinMaxPos = netMinMaxPos_[net];
  int oldValue = oldMinMaxPos.second - oldMinMaxPos.first;
  int newValue = newMinMaxPos.second - newMinMaxPos.first;
  netMinMaxPos_[net] = newMinMaxPos;
  value_ += newValue - oldValue;
}

void IncrNetModel::check() const {
  if ((int) netLimits_.size() != nbNets() + 1) {
    throw std::runtime_error("Net number mismatch");
  }
  if (netLimits_.front() != 0 || netLimits_.back() != (int) netCells_.size()) {
    throw std::runtime_error("Pin number mismatch");
  }
  if (netCells_.size() != netPinOffsets_.size()) {
    throw std::runtime_error("Pin number mismatch");
  }
  for (int c : netCells_) {
    if (c < -1 || c >= nbCells()) {
      throw std::runtime_error("Invalid cell number");
    }
  }
  for (int i = 0; i < nbNets(); ++i) {
    // At least one cell per net
    if (netLimits_[i] + 1 > netLimits_[i + 1]) {
      throw std::runtime_error("Invalid number of pins in nets");
    }
  }
  if ((int) cellLimits_.size() != nbCells() + 1) {
    throw std::runtime_error("Cell number mismatch");
  }
  if (cellLimits_.front() != 0 || cellLimits_.back() != (int) cellNets_.size()) {
    throw std::runtime_error("Pin number mismatch");
  }
  if (cellLimits_.back() != nbPins()) {
    throw std::runtime_error("Pin number mismatch");
  }
  if ((int) cellNets_.size() != nbPins()) {
    throw std::runtime_error("Pin number mismatch");
  }
  assert((int) cellPinOffsets_.size() == nbPins());
  for (int i = 0; i < nbCells(); ++i) {
    if (cellLimits_[i] > cellLimits_[i + 1]) {
      throw std::runtime_error("Invalid number of pins in cell");
    }
  }
  for (int net : cellNets_) {
    if (net < 0 || net >= nbNets()) {
      throw std::runtime_error("Invalid net number");
    }
  }
  auto minMaxPos = computeNetMinMaxPos();
  for (int i = 0; i < nbNets(); ++i) {
    if (minMaxPos[i] != netMinMaxPos_[i]) {
      throw std::runtime_error("Mismatched net bound");
    }
  }
  if (computeValue() != value()) {
    throw std::runtime_error("Wrong incremental value");
  }
}
}  // namespace coloquinte