#include "place_detailed/incr_net_model.hpp"

#include <limits>
#include <numeric>
#include <stdexcept>

namespace coloquinte {

IncrNetModelBuilder::IncrNetModelBuilder(int nbCells) : nbCells_(nbCells) {
  netLimits_.push_back(0);
}

void IncrNetModelBuilder::addNet(const std::vector<int> &cells,
                                 const std::vector<int> &pinOffsets) {
  assert(cells.size() == pinOffsets.size());
  if (cells.size() <= 1) return;
  netLimits_.push_back(netLimits_.back() + cells.size());
  netCells_.insert(netCells_.end(), cells.begin(), cells.end());
  netPinOffsets_.insert(netPinOffsets_.end(), pinOffsets.begin(),
                        pinOffsets.end());
}

IncrNetModel IncrNetModel::xTopology(const Circuit &circuit) {
  IncrNetModelBuilder ret(circuit.nbCells());
  for (int i = 0; i < circuit.nbNets(); ++i) {
    std::vector<int> cells;
    std::vector<int> offsets;
    for (int j = 0; j < circuit.nbPinsNet(i); ++j) {
      int cell = circuit.pinCell(i, j);
      int offset = circuit.pinXOffset(i, j);
      cells.push_back(cell);
      offsets.push_back(offset);
    }
    ret.addNet(cells, offsets);
  }
  return ret.build(circuit.cellX_);
}

IncrNetModel IncrNetModel::yTopology(const Circuit &circuit) {
  IncrNetModelBuilder ret(circuit.nbCells());
  for (int i = 0; i < circuit.nbNets(); ++i) {
    std::vector<int> cells;
    std::vector<int> offsets;
    for (int j = 0; j < circuit.nbPinsNet(i); ++j) {
      int cell = circuit.pinCell(i, j);
      int offset = circuit.pinYOffset(i, j);
      cells.push_back(cell);
      offsets.push_back(offset);
    }
    ret.addNet(cells, offsets);
  }
  return ret.build(circuit.cellY_);
}

void IncrNetModel::exportPlacementX(Circuit &circuit) const {
  assert(circuit.nbCells() == nbCells());
  for (int i = 0; i < circuit.nbCells(); ++i) {
    if (!circuit.fixed(i)) {
      circuit.cellX_[i] = cellPos_[i];
    }
  }
}

void IncrNetModel::exportPlacementY(Circuit &circuit) const {
  assert(circuit.nbCells() == nbCells());
  for (int i = 0; i < circuit.nbCells(); ++i) {
    if (!circuit.fixed(i)) {
      circuit.cellY_[i] = cellPos_[i];
    }
  }
}

IncrNetModel IncrNetModelBuilder::build() const {
  std::vector<int> pos(nbCells(), 0);
  return build(pos);
}

IncrNetModel IncrNetModelBuilder::build(const std::vector<int> &pos) const {
  assert(pos.size() == nbCells());
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
  if (netLimits_.size() != nbNets() + 1) {
    throw std::runtime_error("Net number mismatch");
  }
  if (netLimits_.front() != 0 || netLimits_.back() != netCells_.size()) {
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
  if (cellLimits_.size() != nbCells() + 1) {
    throw std::runtime_error("Cell number mismatch");
  }
  if (cellLimits_.front() != 0 || cellLimits_.back() != cellNets_.size()) {
    throw std::runtime_error("Pin number mismatch");
  }
  if (cellLimits_.back() != nbPins()) {
    throw std::runtime_error("Pin number mismatch");
  }
  if (cellNets_.size() != nbPins()) {
    throw std::runtime_error("Pin number mismatch");
  }
  assert(cellPinOffsets_.size() == nbPins());
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