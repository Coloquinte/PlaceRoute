#include "place_global/net_model.hpp"

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/IterativeLinearSolvers>
#include <eigen3/Eigen/Sparse>
#include <limits>

NetModel NetModel::xTopology(const Circuit &circuit) {
  return fromData(circuit.cellWidths, circuit.cellX, circuit.cellFixed,
                  circuit.netLimits, circuit.pinCells, circuit.pinXOffsets);
}

NetModel NetModel::yTopology(const Circuit &circuit) {
  return fromData(circuit.cellHeights, circuit.cellY, circuit.cellFixed,
                  circuit.netLimits, circuit.pinCells, circuit.pinYOffsets);
}

NetModel NetModel::fromData(const std::vector<int> &cellSizes,
                            const std::vector<int> &pl,
                            const std::vector<char> &cellFixed,
                            const std::vector<int> &netLimits,
                            const std::vector<int> &pinCells,
                            const std::vector<int> &pinOffsets) {
  NetModel ret(cellSizes.size());
  for (int i = 0; i + 1 < netLimits.size(); ++i) {
    float minPos = std::numeric_limits<float>::infinity();
    float maxPos = -std::numeric_limits<float>::infinity();
    int b = netLimits[i];
    int e = netLimits[i + 1];
    std::vector<int> cells;
    std::vector<float> offsets;
    for (int j = b; j < e; ++j) {
      int cell = pinCells[j];
      int offset = pinOffsets[j];
      if (cellFixed[cell]) {
        int pos = pl[cell] + offset;
        minPos = std::min(minPos, (float)pos);
        maxPos = std::max(maxPos, (float)pos);
      } else {
        cells.push_back(cell);
        // Offset to center of cell instead of lower-left
        offsets.push_back(offset - 0.5f * cellSizes[cell]);
      }
    }
    ret.addNet(cells, offsets, minPos, maxPos);
  }
  ret.check();
  return ret;
}

NetModel::NetModel(int nbCells) : nbCells_(nbCells) { netLimits_.push_back(0); }

void NetModel::addNet(const std::vector<int> &cells) {
  std::vector<float> offsets(cells.size(), 0.0f);
  addNet(cells, offsets);
}

void NetModel::addNet(const std::vector<int> &cells,
                      const std::vector<float> &pinOffsets) {
  addNet(cells, pinOffsets, std::numeric_limits<float>::infinity(),
         -std::numeric_limits<float>::infinity());
}

void NetModel::addNet(const std::vector<int> &cells,
                      const std::vector<float> &pinOffsets, float minPin,
                      float maxPin, float weight) {
  assert(cells.size() == pinOffsets.size());
  if (cells.size() <= 1) {
    return;
  }
  netLimits_.push_back(cells.size());
  netWeight_.push_back(weight);
  netCells_.insert(netCells_.end(), cells.begin(), cells.end());
  netPinOffsets_.insert(netPinOffsets_.end(), pinOffsets.begin(),
                        pinOffsets.end());
  netMinPin_.push_back(minPin);
  netMaxPin_.push_back(maxPin);
}

void NetModel::check() const {
  assert(netCells_.size() == pinOffsets_.size());
  assert(netLimits_.size() == nbNets() + 1);
  assert(netWeight_.size() == nbNets());
  assert(netMinPin_.size() == nbNets());
  assert(netMaxPin_.size() == nbNets());
  for (int c : netCells_) {
    assert(c >= 0);
    assert(c < nbCells());
  }
  for (int i = 0; i < nbNets(); ++i) {
    assert(netLimits_[i] < netLimits_[i + 1]);
  }
  assert(netLimits_.front() == 0);
  assert(netLimits.back() == netCells_.size());
}

std::vector<float> NetModel::pinPositions(const std::vector<float> &pl) const {
  std::vector<float> ret(netCells_.size());
  for (int i = 0; i < netCells_.size(); ++i) {
    ret[i] = pl[netCells_[i]] + netPinOffsets_[i];
  }
  return ret;
}