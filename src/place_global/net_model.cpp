#include "place_global/net_model.hpp"

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/IterativeLinearSolvers>
#include <eigen3/Eigen/Sparse>
#include <limits>

NetModel::NetModel(int nbCells) : nbCells_(nbCells) { netLimits_.push_back(0); }

void NetModel::addNet(const std::vector<int> &cells) {
  std::vector<float> offsets(cells.size(), 0.0f);
  addNet(cells, offsets);
}

void NetModel::addNet(const std::vector<int> &cells,
                      const std::vector<float> &pinOffsets) {
  addNet(cells, pinOffsets, std::numeric_limits<float>::infinity(),
         -std::numeric_limits<float>::infinity(), 1.0f);
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