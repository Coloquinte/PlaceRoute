#include "place_global/partitioning.hpp"

#include <ofstream>
#include <unordered_set>

namespace coloquinte {

class PartitioningProblem {
 public:
  PartitioningProblem(int nbNodes);

  int nbNodes() const { return nbPrimaryNodes_ + nbDummyNodes_; }

  int nbEdges() const { return edges_.size(); }

  std::vector<int> solve();

  void writeHgr(const std::string &filename) const;

  void writeFixed(const std::string &filename) const;

 private:
  int nbPrimaryNodes_;
  int nbDummyNodes_;
  int nbEdges_;

  std::vector<std::vector<int> > edges_;
  std::vector<long long> partitionWeights_;
};

void PartitioningProblem::writeHgr(const std::string &filename) const {
    std::ofstream f(filename);

}

Partitioner Partitioner::fromIspdCircuit(const Circuit &circuit,
                                         float sizeFactor) {
  return Partitioner(circuit, HierarchicalDensityPlacement::fromIspdCircuit(
                                  circuit, sizeFactor));
}

void Partitioner::reoptimize(const std::vector<std::pair<int, int> > &bins) {
  std::unordered_set<int> cell_set;
  for (auto [x, y] : bins) {
    for (int c : placement_.binCells(x, y)) {
      cell_set.insert(c);
    }
  }
  std::unordered_set<int> net_set;
  // TODO: pre-build the datastructure here
  for (int n = 0; n < circuit_.nbNets(); ++n) {
    for (int i = 0; i < circuit_.nbPinsNet(n); ++i) {
      int c = circuit_.pinCell(n, i);
      if (cell_set.count(c)) {
        net_set.insert(n);
      }
    }
  }
  std::vector<int> cells(cell_set.begin(), cell_set.end());
  std::sort(cells.begin(), cells.end());
  std::vector<int> nets(net_set.begin(), net_set.end());
  std::sort(nets.begin(), nets.end());
}
}  // namespace coloquinte