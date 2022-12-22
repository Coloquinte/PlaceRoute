#include "place_global/partitioning.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <unordered_set>

namespace coloquinte {

class PartitioningProblem {
 public:
  PartitioningProblem(int nbNodes, int nbPartitions);

  int nbNodes() const { return nodeWeights_.size(); }

  int nbPartitions() const { return partitionWeights_.size(); }

  int nbEdges() const { return edges_.size(); }

  long long totalNodeWeight() const;

  long long totalPartWeight() const;

  void setNodeWeight(int i, long long w) { nodeWeights_[i] = w; }

  void setPartWeight(int i, long long w) { partitionWeights_[i] = w; }

  void addEdge(const std::vector<int> &nodes, const std::vector<int> &fixed);

  std::vector<int> solve();

  void writeHgr(const std::string &filename) const;

  void writeFixed(const std::string &filename) const;

  std::vector<int> readSolution(const std::string &filename);

 private:
  std::vector<std::vector<int> > edges_;
  std::vector<long long> nodeWeights_;
  std::vector<long long> partitionWeights_;
};

PartitioningProblem::PartitioningProblem(int nbNodes, int nbPartitions) {
  nodeWeights_.resize(nbNodes);
  partitionWeights_.resize(nbPartitions);
  nodeWeights_.resize(nbNodes);
  partitionWeights_.resize(nbPartitions);
}

long long PartitioningProblem::totalNodeWeight() const {
  long long ret = 0LL;
  for (long long w : nodeWeights_) {
    ret += w;
  }
  return ret;
}

long long PartitioningProblem::totalPartWeight() const {
  long long ret = 0LL;
  for (long long w : partitionWeights_) {
    ret += w;
  }
  return ret;
}

void PartitioningProblem::addEdge(const std::vector<int> &nodes,
                                  const std::vector<int> &fixed) {
  std::vector<int> uniqueNodes = nodes;
  for (int p : fixed) {
    uniqueNodes.push_back(nbNodes() + p);
  }
  std::sort(uniqueNodes.begin(), uniqueNodes.end());
  uniqueNodes.erase(std::unique(uniqueNodes.begin(), uniqueNodes.end()),
                    uniqueNodes.end());
  if (uniqueNodes.size() <= 1) {
    return;
  }
  edges_.push_back(uniqueNodes);
}

std::vector<int> PartitioningProblem::solve() {
  writeHgr("coloquinte_hypergraph.hgr");
  writeFixed("coloquinte_hypergraph.fixed");
  std::stringstream cmd;
  cmd << "kahypar -h coloquinte_hypergraph.hgr -f coloquinte_hypergraph.fixed ";
  cmd << "-o km1 -m direct -epsilon 0.05 ";
  cmd << "-k " << nbPartitions() << " ";
  cmd << "--use-individual-part-weights --part-weights ";
  for (auto w : partitionWeights_) {
    cmd << w << " ";
  }
  std::cout << "Command: " << cmd.str() << std::endl;
  return readSolution("coloquinte_partitioning.sol");
}

void PartitioningProblem::writeHgr(const std::string &filename) const {
  std::ofstream f(filename);
  f << nbEdges() << " " << nbNodes() + nbPartitions() << " 11\n";
  for (const auto &edge : edges_) {
    f << "1 ";  // Weight
    for (int i = 0; i < edge.size(); ++i) {
      f << " " << edge[i] + 1;
    }
    f << "\n";
  }
  for (int i = 0; i < nbNodes(); ++i) {
    f << nodeWeights_[i] << "\n";
  }
  for (int i = 0; i < nbPartitions(); ++i) {
    // TODO: handle fixed weight
    f << "1\n";
  }
}

void PartitioningProblem::writeFixed(const std::string &filename) const {
  std::ofstream f(filename);
  for (int i = 0; i < nbNodes(); ++i) {
    f << "-1\n";
  }
  for (int i = 0; i < nbPartitions(); ++i) {
    f << i << "\n";
  }
}

std::vector<int> PartitioningProblem::readSolution(
    const std::string &filename) {
  // TODO
  return std::vector<int>();
}

Partitioner Partitioner::fromIspdCircuit(const Circuit &circuit,
                                         float sizeFactor) {
  return Partitioner(circuit, HierarchicalDensityPlacement::fromIspdCircuit(
                                  circuit, sizeFactor));
}

void Partitioner::reoptimize(const std::vector<std::pair<int, int> > &bins) {
  std::vector<std::pair<int, int> > actualBins;
  std::unordered_set<int> cell_set;
  for (auto [x, y] : bins) {
    if (binCapacity(x, y) > 0) {
      actualBins.emplace_back(x, y);
    }
    for (int c : binCells(x, y)) {
      cell_set.insert(c);
    }
  }
  std::unordered_set<int> net_set;
  // TODO: pre-build the reverse datastructure for the circuit
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

  std::unordered_map<int, int> cellToId;
  for (int i = 0; i < cells.size(); ++i) {
    cellToId[cells[i]] = i;
  }

  PartitioningProblem prob(cells.size(), actualBins.size());
  for (int i = 0; i < actualBins.size(); ++i) {
    auto [x, y] = actualBins[i];
    prob.setPartWeight(i, binCapacity(x, y));
  }
  for (int i = 0; i < cells.size(); ++i) {
    prob.setNodeWeight(i, cellDemand(cells[i]));
  }
  for (int n : net_set) {
    std::unordered_set<int> pins;
    std::unordered_set<int> fixedPins;

    for (int i = 0; i < circuit_.nbPinsNet(n); ++i) {
      int c = circuit_.pinCell(n, i);
      if (cell_set.count(c)) {
        pins.insert(cellToId[c]);
      } else if (circuit_.isFixed(c)) {
        // Fixed pins
      } else {
        // Cell outside of the current partition
      }
    }

    std::vector<int> pinsV(pins.begin(), pins.end());
    std::vector<int> fixedPinsV(fixedPins.begin(), fixedPins.end());

    prob.addEdge(pinsV, fixedPinsV);
  }
  prob.solve();
  exit(0);
}

void Partitioner::improveXNeighbours(bool sameParent) {
  for (int i = 0; i + 1 < nbBinsX(); ++i) {
    if ((parentX(i) == parentX(i + 1)) != sameParent) {
      continue;
    }
    for (int j = 0; j < nbBinsY(); ++j) {
      reoptimize({{i, j}, {i + 1, j}});
    }
  }
}

void Partitioner::improveYNeighbours(bool sameParent) {
  for (int j = 0; j + 1 < nbBinsY(); ++j) {
    if ((parentY(j) == parentY(j + 1)) != sameParent) {
      continue;
    }
    for (int i = 0; i < nbBinsX(); ++i) {
      reoptimize({{i, j}, {i, j + 1}});
    }
  }
}

void Partitioner::refine() {
  bool doX = levelX() >= levelY();
  bool doY = levelY() >= levelX();

  if (doX) {
    refineX();
    improveXNeighbours();
    improveXNeighbours(false);
  }
  if (doY) {
    refineY();
    improveYNeighbours();
    improveYNeighbours(false);
  }
}
}  // namespace coloquinte