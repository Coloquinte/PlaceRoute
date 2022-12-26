#include "place_global/partitioning.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <unordered_set>

namespace coloquinte {

class PartitioningProblem {
 public:
  PartitioningProblem(int nbNodes, int nbPartitions, double minMargin = 0.01,
                      double maxMargin = 0.1, bool weightedEdges = false);

  int nbNodes() const { return nodeWeights_.size(); }

  int nbPartitions() const { return partitionWeights_.size(); }

  int nbEdges() const { return edges_.size(); }

  long long totalNodeWeight() const;

  long long totalPartWeight() const;

  std::vector<long long> scaledPartWeight() const;

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
  double minMargin_;
  double maxMargin_;
  bool weightedEdges_;
};

PartitioningProblem::PartitioningProblem(int nbNodes, int nbPartitions,
                                         double minMargin, double maxMargin,
                                         bool weightedEdges) {
  nodeWeights_.resize(nbNodes);
  partitionWeights_.resize(nbPartitions);
  nodeWeights_.resize(nbNodes);
  partitionWeights_.resize(nbPartitions);
  minMargin_ = minMargin;
  maxMargin_ = maxMargin;
  weightedEdges_ = weightedEdges;
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

std::vector<long long> PartitioningProblem::scaledPartWeight() const {
  double totNode = totalNodeWeight();
  double totPart = totalPartWeight();
  double targetAmount = std::clamp(totPart, totNode * (1.0 + minMargin_),
                                   totNode * (1.0 + maxMargin_));
  double scale = targetAmount / totPart;
  std::vector<long long> ret;
  for (long long p : partitionWeights_) {
    ret.push_back((long long)std::round(p * scale));
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
  cmd << "/home/alnurn/Documents/PlaceRoute/thirdparty/kahypar/kahypar/build/"
         "kahypar/application/KaHyPar ";
  cmd << "-p "
         "/home/alnurn/Documents/PlaceRoute/thirdparty/kahypar/kahypar/config/"
         "km1_kKaHyPar_sea20.ini ";
  cmd << "-h coloquinte_hypergraph.hgr -f coloquinte_hypergraph.fixed ";
  cmd << "-o km1 -m direct -e 0.05 -w on ";
  cmd << "-k " << nbPartitions() << " ";
  cmd << "--use-individual-part-weights on --part-weights ";
  for (auto w : scaledPartWeight()) {
    cmd << w << " ";
  }
  std::stringstream outname;
  outname << "coloquinte_hypergraph.hgr.part";
  outname << nbPartitions() << ".epsilon0.05.seed-1.KaHyPar";
  std::cout << "Command: " << cmd.str() << std::endl;
  int res = system(cmd.str().c_str());
  if (res != 0) {
    throw std::runtime_error("Command returned an error value");
  }
  return readSolution(outname.str());
}

void PartitioningProblem::writeHgr(const std::string &filename) const {
  std::ofstream f(filename);
  f << nbEdges() << " " << nbNodes() + nbPartitions() << " ";
  if (weightedEdges_) {
    f << "11\n";
  } else {
    f << "10\n";
  }
  for (const auto &edge : edges_) {
    if (weightedEdges_) {
      f << "1 ";  // Weight
    }
    for (int i = 0; i < edge.size(); ++i) {
      if (i != 0) f << " ";
      f << edge[i] + 1;
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
  std::vector<int> ret;
  std::ifstream f(filename);
  for (int i = 0; i < nbNodes(); ++i) {
    int p;
    f >> p;
    ret.push_back(p);
  }
  return ret;
}

Partitioner Partitioner::fromIspdCircuit(const Circuit &circuit,
                                         float sizeFactor) {
  return Partitioner(circuit, HierarchicalDensityPlacement::fromIspdCircuit(
                                  circuit, sizeFactor));
}

std::vector<int> projectPins(const std::vector<std::pair<int, int> > &bins,
                             int x, int y) {
  int minDist = std::numeric_limits<int>::max();
  for (auto [bx, by] : bins) {
    int dist = std::abs(bx - x) + std::abs(by - y);
    if (dist < minDist) {
      minDist = dist;
    }
  }
  std::vector<int> ret;
  for (int b = 0; b < bins.size(); ++b) {
    auto [bx, by] = bins[b];
    int dist = std::abs(bx - x) + std::abs(by - y);
    if (dist == minDist) {
      ret.emplace_back(b);
    }
  }
  return ret;
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
    binCells_[x][y].clear();
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
      } else {
        int x;
        int y;
        if (circuit_.isFixed(c)) {
          // Fixed pins
          int xpos = circuit_.x(c) + circuit_.pinXOffset(n, i);
          int ypos = circuit_.y(c) + circuit_.pinYOffset(n, i);
          x = findBinByX(xpos);
          y = findBinByY(ypos);
        } else {
          // Cell outside of the current partition
          x = cellBinX(c);
          y = cellBinY(c);
        }
        std::vector<int> newFixed = projectPins(actualBins, x, y);
        for (int p : newFixed) {
          fixedPins.insert(p);
        }
      }
    }

    std::vector<int> pinsV(pins.begin(), pins.end());
    std::vector<int> fixedPinsV(fixedPins.begin(), fixedPins.end());

    prob.addEdge(pinsV, fixedPinsV);
  }
  std::vector<int> assignment = prob.solve();
  setBinCells(actualBins, cells, assignment);
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

void Partitioner::run() {
  while (levelX() > 0 || levelY() > 0) {
    refine();
    improve();
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

void Partitioner::improve() {
  improveXNeighbours();
  improveYNeighbours();
  improveXNeighbours(false);
  improveYNeighbours(false);
}
}  // namespace coloquinte