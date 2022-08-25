#pragma once

#include <cassert>
#include <vector>

#include "coloquinte.hpp"

namespace coloquinte {

class IncrNetModel;

/**
 * @brief Builder to construct IncrNetModel objects
 *
 */
class IncrNetModelBuilder {
 public:
  /**
   * @brief Initialize the datastructure
   */
  explicit IncrNetModelBuilder(int nbCells);

  /**
   * @brief Return the number of cells
   */
  int nbCells() const { return nbCells_; }

  /**
   * @brief Add a new net
   */
  void addNet(const std::vector<int> &cells,
              const std::vector<int> &pinOffsets);

  /**
   * @brief Finalize the object
   */
  IncrNetModel build() const;

  /**
   * @brief Finalize the object with initial positions
   */
  IncrNetModel build(const std::vector<int> &pos) const;

 private:
  int nbCells_;
  std::vector<int> netLimits_;
  std::vector<int> netCells_;
  std::vector<int> netPinOffsets_;
};

/**
 * Incremental computation of the 1D HPWL for detailed placement algorithms
 */
class IncrNetModel {
 public:
  /**
   * @brief Create the horizontal view of the HPWL model
   */
  static IncrNetModel xTopology(const Circuit &);

  /**
   * @brief Create the vertical view of the HPWL model
   */
  static IncrNetModel yTopology(const Circuit &);

  /**
   * @brief Export the horizontal placement to the ISPD circuit
   */
  void exportPlacementX(Circuit &circuit) const;

  /**
   * @brief Export the vertical placement to the ISPD circuit
   */
  void exportPlacementY(Circuit &circuit) const;

  /**
   * @brief Return the number of cells
   */
  int nbCells() const { return cellPos_.size(); }

  /**
   * @brief Return the number of nets
   */
  int nbNets() const { return netLimits_.size() - 1; }

  /**
   * @brief Return the total number of pins
   */
  int nbPins() const { return netLimits_.back(); }

  int cellPos(int cell) const {
    assert(cell >= 0);
    assert(cell < nbCells());
    return cellPos_[cell];
  }

  /**
   * @brief Return the number of pins for a given net
   */
  int nbNetPins(int net) const {
    assert(net >= 0);
    assert(net < nbNets());
    return netLimits_[net + 1] - netLimits_[net];
  }

  /**
   * @brief Return the number of pins for a given net
   */
  int nbCellPins(int cell) const {
    assert(cell >= 0);
    assert(cell < nbCells());
    return cellLimits_[cell + 1] - cellLimits_[cell];
  }

  /**
   * @brief Return the cell for a given net pin
   */
  int pinCell(int net, int pin) const {
    assert(pin >= 0);
    assert(pin < nbNetPins(net));
    return netCells_[netLimits_[net] + pin];
  }

  /**
   * @brief Return the offset for a given net pin
   */
  int netPinOffset(int net, int pin) const {
    assert(pin >= 0);
    assert(pin < nbNetPins(net));
    return netPinOffsets_[netLimits_[net] + pin];
  }

  /**
   * @brief Return the position for a given net pin
   */
  float netPinPosition(int net, int pin) const {
    int c = pinCell(net, pin);
    assert(c >= 0);
    assert(c < nbCells());
    return cellPos_[c] + netPinOffset(net, pin);
  }

  /**
   * @brief Return the net for a given cell pin
   */
  int pinNet(int cell, int pin) const {
    assert(pin >= 0);
    assert(pin < nbCellPins(cell));
    return cellNets_[cellLimits_[cell] + pin];
  }

  /**
   * @brief Return the offset for a given cell pin
   */
  int cellPinOffset(int cell, int pin) const {
    assert(pin >= 0);
    assert(pin < nbCellPins(cell));
    return cellPinOffsets_[cellLimits_[cell] + pin];
  }

  /**
   * @brief Return the position for a given net pin
   */
  float cellPinPosition(int cell, int pin) const {
    assert(cell >= 0);
    assert(cell < nbCells());
    return cellPos_[cell] + cellPinOffset(cell, pin);
  }

  /**
   * @brief Compute the length for a given placement
   */
  long long value() const { return value_; }

  /**
   * @brief Update the cell position and change the value accordingly
   */
  void updateCellPos(int cell, int pos);

  /**
   * @brief Check the consistency of the datastructure
   */
  void check() const;

 private:
  /**
   * @brief Compute the minimum position of a net from scratch
   */
  int computeNetMinPos(int net) const;

  /**
   * @brief Compute the maximum position of a net from scratch
   */
  int computeNetMaxPos(int net) const;

  /**
   * @brief Compute the minimum position of each net from scratch
   */
  std::vector<int> computeNetMinPos() const;

  /**
   * @brief Compute the maximum position of each net from scratch
   */
  std::vector<int> computeNetMaxPos() const;

  /**
   * @brief Compute the total value from scratch
   */
  long long computeValue() const;

  /**
   * @brief Recompute a single net
   */
  void recomputeNet(int net);

  /**
   * @brief Initialize the datastructure; private, you should use
   * IncrNetModelBuilder
   */
  IncrNetModel() {}

  /**
   * @brief Finalize the datastructure during build
   */
  void finalize();

 private:
  // Positions of the cells
  std::vector<int> cellPos_;

  // CSR representation from the nets
  std::vector<int> netLimits_;
  std::vector<int> netCells_;
  std::vector<int> netPinOffsets_;

  // CSR representation from the cells
  std::vector<int> cellLimits_;
  std::vector<int> cellNets_;
  std::vector<int> cellPinOffsets_;

  // Current minimum and maximum for each net
  std::vector<int> netMinPos_;
  std::vector<int> netMaxPos_;

  // Current value
  long long value_;

  friend class IncrNetModelBuilder;
};
}  // namespace coloquinte
