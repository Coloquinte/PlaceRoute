#pragma once

#include <cassert>
#include <vector>

#include "coloquinte.hpp"

/**
 * Representation of the nets as 1D HPWL for optimization algorithms
 */
class NetModel {
 public:
  /**
   * @brief Create the horizontal view of the HPWL model
   */
  static NetModel xTopology(const Circuit &);

  /**
   * @brief Create the vertical view of the HPWL model
   */
  static NetModel yTopology(const Circuit &);

  /**
   * @brief Return the number of cells
   */
  int nbCells() const { return nbCells_; }

  /**
   * @brief Return the number of nets
   */
  int nbNets() const { return netLimits_.size() - 1; }

  /**
   * @brief Return the number of pins for a given net
   */
  int nbPins(int net) const {
    assert(net >= 0);
    assert(net < nbNets());
    return netLimits_[net + 1] - netLimits_[net];
  }

  /**
   * @brief Return the cell for a given net pin
   */
  int pinCell(int net, int pin) const {
    assert(pin >= 0);
    assert(pin < nbPins(net));
    return netCells_[netLimits_[net] + pin];
  }

  /**
   * @brief Return the offset for a given net pin
   */
  float pinOffset(int net, int pin) const {
    assert(pin >= 0);
    assert(pin < nbPins(net));
    return netPinOffsets_[netLimits_[net] + pin];
  }

  /**
   * @brief Return the position for a given net pin with this placement
   */
  float pinPosition(int net, int pin, const std::vector<float> &pl) const {
    int c = pinCell(net, pin);
    assert(c >= -1);
    assert(c < nbCells());
    float pos = c == -1 ? 0.0f : pl[c];
    return pos + pinOffset(net, pin);
  }

  /**
   * @brief Return the weight associated with a given net
   */
  float netWeight(int net) const {
    assert(net < nbNets());
    return netWeight_[net];
  }

  /**
   * @brief Initialize the datastructure
   */
  NetModel(int nbCells);

  /**
   * @brief Add a new net
   */
  void addNet(const std::vector<int> &cells);

  /**
   * @brief Add a new net
   */
  void addNet(const std::vector<int> &cells,
              const std::vector<float> &pinOffsets, float weight = 1.0f);

  /**
   * @brief Add a new net
   */
  void addNet(const std::vector<int> &cells,
              const std::vector<float> &pinOffsets, float minPin, float maxPin,
              float weight = 1.0f);

  /**
   * @brief Compute the length for a given placement
   */
  float value(const std::vector<float> &pl) const;

  /**
   * @brief Solve with a quadratic star model
   *
   * @return Positions for the cells that minimize the wirelength
   */
  std::vector<float> solveStar() const;

  /**
   * @brief Solve using the star model as a local approximation of the actual
   * wirelength
   *
   * @param placement Initial positions to build the approximation
   * @param epsilon Small approximation parameter
   *
   * @return New positions for the cells
   */
  std::vector<float> solveStar(const std::vector<float> &placement,
                               float epsilon) const;

  /**
   * @brief Solve using the star model as a local approximation of the actual
   * wirelength
   *
   * @param netPlacement Initial positions to build the approximation
   * @param epsilon Small approximation parameter
   * @param placementTarget Target positions to pull the cells
   * @param penaltyStrength Magnitude of the penalty towards the target
   * positions
   * @param cutoffDistance Cutoff distance where the strength stops growing
   *
   * @return New positions for the cells
   */
  std::vector<float> solveStar(const std::vector<float> &netPlacement,
                               float epsilon,
                               const std::vector<float> &placementTarget,
                               const std::vector<float> &penaltyStrength,
                               float cutoffDistance) const;

  /**
   * @brief Solve using the bound-to-bound model as a local approximation of the
   * actual wirelength
   *
   * @param placement Initial positions to build the approximation
   * @param epsilon Small approximation parameter
   *
   * @return New positions for the cells
   */
  std::vector<float> solveB2B(const std::vector<float> &placement,
                              float epsilon) const;

  /**
   * @brief Solve the bound-to-bound model as a local approximation of the
   * actual wirelength
   *
   * @param netPlacement Initial positions to build the approximation
   * @param epsilon Small approximation parameter
   * @param placementTarget Target positions to pull the cells
   * @param penaltyStrength Magnitude of the penalty towards the target
   * positions
   * @param cutoffDistance Cutoff distance where the strength stops growing
   *
   * @return New positions for the cells
   */

  std::vector<float> solveB2B(const std::vector<float> &netPlacement,
                              float epsilon,
                              const std::vector<float> &placementTarget,
                              const std::vector<float> &penaltyStrength,
                              float cutoffDistance) const;

  /**
   * @brief Check the consistency of the datastructure
   */
  void check() const;

 private:
  /**
   * @brief Obtain the positions of the pins with a given placement (corresponds
   * to netCells_ and pinOffsets_)
   *
   */
  std::vector<float> pinPositions(const std::vector<float> &pl) const;

  /**
   * @brief Return the pin index, cell, offset and position associated with the
   * leftmost pin for the net. Cell may be -1 for a fixed pin
   */
  std::tuple<int, int, float, float> minPin(int net,
                                            const std::vector<float> &pl) const;

  /**
   * @brief eturn the pin index, cell, offset and position associated with the
   * rightmost pin for the net. Cell may be -1 for a fixed pin
   */
  std::tuple<int, int, float, float> maxPin(int net,
                                            const std::vector<float> &pl) const;

  /**
   * @brief Helper function to create the model from ISPD data
   */
  static NetModel fromData(const std::vector<int> &cellSizes,
                           const std::vector<int> &pl,
                           const std::vector<char> &cellFixed,
                           const std::vector<int> &netLimits,
                           const std::vector<int> &pinCells,
                           const std::vector<int> &pinOffsets);

 private:
  int nbCells_;
  std::vector<int> netLimits_;
  std::vector<int> netWeight_;
  std::vector<int> netCells_;
  std::vector<int> netPinOffsets_;

  friend class MatrixCreator;
};
