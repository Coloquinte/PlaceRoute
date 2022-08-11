#pragma once

#include <vector>

#include "coloquinte.hpp"
#include "place_global/density_legalizer.hpp"
#include "place_global/net_model.hpp"

/**
 * @brief Main class for global placement
 */
class GlobalPlacer {
 public:
  /**
   * @brief Run global placement on the circuit representation
   *
   * @param circuit The circuit to be modified
   * @param effort Effort level, between 0 and 9
   */
  static void place(Circuit &circuit, int effort);

 private:
  /**
   * @brief Initialize the datastructure
   */
  GlobalPlacer(Circuit &circuit);

  /**
   * @brief Return the net model length of the lower-bound placement
   */
  float valueLB() const;

  /**
   * @brief Return the net model length of the upper-bound placement
   */
  float valueUB() const;

  /**
   * @brief Initialize the algorithm's parameters
   */
  void initParameters();

  /**
   * @brief Update the parameters between two placement iterations
   */
  void updateParameters();

  /**
   * @brief Obtain the initial lower-bound placement
   */
  void runInitialLB();

  /**
   * @brief Obtain the lower-bound placement for one iteration
   */
  void runLB();

  /**
   * @brief Obtain the upper-bound placement for one iteration
   */
  void runUB();

  /**
   * @brief Compute the base penalty forces from the area of the cells
   */
  std::vector<float> computeBaseForces() const;

 private:
  Circuit &circuit_;
  DensityLegalizer leg_;
  NetModel xtopo_;
  NetModel ytopo_;
  std::vector<float> baseForces_;

  std::vector<float> xPlacementLB_;
  std::vector<float> yPlacementLB_;
  std::vector<float> xPlacementUB_;
  std::vector<float> yPlacementUB_;

  float epsilon_;
  float cutoffDistance_;
  float updateFactor_;
  int maxNbSteps_;

  float forceFactor_;
  int step_;
};
