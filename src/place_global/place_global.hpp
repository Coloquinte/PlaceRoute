#pragma once

#include <vector>

#include "coloquinte.hpp"
#include "place_global/density_legalizer.hpp"
#include "place_global/net_model.hpp"

namespace coloquinte {
/**
 * @brief Main class for global placement
 */
class GlobalPlacer {
 public:
  struct Parameters {
    /**
     * @brief Maximum number of global placement steps
     */
    int maxNbSteps;

    /**
     * @brief Gap between lower and upper bound placement at which to stop
     * placement early
     */
    float gapTolerance;

    /**
     * @brief Distance at which the full displacement penalty is obtained, as a
     * fraction of the average standard cell length
     */
    float penaltyCutoffDistance;

    /**
     * @brief Initial average strength for the displacement penalty
     */
    float initialPenalty;

    /**
     * @brief Multiplicative factor for the displacement penalty at each
     * iteration
     */
    float penaltyUpdateFactor;

    /**
     * @brief Approximation distance of the continuous model, as a fraction of
     * the average standard cell length
     */
    float approximationDistance;

    /**
     * @brief Maximum number of conjugate gradient steps at each placement
     * iteration
     */
    int maxNbConjugateGradientSteps;

    /**
     * @brief Error tolerance to stop the conjugate gradient solver at each
     * placement iteration
     */
    float conjugateGradientErrorTolerance;

    /**
     * @brief Initialize the parameters with sensible defaults
     *
     * @param effort Placement effort between 1 and 9
     */
    explicit Parameters(int effort = 3);

    /**
     * @brief Check that the parameters make sense
     */
    void check() const;
  };

  /**
   * @brief Run global placement on the circuit representation
   *
   * @param circuit The circuit to be modified
   * @param effort Effort level, between 0 and 9
   */
  static void place(Circuit &circuit, int effort) {
    place(circuit, Parameters(effort));
  }

  /**
   * @brief Run global placement on the circuit representation
   *
   * @param circuit The circuit to be modified
   * @param params Placement parameters
   */
  static void place(Circuit &circuit, const GlobalPlacer::Parameters &params);

 private:
  /**
   * @brief Initialize the datastructure
   */
  explicit GlobalPlacer(Circuit &circuit,
                        const GlobalPlacer::Parameters &params);

  /**
   * @brief Run the whole global placement algorithm
   */
  void run();

  /**
   * @brief Return the net model length of the lower-bound placement
   */
  float valueLB() const;

  /**
   * @brief Return the net model length of the upper-bound placement
   */
  float valueUB() const;

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
   * @brief Compute the average cell size of the circuit
   */
  float computeAverageCellSize() const;

  /**
   * @brief Compute the base penalty forces from the area of the cells
   */
  std::vector<float> computePerCellPenalty() const;

  /**
   * @brief Approximation distance used by the continuous model
   */
  float approximationDistance() const {
    return params_.approximationDistance * averageCellLength_;
  }

  /**
   * @brief Distance at which the full displacement penalty is obtained
   */
  float penaltyCutoffDistance() const {
    return params_.penaltyCutoffDistance * averageCellLength_;
  }

 private:
  Circuit &circuit_;
  DensityLegalizer leg_;
  NetModel xtopo_;
  NetModel ytopo_;

  std::vector<float> xPlacementLB_;
  std::vector<float> yPlacementLB_;
  std::vector<float> xPlacementUB_;
  std::vector<float> yPlacementUB_;

  Parameters params_;

  float averageCellLength_;
  std::vector<float> perCellPenalty_;
  int step_;
  float penalty_;
};
}  // namespace coloquinte