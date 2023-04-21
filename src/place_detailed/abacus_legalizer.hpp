
#include "place_detailed/legalizer.hpp"

namespace coloquinte {
/**
 * Legalization following the "Abacus" algorithm with L1 cost
 *
 * Cells are placed in rows, and are allowed to push other cells to get to a
 * better overall solution.
 */
class AbacusLegalizer : public LegalizerBase {
 public:
  /**
   * @brief Initialization of the datastructure
   */
  AbacusLegalizer(const std::vector<Row> &rows, const std::vector<int> &width,
                  const std::vector<int> &height,
                  const std::vector<CellRowPolarity> &polarities,
                  const std::vector<int> &targetX,
                  const std::vector<int> &targetY,
                  const std::vector<CellOrientation> &targetOrientation);

  /**
   * @brief Run the algorithm
   */
  void run();

  /**
   * @brief Place a single cell optimally
   * Return true if successful
   */
  void placeCell(int cell);

  /**
   * @brief Find the optimal placement in the row
   *
   * Return a pair: true if successful, and the distance with this placement
   */
  std::pair<bool, long long> evaluatePlacement(int cell, int row);

  /**
   * @brief Check consistency
   */
  void check() const;

 private:
  // Algorithm state
  std::vector<std::vector<int> > rowToCells_;
  std::vector<RowLegalizer> rowLegalizers_;
};
}  // namespace coloquinte