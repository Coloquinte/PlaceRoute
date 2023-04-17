

#include "place_detailed/legalizer.hpp"

namespace coloquinte {
class TetrisLegalizer : public LegalizerBase {
  /**
   * @brief Initialization of the datastructure
   */
  TetrisLegalizer(const std::vector<Row> &rows, const std::vector<int> &width,
                  const std::vector<int> &height,
                  const std::vector<CellRowPolarity> &polarities,
                  const std::vector<int> &targetX,
                  const std::vector<int> &targetY);

  /**
   * @brief Run the algorithm
   */
  void run();

  /**
   * @brief Place a single cell
   */
  void placeCell(int c);

  /**
   * @brief Return possible placement intervals for a given width, height and y
   */
  std::vector<std::pair<int, int> > getPossibleIntervals(int w, int h,
                                                         int y) const;

  /// @brief First free position in the row
  std::vector<int> rowFreePos_;
};
}  // namespace coloquinte