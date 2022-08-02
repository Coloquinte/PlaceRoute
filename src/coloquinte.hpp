#pragma once

#include <cassert>
#include <vector>

/**
 * 2D point
 */
struct Point {
  int x;
  int y;
};

/**
 * 2D rectangle
 */
struct Rectangle {
  Rectangle() : minX(0), maxX(0), minY(0), maxY(0) {}
  Rectangle(int minX, int maxX, int minY, int maxY)
      : minX(minX), maxX(maxX), minY(minY), maxY(maxY) {}
  int minX;
  int maxX;
  int minY;
  int maxY;

  int width() const { return maxX - minX; }
  int height() const { return maxY - minY; }
  long long area() const { return (long long)width() * (long long)height(); }
  bool intersects(Rectangle o) const {
    return minX < o.maxX && o.minX < maxX && minY < o.maxY && o.minY < maxY;
  }
  static Rectangle intersection(Rectangle a, Rectangle b) {
    return Rectangle(std::max(a.minX, b.minX), std::min(a.maxX, b.maxX),
                     std::max(a.minY, b.minY), std::min(a.maxY, b.maxY));
  }
};

/**
 * @brief Cost model to use when doing legalization
 */
enum class LegalizationModel {
  /// L1 norm (sum of absolute values)
  L1,
  /// L2 norm (euclidean norm)
  L2,
  /// LInf norm (maximum of absolute values)
  LInf,
  /// Square of euclidean norm
  L2Squared
};

/**
 * @brief Orientation of a cell
 */
enum class CellOrientation {
  /// North (default orientation)
  N,
  /// South
  S,
  /// West
  W,
  /// East
  E,
  /// Flip + North
  FN,
  /// Flip + South
  FS,
  /// Flip + West
  FW,
  /// Flip + East
  FE
};

/**
 * @brief Compute the norm of the 2D vector with the given cost model
 */
float norm(float x, float y, LegalizationModel leg);

/**
 * @brief Compute the norm of the 2D vector with the given cost model
 */
long long norm(int x, int y, LegalizationModel leg);

/**
 * Representation of a flat circuit from an ISPD benchmark
 */

struct Circuit {
  std::vector<int> cellWidths;
  std::vector<int> cellHeights;
  std::vector<char> cellFixed;
  std::vector<int> netLimits;
  std::vector<int> pinCells;
  std::vector<int> pinXOffsets;
  std::vector<int> pinYOffsets;
  std::vector<int> cellX;
  std::vector<int> cellY;
  std::vector<CellOrientation> cellOrientation;
  std::vector<Rectangle> rows;

  Circuit(int nbCells) {}

  int nbCells() const { return cellWidths.size(); }

  int nbNets() const { return netLimits.size() - 1; }

  int nbPins() const { return netLimits.back(); }

  bool isFixed(int cell) const {
    assert(cell < nbCells());
    return cellFixed[cell];
  }

  long long area(int cell) const {
    assert(cell < nbCells());
    return static_cast<long long>(cellWidths[cell]) *
           static_cast<long long>(cellHeights[cell]);
  }

  void setCellX(const std::vector<int> &x) {
    assert(x.size() == nbCells());
    cellX = x;
  }

  void setCellY(const std::vector<int> &y) {
    assert(y.size() == nbCells());
    cellY = y;
  }

  void setCellWidths(const std::vector<int> &widths) {
    assert(widths.size() == nbCells());
    cellWidths = widths;
  }

  void setCellHeights(const std::vector<int> &heights) {
    assert(heights.size() == nbCells());
    cellHeights = heights;
  }

  void setNets(const std::vector<int> &netLimits,
               const std::vector<int> &pinCells,
               const std::vector<int> &pinXOffsets,
               const std::vector<int> &pinYOffsets);

  void setOrientation(const std::vector<CellOrientation> &orient) {
    assert(orient.size() == nbCells());
    cellOrientation = orient;
  }

  /**
   * @brief Return the current orientation of the cell
   */
  CellOrientation orientation(int cell) const {
    assert(cell < nbCells());
    return cellOrientation[cell];
  }

  /**
   * @brief Return the current width of the cell (depends on its orientation)
   */
  int width(int cell) const;

  /**
   * @brief Return the current height of the cell (depends on its orientation)
   */
  int height(int cell) const;

  /**
   * @brief Return the current x position of the cell
   */
  int x(int cell) const {
    assert(cell < nbCells());
    return cellX[cell];
  }

  /**
   * @brief Return the current y position of the cell
   */
  int y(int cell) const {
    assert(cell < nbCells());
    return cellY[cell];
  }

  static Circuit createIspd(int nb_cells, int nb_nets, int *cell_widths,
                            int *cell_heights, char *cell_fixed,
                            int *net_limits, int *pin_cells, int *pin_x_offsets,
                            int *pin_y_offsets, int *cell_x, int *cell_y,
                            int *cell_orientation, int nb_rows, int *row_min_x,
                            int *row_max_x, int *row_min_y, int *row_max_y);

  void check() const;
};

extern "C" {
void place_ispd(int nb_cells, int nb_nets, int *cell_widths, int *cell_heights,
                char *cell_fixed, int *net_limits, int *pin_cells,
                int *pin_x_offsets, int *pin_y_offsets, int *cell_x,
                int *cell_y, int *cell_orientation, int nb_rows, int *row_min_x,
                int *row_max_x, int *row_min_y, int *row_max_y);
void benchmark_quadratic_models(int nb_cells, int nb_nets, int *cell_widths,
                                int *cell_heights, char *cell_fixed,
                                int *net_limits, int *pin_cells,
                                int *pin_x_offsets, int *pin_y_offsets,
                                int *cell_x, int *cell_y, int *cell_orientation,
                                int model_type, int nb_steps, float epsilon,
                                float relaxation);
}
