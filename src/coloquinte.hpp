#pragma once

#include <cassert>
#include <vector>

namespace coloquinte {
/**
 * @brief 2D point
 */
struct Point {
  Point() : x(0), y(0) {}
  Point(int x, int y) : x(x), y(y) {}

  int x;
  int y;
};

/**
 * @brief 2D rectangle
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
 * @brief Flat representation of a circuit
 */
class Circuit {
 public:
  /**
   * @brief Initialize a circuit
   */
  explicit Circuit(int nbCells);

  /**
   * @brief Return the number of cells
   */
  int nbCells() const { return cellWidth_.size(); }

  /**
   * @brief Return the number of nets
   */
  int nbNets() const { return netLimits_.size() - 1; }

  /**
   * @brief Return the number of rows
   */
  int nbRows() const { return rows_.size(); }

  /**
   * @brief Return the total number of pins
   */
  int nbPins() const { return netLimits_.back(); }

  /**
   * @brief Set the x position for all cells
   */
  const std::vector<int> &cellX() const { return cellX_; }

  /**
   * @brief Set the x position for all cells
   */
  void setCellX(const std::vector<int> &x);

  /**
   * @brief Get the y position for all cells
   */
  const std::vector<int> &cellY() const { return cellY_; }

  /**
   * @brief Set the y position for all cells
   */
  void setCellY(const std::vector<int> &y);

  /**
   * @brief Get the fixed status for all cells
   */
  const std::vector<bool> &cellFixed() const { return cellFixed_; }

  /**
   * @brief Set the fixed status for all cells
   */
  void setCellFixed(const std::vector<bool> &f);

  /**
   * @brief Get the orientation for all cells
   */
  const std::vector<CellOrientation> &cellOrientation() const {
    return cellOrientation_;
  }

  /**
   * @brief Set the orientation for all cells
   */
  void setCellOrientation(const std::vector<CellOrientation> &orient);

  /**
   * @brief Get the width for all cells
   */
  const std::vector<int> &cellWidth() const { return cellWidth_; }

  /**
   * @brief Set the width for all cells
   */
  void setCellWidth(const std::vector<int> &widths);

  /**
   * @brief Get the height for all cells
   */
  const std::vector<int> &cellHeight() const { return cellHeight_; }

  /**
   * @brief Set the height for all cells
   */
  void setCellHeight(const std::vector<int> &heights);

  /**
   * @brief Get all rows
   */
  const std::vector<Rectangle> &rows() const { return rows_; }

  /**
   * @brief Set all rows
   */
  void setRows(const std::vector<Rectangle> &r) { rows_ = r; }

  /**
   * @brief Add a single net
   */
  void addNet(const std::vector<int> &cells, const std::vector<int> &xOffsets,
              const std::vector<int> &yOffsets);

  /**
   * @brief Set all nets
   */
  void setNets(const std::vector<int> &limits, const std::vector<int> &cells,
               const std::vector<int> &xOffsets,
               const std::vector<int> &yOffsets);

  /**
   * @brief Return a bounding box of the placement area
   */
  Rectangle computePlacementArea() const;

  /**
   * @brief Return the rows after removing the obstacles
   */
  std::vector<Rectangle> computeRows() const;

  /**
   * @brief Return the current width of the cell (depends on its orientation)
   */
  int placedWidth(int cell) const;

  /**
   * @brief Return the current height of the cell (depends on its orientation)
   */
  int placedHeight(int cell) const;

  /**
   * @brief Return the current x position of the cell
   */
  int x(int cell) const {
    assert(cell < nbCells());
    return cellX_[cell];
  }

  /**
   * @brief Return the current y position of the cell
   */
  int y(int cell) const {
    assert(cell < nbCells());
    return cellY_[cell];
  }

  /**
   * @brief Return the current orientation of the cell
   */
  CellOrientation orientation(int cell) const {
    assert(cell < nbCells());
    return cellOrientation_[cell];
  }

  /**
   * @brief Return true if a cell is fixed
   */
  bool fixed(int cell) const {
    assert(cell < nbCells());
    return cellFixed_[cell];
  }

  /**
   * @brief Get the area for a given cell
   */
  long long area(int cell) const {
    assert(cell < nbCells());
    return static_cast<long long>(cellWidth_[cell]) *
           static_cast<long long>(cellHeight_[cell]);
  }

  /**
   * @brief Return the number of pins for a given net
   *
   */
  int nbPinsNet(int net) const {
    assert(net < nbNets());
    return netLimits_[net + 1] - netLimits_[net];
  }

  /**
   * @brief Return the cell associated with a given pin
   */
  int pinCell(int net, int i) const {
    assert(i < nbPinsNet(net));
    return pinCells_[netLimits_[net] + i];
  }

  /**
   * @brief Return the current x offset of the pin (depends on its orientation)
   */
  int pinXOffset(int net, int i) const;

  /**
   * @brief Return the current y offset of the pin (depends on its orientation)
   */
  int pinYOffset(int net, int i) const;

  /**
   * @brief Return the current half-perimeter wirelength of the circuit
   */
  long long hpwl() const;

  /**
   * @brief Check the consistency of the datastructure
   */
  void check() const;

 public:
  std::vector<int> netLimits_;
  std::vector<int> pinCells_;
  std::vector<int> pinXOffsets_;
  std::vector<int> pinYOffsets_;
  std::vector<int> cellWidth_;
  std::vector<int> cellHeight_;
  std::vector<bool> cellFixed_;
  std::vector<int> cellX_;
  std::vector<int> cellY_;
  std::vector<CellOrientation> cellOrientation_;
  std::vector<Rectangle> rows_;
};

/**
 * @brief Run the whole placement algorithm
 */
void place(Circuit &circuit, int effort);

extern "C" {
/**
 * @brief Ugly interface to easily call Coloquinte from external code
 */
int place_ispd(int nb_cells, int nb_nets, int *cell_widths, int *cell_heights,
               char *cell_fixed, int *net_limits, int *pin_cells,
               int *pin_x_offsets, int *pin_y_offsets, int *cell_x, int *cell_y,
               int *cell_orientation, int nb_rows, int *row_min_x,
               int *row_max_x, int *row_min_y, int *row_max_y, int effort);
}
}  // namespace coloquinte
