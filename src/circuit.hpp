
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

class Circuit {
 public:
  /**
   * @brief Initialize a circuit
   */
  Circuit(int nbCells);

  /**
   * @brief Return the number of cells
   */
  int nbCells() const { return cellWidths.size(); }

  /**
   * @brief Return the number of nets
   */
  int nbNets() const { return netLimits.size() - 1; }

  /**
   * @brief Return the number of rows
   */
  int nbRows() const { return rows.size(); }

  /**
   * @brief Return the total number of pins
   */
  int nbPins() const { return netLimits.back(); }

  /**
   * @brief Return the number of pins for a given net
   *
   */
  int nbPins(int net) const {
    assert(net < nbNets());
    return netLimits[net + 1] - netLimits[net];
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

  /**
   * @brief Return the current orientation of the cell
   */
  CellOrientation orientation(int cell) const {
    assert(cell < nbCells());
    return cellOrientation[cell];
  }

  /**
   * @brief Return true if a cell is fixed
   */
  bool fixed(int cell) const {
    assert(cell < nbCells());
    return cellFixed[cell];
  }

  /**
   * @brief Get the area for a given cell
   */
  long long area(int cell) const {
    assert(cell < nbCells());
    return static_cast<long long>(cellWidths[cell]) *
           static_cast<long long>(cellHeights[cell]);
  }

  /**
   * @brief Return the cell associated with a given pin
   */
  int pinCell(int net, int i) const {
    assert(i < nbPins(net));
    return pinCells[netLimits[net] + i];
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
   * @brief Set the x position for a cell
   */
  void setCellX(int cell, float pos) {
    assert(cell < nbCells());
    cellX[cell] = pos;
  }

  /**
   * @brief Set the y position for a cell
   */
  void setCellY(int cell, float pos) {
    assert(cell < nbCells());
    cellY[cell] = pos;
  }

  /**
   * @brief Set the x position for all cells
   */
  void setCellX(const std::vector<int> &x) {
    assert(x.size() == nbCells());
    cellX = x;
  }

  /**
   * @brief Set the y position for all cells
   */
  void setCellY(const std::vector<int> &y) {
    assert(y.size() == nbCells());
    cellY = y;
  }

  /**
   * @brief Set the fixed status for all cells
   */
  void setCellFixed(const std::vector<char> &f) {
    assert(f.size() == nbCells());
    cellFixed = f;
  }

  /**
   * @brief Set the orientation for all cells
   */
  void setOrientation(const std::vector<CellOrientation> &orient) {
    assert(orient.size() == nbCells());
    cellOrientation = orient;
  }

  /**
   * @brief Set the width for all cells
   */
  void setCellWidths(const std::vector<int> &widths) {
    assert(widths.size() == nbCells());
    cellWidths = widths;
  }

  /**
   * @brief Set the height for all cells
   */
  void setCellHeights(const std::vector<int> &heights) {
    assert(heights.size() == nbCells());
    cellHeights = heights;
  }

  /**
   * @brief Set all nets
   */
  void setNets(const std::vector<int> &limits, const std::vector<int> &cells,
               const std::vector<int> &xOffsets,
               const std::vector<int> &yOffsets);

  /**
   * @brief Set all rows
   */
  void setRows(const std::vector<Rectangle> &r) { rows = r; }

  /**
   * @brief Return a bounding box of the placement area
   */
  Rectangle computePlacementArea() const;

  /**
   * @brief Return the rows after removing the obstacles
   */
  std::vector<Rectangle> computeRows() const;

  /**
   * @brief Direct creation for the C API
   */
  static Circuit createIspd(int nb_cells, int nb_nets, int *cell_widths,
                            int *cell_heights, char *cell_fixed,
                            int *net_limits, int *pin_cells, int *pin_x_offsets,
                            int *pin_y_offsets, int *cell_x, int *cell_y,
                            int *cell_orientation, int nb_rows, int *row_min_x,
                            int *row_max_x, int *row_min_y, int *row_max_y);

  /**
   * @brief Check the consistency of the datastructure
   */
  void check() const;

 public:
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
};