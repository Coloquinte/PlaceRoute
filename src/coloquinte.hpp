#pragma once

#include <cassert>
#include <functional>
#include <iosfwd>
#include <optional>
#include <string>
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

  /**
   * @brief Returns whether the two rectangles intersect
   */
  bool intersects(Rectangle o) const {
    return minX < o.maxX && o.minX < maxX && minY < o.maxY && o.minY < maxY;
  }

  /**
   * @brief Returns whether the rectangle contains another one
   */
  bool contains(Rectangle o) const {
    return minX <= o.minX && maxX >= o.maxX && minY <= o.minY && maxY >= o.maxY;
  }

  /**
   * @brief Return the overlapping area of the rectangles if they intersect
   */
  static Rectangle intersection(Rectangle a, Rectangle b) {
    return Rectangle(std::max(a.minX, b.minX), std::min(a.maxX, b.maxX),
                     std::max(a.minY, b.minY), std::min(a.maxY, b.maxY));
  }

  /**
   * @brief Return a brief description of the rectangle
   */
  std::string toString() const;
};

/**
 * @brief Step of the placement process, for use in callbacks
 */
enum class PlacementStep { LowerBound, UpperBound, Detailed };

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
  /// Square of L1 norm
  L1Squared,
  /// Square of L2 norm
  L2Squared,
  /// Square of LInf norm
  LInfSquared
};

/**
 * @brief Net model to use for continuous optimization
 */
enum class NetModelOption {
  /**
   * @brief Classical bound-to-bound model (connect to the extreme points of the
   * net)
   */
  BoundToBound,
  /**
   * @brief Star model, with a virtual pin at the middle
   */
  Star,
  /**
   * @brief Fully-connected model
   */
  Clique,
  /**
   * @brief Relaxed star model, with weighting more similar to bound-to-bound
   */
  LightStar
};

/**
 * @brief Orientation of a cell
 */
enum class CellOrientation {
  /// @brief North (default orientation)
  N = 0,
  R0 = 0,
  /// @brief South
  S = 1,
  R180 = 1,
  /// @brief West
  W = 2,
  R90 = 2,
  /// @brief East
  E = 3,
  R270 = 3,
  /// @brief Flip + North (vertical mirror, y unchanged)
  FN = 4,
  MY = 4,
  /// @brief Flip + South (horizontal mirror, x unchanged)
  FS = 5,
  MX = 5,
  /// @brief Flip + West (MX then W)
  FW = 6,
  MX90 = 6,
  /// @brief Flip + East (MY then W)
  FE = 7,
  MY90 = 7,
  /// @brief Special value for invalid/forbidden/whatever
  INVALID = 8,
  /// @brief Special value for unknown
  UNKNOWN = 9
};

/**
 * @brief Polarity of a cell with respect to standard cell rows
 *
 * A cell with a polarity can only be placed in rows matching this polarity.
 * If the cell has an even number of rows, it can only be placed in every other
 * row (those with the same polarity). If it has an odd number of rows, it can
 * be placed in any row, but requires a vertical flip for half of them.
 */
enum class CellRowPolarity {
  /// Any row orientation allowed.
  ANY,
  /// Must be in the same orientation as the row at the bottom.
  /// Meaningful for standard cells with odd number of rows.
  SAME,
  /// Must be in the opposite orientation as the row at the bottom.
  /// Meaningful for standard cells with odd number of rows.
  OPPOSITE,
  /// Must be a series of row with alternating orientation starting with N or W.
  /// Flipping is allowed.
  /// Meaningful for standard cells with even number of rows.
  NW,
  /// Must be a series of row with alternating orientation starting with S or E.
  /// Flipping is allowed.
  /// Meaningful for standard cells with even number of rows.
  SE,
};

struct CellPlacement {
  Point position;
  CellOrientation orientation;

  CellPlacement() {}
  CellPlacement(int x, int y, CellOrientation o)
      : position(x, y), orientation(o) {}
};

using PlacementSolution = std::vector<CellPlacement>;

std::string toString(CellOrientation o);
std::string toString(LegalizationModel model);
std::string toString(NetModelOption model);
std::string toString(CellRowPolarity pol);

std::ostream &operator<<(std::ostream &, CellOrientation);
std::ostream &operator<<(std::ostream &, LegalizationModel);
std::ostream &operator<<(std::ostream &, NetModelOption);
std::ostream &operator<<(std::ostream &, CellRowPolarity);

/**
 * @brief Return the "opposite" row orientation (N <--> FS)
 */
CellOrientation oppositeRowOrientation(CellOrientation o);

/**
 * @brief Return the valid orientation for this combination of polarity and row
 * orientation
 */
CellOrientation cellOrientationInRow(CellRowPolarity cellPolarity,
                                     CellOrientation rowOrientation);

/**
 * Return whether the cell view is turned (exchange width and height)
 */
bool isTurn(CellOrientation o);

/**
 * @brief Placement callbacks are called with the current step.
 * The circuit can be accessed (not modified) when the callback is called
 */
using PlacementCallback = std::function<void(PlacementStep)>;

/**
 * @brief Parameters for the rough legalization step
 */
struct RoughLegalizationParameters {
  /**
   * @brief Cost model used for rough legalization
   */
  LegalizationModel costModel;

  /**
   * @brief Number of rough legalization steps at each placement iteration
   */
  int nbSteps;

  /**
   * @brief Size of the rough legalization bin relative to the standard cell
   * height
   */
  double binSize;

  /**
   * @brief Number of x- or y- aligned bins reoptimized together
   */
  int lineReoptSize;

  /**
   * @brief Overlap between two sets of x- or y- reoptimized bins
   */
  int lineReoptOverlap;

  /**
   * @brief Number of diag-aligned bins reoptimized together
   */
  int diagReoptSize;

  /**
   * @brief Overlap between two sets of diag-reoptimized bins
   */
  int diagReoptOverlap;

  /**
   * @brief Size ot the square of bins reoptimized together
   */
  int squareReoptSize;

  /**
   * @brief Overlap between two reoptimized squares bins
   */
  int squareReoptOverlap;

  /**
   * @brief Use unidimensional transportation to reoptimize many bins at once
   */
  bool unidimensionalTransport;

  /**
   * @brief Small quadratic penalty to ensure that closer cells are selected
   * despite the L1 distance, normalized by the placement area
   */
  double quadraticPenalty;

  /**
   * @brief Margin to use on the sides of each row to account for lost space,
   * relative to the standard cell height
   */
  double sideMargin;

  /**
   * @brief Decide how much to coarsen compared based on the current distance to
   * legal
   */
  double coarseningLimit;

  /**
   * @brief Blending between lower-bound and upper-bound placement to decide the
   * target for rough legalization; 0 to use lower bound, 1 to use upper bound
   */
  double targetBlending;

  /**
   * @brief Initialize the parameters
   */
  explicit RoughLegalizationParameters(int effort);

  /**
   * @brief Obtain a string representation
   */
  std::string toString() const;

  /**
   * @brief Check that the parameters make sense
   */
  void check() const;
};

struct PenaltyParameters {
  /**
   * @brief Distance at which the full displacement penalty is obtained,
   * relative to the average standard cell length
   */
  double cutoffDistance;

  /**
   * @brief Per-step update factor for the cutoff distance
   */
  double cutoffDistanceUpdateFactor;

  /**
   * @brief Exponent applied to the cell area to compute the penalty
   */
  double areaExponent;

  /**
   * @brief Initial average strength for the displacement penalty
   */
  double initialValue;

  /**
   * @brief Multiplicative factor for the displacement penalty at each
   * iteration
   */
  double updateFactor;

  /**
   * @brief Blending between lower-bound and upper-bound placement to decide the
   * target for penalty; 0 to use lower bound, 1 to use upper bound
   */
  double targetBlending;

  /**
   * @brief Initialize the parameters
   */
  explicit PenaltyParameters(int effort);

  /**
   * @brief Obtain a string representation
   */
  std::string toString() const;

  /**
   * @brief Check that the parameters make sense
   */
  void check() const;
};

struct ContinuousModelParameters {
  /**
   * @brief Cost model for the continuous optimization
   */
  NetModelOption netModel;

  /**
   * @brief Approximation distance of the continuous model, relative to
   * the average standard cell length
   */
  double approximationDistance;

  /**
   * @brief Per-step update factor for the approximation distance
   */
  double approximationDistanceUpdateFactor;

  /**
   * @brief Maximum number of conjugate gradient steps at each placement
   * iteration
   */
  int maxNbConjugateGradientSteps;

  /**
   * @brief Error tolerance to stop the conjugate gradient solver at each
   * placement iteration
   */
  double conjugateGradientErrorTolerance;

  /**
   * @brief Initialize the parameters
   */
  explicit ContinuousModelParameters(int effort);

  /**
   * @brief Obtain a string representation
   */
  std::string toString() const;

  /**
   * @brief Check that the parameters make sense
   */
  void check() const;
};

/**
 * @brief Parameters for the global placer
 */
struct GlobalPlacerParameters {
  /**
   * @brief Maximum number of global placement steps
   */
  int maxNbSteps;

  /**
   * @brief Number of initial placement steps, without penalization
   */
  int nbInitialSteps;

  /**
   * @brief Number of steps before one rough legalization step
   */
  int nbStepsBeforeRoughLegalization;

  /**
   * @brief Gap between lower and upper bound placement at which to stop
   * global placement early
   */
  double gapTolerance;

  /**
   * @brief Distance between lower and upper bound placement at which to stop
   * global placement early
   */
  double distanceTolerance;

  /**
   * @brief Blending between lower-bound and upper-bound placement at export
   * time; 0 to use lower bound, 1 to use upper bound
   */
  double exportBlending;

  /**
   * @brief Parameters for the continuous model
   */
  ContinuousModelParameters continuousModel;

  /**
   * @brief Parameters for the rough legalization
   */
  RoughLegalizationParameters roughLegalization;

  /**
   * @brief Parameters for the legalization penalty
   */
  PenaltyParameters penalty;

  /**
   * @brief Noise introduced to randomize the algorithm. Note that 0 will erase
   * any effect of the seed
   */
  double noise;

  /**
   * @brief Initialize the parameters with sensible defaults
   */
  explicit GlobalPlacerParameters(int effort);

  /**
   * @brief Obtain a string representation
   */
  std::string toString() const;

  /**
   * @brief Check that the parameters make sense
   */
  void check() const;
};

/**
 * @brief Parameters for the legalization
 */
struct LegalizationParameters {
  /**
   * @brief Cost model used for legalization
   */
  LegalizationModel costModel;

  /**
   * @brief Weight placed on the width when ordering of the legalization
   * heuristic
   *
   * @details This weight decides which cells are legalized first; 0.0 orders
   * the cells by their left side, 1.0 by their right side, 0.5 by their middle.
   */
  double orderingWidth;

  /**
   * @brief Weight placed on the height when ordering of the legalization
   * heuristic
   *
   * @details This weight decides which cells are legalized first; 0.0
   * ignores the height, negative values will legalize high cells first.
   */
  double orderingHeight;

  /**
   * @brief Weight placed on the y position when ordering of the legalization
   * heuristic
   *
   * @details This weight decides which cells are legalized first, and tends to
   * start to legalize up (-1.0) or down (1.0). It should be small in absolute
   * value, as x ordering is preferred.
   */
  double orderingY;

  /**
   * @brief Initialize the parameters with sensible defaults
   */
  explicit LegalizationParameters(int effort);

  /**
   * @brief Obtain a string representation
   */
  std::string toString() const;

  /**
   * @brief Check that the parameters make sense
   */
  void check() const;
};

/**
 * @brief Parameters for the detailed placer
 */
struct DetailedPlacerParameters {
  /**
   * @brief Number of optimization passes
   */
  int nbPasses;

  /**
   * @brief Number of closest neighbours on each side considered during local
   * search
   */
  int localSearchNbNeighbours;

  /**
   * @brief Number of closest rows on each side considered during local
   * search
   */
  int localSearchNbRows;

  /**
   * @brief Number of rows considered together when optimizing shifts
   */
  int shiftNbRows;

  /**
   * @brief Maximum number of cells considered together when optimizing shifts
   */
  int shiftMaxNbCells;

  /**
   * @brief Number of rows considered simultaneously for reordering
   */
  int reorderingNbRows;

  /**
   * @brief Maximum number of cells considered together for reordering
   */
  int reorderingMaxNbCells;

  /**
   * @brief Initialize the parameters with sensible defaults
   */
  explicit DetailedPlacerParameters(int effort);

  /**
   * @brief Obtain a string representation
   */
  std::string toString() const;

  /**
   * @brief Check that the parameters make sense
   */
  void check() const;
};

/**
 * @brief Parameters for the full placer
 */
struct ColoquinteParameters {
  GlobalPlacerParameters global;

  LegalizationParameters legalization;

  DetailedPlacerParameters detailed;

  /**
   * @brief Random seed
   */
  int seed;

  /**
   * @brief Initialize the parameters with sensible defaults
   *
   * @param effort Placement effort between 1 and 9
   * @param seed Random seed
   */
  explicit ColoquinteParameters(int effort = 3, int seed = -1);

  /**
   * @brief Obtain a string representation
   */
  std::string toString() const;

  /**
   * @brief Check that the parameters make sense
   */
  void check() const;
};

/**
 * @brief Representation of a placement row in a circuit
 */
struct Row : public Rectangle {
  /**
   * @brief Placement area for the rows
   */
  CellOrientation orientation;

  /**
   * @brief Construct a row from its area and orientation
   */
  Row(Rectangle a, CellOrientation orient)
      : Rectangle(a), orientation(orient) {}

  /**
   * @brief Construct a row from its coordinates and orientation
   */
  Row(int minX, int maxX, int minY, int maxY, CellOrientation orient)
      : Rectangle(minX, maxX, minY, maxY), orientation(orient) {}
  /**
   * @brief Obtain free placement space after removing the obstacles
   */
  std::vector<Row> freespace(const std::vector<Rectangle> &obstacles) const;
};

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
   * @brief Get the x position for all cells
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
  const std::vector<bool> &cellIsFixed() const { return cellIsFixed_; }

  /**
   * @brief Set the fixed status for all cells
   */
  void setCellIsFixed(const std::vector<bool> &f);

  /**
   * @brief Get the obstruction status for all cells
   */
  const std::vector<bool> &cellIsObstruction() const {
    return cellIsObstruction_;
  }

  /**
   * @brief Set the obstruction status for all cells
   */
  void setCellIsObstruction(const std::vector<bool> &f);

  /**
   * @brief Set the whole placement solution
   */
  void setSolution(const PlacementSolution &sol);

  /**
   * @brief Get the row orientation for all cells
   *
   * This is the orientation that the row at the bottom of the cell should have;
   * for single-row cells, this is always N (north), but multi-row cells may
   * require row with a specific orientation.
   */
  const std::vector<CellRowPolarity> &cellRowPolarity() const {
    return cellRowPolarity_;
  }

  /**
   * @brief Set the row orientation for all cells
   */
  void setCellRowPolarity(const std::vector<CellRowPolarity> &f);

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
   * @brief Get all rows
   */
  const std::vector<Row> &rows() const { return rows_; }

  /**
   * @brief Set all rows
   */
  void setRows(const std::vector<Row> &r);

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
   * @brief Return the height of the rows
   */
  int rowHeight() const;

  /**
   * @brief Setup the rows to fill a given placement area
   */
  void setupRows(Rectangle placementArea, int rowHeight,
                 bool alternatingOrientation = true,
                 bool initialOrientation = true);

  /**
   * @brief Return the rows after removing the obstacles
   */
  std::vector<Row> computeRows(
      const std::vector<Rectangle> &additionalObstacles =
          std::vector<Rectangle>()) const;

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
   * @brief Return the placement occupied by the cell
   */
  Rectangle placement(int cell) const {
    return Rectangle(x(cell), x(cell) + placedWidth(cell), y(cell),
                     y(cell) + placedHeight(cell));
  }

  std::vector<Rectangle> cellPlacement() const {
    std::vector<Rectangle> ret;
    for (int i = 0; i < nbCells(); ++i) {
      ret.push_back(placement(i));
    }
    return ret;
  }

  PlacementSolution solution() const {
    PlacementSolution ret;
    for (int i = 0; i < nbCells(); ++i) {
      ret.emplace_back(x(i), y(i), orientation(i));
    }
    return ret;
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
  bool isFixed(int cell) const {
    assert(cell < nbCells());
    return cellIsFixed_[cell];
  }

  /**
   * @brief Return true if a cell cannot overlap with others
   */
  bool isObstruction(int cell) const {
    assert(cell < nbCells());
    return cellIsObstruction_[cell];
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
   * @brief Run the whole placement algorithm (global + detailed)
   */
  void place(int effort) {
    placeGlobal(effort);
    placeDetailed(effort);
  }

  /**
   * @brief Run the global placement algorithm
   */
  void placeGlobal(int effort) { placeGlobal(ColoquinteParameters(effort)); }

  /**
   * @brief Run the global placement algorithm
   */
  void placeGlobal(const ColoquinteParameters &params,
                   const std::optional<PlacementCallback> &callback = {});

  /**
   * @brief Run the legalization algorithm
   */
  void legalize(int effort) { legalize(ColoquinteParameters(effort)); }

  /**
   * @brief Run the legalization algorithm
   */
  void legalize(const ColoquinteParameters &params,
                const std::optional<PlacementCallback> &callback = {});

  /**
   * @brief Run the detailed placement algorithm
   */
  void placeDetailed(int effort) {
    placeDetailed(ColoquinteParameters(effort));
  }

  /**
   * @brief Run the detailed placement algorithm
   */
  void placeDetailed(const ColoquinteParameters &params,
                     const std::optional<PlacementCallback> &callback = {});

  /**
   * @brief Return a brief description of the circuit
   */
  std::string toString() const;

  /**
   * @brief Obtain a detailed report
   */
  std::string report() const;

  /**
   * @brief Check the consistency of the datastructure
   */
  void check() const;

  /**
   * @brief Compute the mean displacement between two placement solutions
   */
  float meanDisruption(const PlacementSolution &a, const PlacementSolution &b,
                       LegalizationModel costModel);

  /**
   * @brief Compute the root-mean-square displacement between two placement
   * solutions
   */
  float rmsDisruption(const PlacementSolution &a, const PlacementSolution &b,
                      LegalizationModel costModel);

  /**
   * @brief Compute the max displacement between two placement solutions
   */
  float maxDisruption(const PlacementSolution &a, const PlacementSolution &b,
                      LegalizationModel costModel);

  /**
   * @brief Return the flag indicating that the sizes have changes
   */
  bool hasSizeUpdate() const { return hasSizeUpdate_; }

  /**
   * @brief Set the flag indicating that the sizes have changes
   */
  void setSizeUpdate() { hasSizeUpdate_ = true; }

  /**
   * @brief Set the flag indicating that the sizes have changes
   */
  void clearSizeUpdate() { hasSizeUpdate_ = false; }

  /**
   * @brief Indicate that the circuit is currently being used
   */
  void setInUse() { isInUse_ = true; }

  /**
   * @brief Indicate that the circuit is not being used anymore
   */
  void clearInUse() { isInUse_ = false; }

 private:
  std::vector<float> allDistances(const PlacementSolution &a,
                                  const PlacementSolution &b,
                                  LegalizationModel costModel);

  /**
   * Check that the circuit is not being worked on right now
   */
  void checkNotInUse() const;

 public:
  std::vector<int> netLimits_;
  std::vector<int> pinCells_;
  std::vector<int> pinXOffsets_;
  std::vector<int> pinYOffsets_;
  std::vector<int> cellWidth_;
  std::vector<int> cellHeight_;
  std::vector<bool> cellIsFixed_;
  std::vector<bool> cellIsObstruction_;
  std::vector<CellRowPolarity> cellRowPolarity_;
  std::vector<int> cellX_;
  std::vector<int> cellY_;
  std::vector<CellOrientation> cellOrientation_;
  std::vector<Row> rows_;

  // Set whenever the circuit is actively being worked on
  bool isInUse_;

  // Set whenever an update has been made to the circuit
  bool hasSizeUpdate_;
};
}  // namespace coloquinte
