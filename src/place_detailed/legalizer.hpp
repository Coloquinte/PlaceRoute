#pragma once

#include "coloquinte.hpp"


class LegalizerRow {
  public:
    Rectangle area;
    std::vector<int> cellIds;
    std::vector<int> cellX;

    int nbCells() const { return cellIds.size(); }
};

/**
 * Algorithms to obtain a legal placement for standard cells
 *
 */
class Legalizer {
  public:
    /**
     * Run the whole legalization algorithm
     */
    static void legalize(Circuit &circuit);

    /**
     * Initialize the datastructure
     */
    static Legalizer fromIspdCircuit(const Circuit &circuit);
    Legalizer(const std::vector<Rectangle> & rows,
              const std::vector<int> &width,
              const std::vector<int> &targetX,
              const std::vector<int> &targetY);

    int nbRows() const { return rows_.size(); }
    int nbCells() const { return cellWidth_.size(); }

    LegalizationModel costModel() const { return costModel_; }
    void setCostModel(LegalizationModel m) { costModel_ = m; }

    const std::vector<int> &cellWidth() const { return cellWidth_; }
    const std::vector<int> &cellTargetX() const { return cellTargetX_; }
    const std::vector<int> &cellTargetY() const { return cellTargetY_; }

    /**
     * Run the algorithm
     */
    void run();

    /**
     * Compute the coordinates associated with the legalization
     */
    std::vector<int> cellLegalX() const;
    std::vector<int> cellLegalY() const;

    /**
     * Check consistency of the datastructure
     */
    void check() const;

  private:
    /**
     * Compute the distance with the current cost model
     */
    long long distance(int x1, int x2, int y1, int y2) const;

    /**
     * Place a single cell optimally; return true if successful
     */
    bool placeCellOptimally(int cell);

    /**
     * Simulate placing a single cell in a given row
     * Return true if successful and the X coordinate
     */
    std::pair<bool, int> placeCellOptimally(int cell, int row) const;

    /**
     * Materialize the placement of a cell
     */
    void doPlacement(int cell, int row, int x);

    /**
     * Export the placement obtained to the circuit datastructure
     */
    void exportPlacement(Circuit &circuit);

  private:
    // Placement data
    LegalizationModel costModel_;
    std::vector<int> cellWidth_;
    std::vector<int> cellTargetX_;
    std::vector<int> cellTargetY_;

    // Placement status
    std::vector<Rectangle> rows_;
    std::vector<std::vector<int> > rowToCells_;
    std::vector<std::vector<int> > rowToX_;
    std::vector<int> cellToRow_;
    std::vector<int> cellToX_;
    std::vector<int> cellToY_;
};

