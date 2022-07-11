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
              const std::vector<int> &widths,
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
     * Place a single cell to the datastructure
     */
    void placeCell(int cell);

    /**
     * Place a single cell in a given row; return the success and the coordinates
     */
    std::pair<bool, long long> placeCell(int cell, int row);

    /**
     * Place multiple cells in the datastructure with branch-and-bound
     */
    void placeCells(const std::vector<int> &cells, LegalizationModel model);

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

