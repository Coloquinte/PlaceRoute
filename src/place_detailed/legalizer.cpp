
#include "place_detailed/legalizer.hpp"

#include <algorithm>
#include <cassert>
#include <limits>

Legalizer Legalizer::fromIspdCircuit(const Circuit &circuit) {
    // TODO
    std::vector<Rectangle> rows;
    Legalizer ret = Legalizer(rows, circuit.cellWidths, circuit.cellX, circuit.cellY);
    // Represent fixed cells with -1 width so they are not considered
    for (int i = 0; i < circuit.nbCells(); ++i) {
        if (circuit.cellFixed[i]) {
            ret.cellWidth_[i] = -1;
        }
    }
    return ret;
}

Legalizer::Legalizer(const std::vector<Rectangle> &rows,
                     const std::vector<int> &width,
                     const std::vector<int> &targetX,
                     const std::vector<int> &targetY)
    : costModel_(LegalizationModel::L1)
    {
    assert (width.size() == targetX.size());
    assert (width.size() == targetY.size());
    // Sort the rows
    rows_ = rows;
    std::stable_sort(rows_.begin(), rows_.end(), [](Rectangle a, Rectangle b) -> bool {
        return a.minY < b.minY || (a.minY == b.minY && a.minX < b.minX);
    });

    // Sort the cells by target X coordinate
    std::vector<std::pair<int, int> > sortedCells;
    for (int i = 0; i < targetX.size(); ++i) {
        sortedCells.emplace_back(targetX[i], i);
    }
    std::stable_sort(sortedCells.begin(), sortedCells.end());
    for (auto p : sortedCells) {
        int c = p.second;
        cellWidth_.push_back(width[c]);
        cellTargetX_.push_back(targetX[c]);
        cellTargetY_.push_back(targetY[c]);
    }
}

bool Legalizer::placeCellOptimally(int cell) {
    /**
     * Very naive algorithm that tries every possible position
     * TODO: Start from the optimal y position and work from there
     */
    int targetX = cellTargetX_[cell];
    int targetY = cellTargetY_[cell];
    int bestX = 0;
    int bestRow = -1;
    long long bestDist = std::numeric_limits<long long>::max();
    for (int row = 0; row < nbRows(); ++row) {
        int y = rows_[row].minY;
        auto [ok, x] = placeCellOptimally(cell, row);
        if (!ok) {
            continue;
        }
        auto dist = norm(x - targetX, y - targetY, costModel_);
        if (bestRow == -1 || dist < bestDist) {
            bestX = x;
            bestRow = row;
            bestDist = dist;
        }
    }
    if (bestRow == -1) {
        return false;
    }
    doPlacement(cell, bestRow, bestX);
    return true;
}

std::pair<bool, int> Legalizer::placeCellOptimally(int cell, int row) const {
    int minCoord;
    if (rowToCells_[row].empty()) {
        minCoord = rows_[row].minX;
    }
    else {
        int pred = rowToCells_[row].back();
        int x = rowToX_[row].back();
        minCoord = x + cellWidth_[pred];
    }

    int maxCoord = rows_[row].maxX - cellWidth_[cell];

    if (minCoord > maxCoord) {
        return std::make_pair(false, 0);
    }

    // Get the best possible placement
    // OPTIMIZE: with LInf cost we could pack further left without degradation
    int coord = cellTargetX_[cell];
    coord = std::min(coord, maxCoord);
    coord = std::max(coord, minCoord);
    return std::make_pair(true, coord);
}

void Legalizer::doPlacement(int cell, int row, int x) {
    rowToCells_[row].push_back(cell);
    rowToX_[row].push_back(x);
    cellToRow_[cell] = row;
    cellToX_[cell] = x;
    cellToY_[cell] = rows_[row].minY;
}

std::vector<int> Legalizer::cellLegalX() const {
    std::vector<int> ret(nbCells());
    for (int r = 0; r < nbRows(); ++r) {
        for (int i = 0; i < rowToCells_[r].size(); ++i) {
            int c = rowToCells_[r][i];
            ret[c] = rowToX_[r][i];
        }
    }
    return ret;
}

std::vector<int> Legalizer::cellLegalY() const {
    std::vector<int> ret(nbCells());
    for (int r = 0; r < nbRows(); ++r) {
        for (int i = 0; i < rowToCells_[r].size(); ++i) {
            int c = rowToCells_[r][i];
            ret[c] = rows_[r].minY;
        }
    }
    return ret;
}

void Legalizer::exportPlacement(Circuit &circuit) {
    std::vector<int> cellX = cellLegalX();
    std::vector<int> cellY = cellLegalY();
    for (int i = 0; i < circuit.nbCells(); ++i) {
        circuit.cellX[i] = cellX[i];
        circuit.cellY[i] = cellY[i];
    }
}
