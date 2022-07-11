
#include "place_detailed/legalizer.hpp"

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
                     const std::vector<int> &widths,
                     const std::vector<int> &targetX,
                     const std::vector<int> &targetY)
    : rows_(rows)
    , cellWidth_(widths)
    , cellTargetX_(targetX)
    , cellTargetY_(targetY)
    , costModel_(LegalizationModel::L1)
    {}

void Legalizer::placeCell(int cell) {
    // TODO
}

std::pair<bool, long long> Legalizer::placeCell(int cell, int row) {
    // TODO
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
