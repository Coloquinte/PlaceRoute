
#include "place_global/density_legalizer.hpp"

#include <cassert>
#include <cmath>
#include <algorithm>

float DensityLegalizer::distance(float x1, float y1, float x2, float y2, LegalizationModel leg) {
    float dx = x1 - x2;
    float dy = y1 - y2;
    switch (leg) {
        case LegalizationModel::L1:
        return std::abs(dx) + std::abs(dy);
        case LegalizationModel::L2:
        return std::sqrt(dx * dx + dy * dy);
        case LegalizationModel::LInf:
        return std::max(std::abs(dx), std::abs(dy));
        default:
        return dx * dx + dy * dy;
    }
}

DensityLegalizer::DensityLegalizer(Rectangle area, int nbCells) : area_(area) {
    nbCells_ = nbCells;
    cellDemand_.assign(nbCells, 0LL);
    cellTargetX_.assign(nbCells, 0.0f);
    cellTargetY_.assign(nbCells, 0.0f);
    updateBins(1, 1);
    // Dumb initial placement
    std::vector<int> pl;
    for (int i = 0; i < nbCells; ++i) {
        pl.push_back(i);
    }
    binCells_.emplace_back();
    binCells_[0].push_back(pl);
    check();
}

void DensityLegalizer::updateBins(int binsX, int binsY) {
    nbBinsX_ = binsX;
    nbBinsY_ = binsY;
    for (int i = 0; i < binsX; ++i) {
        float n = 2 * i + 1;
        binX_.push_back((n * area_.minX + (2 * binsX - n) * area_.maxX) / (2 * binsX));
    }
    for (int i = 0; i < binsY; ++i) {
        float n = 2 * i + 1;
        binY_.push_back((n * area_.minY + (2 * binsY - n) * area_.maxY) / (2 * binsY));
    }
    // TODO: correctly subdivide the capacity
    long long capacity = area_.area() / (binsX * binsY);
    binCapacity_.assign(binsX, std::vector<long long>(binsY, capacity));
}

float DensityLegalizer::distL1() const {
    float disp = 0.0f;
    for (int i = 0; i < nbBinsX_; ++i) {
        for (int j = 0; j < nbBinsY_; ++j) {
            float x = binX_[i];
            float y = binY_[i];
            for (int c : binCells_[i][j]) {
                float dx = cellTargetX_[c] - x;
                float dy = cellTargetY_[c] - y;
                disp += cellDemand_[c] * (std::abs(dx) + std::abs(dy));
            }
        }
    }
    return disp;
}

float DensityLegalizer::distLInf() const {
    float disp = 0.0f;
    for (int i = 0; i < nbBinsX_; ++i) {
        for (int j = 0; j < nbBinsY_; ++j) {
            float x = binX_[i];
            float y = binY_[i];
            for (int c : binCells_[i][j]) {
                float dx = cellTargetX_[c] - x;
                float dy = cellTargetY_[c] - y;
                disp += cellDemand_[c] * std::max(std::abs(dx), std::abs(dy));
            }
        }
    }
    return disp;
}

float DensityLegalizer::distL2Squared() const {
    float disp = 0.0f;
    for (int i = 0; i < nbBinsX_; ++i) {
        for (int j = 0; j < nbBinsY_; ++j) {
            float x = binX_[i];
            float y = binY_[i];
            for (int c : binCells_[i][j]) {
                float dx = cellTargetX_[c] - x;
                float dy = cellTargetY_[c] - y;
                disp += cellDemand_[c] * (dx * dx + dy * dy);
            }
        }
    }
    return disp;
}

float DensityLegalizer::distL2() const {
    return std::sqrt(distL2Squared());
}

void DensityLegalizer::check() const {
    assert (binCapacity_.size() == nbBinsX());
    assert (binCells_.size() == nbBinsX());
    for (const auto &bc : binCapacity_) {
        assert (bc.size() == nbBinsY());
    }
    for (const auto &bc : binCells_) {
        assert (bc.size() == nbBinsY());
    }
    assert (binX_.size() == nbBinsX());
    assert (binY_.size() == nbBinsY());
    assert (cellDemand_.size() == nbCells());
    assert (cellTargetX_.size() == nbCells());
    assert (cellTargetY_.size() == nbCells());
    // Check that every cell is present exactly once
    std::vector<char> placed(nbCells());
    for (int i = 0; i < nbBinsX_; ++i) {
        for (int j = 0; j < nbBinsY_; ++j) {
            for (int c : binCells_[i][j]) {
                assert (!placed[c]);
                placed[c] = true;
            }
        }
    }
    for (int c = 0; c < nbCells(); ++c) {
        assert (placed[c]);
    }
}

std::vector<std::vector<long long> > DensityLegalizer::binUsage() const {
    std::vector<std::vector<long long> > ret(nbBinsX(), std::vector<long long>(nbBinsY()));
    for (int i = 0; i < nbBinsX_; ++i) {
        for (int j = 0; j < nbBinsY_; ++j) {
            for (int c : binCells_[i][j]) {
                ret[i][j] += cellDemand_[c];
            }
        }
    }
    return ret;
}

long long DensityLegalizer::totalCapacity() const {
    long long ret = 0;
    for (int i = 0; i < nbBinsX_; ++i) {
        for (int j = 0; j < nbBinsY_; ++j) {
            ret += binCapacity_[i][j];
        }
    }
    return ret;
}

long long DensityLegalizer::totalDemand() const {
    long long ret = 0;
    for (int demand : cellDemand_) {
        ret += demand;
    }
    return ret;
}

long long DensityLegalizer::totalOverflow() const {
    long long ret = 0;
    auto usage = binUsage();
    for (int i = 0; i < nbBinsX_; ++i) {
        for (int j = 0; j < nbBinsY_; ++j) {
            ret += std::max(usage[i][j] - binCapacity_[i][j], 0LL);
        }
    }
    return ret;
}

float DensityLegalizer::meanOverflow() const {
    return totalOverflow() / (float) totalDemand();
}

float DensityLegalizer::rmsOverflow() const {
    auto usage = binUsage();
    float fact = 1.0f / totalDemand();
    float ret = 0.0f;
    for (int i = 0; i < nbBinsX_; ++i) {
        for (int j = 0; j < nbBinsY_; ++j) {
            float ovf = std::max(usage[i][j] - binCapacity_[i][j], 0LL) * fact;
            ret += ovf * ovf;
        }
    }
    return std::sqrt(ret);
}

std::vector<float> DensityLegalizer::simpleCoordX() const {
    std::vector<float> ret(nbCells(), 0.0f);
    for (int i = 0; i < nbBinsX_; ++i) {
        for (int j = 0; j < nbBinsY_; ++j) {
            for (int c : binCells_[i][j]) {
                ret[c] = binX_[i];
            }
        }
    }
    return ret;
}

std::vector<float> DensityLegalizer::simpleCoordY() const {
    std::vector<float> ret(nbCells(), 0.0f);
    for (int i = 0; i < nbBinsY_; ++i) {
        for (int j = 0; j < nbBinsX_; ++j) {
            for (int c : binCells_[i][j]) {
                ret[c] = binY_[i];
            }
        }
    }
    return ret;
}

std::vector<std::pair<float, int> > DensityLegalizer::computeCellCosts(float cx1, float cy1, float cx2, float cy2, std::vector<int> cells, LegalizationModel leg) const {
    // TODO: always add a secondary objective using squared distance
    std::vector<std::pair<float, int> > cellCosts;
    for (int c : cells) {
        float x = cellTargetX_[c];
        float y = cellTargetY_[c];
        float cost = distance(x, y, cx1, cy1, leg) - distance(x, y, cx2, cy2, leg);
        cellCosts.emplace_back(cost, c);
    }
    std::sort(cellCosts.begin(), cellCosts.end(), [](std::pair<float, int> a, std::pair<float, int> b) { return a.first < b.first; });
    return cellCosts;
}

std::pair<std::vector<int>, std::vector<int> > DensityLegalizer::doSplit(const std::vector<std::pair<float, int> > &cellCosts, int ind) const {
    std::pair<std::vector<int>, std::vector<int> > ret;
    for (int i = 0; i < ind; ++i) {
        ret.first.push_back(cellCosts[i].second);
    }
    for (int i = ind; i < cellCosts.size(); ++i) {
        ret.second.push_back(cellCosts[i].second);
    }
    return ret;
}

int DensityLegalizer::findIdealSplitPos(const std::vector<std::pair<float, int> > &cellCosts) const {
    // Find the ideal split position
    int splitInd = 0;
    for (; splitInd < cellCosts.size(); ++splitInd) {
        if (cellCosts[splitInd].first > 0.0) {
            break;
        }
    }
    return splitInd;
}

int DensityLegalizer::findConstrainedSplitPos(const std::vector<std::pair<float, int> > &cellCosts, int targetPos, long long capa1, long long capa2) const {
    long long demand1 = 0;
    long long demand2 = 0;
    for (int i = 0; i < targetPos; ++i) {
        demand1 += cellDemand_[cellCosts[i].second];
    }
    for (int i = targetPos; i < cellCosts.size(); ++i) {
        demand2 += cellDemand_[cellCosts[i].second];
    }
    int splitPos = targetPos;
    // Remove from the left if overflowed
    while (splitPos > 0 && demand1 - capa1 > 0 && demand1 - capa1 > demand2 - capa2) {
        int dem = cellDemand_[cellCosts[splitPos-1].second];
        demand1 -= dem;
        demand2 += dem;
        --splitPos;
    }
    // Remove from the right if overflowed
    while (splitPos < cellCosts.size() && demand2 - capa2 > 0 && demand2 - capa2 > demand1 - capa1) {
        int dem = cellDemand_[cellCosts[splitPos].second];
        demand2 -= dem;
        demand1 += dem;
        ++splitPos;
    }
    return splitPos;
}

std::pair<std::vector<int>, std::vector<int> > DensityLegalizer::bisect(float cx1, float cy1, float cx2, float cy2, long long capa1, long long capa2, std::vector<int> cells, LegalizationModel leg) {
    std::vector<std::pair<float, int> > cellCosts = computeCellCosts(cx1, cy1, cx2, cy2, cells, leg);
    int idealSplitPos = findIdealSplitPos(cellCosts);
    int splitPos = findConstrainedSplitPos(cellCosts, idealSplitPos, capa1, capa2);
    return doSplit(cellCosts, splitPos);
}

void DensityLegalizer::rebisect(int x1, int y1, int x2, int y2, LegalizationModel leg) {
    if (x1 == x2 && y1 == y2)
        return;

    // Get all cells
    std::vector<int> cells;
    cells.insert(cells.end(), binCells_[x1][y1].begin(), binCells_[x1][y1].end());
    cells.insert(cells.end(), binCells_[x2][y2].begin(), binCells_[x2][y2].end());

    auto b = bisect(
        binX_[x1], binY_[y1], binX_[x2], binY_[y2],
        binCapacity_[x1][y1], binCapacity_[x2][y2],
        cells, leg
    );
    binCells_[x1][y1] = b.first;
    binCells_[x2][y2] = b.second;
}
