#pragma once

#include "coloquinte.hpp"

/**
 * Representation of an almost legalized placement, where density constraints are met per bin
 */
class DensityLegalizer {
  public:
    DensityLegalizer(Rectangle area, int nbCells);
    DensityLegalizer(const Circuit &circuit);

    // Accessors
    int nbCells() const { return nbCells_; }
    int nbBins() const { return nbBinsX_ * nbBinsY_; }
    int nbBinsX() const { return nbBinsX_; }
    int nbBinsY() const { return nbBinsY_; }

    long long totalCapacity() const;
    long long totalDemand() const;

    const std::vector<int> &cellDemand() const { return cellDemand_; }
    const std::vector<float> &cellTargetX() const { return cellTargetX_; }
    const std::vector<float> &cellTargetY() const { return cellTargetY_; }

    const std::vector<std::vector<long long> > &binCapacity() const { return binCapacity_; }
    std::vector<std::vector<long long> > binUsage() const;

    const std::vector<float> &binX() const { return binX_; }
    const std::vector<float> &binY() const { return binY_; }

    // Problem updates
    void updateBins(int binsX, int binsY);
    void updateCellDemand(std::vector<int> cellDemand) { cellDemand_ = cellDemand; }
    void updateCellTargetX(std::vector<float> cellTargetX) { cellTargetX_ = cellTargetX; }
    void updateCellTargetY(std::vector<float> cellTargetY) { cellTargetY_ = cellTargetY; }

    // How much capacity is overflowing the bins
    long long totalOverflow() const;
    float meanOverflow() const;
    float rmsOverflow() const;

    // Quality metrics for the legalization
    float metrics(LegalizationModel model) const;
    float distL1() const;
    float distL2() const;
    float distLInf() const;
    float distL2Squared() const;

    // Algorithms to create an initial legalization from scratch
    void assign();
    void bisect(LegalizationModel model);

    // Algorithms to improve the current legalization
    void redistributeBisect(LegalizationModel model);
    void redistributeX(LegalizationModel model);
    void redistributeY(LegalizationModel model);

    // Obtain coordinates for the cells from the legalization
    std::vector<float> simpleCoordX() const;
    std::vector<float> simpleCoordY() const;
    // TODO: coordinates with spreading
    std::vector<float> spreadCoordX() const;
    std::vector<float> spreadCoordY() const;

    void check() const;
    void report() const;

    // Utilities
    std::vector<int> allCells() const;

    // Bisection algorithm
    std::pair<std::vector<int>, std::vector<int> > bisect(float cx1, float cy1, float cx2, float cy2, long long capa1, long long capa2, std::vector<int> cells, LegalizationModel leg) const;

  private:
    // Redo the bisection for these bins
    void rebisect(int x1, int y1, int x2, int y2, LegalizationModel leg);

    // Bisection algorithm helpers
    std::vector<std::pair<float, int> > computeCellCosts(float cx1, float cy1, float cx2, float cy2, std::vector<int> cells, LegalizationModel leg) const;
    int findIdealSplitPos(const std::vector<std::pair<float, int> > &cellCosts) const;
    int findConstrainedSplitPos(const std::vector<std::pair<float, int> > &cellCosts, int targetPos, long long capa1, long long capa2) const;
    std::pair<std::vector<int>, std::vector<int> > doSplit(const std::vector<std::pair<float, int> > &cellCosts, int ind) const;

  private:
    // Problem data
    Rectangle area_;
    int nbCells_;
    int nbBinsX_;
    int nbBinsY_;

    // Bin properties
    std::vector<float> binX_;
    std::vector<float> binY_;
    std::vector<int> binLimitX_;
    std::vector<int> binLimitY_;
    std::vector<std::vector<long long> > binCapacity_;

    // Cell properties
    std::vector<int> cellDemand_;
    std::vector<float> cellTargetX_;
    std::vector<float> cellTargetY_;

    // Problem status
    std::vector<std::vector<std::vector<int> > > binCells_;
};



