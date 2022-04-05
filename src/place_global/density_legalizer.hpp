#pragma once

#include "coloquinte.hpp"

struct Rectangle {
    Rectangle(int minX, int maxX, int minY, int maxY) : minX(minX), maxX(maxX), minY(minY), maxY(maxY) {}
    int minX;
    int maxX;
    int minY;
    int maxY;

    int width() const { return maxX - minX; }
    int height() const { return maxY - minY; }
    long long area() const { return (long long) width() * (long long) height(); }
};


enum LegalizationModel { L1, L2, LInf, L2Squared };

/**
 * Representation of an almost legalized placement, where density constraints are met per bin
 */
class DensityLegalizer {
  public:
    DensityLegalizer(Rectangle area, int binsX, int binsY, std::vector<int> cellDemand, std::vector<float> cellTargetX, std::vector<float> cellTargetY);

    int nbCells() const { return cellDemand_.size(); }
    int nbBins() const { return nbBinsX_ * nbBinsY_; }
    int nbBinsX() const { return nbBinsX_; }
    int nbBinsY() const { return nbBinsY_; }
    long long totalCapacity() const;
    long long totalDemand() const;

    void updateTargetX(std::vector<float> cellTargetX) { cellTargetX_ = cellTargetX; }
    void updateTargetY(std::vector<float> cellTargetY) { cellTargetY_ = cellTargetY; }

    // How much capacity is overflowing the bins
    long long totalOverflow() const;
    float meanOverflow() const;
    float rmsOverflow() const;

    // Quality metrics for the legalization
    float distL1() const;
    float distL2() const;
    float distLInf() const;
    float distL2Squared() const;

    // Algorithms to create an initial legalization from scratch
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

  private:
    // Compute the current usage for each bin
    std::vector<std::vector<long long> > binUsage() const;

    // Redo the bisection for these bins
    void rebisect(int x1, int y1, int x2, int y2, LegalizationModel leg);

    // Bisection algorithm and helpers
    std::pair<std::vector<int>, std::vector<int> > bisect(float cx1, float cy1, float cx2, float cy2, long long capa1, long long capa2, std::vector<int> cells, LegalizationModel leg);

    std::vector<std::pair<float, int> > computeCellCosts(float cx1, float cy1, float cx2, float cy2, std::vector<int> cells, LegalizationModel leg) const;
    int findIdealSplitPos(const std::vector<std::pair<float, int> > &cellCosts) const;
    int findConstrainedSplitPos(const std::vector<std::pair<float, int> > &cellCosts, int targetPos, long long capa1, long long capa2) const;
    std::pair<std::vector<int>, std::vector<int> > doSplit(const std::vector<std::pair<float, int> > &cellCosts, int ind) const;

  private:
    // Problem data
    int nbBinsX_;
    int nbBinsY_;
    std::vector<std::vector<long long> > binCapacity_;
    std::vector<float> binX_;
    std::vector<float> binY_;
    std::vector<int> cellDemand_;
    std::vector<float> cellTargetX_;
    std::vector<float> cellTargetY_;

    // Problem status
    std::vector<std::vector<std::vector<int> > > binCells_;
};



