#pragma once

#include "coloquinte.hpp"

#include <vector>

#include <xtensor/xarray.hpp>
#include <xtensor/xtensor.hpp>
#include <eigen3/Eigen/Sparse>


/**
 * Representation of nets of a given size as 1D HPWL for optimization algorithms.
 *
 * This allows us to use more efficient array operations and specialize for fixed terminals
 */
class NetTopologyFixedSize {
  public:
    NetTopologyFixedSize(int nbCells, int netSize, std::vector<int> pinCells, std::vector<float> pinOffsets);
    NetTopologyFixedSize(int nbCells, xt::xtensor<int, 2> pinCells, xt::xtensor<float, 2> pinOffsets);

    int nbCells() const { return nbCells_; }
    int nbNets() const { return nbNets_; }
    int netSize() const { return netSize_; }

    float valueHPWL(xt::xtensor<float, 1> pl) const;
    float valueLSE(xt::xtensor<float, 1> pl, float epsilon) const;
    float valueWA(xt::xtensor<float, 1> pl, float epsilon) const;

    xt::xtensor<float, 1> gradHPWL(xt::xtensor<float, 1> pl) const;
    xt::xtensor<float, 1> gradLSE(xt::xtensor<float, 1> pl, float epsilon) const;
    xt::xtensor<float, 1> gradWA(xt::xtensor<float, 1> pl, float epsilon) const;
    xt::xtensor<float, 1> proximalStep(xt::xtensor<float, 1> pl, float step) const;

    // Cells of the pins for each net
    const xt::xtensor<int, 2> &pinCells() const { return pinCells_; }
    // Offsets of the pins for each net
    const xt::xtensor<float, 2> &pinOffsets() const { return offsets_; }
    // Coordinates of the pins for a given placement
    xt::xtensor<float, 2> pinCoords(xt::xtensor<float, 1> pl) const;

    void check() const;

  protected:
    // Compute the LSE approximation of max on a 2D array
    xt::xtensor<float, 1> lseMax(xt::xtensor<float, 2> pins, float epsilon) const;
    // Compute the WA approximation of max on a 2D array
    xt::xtensor<float, 1> waMax(xt::xtensor<float, 2> pins, float epsilon) const;
    // Compute the gradient of max on a 2D array
    xt::xtensor<float, 2> maxGrad(xt::xtensor<float, 2> pins) const;
    // Compute the gradient of the LSE approximation of max on a 2D array
    xt::xtensor<float, 2> lseMaxGrad(xt::xtensor<float, 2> pins, float epsilon) const;
    // Compute the gradient of the WA approximation of max on a 2D array
    xt::xtensor<float, 2> waMaxGrad(xt::xtensor<float, 2> pins, float epsilon) const;
    // Compute the proximal step of max on a 2D array, like a gradient step
    xt::xtensor<float, 2> maxProximal(xt::xtensor<float, 2> pins, float step) const;
    xt::xtensor<float, 2> maxProximal(xt::xtensor<float, 2> pins, xt::xtensor<float, 1> fixedPins, float step) const;
    // Compute the position after proximal step of max on a 2D array
    xt::xtensor<float, 2> posProximal(xt::xtensor<float, 2> pins, float step) const;

    // Transform a 2D grad tensor (for each pin) to a 1D tensor for each cell
    xt::xtensor<float, 1> toCellGrad(xt::xtensor<float, 2> pinGrad) const;

  protected:
    int nbCells_;
    int nbNets_;
    int netSize_;

    xt::xtensor<int, 2> pinCells_;
    xt::xtensor<float, 2> offsets_;
};

/**
 * Representation of nets of a given size as 1D HPWL for optimization algorithms.
 */
class NetTopologyFixedSizeTerminals : NetTopologyFixedSize {
  public:
    NetTopologyFixedSizeTerminals(int nbCells, int netSize, std::vector<int> pinCells, std::vector<float> pinOffsets, std::vector<float> minPos, std::vector<float> maxPos);
    NetTopologyFixedSizeTerminals(int nbCells, xt::xtensor<int, 2> pinCells, xt::xtensor<float, 2> pinOffsets, xt::xtensor<float, 1> minPos, xt::xtensor<float, 1> maxPos);

    int nbCells() const { return nbCells_; }
    int nbNets() const { return nbNets_; }
    int netSize() const { return netSize_; }
    const xt::xtensor<int, 2> &pinCells() const { return pinCells_; }
    const xt::xtensor<float, 2> &pinOffsets() const { return offsets_; }
    const xt::xtensor<float, 1> &minPins() const { return minPins_; }
    const xt::xtensor<float, 1> &maxPins() const { return maxPins_; }

    float valueHPWL(xt::xtensor<float, 1> pl) const;
    float valueLSE(xt::xtensor<float, 1> pl, float epsilon) const;
    float valueWA(xt::xtensor<float, 1> pl, float epsilon) const;

    xt::xtensor<float, 1> gradHPWL(xt::xtensor<float, 1> pl) const;
    xt::xtensor<float, 1> gradLSE(xt::xtensor<float, 1> pl, float epsilon) const;
    xt::xtensor<float, 1> gradWA(xt::xtensor<float, 1> pl, float epsilon) const;
    xt::xtensor<float, 1> proximalStep(xt::xtensor<float, 1> pl, float step) const;

    // Pair of arrays, with min and max pins included respectively
    std::pair<xt::xtensor<float, 2>, xt::xtensor<float, 2> > pinCoordsMinMax(xt::xtensor<float, 1> pl) const;
    // Single array with min and max pins appended
    xt::xtensor<float, 2> pinCoordsAll(xt::xtensor<float, 1> pl) const;
    // Cells of the pins for each net, with -1 appended for min and max pins
    xt::xtensor<int, 2> pinCellsAll() const;
    // Offsets of the pins for each net, including min and max pins
    xt::xtensor<float, 2> pinOffsetsAll() const;
    // Transform a 2D grad tensor (for each pin) to a 1D tensor for each cell
    xt::xtensor<float, 1> toCellGradMinMax(xt::xtensor<float, 2> pinGrad) const;

    void check() const;

  protected:
    xt::xtensor<float, 1> minPins_;
    xt::xtensor<float, 1> maxPins_;
};

/**
 * Representation of the nets as 1D HPWL for optimization algorithms
 */
class NetTopology {
  public:
    // Constructors
    static NetTopology xTopology(const Circuit &);
    static NetTopology yTopology(const Circuit &);

    static NetTopology fromData(const std::vector<int> &cellSizes, const std::vector<int> & pl, const std::vector<char> &cellFixed, const std::vector<int> &netLimits, const std::vector<int> &pinCells, const std::vector<int> &pinOffsets);

    // Basic accessors
    int nbCells() const { return nbCells_; }
    int nbNets() const { return nbNets_; }

    const std::vector<NetTopologyFixedSize> &nets() const { return nets_; }
    const std::vector<NetTopologyFixedSizeTerminals> &terminalNets() const { return terminalNets_; }

    // Objective value computation
    float valueHPWL(const xt::xtensor<float, 1> &pl) const;
    float valueLSE(const xt::xtensor<float, 1> &pl, float epsilon) const;
    float valueWA(const xt::xtensor<float, 1> &pl, float epsilon) const;

    // Gradient and descent steps
    xt::xtensor<float, 1> gradHPWL(const xt::xtensor<float, 1> &pl) const;
    xt::xtensor<float, 1> gradLSE(const xt::xtensor<float, 1> &pl, float epsilon) const;
    xt::xtensor<float, 1> gradWA(const xt::xtensor<float, 1> &pl, float epsilon) const;
    xt::xtensor<float, 1> proximalStep(const xt::xtensor<float, 1> &pl, float step) const;
    xt::xtensor<float, 1> starStep(const xt::xtensor<float, 1> &pl) const;
    xt::xtensor<float, 1> b2bStep(const xt::xtensor<float, 1> &pl) const;

    // Direct solutions
    xt::xtensor<float, 1> starSolve() const;
    xt::xtensor<float, 1> starSolve(const xt::xtensor<float, 1> &pl) const;
    xt::xtensor<float, 1> b2bSolve(const xt::xtensor<float, 1> &pl) const;

    // Consistency checker
    void check() const;

  private:
    int nbCells_;
    int nbNets_;

    std::vector<NetTopologyFixedSize> nets_;
    std::vector<NetTopologyFixedSizeTerminals> terminalNets_;

    friend class NetTopologyBuilder;
};

class NetTopologyFixedSizeBuilder {
  public:
    NetTopologyFixedSizeBuilder(int nbCells, int netSize) : nbCells_(nbCells), netSize_(netSize) {}
    bool empty() const { return pinCells_.empty(); }
    void push(const std::vector<int> &cells, const std::vector<float> &offsets);
    NetTopologyFixedSize build() const { return NetTopologyFixedSize(nbCells_, netSize_, pinCells_, pinOffsets_); }

  private:
    int nbCells_;
    int netSize_;
    std::vector<int> pinCells_;
    std::vector<float> pinOffsets_;
};

class NetTopologyFixedSizeTerminalsBuilder {
  public:
    NetTopologyFixedSizeTerminalsBuilder(int nbCells, int netSize) : nbCells_(nbCells), netSize_(netSize) {}
    bool empty() const { return pinCells_.empty(); }
    void push(const std::vector<int> &cells, const std::vector<float> &offsets, float minPin, float maxPin);
    NetTopologyFixedSizeTerminals build() const { return NetTopologyFixedSizeTerminals(nbCells_, netSize_, pinCells_, pinOffsets_, minPins_, maxPins_); }

  private:
    int nbCells_;
    int netSize_;
    std::vector<int> pinCells_;
    std::vector<float> pinOffsets_;
    std::vector<float> minPins_;
    std::vector<float> maxPins_;
};

class NetTopologyBuilder {
  public:
    NetTopologyBuilder(int nbCells) : nbCells_(nbCells), nbNets_(0) {}
    void push(const std::vector<int> &cells, const std::vector<float> &offsets);
    void push(const std::vector<int> &cells, const std::vector<float> &offsets, float minPos, float maxPos);
    NetTopology build() const;

  private:
    int nbCells_;
    int nbNets_;
    std::vector<NetTopologyFixedSizeBuilder> netBuilders_;
    std::vector<NetTopologyFixedSizeTerminalsBuilder> terminalNetBuilders_;
};

/**
 * Matrix creation for quadratic approximations
 */
class MatrixBuilder {
  public:
    MatrixBuilder(int nbCells) : nbCells_(nbCells), nbSupps_(0), rhs_(nbCells), initial_(nbCells) {}
    int matSize() const { return nbCells_ + nbSupps_; }

    xt::xtensor<float, 1> solve();
    const std::vector<Eigen::Triplet<float> > &mat() const { return mat_; }
    const std::vector<float> &rhs() const { return rhs_; }
    const std::vector<float> &initial() const { return initial_; }

    // Simple star model, without taking pin positions into account
    static MatrixBuilder createStar(const NetTopology &topo);
    void extendStar(const NetTopologyFixedSize &topo);
    void extendStar(const NetTopologyFixedSizeTerminals &topo);

    // Star model, with local approximation from the current placement
    static MatrixBuilder createStar(const NetTopology &topo, xt::xtensor<float, 1> pl);
    void extendStar(const NetTopologyFixedSize &topo, xt::xtensor<float, 1> pl);
    void extendStar(const NetTopologyFixedSizeTerminals &topo, xt::xtensor<float, 1> pl);

    // Bound to bound model, with local approximation from the current placement
    static MatrixBuilder createB2B(const NetTopology &topo, xt::xtensor<float, 1> pl);
    void extendB2B(const NetTopologyFixedSize &topo, xt::xtensor<float, 1> pl);
    void extendB2B(const NetTopologyFixedSizeTerminals &topo, xt::xtensor<float, 1> pl);

  private:
    void addPin(int c1, int c2, float offs1, float offs2, float weight);
    void addPinOrFixed(int c1, int c2, float offs1, float offs2, float weight);
    void addFixedPin(int c1, float offs1, float pos, float weight);

  private:
    int nbCells_;
    int nbSupps_;
    std::vector<Eigen::Triplet<float> > mat_;
    std::vector<float> rhs_;
    std::vector<float> initial_;
};


