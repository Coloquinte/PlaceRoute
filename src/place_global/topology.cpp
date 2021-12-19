
#include "place_global/topology.hpp"

#include <xtensor/xarray.hpp>
#include <xtensor/xtensor.hpp>
#include <xtensor/xadapt.hpp>
#include <xtensor/xview.hpp>
#include <xtensor/xindex_view.hpp>
#include <xtensor/xstrided_view.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xsort.hpp>

#include <algorithm>
#include <cmath>
#include <cassert>
#include <iostream>

NetTopology NetTopology::xTopology(const Circuit &circuit) {
    return fromData(circuit.cellWidths, circuit.cellX, circuit.cellFixed, circuit.netLimits, circuit.pinCells, circuit.pinXOffsets);
}

NetTopology NetTopology::yTopology(const Circuit &circuit) {
    return fromData(circuit.cellHeights, circuit.cellY, circuit.cellFixed, circuit.netLimits, circuit.pinCells, circuit.pinYOffsets);
}

NetTopology NetTopology::fromData(const std::vector<int> &cellSizes, const std::vector<int> & pl, const std::vector<char> &cellFixed, const std::vector<int> &netLimits, const std::vector<int> &pinCells, const std::vector<int> &pinOffsets) {
    std::vector<NetTopologyFixedSizeBuilder> netBuilders_;
    std::vector<NetTopologyFixedSizeTerminalsBuilder> terminalNetBuilders_;
    int nbCells = cellSizes.size();
    int nbNets = 0;
    for (int i = 0; i + 1 < netLimits.size(); ++i) {
        float minPos = std::numeric_limits<float>::infinity();
        float maxPos = -std::numeric_limits<float>::infinity();
        int b = netLimits[i];
        int e = netLimits[i+1];
        std::vector<int> cells;
        std::vector<float> offsets;
        for (int j = b; j < e; ++j) {
            int cell = pinCells[j];
            int offset = pinOffsets[j];
            if (cellFixed[cell]) {
                int pos = pl[cell] + offset;
                minPos = std::min(minPos, (float) pos);
                maxPos = std::max(maxPos, (float) pos);
            }
            else {
                cells.push_back(cell);
                // Offset to center of cell instead of lower-left
                offsets.push_back(offset - 0.5f * cellSizes[cell]);
            }
        }
        int sz = cells.size();
        if (std::isfinite(minPos)) {
            if (cells.empty()) continue;
            while (terminalNetBuilders_.size() < sz + 1) {
                terminalNetBuilders_.emplace_back(nbCells, terminalNetBuilders_.size());
            }
            terminalNetBuilders_[sz].push(cells, offsets, minPos, maxPos);
        }
        else {
            if (cells.size() <= 1) continue;
            while (netBuilders_.size() < sz + 1) {
                netBuilders_.emplace_back(nbCells, netBuilders_.size());
            }
            netBuilders_[sz].push(cells, offsets);
        }
        ++nbNets;
    }

    NetTopology ret;
    ret.nbCells_ = nbCells;
    ret.nbNets_ = nbNets;
    for (const auto &bd : netBuilders_) {
        if (!bd.empty()) {
            ret.nets_.push_back(bd.build());
        }
    }
    for (const auto &bd : terminalNetBuilders_) {
        if (!bd.empty()) {
            ret.terminalNets_.push_back(bd.build());
        }
    }
    return ret;
}

float NetTopology::valueHPWL(const xt::xtensor<float, 1> &pl) const {
    float val = 0.0;
    for (const auto &bd : nets_) {
        val += bd.valueHPWL(pl);
    }
    for (const auto &bd : terminalNets_) {
        val += bd.valueHPWL(pl);
    }
    return val;
}

float NetTopology::valueLSE(const xt::xtensor<float, 1> &pl, float epsilon) const {
    float val = 0.0;
    for (const auto &bd : nets_) {
        val += bd.valueLSE(pl, epsilon);
    }
    for (const auto &bd : terminalNets_) {
        val += bd.valueLSE(pl, epsilon);
    }
    return val;
}

float NetTopology::valueWA(const xt::xtensor<float, 1> &pl, float epsilon) const {
    float val = 0.0;
    for (const auto &bd : nets_) {
        val += bd.valueWA(pl, epsilon);
    }
    for (const auto &bd : terminalNets_) {
        val += bd.valueWA(pl, epsilon);
    }
    return val;
}

xt::xtensor<float, 1> NetTopology::gradHPWL(const xt::xtensor<float, 1> &pl) const {
    xt::xtensor<float, 1> val = xt::zeros<float>({nbCells()});
    for (const auto &bd : nets_) {
        val += bd.gradHPWL(pl);
    }
    for (const auto &bd : terminalNets_) {
        val += bd.gradHPWL(pl);
    }
    return val;
}

xt::xtensor<float, 1> NetTopology::gradLSE(const xt::xtensor<float, 1> &pl, float epsilon) const {
    xt::xtensor<float, 1> val = xt::zeros<float>({nbCells()});
    for (const auto &bd : nets_) {
        val += bd.gradLSE(pl, epsilon);
    }
    for (const auto &bd : terminalNets_) {
        val += bd.gradLSE(pl, epsilon);
    }
    return val;
}

xt::xtensor<float, 1> NetTopology::gradWA(const xt::xtensor<float, 1> &pl, float epsilon) const {
    xt::xtensor<float, 1> val = xt::zeros<float>({nbCells()});
    for (const auto &bd : nets_) {
        val += bd.gradWA(pl, epsilon);
    }
    for (const auto &bd : terminalNets_) {
        val += bd.gradWA(pl, epsilon);
    }
    return val;
}

xt::xtensor<float, 1> NetTopology::proximalStep(const xt::xtensor<float, 1> &pl, float step) const {
    xt::xtensor<float, 1> val = xt::zeros<float>({nbCells()});
    for (const auto &bd : nets_) {
        val += bd.proximalStep(pl, step);
    }
    for (const auto &bd : terminalNets_) {
        val += bd.proximalStep(pl, step);
    }
    return val;
}

NetTopologyFixedSize::NetTopologyFixedSize(int nbCells, int netSize, std::vector<int> pinCells, std::vector<float> pinOffsets) {
    assert (pinCells.size() == pinOffsets.size());
    assert (pinCells.size() % netSize == 0);
    nbNets_ = pinCells.size() / netSize;
    nbCells_ = nbCells;
    netSize_ = netSize;
    pinCells_ = xt::adapt(pinCells, {nbNets_, netSize_});
    offsets_ = xt::adapt(pinOffsets, {nbNets_, netSize_});
    check();
}

NetTopologyFixedSize::NetTopologyFixedSize(int nbCells, xt::xtensor<int, 2> pinCells, xt::xtensor<float, 2> pinOffsets) {
    assert (pinCells.shape()[0] == pinOffsets.shape()[0]);
    assert (pinCells.shape()[1] == pinOffsets.shape()[1]);
    nbCells_ = nbCells;
    nbNets_ = pinCells.shape()[0];
    netSize_ = pinCells.shape()[1];
    pinCells_ = pinCells;
    offsets_ = pinOffsets;
    check();
}

void NetTopologyFixedSize::check() const {
    assert (pinCells_.shape()[0] == nbNets_);
    assert (offsets_.shape()[0] == nbNets_);
    assert (pinCells_.shape()[1] == netSize_);
    assert (offsets_.shape()[1] == netSize_);
    assert (xt::all(pinCells_ < nbCells()));
    assert (xt::all(pinCells_ >= 0));
}

xt::xtensor<float, 2> NetTopologyFixedSize::pinCoords(xt::xtensor<float, 1> pl) const {
    assert (pl.size() == nbCells());
    xt::xtensor<float, 2> cellCoords = xt::reshape_view(
        xt::index_view(pl, pinCells_),
        {nbNets(), netSize()}
    );
    return cellCoords + offsets_;
}

float NetTopologyFixedSize::valueHPWL(xt::xtensor<float, 1> pl) const {
    auto pins = pinCoords(pl);
    return xt::sum(xt::amax(pins, {1}) - xt::amin(pins, {1}))();
}

float NetTopologyFixedSize::valueLSE(xt::xtensor<float, 1> pl, float epsilon) const {
    auto pins = pinCoords(pl);
    return xt::sum(lseMax(pins, epsilon) + lseMax(-pins, epsilon))();
}

float NetTopologyFixedSize::valueWA(xt::xtensor<float, 1> pl, float epsilon) const {
    auto pins = pinCoords(pl);
    return xt::sum(waMax(pins, epsilon) + waMax(-pins, epsilon))();
}

xt::xtensor<float, 1> NetTopologyFixedSize::gradHPWL(xt::xtensor<float, 1> pl) const {
    auto pins = pinCoords(pl);
    return toCellGrad(maxGrad(pins) - maxGrad(-pins));
}

xt::xtensor<float, 1> NetTopologyFixedSize::gradLSE(xt::xtensor<float, 1> pl, float epsilon) const {
    auto pins = pinCoords(pl);
    return toCellGrad(lseMaxGrad(pins, epsilon) - lseMaxGrad(-pins, epsilon));
}

xt::xtensor<float, 1> NetTopologyFixedSize::gradWA(xt::xtensor<float, 1> pl, float epsilon) const {
    auto pins = pinCoords(pl);
    return toCellGrad(waMaxGrad(pins, epsilon) - waMaxGrad(-pins, epsilon));
}

xt::xtensor<float, 1> NetTopologyFixedSize::proximalStep(xt::xtensor<float, 1> pl, float step) const {
    auto pins = pinCoords(pl);
    return toCellGrad(maxProximal(pins, step) - maxProximal(-pins, step));
}

xt::xtensor<float, 1> NetTopologyFixedSize::lseMax(xt::xtensor<float, 2> pins, float epsilon) const {
    xt::xtensor<float, 1> maxVal = xt::amax(pins, {1});
    xt::xtensor<float, 2> normalized = pins - xt::expand_dims(maxVal, 1);
    return maxVal + epsilon * xt::log(xt::sum(xt::exp(normalized / epsilon), 1));
}

xt::xtensor<float, 1> NetTopologyFixedSize::waMax(xt::xtensor<float, 2> pins, float epsilon) const {
    xt::xtensor<float, 1> maxVal = xt::amax(pins, {1});
    // Normalize by the maximum and remove infinities for the product computation later
    xt::xtensor<float, 2> normalized = xt::maximum(pins - xt::expand_dims(maxVal, 1), -1000.0f*epsilon);
    xt::xtensor<float, 2> expValue = xt::exp(normalized / epsilon);
    return maxVal + xt::sum(normalized * expValue, 1) / xt::sum(expValue, 1);
}

xt::xtensor<float, 2> NetTopologyFixedSize::maxGrad(xt::xtensor<float, 2> pins) const {
    xt::xtensor<int, 1> amx = xt::argmax(pins, 1);
    xt::xtensor<float, 2> ret = xt::zeros<float>(pins.shape());
    // It is a pain to do this in xtensor, so let's go for the simple loop
    for (size_t i = 0; i < nbNets_; ++i) {
        ret(i, amx[i]) += 1.0f;
    }
    return ret;
}

xt::xtensor<float, 2> NetTopologyFixedSize::lseMaxGrad(xt::xtensor<float, 2> pins, float epsilon) const {
    xt::xtensor<float, 1> maxVal = xt::amax(pins, {1});
    xt::xtensor<float, 2> normalized = pins - xt::expand_dims(maxVal, 1);
    xt::xtensor<float, 2> expValue = xt::exp(normalized / epsilon);
    xt::xtensor<float, 1> expSum = xt::sum(expValue, 1);
    return expValue / xt::expand_dims(expSum, 1);
}

xt::xtensor<float, 2> NetTopologyFixedSize::waMaxGrad(xt::xtensor<float, 2> pins, float epsilon) const {
    xt::xtensor<float, 1> maxVal = xt::amax(pins, {1});
    // Normalize by the maximum and remove infinities for the product computation later
    xt::xtensor<float, 2> normalized = xt::maximum(pins - xt::expand_dims(maxVal, 1), -1000.0f*epsilon);
    xt::xtensor<float, 2> expValue = xt::exp(normalized / epsilon);
    xt::xtensor<float, 2> num = normalized * expValue;
    xt::xtensor<float, 1> denum = xt::sum(expValue, 1);
    xt::xtensor<float, 2> numDerivative = expValue + num / epsilon;
    xt::xtensor<float, 2> denumDerivative = expValue / epsilon;
    return (numDerivative * xt::expand_dims(denum, 1) - num * denumDerivative) / xt::square(xt::expand_dims(denum, 1));
}

xt::xtensor<float, 2> NetTopologyFixedSize::posProximal(xt::xtensor<float, 2> pins, float step) const {
    float penalty = 0.5 / step;
    xt::xtensor<float, 1> maxVal = xt::amax(pins, {1});
    xt::xtensor<float, 2> normalized = pins - xt::expand_dims(maxVal, 1);
    xt::xtensor<int, 2> sorted = xt::argsort(normalized, 1);
    xt::xtensor<float, 2> dest = xt::zeros<float>(pins.shape());
    // It is a pain to do this kind of stuff in xtensor, so let's go for the simple loop
    // TODO: use a cumsum function or the like
    for (int i = 0; i < nbNets_; ++i) {
        // Current coefficients of the quadratic penalty function a x^2 + b x + c
        // Start with a linear coefficient of 1 for the maximum
        // Then add the displacement penalty
        float a = 0.0;
        float b = 1.0f;
        for (int j = netSize_-1; j >= 0; --j) {
            // Add penalty * (x-pos)^2 to the function penalty
            float pos = normalized(i, sorted(i, j));
            a += penalty;
            b -= 2 * penalty * pos;
            // Compute the minimum point
            dest(i, j) = -b / (2*a);
        }
    }
    xt::xtensor<float, 1> pos = maxVal + xt::amax(dest, 1);
    return xt::minimum(pins, xt::expand_dims(pos, 1));
}

xt::xtensor<float, 2> NetTopologyFixedSize::maxProximal(xt::xtensor<float, 2> pins, float step) const {
    return posProximal(pins, step) - pins;
}

xt::xtensor<float, 2> NetTopologyFixedSize::maxProximal(xt::xtensor<float, 2> pins, xt::xtensor<float, 1> fixedPins, float step) const {
    xt::xtensor<float, 2> pos = posProximal(pins, step);
    // Stop the movement at the fixed pin
    xt::xtensor<float, 2> satPos = xt::maximum(pos, xt::expand_dims(fixedPins, 1));
    // Pins below the fixed pin do not move at all
    xt::xtensor<float, 2> finalPos = xt::minimum(pins, satPos);
    return finalPos - pins;
}

xt::xtensor<float, 1> NetTopologyFixedSize::toCellGrad(xt::xtensor<float, 2> pinGrad) const {
    xt::xtensor<float, 1> flatPinGrad = xt::flatten(pinGrad);
    xt::xtensor<int, 1> flatCells = xt::flatten(pinCells_);
    xt::xtensor<float, 1> ret = xt::zeros<float>({nbCells()});
    for (size_t i = 0; i < flatPinGrad.size(); ++i) {
        ret[flatCells[i]] += flatPinGrad[i];
    }
    return ret;
}

NetTopologyFixedSizeTerminals::NetTopologyFixedSizeTerminals(int nbCells, int netSize, std::vector<int> pinCells, std::vector<float> pinOffsets, std::vector<float> minPos, std::vector<float> maxPos)
: NetTopologyFixedSize(nbCells, netSize, pinCells, pinOffsets) {
    minPins_ = xt::adapt(minPos, {nbNets_});
    maxPins_ = xt::adapt(maxPos, {nbNets_});
    check();
}

NetTopologyFixedSizeTerminals::NetTopologyFixedSizeTerminals(int nbCells, xt::xtensor<int, 2> pinCells, xt::xtensor<float, 2> pinOffsets, xt::xtensor<float, 1> minPos, xt::xtensor<float, 1> maxPos)
: NetTopologyFixedSize(nbCells, pinCells, pinOffsets) {
    minPins_ = minPos;
    maxPins_ = maxPos;
    check();
}

void NetTopologyFixedSizeTerminals::check() const {
    NetTopologyFixedSize::check();
    assert (minPins_.shape()[0] == nbNets_);
    assert (maxPins_.shape()[0] == nbNets_);
}

float NetTopologyFixedSizeTerminals::valueHPWL(xt::xtensor<float, 1> pl) const {
    auto pinsMinMax = pinCoordsMinMax(pl);
    return xt::sum(xt::amax(pinsMinMax.second, {1}) - xt::amin(pinsMinMax.first, {1}))();
}

float NetTopologyFixedSizeTerminals::valueLSE(xt::xtensor<float, 1> pl, float epsilon) const {
    auto pinsMinMax = pinCoordsMinMax(pl);
    return xt::sum(lseMax(pinsMinMax.second, epsilon) + lseMax(-pinsMinMax.first, epsilon))();
}

float NetTopologyFixedSizeTerminals::valueWA(xt::xtensor<float, 1> pl, float epsilon) const {
    auto pinsMinMax = pinCoordsMinMax(pl);
    return xt::sum(waMax(pinsMinMax.second, epsilon) + waMax(-pinsMinMax.first, epsilon))();
}

xt::xtensor<float, 1> NetTopologyFixedSizeTerminals::gradHPWL(xt::xtensor<float, 1> pl) const {
    auto pinsMinMax = pinCoordsMinMax(pl);
    return toCellGradMinMax(maxGrad(pinsMinMax.second) - maxGrad(-pinsMinMax.first));
}

xt::xtensor<float, 1> NetTopologyFixedSizeTerminals::gradLSE(xt::xtensor<float, 1> pl, float epsilon) const {
    auto pinsMinMax = pinCoordsMinMax(pl);
    return toCellGradMinMax(lseMaxGrad(pinsMinMax.second, epsilon) - lseMaxGrad(-pinsMinMax.first, epsilon));
}

xt::xtensor<float, 1> NetTopologyFixedSizeTerminals::gradWA(xt::xtensor<float, 1> pl, float epsilon) const {
    auto pinsMinMax = pinCoordsMinMax(pl);
    return toCellGradMinMax(waMaxGrad(pinsMinMax.second, epsilon) - waMaxGrad(-pinsMinMax.first, epsilon));
}

xt::xtensor<float, 1> NetTopologyFixedSizeTerminals::proximalStep(xt::xtensor<float, 1> pl, float step) const {
    auto pins = pinCoords(pl);
    return toCellGrad(maxProximal(pins, maxPins_, step) - maxProximal(-pins, -minPins_, step));
}

std::pair<xt::xtensor<float, 2>, xt::xtensor<float, 2> > NetTopologyFixedSizeTerminals::pinCoordsMinMax(xt::xtensor<float, 1> pl) const {
    auto pins = pinCoords(pl);
    xt::xtensor<float, 2> minPins = xt::concatenate(xtuple(pins, xt::expand_dims(minPins_, 1)), 1);
    xt::xtensor<float, 2> maxPins = xt::concatenate(xtuple(pins, xt::expand_dims(maxPins_, 1)), 1);
    return std::make_pair(minPins, maxPins);
}

xt::xtensor<float, 1> NetTopologyFixedSizeTerminals::toCellGradMinMax(xt::xtensor<float, 2> pinGrad) const {
    return toCellGrad(xt::view(pinGrad, xt::all(), xt::range(xt::placeholders::_, -1)));
}

void NetTopologyFixedSizeBuilder::push(const std::vector<int> &cells, const std::vector<float> &offsets) {
    assert (cells.size() == netSize_);
    assert (offsets.size() == netSize_);
    pinCells_.insert(pinCells_.end(), cells.begin(), cells.end());
    pinOffsets_.insert(pinOffsets_.end(), offsets.begin(), offsets.end());
}

void NetTopologyFixedSizeTerminalsBuilder::push(const std::vector<int> &cells, const std::vector<float> &offsets, float minPin, float maxPin) {
    assert (cells.size() == netSize_);
    assert (offsets.size() == netSize_);
    pinCells_.insert(pinCells_.end(), cells.begin(), cells.end());
    pinOffsets_.insert(pinOffsets_.end(), offsets.begin(), offsets.end());
    minPins_.push_back(minPin);
    maxPins_.push_back(maxPin);
}
