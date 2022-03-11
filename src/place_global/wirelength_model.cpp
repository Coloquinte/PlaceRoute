
#include "place_global/wirelength_model.hpp"

#include <xtensor/xarray.hpp>
#include <xtensor/xtensor.hpp>
#include <xtensor/xadapt.hpp>
#include <xtensor/xview.hpp>
#include <xtensor/xindex_view.hpp>
#include <xtensor/xstrided_view.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xsort.hpp>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/IterativeLinearSolvers>

#include <algorithm>
#include <cmath>
#include <cassert>
#include <iostream>

NetWirelength NetWirelength::xTopology(const Circuit &circuit) {
    return fromData(circuit.cellWidths, circuit.cellX, circuit.cellFixed, circuit.netLimits, circuit.pinCells, circuit.pinXOffsets);
}

NetWirelength NetWirelength::yTopology(const Circuit &circuit) {
    return fromData(circuit.cellHeights, circuit.cellY, circuit.cellFixed, circuit.netLimits, circuit.pinCells, circuit.pinYOffsets);
}

NetWirelength NetWirelength::fromData(const std::vector<int> &cellSizes, const std::vector<int> & pl, const std::vector<char> &cellFixed, const std::vector<int> &netLimits, const std::vector<int> &pinCells, const std::vector<int> &pinOffsets) {
    NetWirelengthBuilder builder(cellSizes.size());
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
        builder.push(cells, offsets, minPos, maxPos);
    }
    return builder.build();
}

void NetWirelengthBuilder::push(const std::vector<int> &cells, const std::vector<float> &offsets) {
    assert (cells.size() == offsets.size());
    if (cells.size() <= 1) return;
    int sz = cells.size();
    while (netBuilders_.size() < sz + 1) {
        netBuilders_.emplace_back(nbCells_, netBuilders_.size());
    }
    netBuilders_[sz].push(cells, offsets);
}

void NetWirelengthBuilder::push(const std::vector<int> &cells, const std::vector<float> &offsets, float minPos, float maxPos) {
    assert (cells.size() == offsets.size());
    if (std::isfinite(minPos)) {
        assert (maxPos >= minPos);
        if (cells.empty()) return;
        int sz = cells.size();
        while (terminalNetBuilders_.size() < sz + 1) {
            terminalNetBuilders_.emplace_back(nbCells_, terminalNetBuilders_.size());
        }
        terminalNetBuilders_[sz].push(cells, offsets, minPos, maxPos);
    }
    else {
        push(cells, offsets);
    }
    ++nbNets_;
}

NetWirelength NetWirelengthBuilder::build() const {
    NetWirelength ret;
    ret.nbCells_ = nbCells_;
    ret.nbNets_ = nbNets_;
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

float NetWirelength::valueHPWL(const xt::xtensor<float, 1> &pl) const {
    float val = 0.0;
    for (const auto &bd : nets_) {
        val += bd.valueHPWL(pl);
    }
    for (const auto &bd : terminalNets_) {
        val += bd.valueHPWL(pl);
    }
    return val;
}

float NetWirelength::valueLSE(const xt::xtensor<float, 1> &pl, float epsilon) const {
    float val = 0.0;
    for (const auto &bd : nets_) {
        val += bd.valueLSE(pl, epsilon);
    }
    for (const auto &bd : terminalNets_) {
        val += bd.valueLSE(pl, epsilon);
    }
    return val;
}

float NetWirelength::valueWA(const xt::xtensor<float, 1> &pl, float epsilon) const {
    float val = 0.0;
    for (const auto &bd : nets_) {
        val += bd.valueWA(pl, epsilon);
    }
    for (const auto &bd : terminalNets_) {
        val += bd.valueWA(pl, epsilon);
    }
    return val;
}

xt::xtensor<float, 1> NetWirelength::gradHPWL(const xt::xtensor<float, 1> &pl) const {
    xt::xtensor<float, 1> val = xt::zeros<float>({nbCells()});
    for (const auto &bd : nets_) {
        val += bd.gradHPWL(pl);
    }
    for (const auto &bd : terminalNets_) {
        val += bd.gradHPWL(pl);
    }
    return val;
}

xt::xtensor<float, 1> NetWirelength::gradLSE(const xt::xtensor<float, 1> &pl, float epsilon) const {
    xt::xtensor<float, 1> val = xt::zeros<float>({nbCells()});
    for (const auto &bd : nets_) {
        val += bd.gradLSE(pl, epsilon);
    }
    for (const auto &bd : terminalNets_) {
        val += bd.gradLSE(pl, epsilon);
    }
    return val;
}

xt::xtensor<float, 1> NetWirelength::gradWA(const xt::xtensor<float, 1> &pl, float epsilon) const {
    xt::xtensor<float, 1> val = xt::zeros<float>({nbCells()});
    for (const auto &bd : nets_) {
        val += bd.gradWA(pl, epsilon);
    }
    for (const auto &bd : terminalNets_) {
        val += bd.gradWA(pl, epsilon);
    }
    return val;
}

xt::xtensor<float, 1> NetWirelength::proximalStep(const xt::xtensor<float, 1> &pl, float step) const {
    xt::xtensor<float, 1> val = xt::zeros<float>({nbCells()});
    for (const auto &bd : nets_) {
        val += bd.proximalStep(pl, step);
    }
    for (const auto &bd : terminalNets_) {
        val += bd.proximalStep(pl, step);
    }
    return val;
}

NetWirelengthFixedSize::NetWirelengthFixedSize(int nbCells, int netSize, std::vector<int> pinCells, std::vector<float> pinOffsets) {
    assert (pinCells.size() == pinOffsets.size());
    assert (pinCells.size() % netSize == 0);
    nbNets_ = pinCells.size() / netSize;
    nbCells_ = nbCells;
    netSize_ = netSize;
    pinCells_ = xt::adapt(pinCells, {nbNets_, netSize_});
    offsets_ = xt::adapt(pinOffsets, {nbNets_, netSize_});
    check();
}

NetWirelengthFixedSize::NetWirelengthFixedSize(int nbCells, xt::xtensor<int, 2> pinCells, xt::xtensor<float, 2> pinOffsets) {
    assert (pinCells.shape()[0] == pinOffsets.shape()[0]);
    assert (pinCells.shape()[1] == pinOffsets.shape()[1]);
    nbCells_ = nbCells;
    nbNets_ = pinCells.shape()[0];
    netSize_ = pinCells.shape()[1];
    pinCells_ = pinCells;
    offsets_ = pinOffsets;
    check();
}

void NetWirelengthFixedSize::check() const {
    assert (pinCells_.shape()[0] == nbNets_);
    assert (offsets_.shape()[0] == nbNets_);
    assert (pinCells_.shape()[1] == netSize_);
    assert (offsets_.shape()[1] == netSize_);
    assert (xt::all(pinCells_ < nbCells()));
    assert (xt::all(pinCells_ >= 0));
}

xt::xtensor<float, 2> NetWirelengthFixedSize::pinCoords(xt::xtensor<float, 1> pl) const {
    assert (pl.size() == nbCells());
    xt::xtensor<float, 2> cellCoords = xt::reshape_view(
        xt::index_view(pl, pinCells_),
        {nbNets(), netSize()}
    );
    return cellCoords + offsets_;
}

float NetWirelengthFixedSize::valueHPWL(xt::xtensor<float, 1> pl) const {
    auto pins = pinCoords(pl);
    return xt::sum(xt::amax(pins, {1}) - xt::amin(pins, {1}))();
}

float NetWirelengthFixedSize::valueLSE(xt::xtensor<float, 1> pl, float epsilon) const {
    auto pins = pinCoords(pl);
    return xt::sum(lseMax(pins, epsilon) + lseMax(-pins, epsilon))();
}

float NetWirelengthFixedSize::valueWA(xt::xtensor<float, 1> pl, float epsilon) const {
    auto pins = pinCoords(pl);
    return xt::sum(waMax(pins, epsilon) + waMax(-pins, epsilon))();
}

xt::xtensor<float, 1> NetWirelengthFixedSize::gradHPWL(xt::xtensor<float, 1> pl) const {
    auto pins = pinCoords(pl);
    return toCellGrad(maxGrad(pins) - maxGrad(-pins));
}

xt::xtensor<float, 1> NetWirelengthFixedSize::gradLSE(xt::xtensor<float, 1> pl, float epsilon) const {
    auto pins = pinCoords(pl);
    return toCellGrad(lseMaxGrad(pins, epsilon) - lseMaxGrad(-pins, epsilon));
}

xt::xtensor<float, 1> NetWirelengthFixedSize::gradWA(xt::xtensor<float, 1> pl, float epsilon) const {
    auto pins = pinCoords(pl);
    return toCellGrad(waMaxGrad(pins, epsilon) - waMaxGrad(-pins, epsilon));
}

xt::xtensor<float, 1> NetWirelengthFixedSize::proximalStep(xt::xtensor<float, 1> pl, float step) const {
    auto pins = pinCoords(pl);
    return toCellGrad(maxProximal(pins, step) - maxProximal(-pins, step));
}

xt::xtensor<float, 1> NetWirelengthFixedSize::lseMax(xt::xtensor<float, 2> pins, float epsilon) const {
    xt::xtensor<float, 1> maxVal = xt::amax(pins, {1});
    xt::xtensor<float, 2> normalized = pins - xt::expand_dims(maxVal, 1);
    return maxVal + epsilon * xt::log(xt::sum(xt::exp(normalized / epsilon), 1));
}

xt::xtensor<float, 1> NetWirelengthFixedSize::waMax(xt::xtensor<float, 2> pins, float epsilon) const {
    xt::xtensor<float, 1> maxVal = xt::amax(pins, {1});
    // Normalize by the maximum and remove infinities for the product computation later
    xt::xtensor<float, 2> normalized = xt::maximum(pins - xt::expand_dims(maxVal, 1), -1000.0f*epsilon);
    xt::xtensor<float, 2> expValue = xt::exp(normalized / epsilon);
    return maxVal + xt::sum(normalized * expValue, 1) / xt::sum(expValue, 1);
}

xt::xtensor<float, 2> NetWirelengthFixedSize::maxGrad(xt::xtensor<float, 2> pins) const {
    xt::xtensor<int, 1> amx = xt::argmax(pins, 1);
    xt::xtensor<float, 2> ret = xt::zeros<float>(pins.shape());
    // It is a pain to do this in xtensor, so let's go for the simple loop
    for (size_t i = 0; i < nbNets_; ++i) {
        ret(i, amx[i]) += 1.0f;
    }
    return ret;
}

xt::xtensor<float, 2> NetWirelengthFixedSize::lseMaxGrad(xt::xtensor<float, 2> pins, float epsilon) const {
    xt::xtensor<float, 1> maxVal = xt::amax(pins, {1});
    xt::xtensor<float, 2> normalized = pins - xt::expand_dims(maxVal, 1);
    xt::xtensor<float, 2> expValue = xt::exp(normalized / epsilon);
    xt::xtensor<float, 1> expSum = xt::sum(expValue, 1);
    return expValue / xt::expand_dims(expSum, 1);
}

xt::xtensor<float, 2> NetWirelengthFixedSize::waMaxGrad(xt::xtensor<float, 2> pins, float epsilon) const {
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

xt::xtensor<float, 2> NetWirelengthFixedSize::posProximal(xt::xtensor<float, 2> pins, float step) const {
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

xt::xtensor<float, 2> NetWirelengthFixedSize::maxProximal(xt::xtensor<float, 2> pins, float step) const {
    return posProximal(pins, step) - pins;
}

xt::xtensor<float, 2> NetWirelengthFixedSize::maxProximal(xt::xtensor<float, 2> pins, xt::xtensor<float, 1> fixedPins, float step) const {
    xt::xtensor<float, 2> pos = posProximal(pins, step);
    // Stop the movement at the fixed pin
    xt::xtensor<float, 2> satPos = xt::maximum(pos, xt::expand_dims(fixedPins, 1));
    // Pins below the fixed pin do not move at all
    xt::xtensor<float, 2> finalPos = xt::minimum(pins, satPos);
    return finalPos - pins;
}

xt::xtensor<float, 1> NetWirelengthFixedSize::toCellGrad(xt::xtensor<float, 2> pinGrad) const {
    xt::xtensor<float, 1> flatPinGrad = xt::flatten(pinGrad);
    xt::xtensor<int, 1> flatCells = xt::flatten(pinCells_);
    xt::xtensor<float, 1> ret = xt::zeros<float>({nbCells()});
    for (size_t i = 0; i < flatPinGrad.size(); ++i) {
        ret[flatCells[i]] += flatPinGrad[i];
    }
    return ret;
}

NetWirelengthFixedSizeTerminals::NetWirelengthFixedSizeTerminals(int nbCells, int netSize, std::vector<int> pinCells, std::vector<float> pinOffsets, std::vector<float> minPos, std::vector<float> maxPos)
: NetWirelengthFixedSize(nbCells, netSize, pinCells, pinOffsets) {
    minPins_ = xt::adapt(minPos, {nbNets_});
    maxPins_ = xt::adapt(maxPos, {nbNets_});
    check();
}

NetWirelengthFixedSizeTerminals::NetWirelengthFixedSizeTerminals(int nbCells, xt::xtensor<int, 2> pinCells, xt::xtensor<float, 2> pinOffsets, xt::xtensor<float, 1> minPos, xt::xtensor<float, 1> maxPos)
: NetWirelengthFixedSize(nbCells, pinCells, pinOffsets) {
    minPins_ = minPos;
    maxPins_ = maxPos;
    check();
}

void NetWirelengthFixedSizeTerminals::check() const {
    NetWirelengthFixedSize::check();
    assert (minPins_.shape()[0] == nbNets_);
    assert (maxPins_.shape()[0] == nbNets_);
}

float NetWirelengthFixedSizeTerminals::valueHPWL(xt::xtensor<float, 1> pl) const {
    auto pinsMinMax = pinCoordsMinMax(pl);
    return xt::sum(xt::amax(pinsMinMax.second, {1}) - xt::amin(pinsMinMax.first, {1}))();
}

float NetWirelengthFixedSizeTerminals::valueLSE(xt::xtensor<float, 1> pl, float epsilon) const {
    auto pinsMinMax = pinCoordsMinMax(pl);
    return xt::sum(lseMax(pinsMinMax.second, epsilon) + lseMax(-pinsMinMax.first, epsilon))();
}

float NetWirelengthFixedSizeTerminals::valueWA(xt::xtensor<float, 1> pl, float epsilon) const {
    auto pinsMinMax = pinCoordsMinMax(pl);
    return xt::sum(waMax(pinsMinMax.second, epsilon) + waMax(-pinsMinMax.first, epsilon))();
}

xt::xtensor<float, 1> NetWirelengthFixedSizeTerminals::gradHPWL(xt::xtensor<float, 1> pl) const {
    auto pinsMinMax = pinCoordsMinMax(pl);
    return toCellGradMinMax(maxGrad(pinsMinMax.second) - maxGrad(-pinsMinMax.first));
}

xt::xtensor<float, 1> NetWirelengthFixedSizeTerminals::gradLSE(xt::xtensor<float, 1> pl, float epsilon) const {
    auto pinsMinMax = pinCoordsMinMax(pl);
    return toCellGradMinMax(lseMaxGrad(pinsMinMax.second, epsilon) - lseMaxGrad(-pinsMinMax.first, epsilon));
}

xt::xtensor<float, 1> NetWirelengthFixedSizeTerminals::gradWA(xt::xtensor<float, 1> pl, float epsilon) const {
    auto pinsMinMax = pinCoordsMinMax(pl);
    return toCellGradMinMax(waMaxGrad(pinsMinMax.second, epsilon) - waMaxGrad(-pinsMinMax.first, epsilon));
}

xt::xtensor<float, 1> NetWirelengthFixedSizeTerminals::proximalStep(xt::xtensor<float, 1> pl, float step) const {
    auto pins = pinCoords(pl);
    return toCellGrad(maxProximal(pins, maxPins_, step) - maxProximal(-pins, -minPins_, step));
}

std::pair<xt::xtensor<float, 2>, xt::xtensor<float, 2> > NetWirelengthFixedSizeTerminals::pinCoordsMinMax(xt::xtensor<float, 1> pl) const {
    auto pins = pinCoords(pl);
    xt::xtensor<float, 2> minPins = xt::concatenate(xtuple(pins, xt::expand_dims(minPins_, 1)), 1);
    xt::xtensor<float, 2> maxPins = xt::concatenate(xtuple(pins, xt::expand_dims(maxPins_, 1)), 1);
    return std::make_pair(minPins, maxPins);
}

xt::xtensor<float, 2> NetWirelengthFixedSizeTerminals::pinCoordsAll(xt::xtensor<float, 1> pl) const {
    auto pins = pinCoords(pl);
    return xt::concatenate(xtuple(pins, xt::expand_dims(minPins_, 1), xt::expand_dims(maxPins_, 1)), 1);
}

xt::xtensor<int, 2> NetWirelengthFixedSizeTerminals::pinCellsAll() const {
    return xt::concatenate(xtuple(pinCells(), -xt::ones<int>({nbNets(), 2})), 1);
}

xt::xtensor<float, 2> NetWirelengthFixedSizeTerminals::pinOffsetsAll() const {
    return xt::concatenate(xtuple(pinOffsets(), xt::expand_dims(minPins_, 1), xt::expand_dims(maxPins_, 1)), 1);
}

xt::xtensor<float, 1> NetWirelengthFixedSizeTerminals::toCellGradMinMax(xt::xtensor<float, 2> pinGrad) const {
    return toCellGrad(xt::view(pinGrad, xt::all(), xt::range(xt::placeholders::_, -1)));
}

void NetWirelengthFixedSizeBuilder::push(const std::vector<int> &cells, const std::vector<float> &offsets) {
    assert (cells.size() == netSize_);
    assert (offsets.size() == netSize_);
    pinCells_.insert(pinCells_.end(), cells.begin(), cells.end());
    pinOffsets_.insert(pinOffsets_.end(), offsets.begin(), offsets.end());
}

void NetWirelengthFixedSizeTerminalsBuilder::push(const std::vector<int> &cells, const std::vector<float> &offsets, float minPin, float maxPin) {
    assert (cells.size() == netSize_);
    assert (offsets.size() == netSize_);
    pinCells_.insert(pinCells_.end(), cells.begin(), cells.end());
    pinOffsets_.insert(pinOffsets_.end(), offsets.begin(), offsets.end());
    minPins_.push_back(minPin);
    maxPins_.push_back(maxPin);
}

MatrixBuilder MatrixBuilder::createStar(const NetWirelength &topo) {
    MatrixBuilder bd(topo.nbCells());
    for (const auto &nt : topo.nets()) {
        bd.extendStar(nt);
    }
    for (const auto &nt : topo.terminalNets()) {
        bd.extendStar(nt);
    }
    return bd;
}

void MatrixBuilder::extendStar(const NetWirelengthFixedSize &topo) {
    const auto &cells = topo.pinCells();
    const auto &offsets = topo.pinOffsets();
    for (int i = 0; i < topo.nbNets(); ++i) {
        if (topo.netSize() == 2) {
            addPin(cells(i, 0), cells(i, 1), offsets(i, 0), offsets(i, 1), 0.5);
        }
        else if (topo.netSize() == 3) {
            addPin(cells(i, 0), cells(i, 1), offsets(i, 0), offsets(i, 1), 0.5);
            addPin(cells(i, 0), cells(i, 2), offsets(i, 0), offsets(i, 2), 0.5);
            addPin(cells(i, 2), cells(i, 1), offsets(i, 2), offsets(i, 1), 0.5);
        }
        else {
            float weight = 1.0 / topo.netSize();
            int newCell = nbCells_ + nbSupps_;
            nbSupps_++;
            rhs_.push_back(0.0);
            initial_.push_back(0.0);
            for (int j = 0; j < topo.netSize(); ++j) {
                addPin(cells(i, j), newCell, offsets(i, j), 0.0f, weight);
            }
        }
    }
}

void MatrixBuilder::extendStar(const NetWirelengthFixedSizeTerminals &topo) {
    const auto &cells = topo.pinCells();
    const auto &offsets = topo.pinOffsets();
    for (int i = 0; i < topo.nbNets(); ++i) {
        float pos = 0.5 * (topo.minPins()(i) + topo.maxPins()(i));
        if (topo.netSize() == 1) {
            addFixedPin(cells(i, 0), offsets(i, 0), pos, 0.5);
        }
        else if (topo.netSize() == 2) {
            addFixedPin(cells(i, 1), offsets(i, 1), pos, 0.5);
            addFixedPin(cells(i, 1), offsets(i, 1), pos, 0.5);
            addPin(cells(i, 0), cells(i, 1), offsets(i, 0), offsets(i, 1), 0.5);
        }
        else {
            float weight = 1.0 / (topo.netSize() + 1);
            int newCell = nbCells_ + nbSupps_;
            nbSupps_++;
            rhs_.push_back(0.0);
            initial_.push_back(0.0);
            addFixedPin(newCell, 0.0, pos, weight);
            for (int j = 0; j < topo.netSize(); ++j) {
                addPin(cells(i, j), newCell, offsets(i, j), 0.0f, weight);
            }
        }
    }
}

MatrixBuilder MatrixBuilder::createStar(const NetWirelength &topo, xt::xtensor<float, 1> pl, float epsilon, float relaxation, bool b2bScale) {
    MatrixBuilder bd(topo.nbCells());
    bd.initial_ = std::vector<float>(pl.begin(), pl.end());
    for (const auto &nt : topo.nets()) {
        bd.extendStar(nt, pl, epsilon, relaxation, b2bScale);
    }
    for (const auto &nt : topo.terminalNets()) {
        bd.extendStar(nt, pl, epsilon, relaxation, b2bScale);
    }
    return bd;
}

void MatrixBuilder::extendStar(const NetWirelengthFixedSize &topo, xt::xtensor<float, 1> pl, float epsilon, float relaxation, bool b2bScale) {
    const auto coords = topo.pinCoords(pl);
    const auto &cells = topo.pinCells();
    const auto &offsets = topo.pinOffsets();
    xt::xtensor<float, 1> pinMin = xt::amin(coords, {1});
    xt::xtensor<float, 1> pinMax = xt::amax(coords, {1});
    xt::xtensor<float, 1> pinMed = 0.5 * (pinMin + pinMax);
    xt::xtensor<int, 1> amnt = xt::argmin(coords, 1);
    xt::xtensor<int, 1> amxt = xt::argmax(coords, 1);
    for (int i = 0; i < topo.nbNets(); ++i) {
        if (topo.netSize() == 2) {
            float weight = 1.0 / std::max(pinMax(i) - pinMin(i), epsilon);
            addPin(cells(i, 0), cells(i, 1), offsets(i, 0), offsets(i, 1), weight);
        }
        else {
            int newCell = nbCells_ + nbSupps_;
            nbSupps_++;
            rhs_.push_back(0.0);
            initial_.push_back(pinMed(i));
            float halfDiameter = pinMax(i) - pinMed(i);
            for (int j = 0; j < topo.netSize(); ++j) {
                if (j == amnt(i) || j == amxt(i)) {
                    // Extreme pins are connected directly (gradient of one)
                    float weight = 1.0f / std::max(halfDiameter, epsilon);
                    addPin(cells(i, j), newCell, offsets(i, j), 0.0f, weight);
                }
                else {
                    // Non-extreme pins have a zero force at their current position
                    float pos = coords(i, j);
                    float weight;
                    if (b2bScale) {
                        float w1 = 1.0f / std::max(pinMax(i) - pos, epsilon);
                        float w2 = 1.0f / std::max(pos - pinMin(i), epsilon);
                        // Other possibility, closer to actual B2B but not as good
                        // weight = (w1 + w2) / (topo.netSize() - 1);
                        weight = 2.0 * std::max(w1, w2) / (topo.netSize() - 1);
                    }
                    else {
                        float dist = std::min(pinMax(i) - pos, pos - pinMin(i));
                        weight = 1.0f / std::max(dist, epsilon);
                    }
                    if (relaxation > 0.0f) {
                        float rdist = relaxation * halfDiameter;
                        float cutoff = 1.0f / std::max(rdist, epsilon);
                        weight = std::min(weight, cutoff);
                    }
                    addPin(cells(i, j), newCell, offsets(i, j), pos - pinMed(i), weight);
                }
            }
        }
    }
}

void MatrixBuilder::extendStar(const NetWirelengthFixedSizeTerminals &topo, xt::xtensor<float, 1> pl, float epsilon, float relaxation, bool b2bScale) {
    const auto coords = topo.pinCoordsAll(pl);
    const auto cells = topo.pinCellsAll();
    const auto offsets = topo.pinOffsetsAll();
    xt::xtensor<float, 1> pinMin = xt::amin(coords, {1});
    xt::xtensor<float, 1> pinMax = xt::amax(coords, {1});
    xt::xtensor<float, 1> pinMed = 0.5 * (pinMin + pinMax);
    xt::xtensor<int, 1> amnt = xt::argmin(coords, 1);
    xt::xtensor<int, 1> amxt = xt::argmax(coords, 1);
    for (int i = 0; i < topo.nbNets(); ++i) {
        float minFixed = topo.minPins()(i);
        float maxFixed = topo.maxPins()(i);

        if (topo.netSize() == 1 && minFixed == maxFixed) {
            float weight = 1.0 / std::max(pinMax(i) - pinMin(i), epsilon);
            addFixedPin(cells(i, 0), offsets(i, 0), minFixed, weight);
        }
        else {
            int amn = amnt(i);
            int amx = amxt(i);
            int nbFixed;
            if (topo.minPins()(i) == topo.maxPins()(i)) {
                // Only consider one fixed pin
                nbFixed = 1;
                amn = std::min(amn, topo.netSize());
                amx = std::min(amx, topo.netSize());
            }
            else {
                nbFixed = 2;
            }
            int newCell = nbCells_ + nbSupps_;
            nbSupps_++;
            rhs_.push_back(0.0);
            initial_.push_back(pinMed(i));
            float halfDiameter = pinMax(i) - pinMed(i);
            for (int j = 0; j < topo.netSize() + nbFixed; ++j) {
                if (j == amn || j == amx) {
                    // Extreme pins are connected directly (gradient of one)
                    float weight = 1.0f / std::max(halfDiameter, epsilon);
                    addPinOrFixed(cells(i, j), newCell, offsets(i, j), 0.0f, weight);
                }
                else {
                    // Non-extreme pins have a zero force at their current position
                    float pos = coords(i, j);
                    float weight;
                    if (b2bScale) {
                        float w1 = 1.0f / std::max(pinMax(i) - pos, epsilon);
                        float w2 = 1.0f / std::max(pos - pinMin(i), epsilon);
                        // Other possibility, closer to actual B2B but not as good
                        // weight = (w1 + w2) / (topo.netSize() + nbFixed - 1);
                        weight = 2.0 * std::max(w1, w2) / (topo.netSize() + nbFixed - 1);
                    }
                    else {
                        float dist = std::min(pinMax(i) - pos, pos - pinMin(i));
                        weight = 1.0f / std::max(dist, epsilon);
                    }
                    if (relaxation > 0.0f) {
                        float rdist = relaxation * halfDiameter;
                        float cutoff = 1.0f / std::max(rdist, epsilon);
                        weight = std::min(weight, cutoff);
                    }
                    addPinOrFixed(cells(i, j), newCell, offsets(i, j), pos - pinMed(i), weight);
                }
            }
        }
    }
}

MatrixBuilder MatrixBuilder::createB2B(const NetWirelength &topo, xt::xtensor<float, 1> pl, float epsilon) {
    MatrixBuilder bd(topo.nbCells());
    bd.initial_ = std::vector<float>(pl.begin(), pl.end());
    for (const auto &nt : topo.nets()) {
        bd.extendB2B(nt, pl, epsilon);
    }
    for (const auto &nt : topo.terminalNets()) {
        bd.extendB2B(nt, pl, epsilon);
    }
    return bd;
}

void MatrixBuilder::extendB2B(const NetWirelengthFixedSize &topo, xt::xtensor<float, 1> pl, float epsilon) {
    const auto coords = topo.pinCoords(pl);
    const auto &cells = topo.pinCells();
    const auto &offsets = topo.pinOffsets();
    xt::xtensor<int, 1> amnt = xt::argmin(coords, 1);
    xt::xtensor<int, 1> amxt = xt::argmax(coords, 1);
    for (int i = 0; i < topo.nbNets(); ++i) {
        int amn = amnt(i);
        int amx = amxt(i);
        float factor = 1.0f / (topo.netSize() - 1);
        for (int j = 0; j < topo.netSize(); ++j) {
            float weight;
            if (j == amn) continue;
            weight = factor / std::max(std::abs(coords(i, j) - coords(i, amn)), epsilon);
            addPin(cells(i, j), cells(i, amn), offsets(i, j), offsets(i, amn), weight);
            if (j == amx) continue;
            weight = factor / std::max(std::abs(coords(i, j) - coords(i, amx)), epsilon);
            addPin(cells(i, j), cells(i, amx), offsets(i, j), offsets(i, amx), weight);
        }
    }
}

void MatrixBuilder::extendB2B(const NetWirelengthFixedSizeTerminals &topo, xt::xtensor<float, 1> pl, float epsilon) {
    const auto coords = topo.pinCoordsAll(pl);
    const auto cells = topo.pinCellsAll();
    const auto offsets = topo.pinOffsetsAll();
    xt::xtensor<int, 1> amnt = xt::argmin(coords, 1);
    xt::xtensor<int, 1> amxt = xt::argmax(coords, 1);
    for (int i = 0; i < topo.nbNets(); ++i) {
        int amn = amnt(i);
        int amx = amxt(i);
        int nbFixed;
        if (topo.minPins()(i) == topo.maxPins()(i)) {
            // Only consider one fixed pin
            nbFixed = 1;
            amn = std::min(amn, topo.netSize());
            amx = std::min(amx, topo.netSize());
        }
        else {
            nbFixed = 2;
        }
        float factor = 1.0f / (topo.netSize() + nbFixed - 1);
        for (int j = 0; j < topo.netSize() + nbFixed; ++j) {
            float weight;
            if (j == amn) continue;
            weight = factor / std::max(std::abs(coords(i, j) - coords(i, amn)), epsilon);
            addPinOrFixed(cells(i, j), cells(i, amn), offsets(i, j), offsets(i, amn), weight);
            if (j == amx) continue;
            weight = factor / std::max(std::abs(coords(i, j) - coords(i, amx)), epsilon);
            addPinOrFixed(cells(i, j), cells(i, amx), offsets(i, j), offsets(i, amx), weight);
        }
    }
}

void MatrixBuilder::addPin(int c1, int c2, float offs1, float offs2, float weight) {
    if (c1 == c2) return;
    assert (c1 >= 0);
    assert (c2 >= 0);
    mat_.emplace_back(c1, c2, -weight);
    mat_.emplace_back(c2, c1, -weight);
    mat_.emplace_back(c1, c1, weight);
    mat_.emplace_back(c2, c2, weight);
    rhs_[c1] += weight * (offs2 - offs1);
    rhs_[c2] += weight * (offs1 - offs2);
}

void MatrixBuilder::addFixedPin(int c1, float offs1, float pos, float weight) {
    assert (c1 >= 0);
    mat_.emplace_back(c1, c1, weight);
    rhs_[c1] += weight * (pos - offs1);
}

void MatrixBuilder::addPinOrFixed(int c1, int c2, float offs1, float offs2, float weight) {
    if (c1 == c2) return;
    if (c1 == -1) {
        addFixedPin(c2, offs2, offs1, weight);
        return;
    }
    if (c2 == -1) {
        addFixedPin(c1, offs1, offs2, weight);
        return;
    }
    addPin(c1, c2, offs1, offs2, weight);
}

xt::xtensor<float, 1> MatrixBuilder::solve() {
    Eigen::SparseMatrix<float> mat(matSize(), matSize());
    mat.setFromTriplets(mat_.begin(), mat_.end());
    Eigen::Map<Eigen::Matrix<float, -1, 1> > rhs(rhs_.data(), rhs_.size());
    Eigen::Map<Eigen::Matrix<float, -1, 1> > initial(initial_.data(), initial_.size());
    Eigen::ConjugateGradient<Eigen::SparseMatrix<float> > solver;
    solver.compute(mat);
    solver.setTolerance(1.0e-6f);
    solver.setMaxIterations(1000);
    Eigen::Matrix<float, -1, 1> res = solver.solveWithGuess(rhs, initial);
    // Copy to a std::vector and remove the fake cells
    std::vector<float> ret;
    ret.resize(matSize());
    Eigen::Matrix<float, -1, 1>::Map(ret.data(), ret.size()) = res;
    ret.resize(nbCells_);
    return xt::adapt(ret, {nbCells_});
}

xt::xtensor<float, 1> NetWirelength::starSolve() const {
    MatrixBuilder bd = MatrixBuilder::createStar(*this);
    return bd.solve();
}

xt::xtensor<float, 1> NetWirelength::starSolve(const xt::xtensor<float, 1> &pl, float epsilon, float relaxation, bool b2bScale) const {
    MatrixBuilder bd = MatrixBuilder::createStar(*this, pl, epsilon, relaxation, b2bScale);
    return bd.solve();
}

xt::xtensor<float, 1> NetWirelength::b2bSolve(const xt::xtensor<float, 1> &pl, float epsilon) const {
    MatrixBuilder bd = MatrixBuilder::createB2B(*this, pl, epsilon);
    return bd.solve();
}
