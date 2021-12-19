
#include "place_global/descent.hpp"

#include <iostream>
#include <iomanip>

namespace {
void reportHeader() {
    std::cout << "Model\tSmoothing";
    std::cout << "\tStepSize\tMomentum";
    std::cout << "\tStep#";
    std::cout << "\tHPWL\tValue";
    std::cout << std::endl;
}

void reportFooter() {
    std::cout << std::endl;
}

void report(const NetTopology &topo, xt::xtensor<float, 1> pl, DescentModel model, int stepInd, float stepSize, float momentum, float smoothing) {
    float smoothedValue;
    float hpwlValue = topo.valueHPWL(pl);

    std::cout << std::fixed << std::setprecision(4);

    if (model == DescentModel::HPWL) {
        std::cout << "HPWL\t-";
        smoothedValue = hpwlValue;
    }
    else if (model == DescentModel::LSE) {
        smoothedValue = topo.valueLSE(pl, smoothing);
        std::cout << "LSE\t" << smoothing;
    }
    else if (model == DescentModel::WA) {
        smoothedValue = topo.valueWA(pl, smoothing);
        std::cout << "WA\t" << smoothing;
    }
    if (model == DescentModel::Proximal) {
        std::cout << "Prox\t-";
        smoothedValue = hpwlValue;
    }
    std::cout << "\t" << stepSize;

    if (momentum >= 0.0) {
        std::cout << "\t" << momentum;
    }
    else {
        std::cout << "\tNAG";
    }
    std::cout << "\t" << stepInd;

    std::cout << std::scientific << std::setprecision(4);
    std::cout << "\t" << hpwlValue << "\t" << smoothedValue;
    std::cout << std::endl;
}

xt::xtensor<float, 1> descentStep(const NetTopology &topo, xt::xtensor<float, 1> current, DescentModel model, float stepSize, float smoothing) {
    if (model == DescentModel::HPWL) {
        return -stepSize * topo.gradHPWL(current);
    }
    else if (model == DescentModel::LSE) {
        return -stepSize * topo.gradLSE(current, smoothing);
    }
    else if (model == DescentModel::WA) {
        return -stepSize * topo.gradWA(current, smoothing);
    }
    else { // if (model == DescentModel::Proximal) {
        return topo.proximalStep(current, stepSize);
    }
}
}

xt::xtensor<float, 1> gradientDescent(const NetTopology &topo, xt::xtensor<float, 1> initial, DescentModel model, int nbSteps, float stepSize, float momentum, float smoothing) {
    assert (stepSize >= 0.0);
    assert (momentum >= 0.0 && momentum < 1.0);
    assert (initial.size() == topo.nbCells());
    reportHeader();
    xt::xtensor<float, 1> x = initial;
    xt::xtensor<float, 1> y = initial;
    for (int i = 0; i < nbSteps; ++i) {
        report(topo, x, model, i, stepSize, momentum, smoothing);
        xt::xtensor<float, 1> z = y + descentStep(topo, x, model, stepSize, smoothing);
        y = z + momentum * (z - x);
        x = z;
    }
    report(topo, x, model, nbSteps, stepSize, momentum, smoothing);
    reportFooter();
    return x;
}
