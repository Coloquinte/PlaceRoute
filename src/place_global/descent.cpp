
#include "place_global/descent.hpp"

#include <iostream>
#include <iomanip>

namespace {
void lineSearch(xt::xtensor<float, 1> start, xt::xtensor<float, 1> direction, float &step, float factor, const std::function<float(xt::xtensor<float, 1>)> &eval) {
    /**
     * Simple line search where we look for the minimum of the function; no checking of Wolfe/Armijo conditions
     */
    float value = eval(start + step * direction);
    bool largerStep = false;
    while (true) {
        float nextStep = step / factor;
        float nextValue = eval(start + nextStep * direction);
        if (nextValue < value) {
            step = nextStep;
            value = nextValue;
            largerStep = true;
        }
        else {
            break;
        }
    }
    if (largerStep) return;
    while (true) {
        float nextStep = step * factor;
        float nextValue = eval(start + nextStep * direction);
        if (nextValue < value) {
            step = nextStep;
            value = nextValue;
        }
        else {
            return;
        }
    }
}

float computeValue(const NetTopology &topo, xt::xtensor<float, 1> current, DescentModel model, float smoothing) {
    if (model == DescentModel::HPWL) {
        return topo.valueHPWL(current);
    }
    else if (model == DescentModel::LSE) {
        return topo.valueLSE(current, smoothing);
    }
    else if (model == DescentModel::WA) {
        return topo.valueWA(current, smoothing);
    }
    else { // if (model == DescentModel::Proximal) {
        return topo.valueHPWL(current);
    }
}

xt::xtensor<float, 1> descentDirection(const NetTopology &topo, xt::xtensor<float, 1> current, DescentModel model, float smoothing) {
    if (model == DescentModel::HPWL) {
        return -topo.gradHPWL(current);
    }
    else if (model == DescentModel::LSE) {
        return -topo.gradLSE(current, smoothing);
    }
    else if (model == DescentModel::WA) {
        return -topo.gradWA(current, smoothing);
    }
    else { // if (model == DescentModel::Proximal) {
        return topo.proximalStep(current, smoothing);
    }
}

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

void report(const NetTopology &topo, xt::xtensor<float, 1> pl, DescentModel model, int stepInd, float stepSize, float momentum, float smoothing, bool lineSearch) {
    float hpwlValue = topo.valueHPWL(pl);
    float smoothedValue = computeValue(topo, pl, model, smoothing);

    std::cout << std::fixed << std::setprecision(4);

    if (model == DescentModel::HPWL) {
        std::cout << "HPWL\t";
    }
    else if (model == DescentModel::LSE) {
        std::cout << "LSE\t";
    }
    else if (model == DescentModel::WA) {
        std::cout << "WA\t";
    }
    else if (model == DescentModel::Proximal) {
        std::cout << "Prox\t";
    }
    std::cout << "\t" << smoothing;

    if (lineSearch) {
        std::cout << "\tLS";
    }
    else {
        std::cout << "\t" << stepSize;
    }

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
}

xt::xtensor<float, 1> gradientDescentFixedStep(const NetTopology &topo, xt::xtensor<float, 1> initial, DescentModel model, int nbSteps, float stepSize, float momentum, float smoothing) {
    assert (stepSize > 0.0);
    assert (momentum >= 0.0 && momentum < 1.0);
    assert (initial.size() == topo.nbCells());
    reportHeader();
    xt::xtensor<float, 1> x = initial;
    xt::xtensor<float, 1> y = initial;
    for (int i = 0; i < nbSteps; ++i) {
        report(topo, x, model, i, stepSize, momentum, smoothing, false);
        xt::xtensor<float, 1> z = y + stepSize * descentDirection(topo, x, model, smoothing);
        y = z + momentum * (z - x);
        x = z;
    }
    report(topo, x, model, nbSteps, stepSize, momentum, smoothing, false);
    reportFooter();
    return x;
}

xt::xtensor<float, 1> gradientDescent(const NetTopology &topo, xt::xtensor<float, 1> initial, DescentModel model, int nbSteps, float momentum, float smoothing, float initialStepSize, float stepVariation) {
    assert (initialStepSize > 0.0);
    assert (momentum >= 0.0 && momentum < 1.0);
    assert (stepVariation < 1.0);
    assert (initial.size() == topo.nbCells());
    reportHeader();
    xt::xtensor<float, 1> x = initial;
    xt::xtensor<float, 1> y = initial;
    float currentStep = initialStepSize;
    for (int i = 0; i < nbSteps; ++i) {
        report(topo, x, model, i, initialStepSize, momentum, smoothing, true);
        xt::xtensor<float, 1> direction = descentDirection(topo, x, model, smoothing);
        lineSearch(x, direction, currentStep, stepVariation, [&](xt::xtensor<float, 1> w) -> float { return computeValue(topo, w, model, smoothing); });
        xt::xtensor<float, 1> z = y + currentStep * direction;
        y = z + momentum * (z - x);
        x = z;
    }
    report(topo, x, model, nbSteps, initialStepSize, momentum, smoothing, true);
    reportFooter();
    return x;
}
