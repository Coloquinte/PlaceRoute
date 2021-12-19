#pragma once

#include "place_global/topology.hpp"

enum class DescentModel { HPWL, LSE, WA, Proximal };

xt::xtensor<float, 1> gradientDescent(const NetTopology &topo, xt::xtensor<float, 1> initial, DescentModel model, int nbSteps, float stepSize, float momentum=0.0, float smoothing=1.0);
