#pragma once

#include "place_global/wirelength_model.hpp"

enum class DescentModel { HPWL, LSE, WA, Proximal };
enum class DescentType {
  GradientDescent,
  ConjugateGradientDescent,
  NesterovAcceleratedGradientDescent
};
enum class LineSearchType { FixedStep, Armijo };

xt::xtensor<float, 1> gradientDescentFixedStep(const NetWirelength &topo,
                                               xt::xtensor<float, 1> initial,
                                               DescentModel model, int nbSteps,
                                               float stepSize,
                                               float momentum = 0.0,
                                               float smoothing = 1.0);

xt::xtensor<float, 1> gradientDescent(const NetWirelength &topo,
                                      xt::xtensor<float, 1> initial,
                                      DescentModel model, int nbSteps,
                                      float momentum = 0.0,
                                      float smoothing = 1.0,
                                      float initialStepSize = 1.0,
                                      float stepVariation = 0.9);

xt::xtensor<float, 1> conjugateGradientDescent(const NetWirelength &topo,
                                               xt::xtensor<float, 1> initial,
                                               DescentModel model, int nbSteps,
                                               float smoothing = 1.0,
                                               float initialStepSize = 1.0,
                                               float stepVariation = 0.9);

xt::xtensor<float, 1> nesterovGradientDescent(const NetWirelength &topo,
                                              xt::xtensor<float, 1> initial,
                                              DescentModel model, int nbSteps,
                                              float smoothing = 1.0,
                                              float initialStepSize = 1.0,
                                              float stepVariation = 0.9);
