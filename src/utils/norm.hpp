#pragma once

#include <cmath>
#include <stdexcept>

#include "coloquinte.hpp"

namespace coloquinte {

template <typename T>
inline T computeNorm(T x, T y, LegalizationModel leg) {
  T z;
  switch (leg) {
    case LegalizationModel::L1:
      return std::abs(x) + std::abs(y);
    case LegalizationModel::L2:
      return std::sqrt(x * x + y * y);
    case LegalizationModel::LInf:
      return std::max(std::abs(x), std::abs(y));
    case LegalizationModel::L1Squared:
      z = std::abs(x) + std::abs(y);
      return z * z;
    case LegalizationModel::L2Squared:
      return x * x + y * y;
    case LegalizationModel::LInfSquared:
      z = std::max(std::abs(x), std::abs(y));
      return z * z;
    default:
      throw std::runtime_error("Unknown legalization model");
  }
}

/**
 * @brief Compute the norm of the 2D vector with the given cost model
 */
inline float norm(float x, float y, LegalizationModel leg) {
  return computeNorm<float>(x, y, leg);
}

/**
 * @brief Compute the norm of the 2D vector with the given cost model
 */
inline long long norm(int x, int y, LegalizationModel leg) {
  return computeNorm<long long>(x, y, leg);
}
}  // namespace coloquinte