
#include "coloquinte.hpp"

#include <cmath>

float norm(float x, float y, LegalizationModel leg) {
    switch (leg) {
        case LegalizationModel::L1:
            return std::abs(x) + std::abs(y);
        case LegalizationModel::L2:
            return std::sqrt(x * x + y * y);
        case LegalizationModel::L2Squared:
            return x * x + y * y;
        default:
            return std::max(std::abs(x), std::abs(y));
    }
}

long long norm(int x, int y, LegalizationModel leg) {
    long long lx = x;
    long long ly = y;
    switch (leg) {
        case LegalizationModel::L1:
            return std::abs(lx) + std::abs(ly);
        case LegalizationModel::L2:
            return std::sqrt(lx * lx + ly * ly);
        case LegalizationModel::L2Squared:
            return lx * lx + ly * ly;
        default:
            return std::max(std::abs(lx), std::abs(ly));
    }
}
