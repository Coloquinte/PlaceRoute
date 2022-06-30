#pragma once

#include "coloquinte.hpp"


class GlobalPlacer {
  public:
    static std::pair<std::vector<float>, std::vector<float> > place(const Circuit &circuit);
};



