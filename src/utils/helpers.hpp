#pragma once

#include <cassert>
#include <vector>

namespace coloquinte {
/**
 * @brief Compute the limits when subdividing an interval
 */
inline std::vector<int> computeSubdivisions(int min, int max, int number) {
  assert(number >= 1);
  assert(max >= min);
  std::vector<int> ret;
  for (int i = 0; i < number + 1; ++i) {
    ret.push_back(min + (i * (max - min) / number));
  }
  assert(ret.size() == number + 1);
  assert(ret.front() == min);
  assert(ret.back() == max);
  return ret;
}
}  // namespace coloquinte
