#pragma once

#include "coloquinte.hpp"

/**
 * @brief Main class for detailed placement
 */
class DetailedPlacer {
 public:
  /**
   * @brief Run detailed placement on the circuit representation
   *
   * @param circuit The circuit to be modified
   * @param effort Effort level, between 0 and 9
   */
  static void place(Circuit &circuit, int effort);
};
