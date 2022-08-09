#include "place_detailed/detailed_placement.hpp"

#include <stdexcept>

namespace coloquinte {
void DetailedPlacement::check() const {
  if (rows_.size() != nbRows()) throw std::runtime_error("Row size mismatch");
  if (rowFirstCell_.size() != nbRows())
    throw std::runtime_error("Row size mismatch");
  if (rowLastCell_.size() != nbRows())
    throw std::runtime_error("Row size mismatch");
  if (cellWidth_.size() != nbCells())
    throw std::runtime_error("Cell size mismatch");
  if (cellPred_.size() != nbCells())
    throw std::runtime_error("Cell size mismatch");
  if (cellNext_.size() != nbCells())
    throw std::runtime_error("Cell size mismatch");
  if (cellRow_.size() != nbCells())
    throw std::runtime_error("Cell size mismatch");
  if (cellX_.size() != nbCells())
    throw std::runtime_error("Cell size mismatch");
  if (cellY_.size() != nbCells())
    throw std::runtime_error("Cell size mismatch");
  for (int i = 0; i < nbRows(); ++i) {
    int fc = rowFirstCell(i);
    int lc = rowLastCell(i);
    if ((lc == -1) != (fc == -1)) {
      throw std::runtime_error("Inconcistency between first and last cell");
    }
    if (fc == -1) continue;
    if (cellRow(fc) != i) {
      throw std::runtime_error("Inconsistency in the first row cell");
    }
    if (cellPred(fc) != -1) {
      throw std::runtime_error("Inconsistency in the first row cell");
    }
    if (cellRow(lc) != i) {
      throw std::runtime_error("Inconsistency in the last row cell");
    }
    if (cellNext(lc) != -1) {
      throw std::runtime_error("Inconsistency in the last row cell");
    }
  }
  for (int i = 0; i < nbCells(); ++i) {
    int pc = cellPred(i);
    if (pc != -1) {
      if (cellRow(pc) != cellRow(i)) {
        throw std::runtime_error("Row inconsistency with the predecessor");
      }
    } else {
      if (rowFirstCell(i) != cellRow(i)) {
        throw std::runtime_error("Inconsistent first row cell");
      }
    }
    int nc = cellNext(i);
    if (nc != -1) {
      if (cellRow(nc) != cellRow(i)) {
        throw std::runtime_error("Row inconsistency with the successor");
      }
    } else {
      if (rowLastCell(i) != cellRow(i)) {
        throw std::runtime_error("Inconsistent last row cell");
      }
    }
  }
}
}  // namespace coloquinte