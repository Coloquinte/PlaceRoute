
#include "coloquinte.hpp"

#include <iostream>

#include "place_detailed/place_detailed.hpp"
#include "place_global/density_legalizer.hpp"
#include "place_global/net_model.hpp"
#include "place_global/place_global.hpp"

void place(Circuit &circuit, int effort) {
  std::cout << "Placing circuit with " << circuit.nbCells() << " cells, "
            << circuit.nbNets() << " nets and " << circuit.nbPins() << " pins."
            << std::endl;
  GlobalPlacer::place(circuit, effort);
  std::cout << "Wirelength after global placement: " << circuit.hpwl()
            << std::endl;
  DetailedPlacer::place(circuit, effort);
  std::cout << "Wirelength after detailed placement: " << circuit.hpwl()
            << std::endl;
}

extern "C" {
void place_ispd(int nb_cells, int nb_nets, int *cell_widths, int *cell_heights,
                char *cell_fixed, int *net_limits, int *pin_cells,
                int *pin_x_offsets, int *pin_y_offsets, int *cell_x,
                int *cell_y, int *cell_orientation, int nb_rows, int *row_min_x,
                int *row_max_x, int *row_min_y, int *row_max_y, int effort) {
  Circuit circuit = Circuit::createIspd(
      nb_cells, nb_nets, cell_widths, cell_heights, cell_fixed, net_limits,
      pin_cells, pin_x_offsets, pin_y_offsets, cell_x, cell_y, cell_orientation,
      nb_rows, row_min_x, row_max_x, row_min_y, row_max_y);
  place(circuit, effort);
  for (int i = 0; i < nb_cells; ++i) {
    cell_x[i] = circuit.cellX[i];
    cell_y[i] = circuit.cellY[i];
  }
}
}
