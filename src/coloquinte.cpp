
#include "coloquinte.hpp"

#include <iostream>

#include "place_detailed/place_detailed.hpp"
#include "place_global/density_legalizer.hpp"
#include "place_global/net_model.hpp"
#include "place_global/place_global.hpp"

extern "C" {
void place_ispd(int nb_cells, int nb_nets, int *cell_widths, int *cell_heights,
                char *cell_fixed, int *net_limits, int *pin_cells,
                int *pin_x_offsets, int *pin_y_offsets, int *cell_x,
                int *cell_y, int *cell_orientation, int nb_rows, int *row_min_x,
                int *row_max_x, int *row_min_y, int *row_max_y) {
  Circuit circuit = Circuit::createIspd(
      nb_cells, nb_nets, cell_widths, cell_heights, cell_fixed, net_limits,
      pin_cells, pin_x_offsets, pin_y_offsets, cell_x, cell_y, cell_orientation,
      nb_rows, row_min_x, row_max_x, row_min_y, row_max_y);
  std::cout << "Placing circuit with " << circuit.nbCells() << " cells, "
            << circuit.nbNets() << " nets and " << circuit.nbPins() << " pins."
            << std::endl;
  GlobalPlacer::place(circuit);
  DetailedPlacer::place(circuit);
}

void benchmark_quadratic_models(int nb_cells, int nb_nets, int *cell_widths,
                                int *cell_heights, char *cell_fixed,
                                int *net_limits, int *pin_cells,
                                int *pin_x_offsets, int *pin_y_offsets,
                                int *cell_x, int *cell_y, int *cell_orientation,
                                int model_type, int nb_steps, float epsilon,
                                float relaxation) {
  Circuit circuit = Circuit::createIspd(
      nb_cells, nb_nets, cell_widths, cell_heights, cell_fixed, net_limits,
      pin_cells, pin_x_offsets, pin_y_offsets, cell_x, cell_y, cell_orientation,
      0, NULL, NULL, NULL, NULL);
  auto xtopo = NetModel::xTopology(circuit);

  auto starPlace = xtopo.solveStar();
  std::cout << "INIT\t" << epsilon << "\t" << relaxation << "\t" << 0 << "\t"
            << xtopo.value(starPlace) << std::endl;
  for (int i = 0; i < nb_steps; ++i) {
    if (model_type == 0) {
      starPlace = xtopo.solveStar(starPlace, epsilon);
      std::cout << "STAR\t" << epsilon << "\t" << relaxation << "\t" << i + 1
                << "\t" << xtopo.value(starPlace) << std::endl;
    } else if (model_type == 1) {
      starPlace = xtopo.solveStar(starPlace, epsilon);
      std::cout << "BSTAR\t" << epsilon << "\t" << relaxation << "\t" << i + 1
                << "\t" << xtopo.value(starPlace) << std::endl;
    } else {
      starPlace = xtopo.solveB2B(starPlace, epsilon);
      std::cout << "B2B\t" << epsilon << "\t" << relaxation << "\t" << i + 1
                << "\t" << xtopo.value(starPlace) << std::endl;
    }
  }
}
}
