
#include "coloquinte.hpp"
#include "place_global/wirelength_model.hpp"
#include "place_global/gradient_descent.hpp"
#include "place_global/density_legalizer.hpp"
#include "place_global/place_global.hpp"
#include "place_detailed/place_detailed.hpp"

#include <xtensor/xio.hpp>
#include <xtensor/xrandom.hpp>

#include <iostream>

extern "C" {
void place_ispd(
    int nb_cells,
    int nb_nets,
    int *cell_widths,
    int *cell_heights,
    char *cell_fixed,
    int *net_limits,
    int *pin_cells,
    int *pin_x_offsets,
    int *pin_y_offsets,
    int *cell_x,
    int *cell_y,
    char *cell_flip_x,
    char *cell_flip_y,
    int nb_rows,
    int *row_min_x,
    int *row_max_x,
    int *row_min_y,
    int *row_max_y
) {
    Circuit circuit = Circuit::createIspd(
        nb_cells,
        nb_nets,
        cell_widths,
        cell_heights,
        cell_fixed,
        net_limits,
        pin_cells,
        pin_x_offsets,
        pin_y_offsets,
        cell_x,
        cell_y,
        cell_flip_x,
        cell_flip_y,
        nb_rows,
        row_min_x,
        row_min_y,
        row_max_x,
        row_max_y
    );
    std::cout << "Placing circuit with " << circuit.nbCells() << " cells, "
              << circuit.nbNets() << " nets and "
              << circuit.nbPins() << " pins."
              << std::endl;
    GlobalPlacer::place(circuit);
    DetailedPlacer::place(circuit);
}

void benchmark_quadratic_models(
    int nb_cells,
    int nb_nets,
    int *cell_widths,
    int *cell_heights,
    char *cell_fixed,
    int *net_limits,
    int *pin_cells,
    int *pin_x_offsets,
    int *pin_y_offsets,
    int *cell_x,
    int *cell_y,
    char *cell_flip_x,
    char *cell_flip_y,
    int model_type,
    int nb_steps,
    float epsilon,
    float relaxation
) {
    Circuit circuit = Circuit::createIspd(
        nb_cells,
        nb_nets,
        cell_widths,
        cell_heights,
        cell_fixed,
        net_limits,
        pin_cells,
        pin_x_offsets,
        pin_y_offsets,
        cell_x,
        cell_y,
        cell_flip_x,
        cell_flip_y,
        0, NULL, NULL, NULL, NULL
    );
    auto xtopo = NetWirelength::xTopology(circuit);

    auto starPlace = xtopo.starSolve();
    std::cout << "INIT\t" << epsilon << "\t" << relaxation << "\t" << 0 << "\t" << xtopo.valueHPWL(starPlace) << std::endl;
    for (int i = 0; i < nb_steps; ++i) {
        if (model_type == 0) {
            starPlace = xtopo.starSolve(starPlace, epsilon, relaxation, false);
            std::cout << "STAR\t" << epsilon << "\t" << relaxation << "\t" << i + 1 << "\t" << xtopo.valueHPWL(starPlace) << std::endl;
        }
        else if (model_type == 1) {
            starPlace = xtopo.starSolve(starPlace, epsilon, relaxation, true);
            std::cout << "BSTAR\t" << epsilon << "\t" << relaxation << "\t" << i + 1 << "\t" << xtopo.valueHPWL(starPlace) << std::endl;
        }
        else {
            starPlace = xtopo.b2bSolve(starPlace, epsilon);
            std::cout << "B2B\t" << epsilon << "\t" << relaxation << "\t" << i + 1 << "\t" << xtopo.valueHPWL(starPlace) << std::endl;
        }
    }
}
}

Circuit Circuit::createIspd(
    int nb_cells,
    int nb_nets,
    int *cell_widths,
    int *cell_heights,
    char *cell_fixed,
    int *net_limits,
    int *pin_cells,
    int *pin_x_offsets,
    int *pin_y_offsets,
    int *cell_x,
    int *cell_y,
    char *cell_flip_x,
    char *cell_flip_y,
    int nb_rows,
    int *row_min_x,
    int *row_max_x,
    int *row_min_y,
    int *row_max_y
) {
    Circuit ret;
    ret.cellWidths.assign(cell_widths, cell_widths + nb_cells);
    ret.cellHeights.assign(cell_heights, cell_heights + nb_cells);
    ret.cellFixed.assign(cell_fixed, cell_fixed + nb_cells);
    ret.netLimits.assign(net_limits, net_limits + nb_nets + 1);
    int nbPins = ret.netLimits.back();
    ret.pinCells.assign(pin_cells, pin_cells + nbPins);
    ret.pinXOffsets.assign(pin_x_offsets, pin_x_offsets + nbPins);
    ret.pinYOffsets.assign(pin_y_offsets, pin_y_offsets + nbPins);
    ret.cellX.assign(cell_x, cell_x + nb_cells);
    ret.cellY.assign(cell_y, cell_y + nb_cells);
    ret.cellFlipX.assign(cell_flip_x, cell_flip_x + nb_cells);
    ret.cellFlipY.assign(cell_flip_y, cell_flip_y + nb_cells);
    for (int i = 0; i < nb_rows; ++i) {
        ret.rows.emplace_back(row_min_x[i], row_max_x[i], row_min_y[i], row_max_y[i]);
    }
    ret.check();
    return ret;
}

void Circuit::check() const {
}


