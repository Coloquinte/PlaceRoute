
#include "coloquinte.hpp"
#include "place_global/wirelength_model.hpp"
#include "place_global/gradient_descent.hpp"
#include "place_global/density_legalizer.hpp"

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
    char *cell_flip_y
) {
    Circuit circuit = Circuit::create_ispd(
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
        cell_flip_y
    );
    std::cout << "Placing circuit with " << circuit.nbCells() << " cells, "
              << circuit.nbNets() << " nets and "
              << circuit.nbPins() << " pins."
              << std::endl;

    auto xtopo = NetWirelength::xTopology(circuit);
    auto ytopo = NetWirelength::yTopology(circuit);

    auto xplace = xtopo.starSolve();
    auto yplace = ytopo.starSolve();
    std::cout << "Initial HPWL: " << xtopo.valueHPWL(xplace) + ytopo.valueHPWL(yplace) << std::endl;

    std::vector<int> cellDemand;
    std::vector<float> cellTargetX;
    std::vector<float> cellTargetY;
    for (int i = 0; i < circuit.nbCells(); ++i) {
        cellDemand.push_back(circuit.isFixed(i) ? 0 : circuit.getArea(i));
        cellTargetX.push_back(xplace[i]);
        cellTargetY.push_back(yplace[i]);
    }
    DensityLegalizer leg(Rectangle(0, 10000, 0, 10000), xtopo.nbCells());
    leg.updateBins(10, 10);
    leg.updateCellDemand(cellDemand);
    leg.updateCellTargetX(cellTargetX);
    leg.updateCellTargetY(cellTargetY);
    leg.assign();
    std::cout << "Overflow: " << leg.meanOverflow() << std::endl;
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
    Circuit circuit = Circuit::create_ispd(
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
        cell_flip_y
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

Circuit Circuit::create_ispd(
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
    char *cell_flip_y
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
    ret.check();
    return ret;
}

void Circuit::check() const {
}


