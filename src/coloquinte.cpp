
#include "coloquinte.hpp"
#include "place_global/topology.hpp"
#include "place_global/descent.hpp"

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
    auto xtopo = NetTopology::xTopology(circuit);

    xt::random::seed(0);
    //std::cout << "Initial X HPWL: " << xtopo.valueHPWL(xplace) << std::endl;
    //std::cout << "Initial X LSE: " << xtopo.valueLSE(xplace, 1.0) << std::endl;
    //std::cout << "Initial X WA: " << xtopo.valueWA(xplace, 1.0) << std::endl;
    //std::cout << "Gradient X HPWL: " << xtopo.gradHPWL(xplace) << std::endl;
    //std::cout << "Gradient X LSE: " << xtopo.gradLSE(xplace, 1.0) << std::endl;
    //std::cout << "Gradient X WA: " << xtopo.gradWA(xplace, 1.0) << std::endl;
    //std::cout << "Proximal step X: " << xtopo.proximalStep(xplace, 1.0) << std::endl;
    /*
    int nbSteps = 20;
    float epsilon = 1.0;
    for (float relaxation : {1.0, 0.5, 0.1, 0.0}) {
        auto starPlace = xtopo.starSolve();
        std::cout << "INIT\t" << epsilon << "\t" << relaxation << "\t" << 0 << "\t" << xtopo.valueHPWL(starPlace) << std::endl;
        for (int i = 0; i < nbSteps; ++i) {
            starPlace = xtopo.starSolve(starPlace, epsilon, relaxation, true);
            std::cout << "STAR\t" << epsilon << "\t" << relaxation << "\t" << i + 1 << "\t" << xtopo.valueHPWL(starPlace) << std::endl;
        }
    }
    */
    auto starPlace = xtopo.starSolve();
    int nbSteps = 100;
    float stepSize = 1.0;
    float momentum = 0.9;
    //DescentModel model = DescentModel::Proximal;
    //for (float momentum : {0.0, 0.9, 0.95, 0.975, 0.9875}) {
    for (float smoothing : {1.0, 2.0, 4.0}) {
        for (DescentModel model : {DescentModel::WA, DescentModel::Proximal}) {
            conjugateGradientDescent(xtopo, starPlace, model, nbSteps, smoothing, stepSize);
            gradientDescent(xtopo, starPlace, model, nbSteps, momentum, smoothing, stepSize);
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


