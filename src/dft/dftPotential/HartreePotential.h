#pragma once

#include <vector>
#include <cstddef>

#include <glm/glm.hpp>

#include "DFTGrid.h"
#include "DFTDensity.h"
#include "MultigridSolver.h"

class HartreePotential
{
public:

    HartreePotential();

    explicit HartreePotential(
        const DFTGrid& grid
    );

    void initialize(
        const DFTGrid& grid
    );

    void calculate(
        const DFTDensity& density,
        int iterations = 20
    );

    double get(
        std::size_t index
    ) const;

    const std::vector<double>&
    getPotential() const;

    double getResidual() const;

    int getCycles() const;

    bool hasConverged() const;

private:

    const DFTGrid* grid;

    std::vector<double> potential;

    MultigridSolver multigridSolver;

    void applyBoundaryConditions(
        const DFTDensity& density
    );

    bool isBoundaryPoint(
        int x,
        int y,
        int z
    ) const;
};