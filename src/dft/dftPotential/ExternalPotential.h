#pragma once

#include <vector>
#include <cstddef>

#include <glm/glm.hpp>

#include "DFTGrid.h"

class ExternalPotential
{
public:

    ExternalPotential();

    explicit ExternalPotential(
        const DFTGrid& grid
    );

    void initialize(
        const DFTGrid& grid
    );

    void calculate(
        const std::vector<int>& nuclearCharges,
        const std::vector<glm::dvec3>& nuclearPositions
    );

    double get(
        std::size_t index
    ) const;

    const std::vector<double>& getPotential() const;

    void setSoftening(
        double value
    );

    double getSoftening() const;

private:

    const DFTGrid* grid;

    std::vector<double> potential;

    double softening;
};