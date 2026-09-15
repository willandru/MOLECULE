#pragma once

#include <vector>
#include <cstddef>

#include <glm/glm.hpp>

#include "DFTGrid.h"
#include "DFTDensity.h"

class NuclearForces
{
public:

    NuclearForces();

    explicit NuclearForces(
        const DFTGrid& grid
    );

    void initialize(
        const DFTGrid& grid
    );

    std::vector<glm::dvec3> calculate(
        const DFTDensity& density,
        const std::vector<int>& nuclearCharges,
        const std::vector<glm::dvec3>& nuclearPositions,
        double softening
    ) const;

private:

    const DFTGrid* grid;
};