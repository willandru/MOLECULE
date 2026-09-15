#pragma once

#include <vector>
#include <cstddef>

#include <glm/glm.hpp>

#include "DFTGrid.h"
#include "DFTDensity.h"

class PotentialEnergy
{
public:

    PotentialEnergy();

    explicit PotentialEnergy(
        const DFTGrid& grid
    );

    void initialize(
        const DFTGrid& grid
    );

    double calculateElectronNuclearEnergy(
        const DFTDensity& density,
        const std::vector<double>& externalPotential
    ) const;

    double calculateHartreeEnergy(
        const DFTDensity& density,
        const std::vector<double>& hartreePotential
    ) const;

    double calculateExchangeEnergy(
        const DFTDensity& density
    ) const;

    double calculateNuclearRepulsionEnergy(
        const std::vector<int>& nuclearCharges,
        const std::vector<glm::dvec3>& nuclearPositions
    ) const;

private:

    const DFTGrid* grid;
};